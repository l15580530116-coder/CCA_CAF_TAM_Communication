# ============================================================
# 分析目的: 整理 TCGA-CHOL 临床数据，提取生存信息并清洗
# 输入文件: data/clinical/TCGA_CHOL/TCGA_CHOL_clinical.csv
# 输出文件:
#   data/clinical/TCGA_CHOL/TCGA_CHOL_survival_clean.csv
#   tables/TCGA_CHOL_clinical_summary.csv
# 主要方法: 基于 days_to_death / days_to_last_follow_up 构建 OS
#          自动识别列名（兼容不同 TCGA 数据版本）
# ============================================================

# ---- 加载包 ----
library(data.table)

# ---- 读取原始临床数据 ----
clin_raw <- fread("data/clinical/TCGA_CHOL/TCGA_CHOL_clinical.csv", data.table = FALSE)
n_raw <- nrow(clin_raw)

cat(sprintf("原始临床记录数: %d\n", n_raw))

# ---- Step 1: 自动识别 patient ID 列 ----
id_col <- NULL
for (col in c("submitter_id", "bcr_patient_barcode", "patient_id", "case_submitter_id")) {
  if (col %in% colnames(clin_raw) && sum(!is.na(clin_raw[[col]])) > 0) {
    id_col <- col
    break
  }
}
if (is.null(id_col)) stop("未找到 patient ID 列")
cat(sprintf("Patient ID 列: %s\n", id_col))

# ---- Step 2: 提取核心生存变量 ----
# 自动识别变量名
find_col <- function(candidates, data) {
  for (c in candidates) {
    if (c %in% colnames(data)) return(c)
  }
  return(NULL)
}

death_col  <- find_col(c("days_to_death"), clin_raw)
fup_col    <- find_col(c("days_to_last_follow_up"), clin_raw)
vital_col  <- find_col(c("vital_status"), clin_raw)

if (is.null(death_col) || is.null(fup_col)) {
  stop("缺少 days_to_death 或 days_to_last_follow_up 列")
}

cat(sprintf("生存变量: death=%s, follow_up=%s, vital=%s\n", death_col, fup_col, vital_col %||% "N/A"))

# ---- Step 3: 构建 OS_time 和 OS_status ----
# OS_status: 1=死亡事件, 0=删失
clin_os <- data.frame(
  patient_id = clin_raw[[id_col]],
  OS_time    = NA_real_,
  OS_status  = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(clin_raw))) {
  dtd <- clin_raw[[death_col]][i]
  flf <- clin_raw[[fup_col]][i]

  if (!is.na(dtd) && dtd > 0) {
    clin_os$OS_time[i]   <- dtd
    clin_os$OS_status[i] <- 1
  } else if (!is.na(flf) && flf > 0) {
    clin_os$OS_time[i]   <- flf
    clin_os$OS_status[i] <- 0
  }
}

n_dead  <- sum(clin_os$OS_status == 1, na.rm = TRUE)
n_alive <- sum(clin_os$OS_status == 0, na.rm = TRUE)
n_miss  <- sum(is.na(clin_os$OS_time))
n_remove <- n_raw - (n_dead + n_alive)

cat(sprintf("\n=== 生存信息构建 ===\n"))
cat(sprintf("  死亡事件 (OS_status=1): %d\n", n_dead))
cat(sprintf("  删失     (OS_status=0): %d\n", n_alive))
cat(sprintf("  无法构建 OS:           %d\n", n_miss))
cat(sprintf("  OS_time 范围:          %.0f - %.0f 天\n",
            min(clin_os$OS_time, na.rm = TRUE), max(clin_os$OS_time, na.rm = TRUE)))

# ---- Step 4: 附加临床协变量 ----
# 定义需要的变量（候选名列表）
var_map <- list(
  age               = c("age_at_index", "age_at_diagnosis", "age"),
  gender            = c("gender"),
  race              = c("race"),
  vital_status      = c("vital_status"),
  stage             = c("ajcc_pathologic_stage", "ajcc_pathologic_t_stage"),
  T_stage           = c("ajcc_pathologic_t"),
  N_stage           = c("ajcc_pathologic_n"),
  M_stage           = c("ajcc_pathologic_m"),
  grade             = c("tumor_grade", "grade"),
  residual_disease  = c("residual_disease", "residual_tumor")
)

# 输出列名（友好名称）
out_names <- c("patient_id", "OS_time", "OS_status")

for (v in names(var_map)) {
  col_found <- find_col(var_map[[v]], clin_raw)
  if (!is.null(col_found)) {
    clin_os[[col_found]] <- clin_raw[[col_found]]
    out_names <- c(out_names, col_found)
    cat(sprintf("  临床变量 %-12s -> %s\n", v, col_found))
  } else {
    cat(sprintf("  WARNING: 未找到变量 %s，已跳过\n", v))
  }
}

# 确保数据框包含所有要输出的列
clin_os <- clin_os[, out_names, drop = FALSE]

# ---- Step 5: 过滤无效记录 ----
# 删除 OS_time 缺失 或 <= 0
before_filter <- nrow(clin_os)
valid_os <- !is.na(clin_os$OS_time) & clin_os$OS_time > 0
clin_clean <- clin_os[valid_os, , drop = FALSE]
after_filter <- nrow(clin_clean)
removed_os <- before_filter - after_filter

cat(sprintf("\n=== 过滤统计 ===\n"))
cat(sprintf("  过滤前: %d\n", before_filter))
cat(sprintf("  过滤后: %d\n", after_filter))
cat(sprintf("  删除:   %d (OS_time 缺失或 <= 0)\n", removed_os))

# ---- Step 6: 保存清洗后数据 ----
dir.create("data/clinical/TCGA_CHOL", recursive = TRUE, showWarnings = FALSE)

write.csv(clin_clean,
          "data/clinical/TCGA_CHOL/TCGA_CHOL_survival_clean.csv",
          row.names = FALSE)
cat(sprintf("\n已保存: data/clinical/TCGA_CHOL/TCGA_CHOL_survival_clean.csv\n"))

# ---- Step 7: 生成临床概览表 ----
dir.create("tables", recursive = TRUE, showWarnings = FALSE)

summary_rows <- list(
  c("原始临床记录数",            as.character(n_raw)),
  c("清洗后可用生存样本数",       as.character(after_filter)),
  c("死亡事件数 (OS_status=1)",  as.character(n_dead)),
  c("删失数 (OS_status=0)",      as.character(n_alive)),
  c("OS_time 范围 (天)",         sprintf("%.0f - %.0f",
    min(clin_clean$OS_time, na.rm = TRUE),
    max(clin_clean$OS_time, na.rm = TRUE))),
  c("删除样本数 (OS无法构建)",   as.character(removed_os)),
  c("Patient ID 列",             id_col),
  c("输出列名",                  paste(out_names, collapse = ", "))
)

# 添加各临床变量的非缺失样本数
for (v in out_names) {
  if (v != "patient_id" && v != "OS_time" && v != "OS_status") {
    n <- sum(!is.na(clin_clean[[v]]))
    summary_rows <- c(summary_rows, list(c(sprintf("  %s (有效值)", v), sprintf("%d / %d", n, after_filter))))
  }
}

summary_df <- do.call(rbind, summary_rows)
colnames(summary_df) <- c("指标", "值")

write.csv(summary_df,
          "tables/TCGA_CHOL_clinical_summary.csv",
          row.names = FALSE)
cat(sprintf("已保存: tables/TCGA_CHOL_clinical_summary.csv\n"))

# ---- 最终报告 ----
cat(sprintf("\n========================================\n"))
cat(sprintf("TCGA-CHOL 临床数据清洗完成\n"))
cat(sprintf("========================================\n"))
cat(sprintf("可用于生存分析的样本数: %d\n", after_filter))
cat(sprintf("  其中死亡事件: %d\n", n_dead))
cat(sprintf("  其中删失:     %d\n", n_alive))
cat(sprintf("========================================\n"))
