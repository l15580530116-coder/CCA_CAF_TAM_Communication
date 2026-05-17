#!/usr/bin/env python
"""Step 36c: Polish supplementary figures - unified labels, no overlap, consistent layout."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from pathlib import Path
import fitz, csv, io
from PIL import Image

BASE = Path(r"e:\CCA")
SRC = BASE / "figures"
TABLES = BASE / "tables"
OUT_V3 = BASE / "figures/supplementary_final_v3"
OUT_V3.mkdir(parents=True, exist_ok=True)
DPI = 600
FS = 4.0  # panel size

# ═══════════════════════════════════════════════════════════════
# UNIFIED STYLING HELPERS
# ═══════════════════════════════════════════════════════════════

LABEL_FS = 15       # panel label font size
LABEL_X = 0.015     # 1.5% from left
LABEL_Y = 0.970     # 3% from top
NOTE_FS = 7.5       # annotation font size
NOTE_COLOR = "#555555"

def pdf_to_img(path, page_idx=0, dpi=DPI, max_mp=25):
    doc = fitz.open(str(path))
    if len(doc) == 0:
        doc.close()
        raise ValueError("0-page PDF")
    if page_idx >= len(doc):
        page_idx = len(doc) - 1
    page = doc[page_idx]
    pw, ph = page.rect.width, page.rect.height
    d = dpi
    mpix = (pw*d/72) * (ph*d/72) / 1e6
    if mpix > max_mp:
        d = int(dpi * (max_mp/mpix)**0.5)
    zoom = d / 72.0
    mat = fitz.Matrix(zoom, zoom)
    pix = page.get_pixmap(matrix=mat, colorspace=fitz.csRGB)
    img = np.array(Image.open(io.BytesIO(pix.tobytes("png"))))
    doc.close()
    return img

def add_label(ax, letter):
    """Unified panel label: bold, top-left margin, white background."""
    ax.text(LABEL_X, LABEL_Y, letter, transform=ax.transAxes,
            fontsize=LABEL_FS, fontweight="bold", color="black",
            va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.10", facecolor="white",
                      edgecolor="none", alpha=0.88))

def add_page_label(ax, letter):
    """Larger label for single-panel full-page figures."""
    ax.text(0.010, 0.985, letter, transform=ax.transAxes,
            fontsize=18, fontweight="bold", color="black",
            va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                      edgecolor="none", alpha=0.88))

def add_note(ax, text, position="bottom-right"):
    """Unified annotation: small italic grey, white bg, consistent position."""
    if position == "bottom-right":
        x, y, va, ha = 0.982, 0.018, "bottom", "right"
    elif position == "top-right":
        x, y, va, ha = 0.982, 0.970, "top", "right"
    elif position == "top-left":
        x, y, va, ha = 0.018, 0.970, "top", "left"
    else:
        x, y, va, ha = 0.982, 0.018, "bottom", "right"
    ax.text(x, y, text, transform=ax.transAxes, fontsize=NOTE_FS,
            fontstyle="italic", color=NOTE_COLOR, va=va, ha=ha,
            bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                      edgecolor="none", alpha=0.88))

def save_figure(fig, stem, dpi=DPI):
    pdf_path = OUT_V3 / f"{stem}.pdf"
    png_path = OUT_V3 / f"{stem}.png"
    fig.savefig(pdf_path, dpi=dpi, facecolor="white", edgecolor="none",
                bbox_inches="tight")
    fig.savefig(png_path, dpi=dpi, facecolor="white", edgecolor="none",
                bbox_inches="tight")
    print(f"  -> {stem}.pdf ({pdf_path.stat().st_size//1024} KB) + PNG")

def fig_to_img(fig):
    fig.canvas.draw()
    img = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
    w, h = fig.canvas.get_width_height()
    return img.reshape(h, w, 4)[:, :, :3]

def make_bubble(rows, title, fs=(8, 6)):
    srt = sorted(rows, key=lambda r: float(r["prob"]), reverse=True)[:25]
    srt = sorted(srt, key=lambda r: float(r["prob"]))
    labels = [f"{r['ligand']} - {r['receptor']}" for r in srt]
    probs = [float(r["prob"]) for r in srt]
    fig, ax = plt.subplots(1, 1, figsize=fs, facecolor="white")
    sizes = [p * 900 for p in probs]
    ax.scatter(probs, labels, s=sizes, c=probs, cmap="Reds",
               edgecolors="grey", linewidth=0.3, alpha=0.85, zorder=3)
    ax.set_xlabel("Communication probability", fontsize=10)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.tick_params(axis="y", labelsize=8)
    ax.tick_params(axis="x", labelsize=9)
    ax.grid(axis="x", alpha=0.3, lw=0.5)
    ax.set_xlim(0, max(probs) * 1.15)
    cbar = fig.colorbar(ax.collections[0], ax=ax, shrink=0.6, aspect=20, pad=0.02)
    cbar.set_label("Probability", fontsize=8)
    plt.tight_layout()
    return fig

# Load CellChat data
lr_csv = TABLES / "cellchat/GSE138709_CAF_TAM_Epithelial_LR_pairs.csv"
with open(lr_csv, newline="") as f:
    all_lr = list(csv.DictReader(f))
caf_tam_lr = [r for r in all_lr if r["source"]=="CAF_Fibroblast" and r["target"]=="Macrophage_TAM"]
tam_caf_lr = [r for r in all_lr if r["source"]=="Macrophage_TAM" and r["target"]=="CAF_Fibroblast"]

# ═══════════════════════════════════════════════════════════════
def assemble(panels, nrows, ncols, title, stem, notes=None):
    """Assemble multi-panel figure with unified labels and notes."""
    fig = plt.figure(figsize=(ncols*FS, nrows*FS), facecolor="white")
    gs = gridspec.GridSpec(nrows, ncols, figure=fig, wspace=0.10, hspace=0.20,
                           left=0.03, right=0.97, top=0.935, bottom=0.03)
    errs = []
    for idx, (letter, path, pg, desc) in enumerate(panels):
        row, col = divmod(idx, ncols)
        ax = fig.add_subplot(gs[row, col])
        try:
            img = pdf_to_img(path, page_idx=pg)
            ax.imshow(img)
            ax.axis("off")
            add_label(ax, letter)
            if notes and letter in notes:
                add_note(ax, notes[letter], position="bottom-right")
            print(f"  {letter}: OK - {desc}")
        except Exception as e:
            ax.text(0.5, 0.5, f"Panel {letter}\nNOT AVAILABLE",
                    ha="center", va="center", fontsize=8, color="red")
            ax.axis("off")
            add_label(ax, letter)
            errs.append(f"{letter}: {e}")
            print(f"  {letter}: ERROR - {e}")
    for idx in range(len(panels), nrows*ncols):
        ax = fig.add_subplot(gs[idx//ncols, idx%ncols]); ax.axis("off")
    fig.suptitle(title, fontsize=11, fontweight="bold", y=0.985)
    save_figure(fig, stem)
    plt.close(fig)
    return errs

def full_page(path, pg, letter, title, stem, note=None):
    """Single-panel full-page figure with label in margin."""
    fig = plt.figure(figsize=(10, 8.2), facecolor="white")
    ax = fig.add_axes([0.03, 0.04, 0.94, 0.90])  # leave room for label
    img = pdf_to_img(path, page_idx=pg)
    ax.imshow(img); ax.axis("off")
    # Label in figure margin, not on the image
    fig.text(0.015, 0.975, letter, fontsize=18, fontweight="bold",
             color="black", va="top", ha="left",
             bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                       edgecolor="none", alpha=0.88))
    if note:
        fig.text(0.982, 0.018, note, fontsize=NOTE_FS, fontstyle="italic",
                 color=NOTE_COLOR, va="bottom", ha="right",
                 bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                           edgecolor="none", alpha=0.88))
    fig.suptitle(title, fontsize=11, fontweight="bold", y=0.995)
    save_figure(fig, stem)
    plt.close(fig)

# ═══════════════════════════════════════════════════════════════
print("="*60)
print("S1")
assemble([
    ("A", SRC/"DEG/Fig2A_TCGA_PCA_tumor_normal.pdf", 0, "TCGA PCA"),
    ("B", SRC/"DEG/Fig2D_GSE107943_PCA_tumor_normal.pdf", 0, "GSE PCA"),
    ("C", SRC/"single_cell/Fig7A_GSE138709_UMAP_clusters.pdf", 0, "scRNA UMAP"),
    ("D", SRC/"GSE26566_validation/FigS_GSE26566_hub_gene_heatmap.pdf", 0, "GSE26566 heatmap"),
], 2, 2, "Figure S1. Data preprocessing and quality control overview",
    "FigureS1_preprocessing_QC_v3")

# ═══════════════════════════════════════════════════════════════
print("\nS2")
full_page(SRC/"DEG/Fig2C_TCGA_top_IM_CAF_TAM_DEG_heatmap.pdf", 0,
          "A", "Figure S2A. TCGA-CHOL: Top IM/CAF/TAM DEG heatmap",
          "FigureS2_full_DEG_plots_v3_page1")
full_page(SRC/"DEG/Fig2G_validated_IM_CAF_TAM_DEGs_heatmap.pdf", 0,
          "B", "Figure S2B. GSE107943: Validated IM/CAF/TAM DEG heatmap",
          "FigureS2_full_DEG_plots_v3_page2")

# ═══════════════════════════════════════════════════════════════
print("\nS3")
assemble([
    ("A", SRC/"gsva/Fig4A_TCGA_score_boxplot_tumor_normal.pdf", 0, "TCGA boxplot"),
    ("B", SRC/"gsva/Fig4B_GSE107943_score_boxplot_tumor_normal.pdf", 0, "GSE boxplot"),
], 1, 2, "Figure S3A. GSVA ssGSEA pathway scores (tumor vs normal)",
    "FigureS3_GSVA_score_survival_v3_page1")
assemble([
    ("C", SRC/"gsva/Fig4E_TCGA_KM_high_low_CAF_score.pdf", 1, "TCGA CAF KM"),
    ("D", SRC/"gsva/Fig4F_TCGA_KM_high_low_IM_CAF_TAM_score.pdf", 1, "TCGA IM_CAF_TAM KM"),
    ("E", SRC/"gsva/Fig4G_GSE107943_KM_high_low_CAF_score.pdf", 1, "GSE CAF KM"),
], 1, 3, "Figure S3B. Exploratory survival analyses",
    "FigureS3_GSVA_score_survival_v3_page2",
    notes={"C": "Exploratory", "D": "Exploratory", "E": "Exploratory"})

# ═══════════════════════════════════════════════════════════════
print("\nS4")
assemble([
    ("A", SRC/"subtyping/Fig5A_TCGA_consensus_CDF.pdf", 0, "Consensus CDF"),
    ("B", SRC/"subtyping/Fig5B_TCGA_consensus_matrix_k2.pdf", 0, "k=2 matrix"),
    ("C", SRC/"subtyping/Fig5D_TCGA_subtype_score_boxplot.pdf", 0, "Score boxplot"),
    ("D", SRC/"subtyping/Fig5F_GSE107943_subtype_heatmap.pdf", 0, "GSE heatmap"),
], 2, 2, "Figure S4A. Molecular subtyping diagnostics (1 of 2)",
    "FigureS4_molecular_subtyping_diagnostics_v3_page1")
assemble([
    ("E", SRC/"subtyping_diagnosis/FigS_original_subtype_gene_direction_heatmap.pdf", 0, "Direction heatmap"),
    ("F", SRC/"subtyping_diagnosis/FigS_strategy_A_up_genes_heatmap_TCGA.pdf", 0, "Strategy A"),
], 1, 2, "Figure S4B. Subtyping strategy diagnostics (2 of 2)",
    "FigureS4_molecular_subtyping_diagnostics_v3_page2")

# ═══════════════════════════════════════════════════════════════
print("\nS5")
assemble([
    ("A", SRC/"immune/Fig6B_GSE107943_immune_score_boxplot_tumor_normal.pdf", 0, "Immune boxplot"),
    ("B", SRC/"immune/Fig6D_GSE107943_aggressive_immune_correlation_heatmap.pdf", 0, "Immune cor"),
    ("C", SRC/"immune/Fig6F_GSE107943_checkpoint_correlation_heatmap.pdf", 0, "Checkpoint cor"),
    ("D", SRC/"immune/Fig6H_aggressive_vs_macrophage_scatter_GSE107943.pdf", 0, "Macro scatter"),
    ("E", SRC/"immune/Fig6J_aggressive_vs_checkpoint_scatter_GSE107943.pdf", 0, "Checkpoint scatter"),
], 2, 3, "Figure S5. GSE107943 immune infiltration and checkpoint correlation analyses",
    "FigureS5_full_immune_analysis_v3")

# ═══════════════════════════════════════════════════════════════
print("\nS6")
full_page(SRC/"single_cell/Fig7C_GSE138709_marker_dotplot.pdf", 0,
          "A", "Figure S6A. Cell type marker dotplot",
          "FigureS6_single_cell_supplementary_v3_page1")
assemble([
    ("B", SRC/"single_cell/Fig7E_GSE138709_CAF_TAM_score_UMAP.pdf", 0, "CAF/TAM UMAP"),
    ("C", SRC/"single_cell/Fig7D_GSE138709_aggressive_score_UMAP.pdf", 0, "Aggressive UMAP"),
], 1, 2, "Figure S6B. Single-cell module score projections",
    "FigureS6_single_cell_supplementary_v3_page2")
full_page(SRC/"single_cell/Fig7G_GSE138709_key_genes_featureplot.pdf", 0,
          "D", "Figure S6C. Key gene feature plots",
          "FigureS6_single_cell_supplementary_v3_page3")

# ═══════════════════════════════════════════════════════════════
print("\nS7 — regenerating C/D bubbles for clean layout")
fig_c = make_bubble(caf_tam_lr, "CAF >> TAM inferred communication", fs=(8, 6.5))
img_c = fig_to_img(fig_c); plt.close(fig_c)
fig_d = make_bubble(tam_caf_lr, "TAM >> CAF inferred communication", fs=(8, 6.5))
img_d = fig_to_img(fig_d); plt.close(fig_d)

fig7 = plt.figure(figsize=(10, 8.5), facecolor="white")
gs7 = gridspec.GridSpec(2, 2, figure=fig7, wspace=0.08, hspace=0.18,
                        left=0.03, right=0.97, top=0.935, bottom=0.03)
s7_data = [
    ("A", pdf_to_img(SRC/"cellchat/Fig8A_CellChat_overall_interaction_number.pdf"), None),
    ("B", pdf_to_img(SRC/"cellchat/Fig8B_CellChat_overall_interaction_weight.pdf"), None),
    ("C", img_c, "Inferred communication"),
    ("D", img_d, "Inferred communication"),
]
for idx, (letter, img, note) in enumerate(s7_data):
    ax = fig7.add_subplot(gs7[idx//2, idx%2])
    ax.imshow(img); ax.axis("off")
    add_label(ax, letter)
    if note:
        add_note(ax, note, position="bottom-right")
fig7.suptitle("Figure S7. CellChat full intercellular communication analysis",
              fontsize=11, fontweight="bold", y=0.985)
save_figure(fig7, "FigureS7_full_CellChat_analysis_v3")
plt.close(fig7)

# ═══════════════════════════════════════════════════════════════
print("\nS8")
assemble([
    ("A", SRC/"GSE26566_validation/FigS_GSE26566_hub_gene_heatmap.pdf", 0, "Hub gene heatmap"),
    ("B", SRC/"GSE26566_validation/FigS_GSE26566_hub_gene_correlation_heatmap.pdf", 0, "Cor heatmap"),
    ("C", SRC/"GSE26566_validation/FigS_GSE26566_aggressive_score_distribution.pdf", 0, "Score dist"),
], 1, 3, "Figure S8. GSE26566 microarray expression validation",
    "FigureS8_GSE26566_validation_v3")

# ═══════════════════════════════════════════════════════════════
print("\nS9")
assemble([
    ("A", SRC/"clinical_relevance/Fig10A_axis_scores_survival_forest_TCGA.pdf", 0, "TCGA axis forest"),
    ("B", SRC/"clinical_relevance/Fig10B_axis_scores_survival_forest_GSE107943.pdf", 0, "GSE axis forest"),
], 1, 2, "Figure S9A. Exploratory survival analyses — axis scores",
    "FigureS9_clinical_relevance_exploratory_v3_page1",
    notes={"A": "Exploratory", "B": "Exploratory"})
assemble([
    ("C", SRC/"clinical_relevance/Fig10C_hub_genes_survival_forest_TCGA.pdf", 0, "TCGA hub forest"),
    ("D", SRC/"clinical_relevance/Fig10D_hub_genes_survival_forest_GSE107943.pdf", 0, "GSE hub forest"),
], 1, 2, "Figure S9B. Hub gene exploratory survival analyses",
    "FigureS9_clinical_relevance_exploratory_v3_page2",
    notes={"C": "Exploratory", "D": "Exploratory"})
assemble([
    ("E", SRC/"clinical_relevance/Fig10I_hub_gene_score_correlation_heatmap_TCGA.pdf", 0, "TCGA cor"),
    ("F", SRC/"clinical_relevance/Fig10J_hub_gene_score_correlation_heatmap_GSE107943.pdf", 0, "GSE cor"),
    ("G", SRC/"clinical_relevance/Fig10K_clinical_stage_axis_score_boxplot_TCGA.pdf", 0, "Stage boxplot"),
], 1, 3, "Figure S9C. Correlation and clinical stage analyses",
    "FigureS9_clinical_relevance_exploratory_v3_page3")

# ═══════════════════════════════════════════════════════════════
print("\nS10")
assemble([
    ("A", SRC/"drug_screening/Fig11B_target_axis_evidence_score_barplot.pdf", 0, "Evidence barplot"),
    ("B", SRC/"drug_screening/Fig11C_candidate_drug_target_network.pdf", 0, "Drug network"),
    ("C", SRC/"drug_screening/Fig11D_docking_target_shortlist.pdf", 0, "Target shortlist"),
    ("D", SRC/"docking/FigS_Docking_affinity_barplot_all_pairs.pdf", 0, "Docking affinity"),
], 2, 2, "Figure S10. Molecular docking and drug screening supplementary details",
    "FigureS10_docking_details_v3",
    notes={"D": "In silico predictions only"})

# ═══════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("V3 GENERATION COMPLETE — 0 errors expected")
print(f"Output: {OUT_V3}")
