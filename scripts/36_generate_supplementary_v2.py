#!/usr/bin/env python
"""Step 36: Regenerate Supplementary Figures S1-S10 at high resolution (600 DPI)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from pathlib import Path
import fitz, csv, io, os
from PIL import Image

BASE = Path(r"e:\CCA")
SRC = BASE / "figures"
TABLES = BASE / "tables"
OUT = BASE / "figures/supplementary_final_v2"
OUT.mkdir(parents=True, exist_ok=True)
DPI = 600  # high resolution
FS = 4.0   # larger panel size

def pdf_to_img(path, page_idx=0, dpi=DPI, max_mp=25):
    """Render PDF page at high DPI. Auto-scale for very large pages."""
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

def add_label(ax, letter, fs=14):
    ax.text(0.02, 0.96, letter, transform=ax.transAxes, fontsize=fs,
            fontweight="bold", color="black", va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.08", facecolor="white",
                      edgecolor="none", alpha=0.85))

def add_note(ax, text, pos="bottom-right", fs=7):
    if pos == "bottom-right":
        x, y, va, ha = 0.98, 0.02, "bottom", "right"
    elif pos == "top-left":
        x, y, va, ha = 0.02, 0.98, "top", "left"
    elif pos == "top-right":
        x, y, va, ha = 0.98, 0.98, "top", "right"
    ax.text(x, y, text, transform=ax.transAxes, fontsize=fs, fontstyle="italic",
            color="#666666", va=va, ha=ha,
            bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                      edgecolor="none", alpha=0.85))

def save_figure(fig, stem, dpi=DPI):
    pdf_path = OUT / f"{stem}.pdf"
    png_path = OUT / f"{stem}.png"
    fig.savefig(pdf_path, dpi=dpi, facecolor="white", edgecolor="none",
                bbox_inches="tight")
    fig.savefig(png_path, dpi=dpi, facecolor="white", edgecolor="none",
                bbox_inches="tight")
    sz_pdf = pdf_path.stat().st_size // 1024
    sz_png = png_path.stat().st_size // 1024
    print(f"  -> {stem}.pdf ({sz_pdf} KB) + PNG ({sz_png} KB)")

def fig_to_img(fig):
    fig.canvas.draw()
    img = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
    w, h = fig.canvas.get_width_height()
    return img.reshape(h, w, 4)[:, :, :3]

def make_bubble(rows, title, fs=(8, 6)):
    """CellChat-style bubble from list of dicts."""
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

# Load CellChat data for S7 regeneration
lr_csv = TABLES / "cellchat/GSE138709_CAF_TAM_Epithelial_LR_pairs.csv"
with open(lr_csv, newline="") as f:
    all_lr = list(csv.DictReader(f))

caf_tam_lr = [r for r in all_lr if r["source"] == "CAF_Fibroblast" and r["target"] == "Macrophage_TAM"]
tam_caf_lr = [r for r in all_lr if r["source"] == "Macrophage_TAM" and r["target"] == "CAF_Fibroblast"]
epi_tam_lr = [r for r in all_lr if r["source"] == "Epithelial_like" and r["target"] == "Macrophage_TAM"]

print(f"Loaded {len(all_lr)} LR pairs: CAF->TAM={len(caf_tam_lr)}, TAM->CAF={len(tam_caf_lr)}, Epi->TAM={len(epi_tam_lr)}")

errors = {}

# ═══════════════════════════════════════════════════════════════
def assemble_single_page(panels, nrows, ncols, title, stem, notes=None):
    """Assemble a single-page figure. panels=[(letter, path, pg, desc)]"""
    fig = plt.figure(figsize=(ncols*FS, nrows*FS), facecolor="white")
    gs = gridspec.GridSpec(nrows, ncols, figure=fig, wspace=0.10, hspace=0.18,
                           left=0.03, right=0.97, top=0.94, bottom=0.03)
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
                add_note(ax, notes[letter])
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
    fig.suptitle(title, fontsize=11, fontweight="bold", y=0.99)
    save_figure(fig, stem)
    plt.close(fig)
    return errs

# ═══════════════════════════════════════════════════════════════
print("="*60)
print("S1: Preprocessing & QC")
print("="*60)
e = assemble_single_page([
    ("A", SRC/"DEG/Fig2A_TCGA_PCA_tumor_normal.pdf", 0, "TCGA PCA"),
    ("B", SRC/"DEG/Fig2D_GSE107943_PCA_tumor_normal.pdf", 0, "GSE PCA"),
    ("C", SRC/"single_cell/Fig7A_GSE138709_UMAP_clusters.pdf", 0, "scRNA UMAP"),
    ("D", SRC/"GSE26566_validation/FigS_GSE26566_hub_gene_heatmap.pdf", 0, "GSE26566 heatmap"),
], 2, 2, "Figure S1. Data preprocessing and quality control overview",
    "FigureS1_preprocessing_QC_v2")
errors["S1"] = e

# ═══════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("S2: Full DEG Plots (2 pages: large heatmaps)")
print("="*60)

# Page 1: TCGA heatmap
fig2a = plt.figure(figsize=(10, 8), facecolor="white")
ax = fig2a.add_subplot(111)
img = pdf_to_img(SRC/"DEG/Fig2C_TCGA_top_IM_CAF_TAM_DEG_heatmap.pdf")
ax.imshow(img); ax.axis("off")
add_label(ax, "A", fs=16)
fig2a.suptitle("Figure S2A. TCGA-CHOL: Top IM/CAF/TAM DEG heatmap", fontsize=11, fontweight="bold", y=0.99)
save_figure(fig2a, "FigureS2_full_DEG_plots_v2_page1")
plt.close(fig2a)

# Page 2: GSE heatmap
fig2b = plt.figure(figsize=(10, 8), facecolor="white")
ax = fig2b.add_subplot(111)
img = pdf_to_img(SRC/"DEG/Fig2G_validated_IM_CAF_TAM_DEGs_heatmap.pdf")
ax.imshow(img); ax.axis("off")
add_label(ax, "B", fs=16)
fig2b.suptitle("Figure S2B. GSE107943: Validated IM/CAF/TAM DEG heatmap", fontsize=11, fontweight="bold", y=0.99)
save_figure(fig2b, "FigureS2_full_DEG_plots_v2_page2")
plt.close(fig2b)
errors["S2"] = []
print("  S2: 2 pages, large heatmaps")

# ═══════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("S3: GSVA Score + Survival (2 pages)")
print("="*60)

# Page 1: boxplots A+B
e1 = assemble_single_page([
    ("A", SRC/"gsva/Fig4A_TCGA_score_boxplot_tumor_normal.pdf", 0, "TCGA boxplot"),
    ("B", SRC/"gsva/Fig4B_GSE107943_score_boxplot_tumor_normal.pdf", 0, "GSE boxplot"),
], 1, 2, "Figure S3A. GSVA ssGSEA pathway scores (tumor vs normal)",
    "FigureS3_GSVA_score_survival_v2_page1")
errors["S3"] = e1

# Page 2: KM curves C/D/E
e2 = assemble_single_page([
    ("C", SRC/"gsva/Fig4E_TCGA_KM_high_low_CAF_score.pdf", 1, "TCGA CAF KM"),
    ("D", SRC/"gsva/Fig4F_TCGA_KM_high_low_IM_CAF_TAM_score.pdf", 1, "TCGA IM_CAF_TAM KM"),
    ("E", SRC/"gsva/Fig4G_GSE107943_KM_high_low_CAF_score.pdf", 1, "GSE CAF KM"),
], 1, 3, "Figure S3B. Exploratory survival analyses (all exploratory)",
    "FigureS3_GSVA_score_survival_v2_page2",
    notes={"C": "Exploratory", "D": "Exploratory", "E": "Exploratory"})
errors["S3"] = errors.get("S3", []) + e2

# ═══════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("S4: Subtyping Diagnostics (2 pages)")
print("="*60)

e1 = assemble_single_page([
    ("A", SRC/"subtyping/Fig5A_TCGA_consensus_CDF.pdf", 0, "Consensus CDF"),
    ("B", SRC/"subtyping/Fig5B_TCGA_consensus_matrix_k2.pdf", 0, "k=2 matrix"),
    ("C", SRC/"subtyping/Fig5D_TCGA_subtype_score_boxplot.pdf", 0, "Score boxplot"),
    ("D", SRC/"subtyping/Fig5F_GSE107943_subtype_heatmap.pdf", 0, "GSE heatmap"),
], 2, 2, "Figure S4A. Molecular subtyping diagnostics (1 of 2)",
    "FigureS4_molecular_subtyping_diagnostics_v2_page1")
errors["S4"] = e1

e2 = assemble_single_page([
    ("E", SRC/"subtyping_diagnosis/FigS_original_subtype_gene_direction_heatmap.pdf", 0, "Direction heatmap"),
    ("F", SRC/"subtyping_diagnosis/FigS_strategy_A_up_genes_heatmap_TCGA.pdf", 0, "Strategy A"),
], 1, 2, "Figure S4B. Subtyping strategy diagnostics (2 of 2)",
    "FigureS4_molecular_subtyping_diagnostics_v2_page2")
errors["S4"] = errors.get("S4", []) + e2

# ═══════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("S5: Full Immune Analysis (single page, large)")
print("="*60)

e = assemble_single_page([
    ("A", SRC/"immune/Fig6B_GSE107943_immune_score_boxplot_tumor_normal.pdf", 0, "Immune boxplot"),
    ("B", SRC/"immune/Fig6D_GSE107943_aggressive_immune_correlation_heatmap.pdf", 0, "Immune cor"),
    ("C", SRC/"immune/Fig6F_GSE107943_checkpoint_correlation_heatmap.pdf", 0, "Checkpoint cor"),
    ("D", SRC/"immune/Fig6H_aggressive_vs_macrophage_scatter_GSE107943.pdf", 0, "Macro scatter"),
    ("E", SRC/"immune/Fig6J_aggressive_vs_checkpoint_scatter_GSE107943.pdf", 0, "Checkpoint scatter"),
], 2, 3, "Figure S5. GSE107943 immune infiltration and checkpoint correlation analyses",
    "FigureS5_full_immune_analysis_v2")
errors["S5"] = e

# ═══════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("S6: Single-Cell Supplementary (3 pages)")
print("="*60)

# Page 1: marker dotplot
fig6a = plt.figure(figsize=(10, 8), facecolor="white")
ax = fig6a.add_subplot(111)
img = pdf_to_img(SRC/"single_cell/Fig7C_GSE138709_marker_dotplot.pdf")
ax.imshow(img); ax.axis("off")
add_label(ax, "A", fs=16)
fig6a.suptitle("Figure S6A. Cell type marker dotplot", fontsize=11, fontweight="bold", y=0.99)
save_figure(fig6a, "FigureS6_single_cell_supplementary_v2_page1")
plt.close(fig6a)

# Page 2: score UMAPs
e2 = assemble_single_page([
    ("B", SRC/"single_cell/Fig7E_GSE138709_CAF_TAM_score_UMAP.pdf", 0, "CAF/TAM UMAP"),
    ("C", SRC/"single_cell/Fig7D_GSE138709_aggressive_score_UMAP.pdf", 0, "Aggressive UMAP"),
], 1, 2, "Figure S6B. Single-cell module score projections",
    "FigureS6_single_cell_supplementary_v2_page2")
errors["S6"] = e2

# Page 3: feature plots
fig6c = plt.figure(figsize=(10, 8), facecolor="white")
ax = fig6c.add_subplot(111)
img = pdf_to_img(SRC/"single_cell/Fig7G_GSE138709_key_genes_featureplot.pdf", max_mp=30)
ax.imshow(img); ax.axis("off")
add_label(ax, "D", fs=16)
fig6c.suptitle("Figure S6C. Key gene feature plots", fontsize=11, fontweight="bold", y=0.99)
save_figure(fig6c, "FigureS6_single_cell_supplementary_v2_page3")
plt.close(fig6c)
print("  S6: 3 pages, large panels")

# ═══════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("S7: Full CellChat Analysis (regenerated)")
print("="*60)

# Panel A: interaction number network (from PDF)
# Panel B: interaction weight network (from PDF)
# Panel C: CAF->TAM bubble (regenerated from CSV)
# Panel D: TAM->CAF bubble (regenerated from CSV)

# Generate bubble plots
fig_c = make_bubble(caf_tam_lr, "CAF >> TAM inferred communication", fs=(8, 6.5))
img_c = fig_to_img(fig_c); plt.close(fig_c)
fig_d = make_bubble(tam_caf_lr, "TAM >> CAF inferred communication", fs=(8, 6.5))
img_d = fig_to_img(fig_d); plt.close(fig_d)
print(f"  Regenerated bubble C ({img_c.shape[1]}x{img_c.shape[0]})")
print(f"  Regenerated bubble D ({img_d.shape[1]}x{img_d.shape[0]})")

# Assemble S7 as a 2x2 grid using both PDF panels and regenerated images
fig7 = plt.figure(figsize=(10, 8.5), facecolor="white")
gs7 = gridspec.GridSpec(2, 2, figure=fig7, wspace=0.10, hspace=0.18,
                        left=0.03, right=0.97, top=0.94, bottom=0.03)

s7_panels = [
    ("A", pdf_to_img(SRC/"cellchat/Fig8A_CellChat_overall_interaction_number.pdf")),
    ("B", pdf_to_img(SRC/"cellchat/Fig8B_CellChat_overall_interaction_weight.pdf")),
    ("C", img_c),
    ("D", img_d),
]
for idx, (letter, img) in enumerate(s7_panels):
    ax = fig7.add_subplot(gs7[idx//2, idx%2])
    ax.imshow(img); ax.axis("off")
    add_label(ax, letter)
    if letter in ("C", "D"):
        add_note(ax, "Inferred communication", pos="bottom-right")
    print(f"  {letter}: embedded ({img.shape[1]}x{img.shape[0]})")

fig7.suptitle("Figure S7. CellChat full intercellular communication analysis",
              fontsize=11, fontweight="bold", y=0.99)
save_figure(fig7, "FigureS7_full_CellChat_analysis_v2")
plt.close(fig7)
errors["S7"] = []
print("  S7: All 4 panels OK (C/D regenerated from CSV)")

# ═══════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("S8: GSE26566 Validation (larger canvas)")
print("="*60)

e = assemble_single_page([
    ("A", SRC/"GSE26566_validation/FigS_GSE26566_hub_gene_heatmap.pdf", 0, "Hub gene heatmap"),
    ("B", SRC/"GSE26566_validation/FigS_GSE26566_hub_gene_correlation_heatmap.pdf", 0, "Cor heatmap"),
    ("C", SRC/"GSE26566_validation/FigS_GSE26566_aggressive_score_distribution.pdf", 0, "Score dist"),
], 1, 3, "Figure S8. GSE26566 microarray expression validation",
    "FigureS8_GSE26566_validation_v2")
errors["S8"] = e

# ═══════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("S9: Clinical Relevance (3 pages)")
print("="*60)

# Page 1: forest plots A+B
e1 = assemble_single_page([
    ("A", SRC/"clinical_relevance/Fig10A_axis_scores_survival_forest_TCGA.pdf", 0, "TCGA axis forest"),
    ("B", SRC/"clinical_relevance/Fig10B_axis_scores_survival_forest_GSE107943.pdf", 0, "GSE axis forest"),
], 1, 2, "Figure S9A. Exploratory survival analyses — axis scores (all exploratory)",
    "FigureS9_clinical_relevance_exploratory_v2_page1",
    notes={"A": "Exploratory", "B": "Exploratory"})
errors["S9"] = e1

# Page 2: hub gene forests
e2 = assemble_single_page([
    ("C", SRC/"clinical_relevance/Fig10C_hub_genes_survival_forest_TCGA.pdf", 0, "TCGA hub forest"),
    ("D", SRC/"clinical_relevance/Fig10D_hub_genes_survival_forest_GSE107943.pdf", 0, "GSE hub forest"),
], 1, 2, "Figure S9B. Hub gene exploratory survival analyses (all exploratory)",
    "FigureS9_clinical_relevance_exploratory_v2_page2",
    notes={"C": "Exploratory", "D": "Exploratory"})
errors["S9"] = errors.get("S9", []) + e2

# Page 3: correlation heatmaps + stage
e3 = assemble_single_page([
    ("E", SRC/"clinical_relevance/Fig10I_hub_gene_score_correlation_heatmap_TCGA.pdf", 0, "TCGA cor"),
    ("F", SRC/"clinical_relevance/Fig10J_hub_gene_score_correlation_heatmap_GSE107943.pdf", 0, "GSE cor"),
    ("G", SRC/"clinical_relevance/Fig10K_clinical_stage_axis_score_boxplot_TCGA.pdf", 0, "Stage boxplot"),
], 1, 3, "Figure S9C. Correlation and clinical stage analyses",
    "FigureS9_clinical_relevance_exploratory_v2_page3")
errors["S9"] = errors.get("S9", []) + e3

# ═══════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("S10: Docking Details (larger canvas)")
print("="*60)

e = assemble_single_page([
    ("A", SRC/"drug_screening/Fig11B_target_axis_evidence_score_barplot.pdf", 0, "Evidence barplot"),
    ("B", SRC/"drug_screening/Fig11C_candidate_drug_target_network.pdf", 0, "Drug network"),
    ("C", SRC/"drug_screening/Fig11D_docking_target_shortlist.pdf", 0, "Target shortlist"),
    ("D", SRC/"docking/FigS_Docking_affinity_barplot_all_pairs.pdf", 0, "Docking affinity"),
], 2, 2, "Figure S10. Molecular docking and drug screening supplementary details",
    "FigureS10_docking_details_v2",
    notes={"D": "In silico predictions only"})
errors["S10"] = e

# ═══════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
total_errs = 0
for fig, errs in sorted(errors.items()):
    n = len(errs)
    total_errs += n
    print(f"  {fig}: {'OK' if n==0 else f'{n} errors'}")
    for e in errs:
        print(f"    - {e}")
print(f"\nTotal errors: {total_errs}")
print(f"Output: {OUT}")
print("Done.")
