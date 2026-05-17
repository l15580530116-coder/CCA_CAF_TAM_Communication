#!/usr/bin/env python
"""Assemble Supplementary Figures S2–S10 from existing source panel PDFs."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from pathlib import Path
import fitz
from PIL import Image
import io

BASE = Path(r"e:\CCA")
SRC = BASE / "figures"
OUT = BASE / "figures/supplementary_final"
OUT.mkdir(parents=True, exist_ok=True)
DPI = 300
FS = 3.2  # inches per panel

def pdf_to_img(path, page_idx=0):
    """Render PDF page to numpy array."""
    doc = fitz.open(str(path))
    if page_idx >= len(doc):
        page_idx = len(doc) - 1
    page = doc[page_idx]
    zoom = DPI / 72.0
    mat = fitz.Matrix(zoom, zoom)
    pix = page.get_pixmap(matrix=mat, colorspace=fitz.csRGB)
    img = np.array(Image.open(io.BytesIO(pix.tobytes("png"))))
    doc.close()
    return img

def add_label(ax, letter, fontsize=13):
    ax.text(0.02, 0.96, letter, transform=ax.transAxes, fontsize=fontsize,
            fontweight="bold", color="black", va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.1", facecolor="white",
                      edgecolor="none", alpha=0.85))

def add_note(ax, text):
    ax.text(0.98, 0.02, text, transform=ax.transAxes, fontsize=6.5,
            fontstyle="italic", color="#888888", va="bottom", ha="right")

def assemble_figure(panels, nrows, ncols, title, stem, caveats=None):
    """
    panels: list of (letter, path, page_idx, desc, note) tuples
    caveats: dict {letter: str} for small annotation text
    """
    fig = plt.figure(figsize=(ncols*FS, nrows*FS), facecolor="white")
    gs = gridspec.GridSpec(nrows, ncols, figure=fig,
                           wspace=0.10, hspace=0.18,
                           left=0.03, right=0.97, top=0.94, bottom=0.03)

    errors = []
    for idx, item in enumerate(panels):
        letter, path, pg, desc, note = item
        row, col = divmod(idx, ncols)
        ax = fig.add_subplot(gs[row, col])
        try:
            img = pdf_to_img(path, page_idx=pg)
            ax.imshow(img)
            ax.axis("off")
            add_label(ax, letter)
            if caveats and letter in caveats:
                add_note(ax, caveats[letter])
            print(f"  {letter}: {desc} — OK")
        except Exception as e:
            ax.text(0.5, 0.5, f"Panel {letter}\nNOT AVAILABLE\n{str(e)[:60]}",
                    ha="center", va="center", fontsize=7, color="red")
            ax.axis("off")
            add_label(ax, letter)
            errors.append(f"{letter}: {e}")
            print(f"  {letter}: ERROR — {e}")

    # Hide unused slots
    for idx in range(len(panels), nrows*ncols):
        row, col = divmod(idx, ncols)
        ax = fig.add_subplot(gs[row, col])
        ax.axis("off")

    fig.suptitle(title, fontsize=10, fontweight="bold", y=0.99)

    pdf_path = OUT / f"{stem}.pdf"
    png_path = OUT / f"{stem}.png"
    fig.savefig(pdf_path, dpi=DPI, facecolor="white", edgecolor="none",
                bbox_inches="tight")
    fig.savefig(png_path, dpi=DPI, facecolor="white", edgecolor="none",
                bbox_inches="tight")
    plt.close(fig)
    sz_pdf = pdf_path.stat().st_size // 1024
    sz_png = png_path.stat().st_size // 1024
    print(f"  Saved: {stem}.pdf ({sz_pdf} KB) + PNG ({sz_png} KB)")
    return errors

# ═══════════════════════════════════════════════════════════════
print("Assembling Supplementary Figures S2–S10\n")

all_errors = {}

# ── S2: Full DEG plots ──
print("=== S2: Full DEG Plots ===")
s2_panels = [
    ("A", SRC/"DEG/Fig2C_TCGA_top_IM_CAF_TAM_DEG_heatmap.pdf", 0,
     "TCGA-CHOL: Top IM/CAF/TAM DEG heatmap", None),
    ("B", SRC/"DEG/Fig2G_validated_IM_CAF_TAM_DEGs_heatmap.pdf", 0,
     "GSE107943: Validated IM/CAF/TAM DEG heatmap", None),
]
e = assemble_figure(s2_panels, 1, 2,
    "Figure S2. Full differential expression heatmaps",
    "FigureS2_full_DEG_plots")
all_errors["S2"] = e

# ── S3: GSVA score + survival ──
print("\n=== S3: GSVA Score + Survival ===")
s3_panels = [
    ("A", SRC/"gsva/Fig4A_TCGA_score_boxplot_tumor_normal.pdf", 0,
     "TCGA-CHOL: Score boxplot (T vs N)", None),
    ("B", SRC/"gsva/Fig4B_GSE107943_score_boxplot_tumor_normal.pdf", 0,
     "GSE107943: Score boxplot (T vs N)", None),
    ("C", SRC/"gsva/Fig4E_TCGA_KM_high_low_CAF_score.pdf", 1,
     "TCGA-CHOL: CAF score KM", "2-page PDF, page 1 used"),
    ("D", SRC/"gsva/Fig4F_TCGA_KM_high_low_IM_CAF_TAM_score.pdf", 1,
     "TCGA-CHOL: IM_CAF_TAM score KM", "2-page PDF, page 1 used"),
    ("E", SRC/"gsva/Fig4G_GSE107943_KM_high_low_CAF_score.pdf", 1,
     "GSE107943: CAF score KM", "2-page PDF, page 1 used"),
]
e = assemble_figure(s3_panels, 2, 3,
    "Figure S3. GSVA ssGSEA pathway scores and exploratory survival analyses",
    "FigureS3_GSVA_score_survival")
all_errors["S3"] = e

# ── S4: Molecular subtyping diagnostics ──
print("\n=== S4: Molecular Subtyping Diagnostics ===")
s4_panels = [
    ("A", SRC/"subtyping/Fig5A_TCGA_consensus_CDF.pdf", 0,
     "Consensus CDF (k=2–6)", None),
    ("B", SRC/"subtyping/Fig5B_TCGA_consensus_matrix_k2.pdf", 0,
     "Consensus matrix (k=2)", None),
    ("C", SRC/"subtyping/Fig5C_TCGA_subtype_heatmap.pdf", 0,
     "TCGA subtype heatmap", None),
    ("D", SRC/"subtyping/Fig5F_GSE107943_subtype_heatmap.pdf", 0,
     "GSE107943 subtype heatmap", None),
    ("E", SRC/"subtyping_diagnosis/FigS_original_subtype_gene_direction_heatmap.pdf", 0,
     "Gene direction heatmap (diagnostic)", None),
    ("F", SRC/"subtyping_diagnosis/FigS_strategy_A_up_genes_heatmap_TCGA.pdf", 0,
     "Strategy A: up-genes heatmap", None),
]
e = assemble_figure(s4_panels, 3, 2,
    "Figure S4. Molecular subtyping diagnostic analyses",
    "FigureS4_molecular_subtyping_diagnostics")
all_errors["S4"] = e

# ── S5: Full immune analysis ──
print("\n=== S5: Full Immune Analysis (GSE107943) ===")
s5_panels = [
    ("A", SRC/"immune/Fig6B_GSE107943_immune_score_boxplot_tumor_normal.pdf", 0,
     "GSE107943: Immune scores boxplot", None),
    ("B", SRC/"immune/Fig6D_GSE107943_aggressive_immune_correlation_heatmap.pdf", 0,
     "GSE107943: Immune correlation heatmap", None),
    ("C", SRC/"immune/Fig6F_GSE107943_checkpoint_correlation_heatmap.pdf", 0,
     "GSE107943: Checkpoint correlation heatmap", None),
    ("D", SRC/"immune/Fig6H_aggressive_vs_macrophage_scatter_GSE107943.pdf", 0,
     "GSE107943: Aggressive vs Macrophage scatter", None),
    ("E", SRC/"immune/Fig6J_aggressive_vs_checkpoint_scatter_GSE107943.pdf", 0,
     "GSE107943: Aggressive vs Checkpoint scatter", None),
]
e = assemble_figure(s5_panels, 2, 3,
    "Figure S5. GSE107943 immune infiltration and checkpoint correlation analyses",
    "FigureS5_full_immune_analysis")
all_errors["S5"] = e

# ── S6: Single-cell supplementary ──
print("\n=== S6: Single-Cell Supplementary ===")
s6_panels = [
    ("A", SRC/"single_cell/Fig7C_GSE138709_marker_dotplot.pdf", 0,
     "Marker dotplot by cell type", None),
    ("B", SRC/"single_cell/Fig7E_GSE138709_CAF_TAM_score_UMAP.pdf", 0,
     "CAF and TAM module score UMAP", None),
    ("C", SRC/"single_cell/Fig7G_GSE138709_key_genes_featureplot.pdf", 0,
     "Key gene feature plots", "Large file (21 MB)"),
]
e = assemble_figure(s6_panels, 1, 3,
    "Figure S6. Single-cell RNA-seq supplementary analyses",
    "FigureS6_single_cell_supplementary")
all_errors["S6"] = e

# ── S7: Full CellChat analysis ──
print("\n=== S7: Full CellChat Analysis ===")
s7_panels = [
    ("A", SRC/"cellchat/Fig8A_CellChat_overall_interaction_number.pdf", 0,
     "Overall interaction numbers", None),
    ("B", SRC/"cellchat/Fig8B_CellChat_overall_interaction_weight.pdf", 0,
     "Overall interaction weights", None),
    ("C", SRC/"cellchat/Fig8E_CellChat_top_pathways_bubble.pdf", 0,
     "Top pathways bubble plot", None),
    ("D", SRC/"cellchat/Fig8H_CellChat_CAF_TAM_Epithelial_key_pathways.pdf", 0,
     "CAF–TAM–Epithelial key pathways", None),
]
e = assemble_figure(s7_panels, 2, 2,
    "Figure S7. CellChat full intercellular communication analysis",
    "FigureS7_full_CellChat_analysis")
all_errors["S7"] = e

# ── S8: GSE26566 validation ──
print("\n=== S8: GSE26566 Validation ===")
s8_panels = [
    ("A", SRC/"GSE26566_validation/FigS_GSE26566_hub_gene_heatmap.pdf", 0,
     "Hub gene expression heatmap", None),
    ("B", SRC/"GSE26566_validation/FigS_GSE26566_hub_gene_correlation_heatmap.pdf", 0,
     "Hub gene correlation heatmap", None),
    ("C", SRC/"GSE26566_validation/FigS_GSE26566_aggressive_score_distribution.pdf", 0,
     "Aggressive score distribution", None),
]
e = assemble_figure(s8_panels, 1, 3,
    "Figure S8. GSE26566 microarray expression validation",
    "FigureS8_GSE26566_validation")
all_errors["S8"] = e

# ── S9: Clinical relevance exploratory ──
print("\n=== S9: Clinical Relevance Exploratory ===")
s9_panels = [
    ("A", SRC/"clinical_relevance/Fig10A_axis_scores_survival_forest_TCGA.pdf", 0,
     "TCGA: Axis scores forest plot", None),
    ("B", SRC/"clinical_relevance/Fig10C_hub_genes_survival_forest_TCGA.pdf", 0,
     "TCGA: Hub genes forest plot", None),
    ("C", SRC/"clinical_relevance/Fig10D_hub_genes_survival_forest_GSE107943.pdf", 0,
     "GSE107943: Hub genes forest plot", None),
    ("D", SRC/"clinical_relevance/Fig10I_hub_gene_score_correlation_heatmap_TCGA.pdf", 0,
     "TCGA: Correlation heatmap", None),
    ("E", SRC/"clinical_relevance/Fig10J_hub_gene_score_correlation_heatmap_GSE107943.pdf", 0,
     "GSE107943: Correlation heatmap", None),
    ("F", SRC/"clinical_relevance/Fig10K_clinical_stage_axis_score_boxplot_TCGA.pdf", 0,
     "TCGA: Clinical stage boxplot", None),
]
caveats_s9 = {"A": "All survival analyses are exploratory",
              "B": "All survival analyses are exploratory",
              "C": "All survival analyses are exploratory"}
e = assemble_figure(s9_panels, 3, 2,
    "Figure S9. Exploratory clinical relevance analyses",
    "FigureS9_clinical_relevance_exploratory",
    caveats=caveats_s9)
all_errors["S9"] = e

# ── S10: Docking details ──
print("\n=== S10: Docking Details ===")
s10_panels = [
    ("A", SRC/"drug_screening/Fig11B_target_axis_evidence_score_barplot.pdf", 0,
     "Target-axis evidence barplot", None),
    ("B", SRC/"drug_screening/Fig11C_candidate_drug_target_network.pdf", 0,
     "Candidate drug-target network", None),
    ("C", SRC/"drug_screening/Fig11D_docking_target_shortlist.pdf", 0,
     "Docking target shortlist", None),
    ("D", SRC/"docking/FigS_Docking_affinity_barplot_all_pairs.pdf", 0,
     "Docking affinity barplot (all pairs)", "See main Figure 7B for 4-pair version"),
]
e = assemble_figure(s10_panels, 2, 2,
    "Figure S10. Molecular docking and drug screening supplementary details",
    "FigureS10_docking_details",
    caveats={"D": "In silico predictions only"})
all_errors["S10"] = e

# ═══════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("ASSEMBLY COMPLETE")
print("="*60)
print(f"\nOutput directory: {OUT}")
for stem in [f"FigureS{i}" for i in range(2, 11)]:
    for ext in [".pdf", ".png"]:
        f = OUT / f"{stem}*{ext}"
        matches = list(OUT.glob(f"{stem}*{ext}"))
        if matches:
            sz = matches[0].stat().st_size // 1024
            print(f"  {matches[0].name} ({sz} KB)")

print(f"\nError summary:")
for fig, errs in all_errors.items():
    if errs:
        print(f"  {fig}: {len(errs)} error(s) — {errs}")
    else:
        print(f"  {fig}: 0 errors ✓")

print("\nDone.")
