#!/usr/bin/env python
"""Assemble S4–S10 with robust handling of 0-page PDFs and large renders."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from pathlib import Path
import fitz
from PIL import Image
import io, sys

BASE = Path(r"e:\CCA")
SRC = BASE / "figures"
OUT = BASE / "figures/supplementary_final"
OUT.mkdir(parents=True, exist_ok=True)
DPI = 200  # lower DPI for supplementary (adequate, faster)
FS = 3.0

def pdf_to_img(path, page_idx=0, max_megapixels=12):
    doc = fitz.open(str(path))
    if len(doc) == 0:
        doc.close()
        raise ValueError("0-page PDF")
    if page_idx >= len(doc):
        page_idx = len(doc) - 1
    page = doc[page_idx]
    # Auto-scale DPI for large pages
    pw, ph = page.rect.width, page.rect.height
    dpi = DPI
    mpix = (pw*dpi/72) * (ph*dpi/72) / 1e6
    if mpix > max_megapixels:
        dpi = int(DPI * (max_megapixels/mpix)**0.5)
    zoom = dpi / 72.0
    mat = fitz.Matrix(zoom, zoom)
    pix = page.get_pixmap(matrix=mat, colorspace=fitz.csRGB)
    img = np.array(Image.open(io.BytesIO(pix.tobytes("png"))))
    doc.close()
    return img

def add_label(ax, letter):
    ax.text(0.02, 0.96, letter, transform=ax.transAxes, fontsize=12,
            fontweight="bold", color="black", va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.08", facecolor="white",
                      edgecolor="none", alpha=0.85))

def assemble_figure(panels, nrows, ncols, title, stem):
    fig = plt.figure(figsize=(ncols*FS, nrows*FS), facecolor="white")
    gs = gridspec.GridSpec(nrows, ncols, figure=fig,
                           wspace=0.10, hspace=0.20,
                           left=0.03, right=0.97, top=0.93, bottom=0.03)
    errors = []
    for idx, (letter, path, pg, desc) in enumerate(panels):
        row, col = divmod(idx, ncols)
        ax = fig.add_subplot(gs[row, col])
        try:
            img = pdf_to_img(path, page_idx=pg)
            ax.imshow(img)
            ax.axis("off")
            add_label(ax, letter)
            print(f"  {letter}: OK — {desc}")
        except Exception as e:
            emsg = str(e)[:80]
            ax.text(0.5, 0.5, f"Panel {letter}\nNOT AVAILABLE\n{emsg}",
                    ha="center", va="center", fontsize=7, color="red")
            ax.axis("off")
            add_label(ax, letter)
            errors.append(f"{letter}: {emsg}")
            print(f"  {letter}: MISSING — {emsg}")

    for idx in range(len(panels), nrows*ncols):
        ax = fig.add_subplot(gs[idx//ncols, idx%ncols])
        ax.axis("off")

    fig.suptitle(title, fontsize=9, fontweight="bold", y=0.99)

    for fmt in ["pdf", "png"]:
        path = OUT / f"{stem}.{fmt}"
        fig.savefig(path, dpi=DPI, facecolor="white", edgecolor="none",
                    bbox_inches="tight")
        sz = path.stat().st_size // 1024
        print(f"  -> {stem}.{fmt} ({sz} KB)")
    plt.close(fig)
    return errors

all_errors = {}

# ── S4: Subtyping diagnostics ──
print("=== S4: Molecular Subtyping Diagnostics ===\n")
s4 = [
    ("A", SRC/"subtyping/Fig5A_TCGA_consensus_CDF.pdf", 0, "Consensus CDF"),
    ("B", SRC/"subtyping/Fig5B_TCGA_consensus_matrix_k2.pdf", 0, "Consensus matrix k=2"),
    ("C", SRC/"subtyping/Fig5D_TCGA_subtype_score_boxplot.pdf", 0, "TCGA score boxplot"),
    ("D", SRC/"subtyping/Fig5F_GSE107943_subtype_heatmap.pdf", 0, "GSE subtype heatmap"),
    ("E", SRC/"subtyping_diagnosis/FigS_original_subtype_gene_direction_heatmap.pdf", 0, "Gene direction heatmap"),
    ("F", SRC/"subtyping_diagnosis/FigS_strategy_A_up_genes_heatmap_TCGA.pdf", 0, "Strategy A heatmap"),
]
e = assemble_figure(s4, 3, 2,
    "Figure S4. Molecular subtyping diagnostic analyses", "FigureS4_molecular_subtyping_diagnostics")
all_errors["S4"] = e

# ── S5: Immune GSE ──
print("\n=== S5: Full Immune Analysis ===\n")
s5 = [
    ("A", SRC/"immune/Fig6B_GSE107943_immune_score_boxplot_tumor_normal.pdf", 0, "GSE immune boxplot"),
    ("B", SRC/"immune/Fig6D_GSE107943_aggressive_immune_correlation_heatmap.pdf", 0, "GSE immune cor heatmap"),
    ("C", SRC/"immune/Fig6F_GSE107943_checkpoint_correlation_heatmap.pdf", 0, "GSE checkpoint cor"),
    ("D", SRC/"immune/Fig6H_aggressive_vs_macrophage_scatter_GSE107943.pdf", 0, "GSE macro scatter"),
    ("E", SRC/"immune/Fig6J_aggressive_vs_checkpoint_scatter_GSE107943.pdf", 0, "GSE checkpoint scatter"),
]
e = assemble_figure(s5, 2, 3,
    "Figure S5. GSE107943 immune infiltration and checkpoint correlation", "FigureS5_full_immune_analysis")
all_errors["S5"] = e

# ── S6: Single-cell ──
print("\n=== S6: Single-Cell Supplementary ===\n")
s6 = [
    ("A", SRC/"single_cell/Fig7C_GSE138709_marker_dotplot.pdf", 0, "Marker dotplot"),
    ("B", SRC/"single_cell/Fig7E_GSE138709_CAF_TAM_score_UMAP.pdf", 0, "CAF/TAM score UMAP"),
    ("C", SRC/"single_cell/Fig7G_GSE138709_key_genes_featureplot.pdf", 0, "Feature plots"),
]
e = assemble_figure(s6, 1, 3,
    "Figure S6. Single-cell RNA-seq supplementary analyses", "FigureS6_single_cell_supplementary")
all_errors["S6"] = e

# ── S7: CellChat ──
print("\n=== S7: Full CellChat Analysis ===\n")
s7 = [
    ("A", SRC/"cellchat/Fig8A_CellChat_overall_interaction_number.pdf", 0, "Interaction numbers"),
    ("B", SRC/"cellchat/Fig8B_CellChat_overall_interaction_weight.pdf", 0, "Interaction weights"),
    ("C", SRC/"cellchat/Fig8E_CellChat_top_pathways_bubble.pdf", 0, "Top pathways bubble"),
    ("D", SRC/"cellchat/Fig8H_CellChat_CAF_TAM_Epithelial_key_pathways.pdf", 0, "Key pathways"),
]
e = assemble_figure(s7, 2, 2,
    "Figure S7. CellChat full intercellular communication analysis", "FigureS7_full_CellChat_analysis")
all_errors["S7"] = e

# ── S8: GSE26566 ──
print("\n=== S8: GSE26566 Validation ===\n")
s8 = [
    ("A", SRC/"GSE26566_validation/FigS_GSE26566_hub_gene_heatmap.pdf", 0, "Hub gene heatmap"),
    ("B", SRC/"GSE26566_validation/FigS_GSE26566_hub_gene_correlation_heatmap.pdf", 0, "Correlation heatmap"),
    ("C", SRC/"GSE26566_validation/FigS_GSE26566_aggressive_score_distribution.pdf", 0, "Score distribution"),
]
e = assemble_figure(s8, 1, 3,
    "Figure S8. GSE26566 microarray expression validation", "FigureS8_GSE26566_validation")
all_errors["S8"] = e

# ── S9: Clinical ──
print("\n=== S9: Clinical Relevance ===\n")
s9 = [
    ("A", SRC/"clinical_relevance/Fig10A_axis_scores_survival_forest_TCGA.pdf", 0, "TCGA axis forest"),
    ("B", SRC/"clinical_relevance/Fig10C_hub_genes_survival_forest_TCGA.pdf", 0, "TCGA hub forest"),
    ("C", SRC/"clinical_relevance/Fig10D_hub_genes_survival_forest_GSE107943.pdf", 0, "GSE hub forest"),
    ("D", SRC/"clinical_relevance/Fig10I_hub_gene_score_correlation_heatmap_TCGA.pdf", 0, "TCGA cor heatmap"),
    ("E", SRC/"clinical_relevance/Fig10J_hub_gene_score_correlation_heatmap_GSE107943.pdf", 0, "GSE cor heatmap"),
    ("F", SRC/"clinical_relevance/Fig10K_clinical_stage_axis_score_boxplot_TCGA.pdf", 0, "Stage boxplot"),
]
e = assemble_figure(s9, 3, 2,
    "Figure S9. Exploratory clinical relevance analyses", "FigureS9_clinical_relevance_exploratory")
all_errors["S9"] = e

# ── S10: Docking ──
print("\n=== S10: Docking Details ===\n")
s10 = [
    ("A", SRC/"drug_screening/Fig11B_target_axis_evidence_score_barplot.pdf", 0, "Evidence barplot"),
    ("B", SRC/"drug_screening/Fig11C_candidate_drug_target_network.pdf", 0, "Drug-target network"),
    ("C", SRC/"drug_screening/Fig11D_docking_target_shortlist.pdf", 0, "Target shortlist"),
    ("D", SRC/"docking/FigS_Docking_affinity_barplot_all_pairs.pdf", 0, "Docking affinity (all pairs)"),
]
e = assemble_figure(s10, 2, 2,
    "Figure S10. Molecular docking and drug screening supplementary", "FigureS10_docking_details")
all_errors["S10"] = e

# ── Summary ──
print("\n" + "="*60)
print("SUMMARY")
print("="*60)
total_errs = 0
for fig, errs in sorted(all_errors.items()):
    n = len(errs)
    total_errs += n
    status = f"{n} error(s)" if n else "OK"
    print(f"  {fig}: {status}")
    for e in errs:
        print(f"    - {e}")
print(f"\nTotal errors: {total_errs}")
print(f"Output: {OUT}")
