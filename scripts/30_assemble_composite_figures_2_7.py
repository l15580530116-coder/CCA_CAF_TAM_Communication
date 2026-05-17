#!/usr/bin/env python3
"""Assemble composite Figures 2–7 from individual panel PDFs using matplotlib + pymupdf."""
import fitz  # pymupdf
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
from pathlib import Path
import numpy as np
from PIL import Image
import io

BASE = Path(r"e:\CCA")
OUT = BASE / "figures" / "main"
OUT.mkdir(parents=True, exist_ok=True)

DPI_RENDER = 300  # resolution for PDF→image conversion
FIGSIZE_SCALE = 3.5  # per-panel size multiplier

# ── Panel mapping ──────────────────────────────────────────────
# From final_figure_panel_mapping_v3.csv

FIG2_PANELS = {
    "A": BASE / "figures/DEG/Fig2A_TCGA_PCA_tumor_normal.pdf",
    "B": BASE / "figures/DEG/Fig2B_TCGA_volcano.pdf",
    "C": BASE / "figures/DEG/Fig2D_GSE107943_PCA_tumor_normal.pdf",
    "D": BASE / "figures/DEG/Fig2E_GSE107943_paired_volcano.pdf",
    "E": BASE / "figures/DEG/Fig2F_TCGA_GSE107943_overlap_upset_or_venn.pdf",
    "F": BASE / "figures/enrichment/Fig3A_GO_BP_up_dotplot.pdf",
    "G": BASE / "figures/enrichment/Fig3B_GO_BP_down_dotplot.pdf",
    "H": BASE / "figures/enrichment/Fig3E_up_down_pathway_barplot.pdf",
}

FIG3_PANELS = {
    "A": BASE / "figures/gsva/Fig4C_TCGA_score_heatmap.pdf",
    "B": BASE / "figures/gsva/Fig4D_GSE107943_score_heatmap.pdf",
    "C": BASE / "figures/gsva/Fig4H_GSE107943_KM_high_low_IM_CAF_TAM_score.pdf",
    "D": BASE / "figures/gsva/Fig4I_score_correlation_heatmap.pdf",
    "E": BASE / "figures/immune/Fig6C_TCGA_aggressive_immune_correlation_heatmap.pdf",
    "F": BASE / "figures/immune/Fig6E_TCGA_checkpoint_correlation_heatmap.pdf",
    "G": BASE / "figures/immune/Fig6G_aggressive_vs_macrophage_scatter_TCGA.pdf",
    "H": BASE / "figures/immune/Fig6I_aggressive_vs_checkpoint_scatter_TCGA.pdf",
}

FIG4_PANELS = {
    "A": BASE / "figures/single_cell/Fig7A_GSE138709_UMAP_clusters.pdf",
    "B": BASE / "figures/single_cell/Fig7B_GSE138709_UMAP_celltypes.pdf",
    "C": BASE / "figures/single_cell/Fig7D_GSE138709_aggressive_score_UMAP.pdf",
    "D": BASE / "figures/single_cell/Fig7F_GSE138709_score_by_celltype_violin.pdf",
    "E": BASE / "figures/single_cell/Fig7H_GSE138709_key_genes_dotplot_by_celltype.pdf",
}

FIG5_PANELS = {
    "A": BASE / "figures/cellchat/Fig8C_CellChat_CAF_TAM_Epithelial_network.pdf",
    "B": BASE / "figures/cellchat/Fig8D_CellChat_outgoing_incoming_heatmap.pdf",
    "C": BASE / "figures/cellchat/Fig8F_CellChat_CAF_to_TAM_bubble.pdf",
    "D": BASE / "figures/cellchat/Fig8G_CellChat_TAM_to_CAF_bubble.pdf",
}

FIG6_PANELS = {
    "A": BASE / "figures/integrated/Fig9A_integrated_evidence_heatmap.pdf",
    "B": BASE / "figures/integrated/Fig9B_top_hub_genes_barplot.pdf",
    "C": BASE / "figures/integrated/Fig9C_prioritized_LR_axis_dotplot.pdf",
    "D": BASE / "figures/clinical_relevance/Fig10B_axis_scores_survival_forest_GSE107943.pdf",
}

FIG7_PANELS = {
    "A": BASE / "figures/drug_screening/Fig11A_target_druggability_heatmap.pdf",
    "B": BASE / "figures/docking/FigS_Docking_affinity_barplot_all_pairs.pdf",
    "C": BASE / "figures/integrated/Fig9D_CAF_TAM_Epithelial_mechanism_schematic.pdf",
}

# ── Labels/annotations per figure ──
FIG_LABELS = {
    2: {
        "E": "Intersection size = 4,560\n4,534 (99.4%) concordant direction",
        "F": "GO BP Upregulated DEGs",
        "G": "GO BP Downregulated DEGs",
        "H": "Up vs Down Pathway Comparison",
        "annotations": {
            "C": "GSE107943 KM (exploratory)",
        }
    },
    3: {
        "C": "GSE107943 OS\nExploratory",
        "annotations": {
            "C": "All survival analyses are exploratory",
        }
    },
    4: {},
    5: {
        "annotations": {
            "title": "Inferred communication (CellChat v2.2)"
        }
    },
    6: {
        "D": "GSE107943 — Exploratory",
        "annotations": {
            "D": "Exploratory — not validated"
        }
    },
    7: {
        "B": "In silico docking predictions",
        "annotations": {
            "B": "Computational prediction; not experimental binding affinity"
        }
    },
}

# ── Helper: render PDF page to numpy array ──
def pdf_to_image(pdf_path, dpi=DPI_RENDER):
    """Render first page of PDF to RGB numpy array at given DPI."""
    doc = fitz.open(str(pdf_path))
    page = doc[0]
    # Calculate matrix for desired DPI
    zoom = dpi / 72.0
    mat = fitz.Matrix(zoom, zoom)
    pix = page.get_pixmap(matrix=mat, colorspace=fitz.csRGB)
    img_data = pix.tobytes("png")
    doc.close()
    img = Image.open(io.BytesIO(img_data))
    return np.array(img)


def add_panel_label(ax, letter, x=0.02, y=0.96, fontsize=14, color="black",
                    weight="bold"):
    """Add bold uppercase panel letter in top-left corner."""
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=fontsize,
            fontweight=weight, color=color, va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.1", facecolor="white",
                      edgecolor="none", alpha=0.85))


def assemble_figure(panels_dict, fig_id, nrows, ncols, panel_order,
                    figsize_per_panel=FIGSIZE_SCALE, title=None,
                    save_pdf=True, save_png=True):
    """Generic composite figure assembler.

    Parameters
    ----------
    panels_dict : dict
        {letter: Path} mapping for all available panel PDFs
    fig_id : int
        Figure number (2–7)
    nrows, ncols : int
        Grid dimensions
    panel_order : list of str or None
        Ordered list of panel letters in grid positions; None → use
        panels_dict keys in order. Use None for empty slots.
    figsize_per_panel : float
        Inches per panel dimension
    title : str or None
        Optional super-title
    """
    letters = [k for k in panel_order if k is not None]
    n_panels = len(letters)

    fig = plt.figure(
        fig_id, figsize=(ncols * figsize_per_panel, nrows * figsize_per_panel),
        facecolor="white"
    )
    gs = gridspec.GridSpec(nrows, ncols, figure=fig,
                           wspace=0.12, hspace=0.15,
                           left=0.03, right=0.97, top=0.95, bottom=0.03)

    for idx, letter in enumerate(panel_order):
        if letter is None:
            continue  # empty slot
        pdf_path = panels_dict[letter]
        if not pdf_path.exists():
            print(f"  WARNING: {pdf_path} not found — skipping panel {letter}")
            continue

        row, col = divmod(idx, ncols)
        ax = fig.add_subplot(gs[row, col])

        # Render
        try:
            img = pdf_to_image(pdf_path)
        except Exception as e:
            print(f"  ERROR rendering panel {letter}: {e}")
            ax.text(0.5, 0.5, f"Panel {letter}\n(render error)",
                    ha="center", va="center", fontsize=8, color="red")
            ax.axis("off")
            continue

        ax.imshow(img)
        ax.axis("off")

        # Panel letter label
        add_panel_label(ax, letter)

        # Annotation overlay (from FIG_LABELS)
        label_info = FIG_LABELS.get(fig_id, {})
        annot = label_info.get("annotations", {})
        if letter in annot:
            ax.text(0.98, 0.02, annot[letter], transform=ax.transAxes,
                    fontsize=7, fontstyle="italic", color="#666666",
                    va="bottom", ha="right")

    # Remove unused subplots
    total_slots = nrows * ncols
    for idx in range(len(panel_order), total_slots):
        row, col = divmod(idx, ncols)
        ax = fig.add_subplot(gs[row, col])
        ax.axis("off")
        ax.set_visible(False)

    # Titles
    fig_titles = {
        2: "Figure 2. Cross-cohort DEG validation and functional enrichment",
        3: "Figure 3. Aggressive microenvironment score and immune landscape",
        4: "Figure 4. Single-cell transcriptomic localization",
        5: "Figure 5. CellChat-inferred CAF–TAM–Epithelial communication network",
        6: "Figure 6. Integrated hub gene prioritization and exploratory clinical associations",
        7: "Figure 7. Therapeutic target landscape and in silico docking",
    }
    suptitle = fig_titles.get(fig_id, "")
    if suptitle:
        fig.suptitle(suptitle, fontsize=11, fontweight="bold", y=0.99)

    # Save
    stem = OUT / f"Fig{fig_id}"
    suffixes = []
    if save_pdf:
        fig.savefig(f"{stem}.pdf", dpi=DPI_RENDER, facecolor="white",
                    edgecolor="none", bbox_inches="tight")
        suffixes.append("PDF")
    if save_png:
        fig.savefig(f"{stem}.png", dpi=DPI_RENDER, facecolor="white",
                    edgecolor="none", bbox_inches="tight")
        suffixes.append("PNG")

    plt.close(fig)
    print(f"  Figure {fig_id} saved: {', '.join(suffixes)}")
    return True


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("=" * 60)
    print("Assembling composite Figures 2–7")
    print("=" * 60)

    # ── Figure 2: 4 rows × 2 columns ──
    print("\nFigure 2 — DEG cross-validation (A–H)")
    fig2_order = ["A", "B", "C", "D", "E", "F", "G", "H"]
    assemble_figure(FIG2_PANELS, 2, nrows=4, ncols=2,
                    panel_order=fig2_order)

    # ── Figure 3: 4 rows × 2 columns ──
    print("\nFigure 3 — Aggressive score + immune (A–H)")
    fig3_order = ["A", "B", "C", "D", "E", "F", "G", "H"]
    assemble_figure(FIG3_PANELS, 3, nrows=4, ncols=2,
                    panel_order=fig3_order)

    # ── Figure 4: 2 rows × 3 columns, last slot empty ──
    print("\nFigure 4 — Single-cell localization (A–E)")
    fig4_order = ["A", "B", "C", "D", "E", None]
    assemble_figure(FIG4_PANELS, 4, nrows=2, ncols=3,
                    panel_order=fig4_order)

    # ── Figure 5: 2 rows × 2 columns ──
    print("\nFigure 5 — CellChat communication (A–D)")
    fig5_order = ["A", "B", "C", "D"]
    assemble_figure(FIG5_PANELS, 5, nrows=2, ncols=2,
                    panel_order=fig5_order)

    # ── Figure 6: 2 rows × 2 columns ──
    print("\nFigure 6 — Hub genes + clinical relevance (A–D)")
    fig6_order = ["A", "B", "C", "D"]
    assemble_figure(FIG6_PANELS, 6, nrows=2, ncols=2,
                    panel_order=fig6_order)

    # ── Figure 7: 1 row × 3 columns ──
    print("\nFigure 7 — Therapeutic implications (A–C)")
    fig7_order = ["A", "B", "C"]
    assemble_figure(FIG7_PANELS, 7, nrows=1, ncols=3,
                    panel_order=fig7_order,
                    figsize_per_panel=3.0)

    print(f"\n{'=' * 60}")
    print("All figures assembled.")
    print(f"Output directory: {OUT}")
    for i in range(2, 8):
        print(f"  Fig{i}_*.pdf / Fig{i}_*.png")
    print("=" * 60)
