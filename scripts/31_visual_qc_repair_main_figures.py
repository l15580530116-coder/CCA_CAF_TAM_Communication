#!/usr/bin/env python
"""Step 31: Visual QC and repair of composite Figures 2–7.

Issues found in v7 composites:
  - Fig3C: source PDF is 2-page; page 1 (not 0) contains the KM curve
  - Fig5 C/D: source PDFs are 0-page (empty) — regenerate bubble plots from CSV
  - Fig6D: annotation text overlaps forest plot axis — reposition
  - Fig7B: annotation text overlaps x-axis labels — reposition
  - Fig2: check for stray KM/exploratory text in non-KM panels
  - Fig7C: retain existing; note BioRender recommendation
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch
import numpy as np
import csv
from pathlib import Path
import fitz
from PIL import Image
import io

BASE = Path(r"e:\CCA")
SRC = BASE / "figures"
OUT = BASE / "figures" / "main_repaired"
OUT.mkdir(parents=True, exist_ok=True)

DPI = 300
FIGSIZE = 3.4  # inches per panel

# ═══════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════

def pdf_page_to_image(pdf_path, page_idx=0, dpi=DPI):
    """Render a specific page of a PDF to RGB numpy array."""
    doc = fitz.open(str(pdf_path))
    if page_idx >= len(doc):
        print(f"  WARNING: {pdf_path.name} has {len(doc)} pages, "
              f"requested page {page_idx}")
        page_idx = len(doc) - 1
    page = doc[page_idx]
    zoom = dpi / 72.0
    mat = fitz.Matrix(zoom, zoom)
    pix = page.get_pixmap(matrix=mat, colorspace=fitz.csRGB)
    img_data = pix.tobytes("png")
    doc.close()
    return np.array(Image.open(io.BytesIO(img_data)))


def add_panel_label(ax, letter, x=0.02, y=0.96, fontsize=14, weight="bold"):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=fontsize,
            fontweight=weight, color="black", va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.1", facecolor="white",
                      edgecolor="none", alpha=0.85))


def add_caveat_annotation(ax, text, position="top-right", fontsize=7):
    """Add a small safe-wording annotation."""
    if position == "top-right":
        x, y, ha, va = 0.98, 0.98, "right", "top"
    elif position == "top-left":
        x, y, ha, va = 0.02, 0.98, "left", "top"
    else:
        x, y, ha, va = 0.98, 0.02, "right", "bottom"
    ax.text(x, y, text, transform=ax.transAxes, fontsize=fontsize,
            fontstyle="italic", color="#555555", va=va, ha=ha,
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                      edgecolor="#cccccc", alpha=0.90))


def panel_from_pdf(pdf_path, page_idx=0):
    """Load panel image from PDF. Returns (img_array, error_msg)."""
    try:
        if not pdf_path.exists():
            return None, "file missing"
        doc = fitz.open(str(pdf_path))
        if len(doc) == 0:
            doc.close()
            return None, "0-page PDF (empty)"
        if page_idx >= len(doc):
            doc.close()
            return None, f"page {page_idx} not in {len(doc)}-page doc"
        page = doc[page_idx]
        zoom = DPI / 72.0
        mat = fitz.Matrix(zoom, zoom)
        pix = page.get_pixmap(matrix=mat, colorspace=fitz.csRGB)
        img_data = pix.tobytes("png")
        doc.close()
        return np.array(Image.open(io.BytesIO(img_data))), None
    except Exception as e:
        return None, str(e)


def make_bubble_plot(rows, title, figsize=(7, 5)):
    """Generate a CellChat-style bubble plot from list of dict rows."""
    # Sort by prob descending, take top 25
    rows_sorted = sorted(rows, key=lambda r: float(r["prob"]), reverse=True)[:25]
    # Then ascending for horizontal display
    rows_sorted = sorted(rows_sorted, key=lambda r: float(r["prob"]))
    pair_labels = [f"{r['ligand']} – {r['receptor']}" for r in rows_sorted]
    probs = [float(r["prob"]) for r in rows_sorted]

    fig, ax = plt.subplots(1, 1, figsize=figsize, facecolor="white")
    sizes = [p * 800 for p in probs]

    scatter = ax.scatter(probs, pair_labels, s=sizes,
                         c=probs, cmap="Reds", edgecolors="grey",
                         linewidth=0.3, alpha=0.85, zorder=3)

    ax.set_xlabel("Communication probability", fontsize=9)
    ax.set_ylabel("")
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.tick_params(axis="y", labelsize=7)
    ax.tick_params(axis="x", labelsize=8)
    ax.grid(axis="x", alpha=0.3, lw=0.5)
    ax.set_xlim(0, max(probs) * 1.15)

    cbar = fig.colorbar(scatter, ax=ax, shrink=0.6, aspect=20, pad=0.02)
    cbar.set_label("Probability", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    plt.tight_layout()
    return fig


def fig_to_image(fig):
    """Convert matplotlib figure to numpy array."""
    fig.canvas.draw()
    img_data = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
    w, h = fig.canvas.get_width_height()
    return img_data.reshape(h, w, 4)[:, :, :3]  # drop alpha


def assemble_and_save(panels, panel_order, nrows, ncols, fig_id,
                      fig_title, output_stem, caveats=None):
    """Generic figure assembler. panels = {letter: (img_or_None, error_msg)}"""
    fig = plt.figure(figsize=(ncols * FIGSIZE, nrows * FIGSIZE),
                     facecolor="white")
    gs = gridspec.GridSpec(nrows, ncols, figure=fig,
                           wspace=0.12, hspace=0.18,
                           left=0.03, right=0.97, top=0.94, bottom=0.03)

    for idx, letter in enumerate(panel_order):
        if letter is None:
            continue
        row, col = divmod(idx, ncols)
        ax = fig.add_subplot(gs[row, col])
        img, err = panels.get(letter, (None, "not found"))

        if img is not None:
            ax.imshow(img)
            ax.axis("off")
            add_panel_label(ax, letter)

            # Add caveat annotations per figure
            if caveats and letter in caveats:
                for pos, txt in caveats[letter]:
                    add_caveat_annotation(ax, txt, position=pos, fontsize=7)
        else:
            # Error placeholder
            ax.text(0.5, 0.5, f"Panel {letter}\n{err}",
                    ha="center", va="center", fontsize=8, color="red")
            ax.axis("off")
            add_panel_label(ax, letter)

    # Hide unused slots
    total = nrows * ncols
    for idx in range(len(panel_order), total):
        row, col = divmod(idx, ncols)
        ax = fig.add_subplot(gs[row, col])
        ax.axis("off")

    # Common caveats as list at bottom if needed
    fig.suptitle(fig_title, fontsize=11, fontweight="bold", y=0.99)

    pdf_path = f"{output_stem}.pdf"
    png_path = f"{output_stem}.png"
    fig.savefig(pdf_path, dpi=DPI, facecolor="white", edgecolor="none",
                bbox_inches="tight")
    fig.savefig(png_path, dpi=DPI, facecolor="white", edgecolor="none",
                bbox_inches="tight")
    plt.close(fig)
    ks_pdf = Path(pdf_path).stat().st_size // 1024
    ks_png = Path(png_path).stat().st_size // 1024
    print(f"  Saved: {Path(pdf_path).name} ({ks_pdf} KB) + PNG ({ks_png} KB)")


# ═══════════════════════════════════════════════════════════════
# QC: Check all source PDFs
# ═══════════════════════════════════════════════════════════════

def qc_all_sources():
    """Return {fig_id: {letter: (valid, n_pages, note)}}"""
    mappings = {
        "Fig2": {
            "A": SRC / "DEG/Fig2A_TCGA_PCA_tumor_normal.pdf",
            "B": SRC / "DEG/Fig2B_TCGA_volcano.pdf",
            "C": SRC / "DEG/Fig2D_GSE107943_PCA_tumor_normal.pdf",
            "D": SRC / "DEG/Fig2E_GSE107943_paired_volcano.pdf",
            "E": SRC / "DEG/Fig2F_TCGA_GSE107943_overlap_upset_or_venn.pdf",
            "F": SRC / "enrichment/Fig3A_GO_BP_up_dotplot.pdf",
            "G": SRC / "enrichment/Fig3B_GO_BP_down_dotplot.pdf",
            "H": SRC / "enrichment/Fig3E_up_down_pathway_barplot.pdf",
        },
        "Fig3": {
            "A": SRC / "gsva/Fig4C_TCGA_score_heatmap.pdf",
            "B": SRC / "gsva/Fig4D_GSE107943_score_heatmap.pdf",
            "C": SRC / "gsva/Fig4H_GSE107943_KM_high_low_IM_CAF_TAM_score.pdf",
            "D": SRC / "gsva/Fig4I_score_correlation_heatmap.pdf",
            "E": SRC / "immune/Fig6C_TCGA_aggressive_immune_correlation_heatmap.pdf",
            "F": SRC / "immune/Fig6E_TCGA_checkpoint_correlation_heatmap.pdf",
            "G": SRC / "immune/Fig6G_aggressive_vs_macrophage_scatter_TCGA.pdf",
            "H": SRC / "immune/Fig6I_aggressive_vs_checkpoint_scatter_TCGA.pdf",
        },
        "Fig4": {
            "A": SRC / "single_cell/Fig7A_GSE138709_UMAP_clusters.pdf",
            "B": SRC / "single_cell/Fig7B_GSE138709_UMAP_celltypes.pdf",
            "C": SRC / "single_cell/Fig7D_GSE138709_aggressive_score_UMAP.pdf",
            "D": SRC / "single_cell/Fig7F_GSE138709_score_by_celltype_violin.pdf",
            "E": SRC / "single_cell/Fig7H_GSE138709_key_genes_dotplot_by_celltype.pdf",
        },
        "Fig5": {
            "A": SRC / "cellchat/Fig8C_CellChat_CAF_TAM_Epithelial_network.pdf",
            "B": SRC / "cellchat/Fig8D_CellChat_outgoing_incoming_heatmap.pdf",
            "C": SRC / "cellchat/Fig8F_CellChat_CAF_to_TAM_bubble.pdf",
            "D": SRC / "cellchat/Fig8G_CellChat_TAM_to_CAF_bubble.pdf",
        },
        "Fig6": {
            "A": SRC / "integrated/Fig9A_integrated_evidence_heatmap.pdf",
            "B": SRC / "integrated/Fig9B_top_hub_genes_barplot.pdf",
            "C": SRC / "integrated/Fig9C_prioritized_LR_axis_dotplot.pdf",
            "D": SRC / "clinical_relevance/Fig10B_axis_scores_survival_forest_GSE107943.pdf",
        },
        "Fig7": {
            "A": SRC / "drug_screening/Fig11A_target_druggability_heatmap.pdf",
            "B": SRC / "docking/FigS_Docking_affinity_barplot_all_pairs.pdf",
            "C": SRC / "integrated/Fig9D_CAF_TAM_Epithelial_mechanism_schematic.pdf",
        },
    }
    results = {}
    for fig_id, panels in mappings.items():
        results[fig_id] = {}
        for letter, path in panels.items():
            if not path.exists():
                results[fig_id][letter] = (False, 0, "MISSING")
                continue
            try:
                doc = fitz.open(str(path))
                n = len(doc)
                if n == 0:
                    results[fig_id][letter] = (False, 0, "0 pages (empty)")
                else:
                    # Check content
                    txt = doc[0].get_text()
                    has_txt = len(txt.strip()) > 20
                    results[fig_id][letter] = (True, n,
                        f"OK ({n}p, {len(txt)} chars)" if has_txt
                        else f"OK ({n}p, but text sparse)")
                doc.close()
            except Exception as e:
                results[fig_id][letter] = (False, 0, f"CORRUPT: {e}")
    return results


# ═══════════════════════════════════════════════════════════════
# FIGURE 2: Check for stray text, re-assemble clean
# ═══════════════════════════════════════════════════════════════

def repair_figure2():
    print("\n--- Figure 2: QC + re-assemble ---")
    panels = {}
    for letter, path in [
        ("A", SRC / "DEG/Fig2A_TCGA_PCA_tumor_normal.pdf"),
        ("B", SRC / "DEG/Fig2B_TCGA_volcano.pdf"),
        ("C", SRC / "DEG/Fig2D_GSE107943_PCA_tumor_normal.pdf"),
        ("D", SRC / "DEG/Fig2E_GSE107943_paired_volcano.pdf"),
        ("E", SRC / "DEG/Fig2F_TCGA_GSE107943_overlap_upset_or_venn.pdf"),
        ("F", SRC / "enrichment/Fig3A_GO_BP_up_dotplot.pdf"),
        ("G", SRC / "enrichment/Fig3B_GO_BP_down_dotplot.pdf"),
        ("H", SRC / "enrichment/Fig3E_up_down_pathway_barplot.pdf"),
    ]:
        img, err = panel_from_pdf(path)
        if img is not None:
            # Check for stray "KM" or "exploratory" text in non-KM panels
            doc = fitz.open(str(path))
            txt = doc[0].get_text()
            doc.close()
            if letter not in ("C",) and ("KM" in txt or "survival" in txt.lower()):
                print(f"  Panel {letter}: contains 'KM'/'survival' text — may be stray")
            else:
                print(f"  Panel {letter}: OK")
        else:
            print(f"  Panel {letter}: ERROR — {err}")
        panels[letter] = (img, err)

    order = ["A", "B", "C", "D", "E", "F", "G", "H"]
    # Fig2 has no KM/exploratory annotations to add for panels A-H
    assemble_and_save(panels, order, nrows=4, ncols=2, fig_id=2,
                      fig_title="Figure 2. Cross-cohort DEG validation and functional enrichment",
                      output_stem=str(OUT / "Fig2_DEG_cross_validation_repaired"))


# ═══════════════════════════════════════════════════════════════
# FIGURE 3: Fix Panel C (use page 1 of 2-page KM PDF)
# ═══════════════════════════════════════════════════════════════

def repair_figure3():
    print("\n--- Figure 3: Fix Panel C (KM curve) ---")
    panels = {}
    for letter, path, pg in [
        ("A", SRC / "gsva/Fig4C_TCGA_score_heatmap.pdf", 0),
        ("B", SRC / "gsva/Fig4D_GSE107943_score_heatmap.pdf", 0),
        ("C", SRC / "gsva/Fig4H_GSE107943_KM_high_low_IM_CAF_TAM_score.pdf", 1),
        ("D", SRC / "gsva/Fig4I_score_correlation_heatmap.pdf", 0),
        ("E", SRC / "immune/Fig6C_TCGA_aggressive_immune_correlation_heatmap.pdf", 0),
        ("F", SRC / "immune/Fig6E_TCGA_checkpoint_correlation_heatmap.pdf", 0),
        ("G", SRC / "immune/Fig6G_aggressive_vs_macrophage_scatter_TCGA.pdf", 0),
        ("H", SRC / "immune/Fig6I_aggressive_vs_checkpoint_scatter_TCGA.pdf", 0),
    ]:
        img, err = panel_from_pdf(path, page_idx=pg)
        status = "OK" if img is not None else f"ERROR: {err}"
        if letter == "C" and pg == 1:
            status += " (used page 1 of 2-page KM PDF)"
        print(f"  Panel {letter}: {status}")
        panels[letter] = (img, err)

    caveats = {
        "C": [("top-right", "All survival analyses\nare exploratory")],
    }

    order = ["A", "B", "C", "D", "E", "F", "G", "H"]
    assemble_and_save(panels, order, nrows=4, ncols=2, fig_id=3,
                      fig_title="Figure 3. Aggressive microenvironment score and immune landscape",
                      output_stem=str(OUT / "Fig3_aggressive_score_immune_repaired"),
                      caveats=caveats)


# ═══════════════════════════════════════════════════════════════
# FIGURE 4: Re-assemble (no issues found, but include for completeness)
# ═══════════════════════════════════════════════════════════════

def repair_figure4():
    print("\n--- Figure 4: Re-assemble (no known issues) ---")
    panels = {}
    for letter, path in [
        ("A", SRC / "single_cell/Fig7A_GSE138709_UMAP_clusters.pdf"),
        ("B", SRC / "single_cell/Fig7B_GSE138709_UMAP_celltypes.pdf"),
        ("C", SRC / "single_cell/Fig7D_GSE138709_aggressive_score_UMAP.pdf"),
        ("D", SRC / "single_cell/Fig7F_GSE138709_score_by_celltype_violin.pdf"),
        ("E", SRC / "single_cell/Fig7H_GSE138709_key_genes_dotplot_by_celltype.pdf"),
    ]:
        img, err = panel_from_pdf(path)
        status = "OK" if img is not None else f"ERROR: {err}"
        print(f"  Panel {letter}: {status}")
        panels[letter] = (img, err)

    order = ["A", "B", "C", "D", "E", None]
    assemble_and_save(panels, order, nrows=2, ncols=3, fig_id=4,
                      fig_title="Figure 4. Single-cell transcriptomic localization",
                      output_stem=str(OUT / "Fig4_single_cell_localization_repaired"))


# ═══════════════════════════════════════════════════════════════
# FIGURE 5: Regenerate C/D bubble plots from CSV data
# ═══════════════════════════════════════════════════════════════

def repair_figure5():
    print("\n--- Figure 5: Regenerate C/D bubble plots from CSV ---")

    # Load CellChat data using csv module
    lr_csv = BASE / "tables/cellchat/GSE138709_CAF_TAM_Epithelial_LR_pairs.csv"
    if not lr_csv.exists():
        print(f"  ERROR: {lr_csv} not found!")
        return

    with open(lr_csv, newline="") as f:
        reader = csv.DictReader(f)
        lr_all = list(reader)
    print(f"  Loaded {len(lr_all)} LR pairs from CSV")

    # Panel C: CAF → TAM
    caf_to_tam = [r for r in lr_all
                  if r["source"] == "CAF_Fibroblast"
                  and r["target"] == "Macrophage_TAM"]
    print(f"  CAF→TAM pairs: {len(caf_to_tam)}")
    fig_c = make_bubble_plot(caf_to_tam, "CAF → TAM (inferred communication)")
    img_c = fig_to_image(fig_c)
    plt.close(fig_c)

    # Panel D: TAM → CAF
    tam_to_caf = [r for r in lr_all
                  if r["source"] == "Macrophage_TAM"
                  and r["target"] == "CAF_Fibroblast"]
    print(f"  TAM→CAF pairs: {len(tam_to_caf)}")
    fig_d = make_bubble_plot(tam_to_caf, "TAM → CAF (inferred communication)")
    img_d = fig_to_image(fig_d)
    plt.close(fig_d)

    # Panel A and B from existing PDFs
    panels = {}
    for letter, path in [
        ("A", SRC / "cellchat/Fig8C_CellChat_CAF_TAM_Epithelial_network.pdf"),
        ("B", SRC / "cellchat/Fig8D_CellChat_outgoing_incoming_heatmap.pdf"),
    ]:
        img, err = panel_from_pdf(path)
        status = "OK" if img is not None else f"ERROR: {err}"
        print(f"  Panel {letter}: {status}")
        panels[letter] = (img, err)

    panels["C"] = (img_c, None)
    panels["D"] = (img_d, None)
    print(f"  Panel C: regenerated from CSV ({img_c.shape[1]}x{img_c.shape[0]})")
    print(f"  Panel D: regenerated from CSV ({img_d.shape[1]}x{img_d.shape[0]})")

    order = ["A", "B", "C", "D"]
    assemble_and_save(panels, order, nrows=2, ncols=2, fig_id=5,
                      fig_title="Figure 5. CellChat-inferred CAF–TAM–Epithelial communication network",
                      output_stem=str(OUT / "Fig5_CellChat_communication_repaired"))


# ═══════════════════════════════════════════════════════════════
# FIGURE 6: Fix Panel D annotation overlap
# ═══════════════════════════════════════════════════════════════

def repair_figure6():
    print("\n--- Figure 6: Fix Panel D annotation overlap ---")
    panels = {}
    for letter, path in [
        ("A", SRC / "integrated/Fig9A_integrated_evidence_heatmap.pdf"),
        ("B", SRC / "integrated/Fig9B_top_hub_genes_barplot.pdf"),
        ("C", SRC / "integrated/Fig9C_prioritized_LR_axis_dotplot.pdf"),
        ("D", SRC / "clinical_relevance/Fig10B_axis_scores_survival_forest_GSE107943.pdf"),
    ]:
        img, err = panel_from_pdf(path)
        status = "OK" if img is not None else f"ERROR: {err}"
        print(f"  Panel {letter}: {status}")
        panels[letter] = (img, err)

    caveats = {
        "D": [("top-left", "Exploratory — not validated\nGSE107943 only")],
    }

    order = ["A", "B", "C", "D"]
    assemble_and_save(panels, order, nrows=2, ncols=2, fig_id=6,
                      fig_title="Figure 6. Integrated hub gene prioritization and exploratory clinical associations",
                      output_stem=str(OUT / "Fig6_hub_genes_clinical_relevance_repaired"),
                      caveats=caveats)


# ═══════════════════════════════════════════════════════════════
# FIGURE 7: Fix Panel B annotation, flag Panel C for BioRender
# ═══════════════════════════════════════════════════════════════

def repair_figure7():
    print("\n--- Figure 7: Fix Panel B annotation ---")
    panels = {}
    for letter, path in [
        ("A", SRC / "drug_screening/Fig11A_target_druggability_heatmap.pdf"),
        ("B", SRC / "docking/FigS_Docking_affinity_barplot_all_pairs.pdf"),
        ("C", SRC / "integrated/Fig9D_CAF_TAM_Epithelial_mechanism_schematic.pdf"),
    ]:
        img, err = panel_from_pdf(path)
        status = "OK" if img is not None else f"ERROR: {err}"
        print(f"  Panel {letter}: {status}")
        panels[letter] = (img, err)

    caveats = {
        "B": [("top-left", "In silico docking predictions\nNot experimental binding data")],
        "C": [("top-right", "Recommend BioRender\nrefinement before publication")],
    }

    order = ["A", "B", "C"]
    assemble_and_save(panels, order, nrows=1, ncols=3, fig_id=7,
                      fig_title="Figure 7. Therapeutic target landscape and in silico docking",
                      output_stem=str(OUT / "Fig7_therapeutic_implications_repaired"))


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("=" * 60)
    print("Step 31: Visual QC and Repair of Main Figures 2–7")
    print("=" * 60)

    # --- QC Phase ---
    print("\n>>> QC PHASE: Checking all source PDFs...")
    qc = qc_all_sources()
    issues_found = 0
    for fig_id, panels in qc.items():
        for letter, (valid, npages, note) in panels.items():
            if not valid:
                issues_found += 1
                print(f"  ISSUE: {fig_id} Panel {letter}: {note}")
    print(f"  Total issues found: {issues_found}")

    # --- Repair Phase ---
    print("\n>>> REPAIR PHASE: Fixing and re-assembling...")

    repair_figure2()   # Check stray text, re-assemble
    repair_figure3()   # Fix Panel C (page 1 of 2-page KM)
    repair_figure4()   # Re-assemble (clean)
    repair_figure5()   # Regenerate C/D bubble plots
    repair_figure6()   # Fix Panel D annotation
    repair_figure7()   # Fix Panel B annotation

    # --- Final listing ---
    print(f"\n{'=' * 60}")
    print("Repaired figures in: figures/main_repaired/")
    for f in sorted(OUT.glob("*")):
        sz_kb = f.stat().st_size // 1024
        print(f"  {f.name} ({sz_kb} KB)")
    print("=" * 60)
    print("Done. QC report: manuscript/main_figure_visual_qc_report_v8.md")
    print("Repair log: manuscript/main_figure_repair_log_v8.md")
