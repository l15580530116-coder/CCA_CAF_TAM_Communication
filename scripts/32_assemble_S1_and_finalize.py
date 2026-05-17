#!/usr/bin/env python
"""Step 32: Assemble Supplementary Figure S1 + finalize main_final directory."""
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
OUT_SUPP = BASE / "figures/supplementary_final"
OUT_SUPP.mkdir(parents=True, exist_ok=True)
DPI = 300
FS = 3.2  # panel size for supplementary

def pdf_to_img(path, page_idx=0):
    doc = fitz.open(str(path))
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

# ═══════════════════════════════════════════════════════════════
# Figure S1: Preprocessing and Quality Control Overview
# Panels:
#   A: TCGA PCA (T vs N) — discovery cohort overview
#   B: GSE107943 PCA (T vs N) — validation cohort overview
#   C: GSE138709 UMAP clusters — single-cell data overview
#   D: GSE26566 hub gene heatmap — microarray validation overview
# ═══════════════════════════════════════════════════════════════

print("Assembling Supplementary Figure S1...")

panels_s1 = [
    ("A", BASE / "figures/DEG/Fig2A_TCGA_PCA_tumor_normal.pdf", 0,
     "TCGA-CHOL: PCA (tumor vs normal)"),
    ("B", BASE / "figures/DEG/Fig2D_GSE107943_PCA_tumor_normal.pdf", 0,
     "GSE107943: PCA (tumor vs normal)"),
    ("C", BASE / "figures/single_cell/Fig7A_GSE138709_UMAP_clusters.pdf", 0,
     "GSE138709: UMAP (22 clusters)"),
    ("D", BASE / "figures/GSE26566_validation/FigS_GSE26566_hub_gene_heatmap.pdf", 0,
     "GSE26566: Hub gene expression"),
]

fig, axes = plt.subplots(2, 2, figsize=(8.5, 7.5), facecolor="white")
for ax_flat, (letter, path, pg, desc) in zip(axes.flat, panels_s1):
    try:
        img = pdf_to_img(path, page_idx=pg)
        ax_flat.imshow(img)
        print(f"  Panel {letter}: {desc} ({img.shape[1]}x{img.shape[0]})")
    except Exception as e:
        ax_flat.text(0.5, 0.5, f"Panel {letter}\nNot available\n{e}",
                     ha="center", va="center", fontsize=7, color="red")
        print(f"  Panel {letter}: ERROR — {e}")
    ax_flat.axis("off")
    add_label(ax_flat, letter)
    # Small description below each panel
    ax_flat.set_title(desc, fontsize=7, color="#555555", pad=3)

fig.suptitle("Figure S1. Data preprocessing and quality control overview",
             fontsize=10, fontweight="bold", y=0.99)

# Save
s1_pdf = OUT_SUPP / "FigureS1_preprocessing_QC.pdf"
s1_png = OUT_SUPP / "FigureS1_preprocessing_QC.png"
fig.savefig(s1_pdf, dpi=DPI, facecolor="white", edgecolor="none",
            bbox_inches="tight")
fig.savefig(s1_png, dpi=DPI, facecolor="white", edgecolor="none",
            bbox_inches="tight")
plt.close(fig)
print(f"  Saved: FigureS1_preprocessing_QC.pdf ({s1_pdf.stat().st_size//1024} KB)")
print(f"  Saved: FigureS1_preprocessing_QC.png ({s1_png.stat().st_size//1024} KB)")

# ═══════════════════════════════════════════════════════════════
# Final QC: Verify all 7 main figures in main_final/
# ═══════════════════════════════════════════════════════════════

print("\n=== Final Main Figure Verification ===")
main_final = BASE / "figures/main_final"
expected = [
    "Figure1_workflow_schematic",
    "Figure2_DEG_cross_validation",
    "Figure3_aggressive_score_immune",
    "Figure4_single_cell_localization",
    "Figure5_CellChat_communication",
    "Figure6_hub_genes_clinical_relevance",
    "Figure7_therapeutic_implications",
]

all_ok = True
for name in expected:
    pdf = main_final / f"{name}.pdf"
    png = main_final / f"{name}.png"
    pdf_ok = pdf.exists()
    png_ok = png.exists()
    # Quick content check
    if pdf_ok:
        doc = fitz.open(str(pdf))
        npages = len(doc)
        has_content = npages > 0
        doc.close()
    else:
        has_content = False
    status = "OK" if (pdf_ok and png_ok and has_content) else "MISSING/EMPTY"
    if status != "OK":
        all_ok = False
    print(f"  {name}: PDF={pdf_ok}, PNG={png_ok}, content={has_content} → {status}")

print(f"\nAll figures present: {all_ok}")

# Specific check: Figure 5 must be v2 (no render error, has network)
fig5 = main_final / "Figure5_CellChat_communication.pdf"
doc5 = fitz.open(str(fig5))
txt5 = doc5[0].get_text()
doc5.close()
is_v2 = "Network viz failed" not in txt5 and "inferred" in txt5.lower()
print(f"Figure 5 is v2 (network, no error): {is_v2}")

print("\nDone.")
