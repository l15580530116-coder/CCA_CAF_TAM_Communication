#!/usr/bin/env python
"""Step 38: Redesign and polish Main Figures 1-7 for submission quality."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path
import fitz, io, os
from PIL import Image

BASE = Path(r"e:\CCA")
SRC = BASE / "figures"
OUT = BASE / "figures/main_final_v2"
OUT.mkdir(parents=True, exist_ok=True)
DPI = 600
LABEL_FS = 17
NOTE_FS = 8

def pdf_img(path, pg=0, dpi=DPI, max_mp=20):
    doc = fitz.open(str(path))
    if len(doc)==0: doc.close(); raise ValueError("0-page")
    if pg>=len(doc): pg=len(doc)-1
    page = doc[pg]
    pw, ph = page.rect.width, page.rect.height
    d = dpi
    if (pw*d/72)*(ph*d/72)/1e6 > max_mp:
        d = int(dpi*(max_mp/((pw*dpi/72)*(ph*dpi/72)/1e6))**0.5)
    zoom = d/72.0
    mat = fitz.Matrix(zoom, zoom)
    pix = page.get_pixmap(matrix=mat, colorspace=fitz.csRGB)
    img = np.array(Image.open(io.BytesIO(pix.tobytes("png"))))
    doc.close()
    return img

def add_label(ax, letter, fs=LABEL_FS):
    ax.text(0.015, 0.975, letter, transform=ax.transAxes, fontsize=fs,
            fontweight="bold", color="black", va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.10", facecolor="white",
                      edgecolor="none", alpha=0.88))

def add_note(ax, text, pos="bottom-right", fs=NOTE_FS):
    if pos=="bottom-right": x,y,va,ha = 0.982, 0.018, "bottom","right"
    elif pos=="top-right": x,y,va,ha = 0.982, 0.975, "top","right"
    else: x,y,va,ha = 0.982, 0.018, "bottom","right"
    ax.text(x,y,text,transform=ax.transAxes,fontsize=fs,fontstyle="italic",
            color="#555555",va=va,ha=ha,
            bbox=dict(boxstyle="round,pad=0.10",facecolor="white",
                      edgecolor="none",alpha=0.88))

def save_fig(fig, stem):
    for ext in ["pdf","png"]:
        p = OUT / f"{stem}.{ext}"
        fig.savefig(p, dpi=DPI, facecolor="white", edgecolor="none", bbox_inches="tight")
        print(f"  {stem}.{ext} ({p.stat().st_size//1024} KB)")

def assemble(panels, nrows, ncols, title, stem, notes=None, figsize_fs=4.0):
    """panels=[(letter, path, pg, desc, img_override_or_None)]"""
    fig = plt.figure(figsize=(ncols*figsize_fs, nrows*figsize_fs), facecolor="white")
    gs = gridspec.GridSpec(nrows, ncols, figure=fig, wspace=0.10, hspace=0.20,
                           left=0.03, right=0.97, top=0.935, bottom=0.03)
    for idx, item in enumerate(panels):
        letter, path, pg, desc = item[:4]
        override = item[4] if len(item)>4 else None
        row, col = divmod(idx, ncols)
        ax = fig.add_subplot(gs[row, col])
        try:
            img = override if override is not None else pdf_img(path, pg=pg)
            ax.imshow(img); ax.axis("off")
            add_label(ax, letter)
            if notes and letter in notes:
                add_note(ax, notes[letter])
        except Exception as e:
            ax.text(0.5,0.5,f"Panel {letter}\n{e}",ha="center",va="center",fontsize=7,color="red")
            ax.axis("off")
    for idx in range(len(panels), nrows*ncols):
        ax=fig.add_subplot(gs[idx//ncols, idx%ncols]); ax.axis("off")
    fig.suptitle(title, fontsize=12, fontweight="bold", y=0.99)
    save_fig(fig, stem)
    plt.close(fig)

def fig_to_img(fig):
    fig.canvas.draw()
    img = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
    w, h = fig.canvas.get_width_height()
    return img.reshape(h, w, 4)[:,:,:3]

# ═══════════════════════════════════════════════════════════════
# FIGURE 1: Clean 6-module horizontal workflow
# ═══════════════════════════════════════════════════════════════
print("=== Figure 1: Redesigned workflow ===")
C = {"data":"#4472C4","deg":"#2F9B8E","score":"#ED7D31","single":"#9B59B6",
     "cellchat":"#C0392B","thera":"#27AE60","arrow":"#666666"}

fig1, ax1 = plt.subplots(1,1,figsize=(24,7),facecolor="white")
ax1.set_xlim(0,24); ax1.set_ylim(0,7); ax1.axis("off")

modules = [
    ("A","Data\nCollection",(1.8,4.0,3.2,5.6),C["data"],
     ["4 Public Cohorts","TCGA-CHOL (n=44)","GSE107943 (n=57)","GSE138709 (32,626 cells)","GSE26566 (n=169)"]),
    ("B","DEG\nValidation",(5.3,4.0,3.2,5.6),C["deg"],
     ["DESeq2 + edgeR","4,534 same-direction","344 validated","IM/CAF/TAM DEGs","99.4% concordance"]),
    ("C","Score &\nImmune Analysis",(8.8,4.0,3.2,5.6),C["score"],
     ["Aggressive score","= ECM+CAF+TAM-Metab","CAF rho=0.77","HAVCR2 rho=0.54","Survival: exploratory"]),
    ("D","Single-Cell\nLocalization",(12.3,4.0,3.2,5.6),C["single"],
     ["32,626 cells, 8 types","CAF: highest score","TAM: checkpoint expr.","HAVCR2,IDO1,CD86","mRNA-level only"]),
    ("E","CellChat\nInference",(15.8,4.0,3.2,5.6),C["cellchat"],
     ["CAF-TAM-Epi network","COL1A1/2-CD44","MIF-CD74/CXCR4","TGFB1-TGFBR","Inferred only"]),
    ("F","Candidate\nTargets",(19.3,4.0,3.2,5.6),C["thera"],
     ["Hub: COL1A1,TREM2","14 druggable targets","CXCR4/IDO1 docking","Plerixafor -8.01","In silico only"]),
]

for letter, title, (cx,cy,w,h), color, lines in modules:
    x0, y0 = cx-w/2, cy-h/2
    box = FancyBboxPatch((x0,y0),w,h, boxstyle="round,pad=0.1,rounding_size=0.2",
                          facecolor=color, edgecolor="white", lw=2, alpha=0.16)
    ax1.add_patch(box)
    circle = Circle((cx, cy+h/2-0.5), 0.32, color=color, ec="white", lw=1.5, zorder=5)
    ax1.add_patch(circle)
    ax1.text(cx, cy+h/2-0.5, letter, ha="center", va="center", fontsize=13,
             fontweight="bold", color="white")
    ax1.text(cx, cy+h/2-1.05, title, ha="center", va="top", fontsize=9.5,
             fontweight="bold", color=color)
    for i, line in enumerate(lines):
        ax1.text(cx, cy+h/2-1.55-i*0.32, line, ha="center", va="top",
                 fontsize=8, color="#333333")

# Arrows
for i in range(len(modules)-1):
    x1 = modules[i][2][0]+modules[i][2][2]/2+0.05
    x2 = modules[i+1][2][0]-modules[i+1][2][2]/2-0.05
    ax1.annotate("", xy=(x2,4.0), xytext=(x1,4.0),
                 arrowprops=dict(arrowstyle="->",color=C["arrow"],lw=2.5))

ax1.text(12, 1.65, "CellChat and docking results are computational predictions; survival analyses are exploratory.",
         ha="center", fontsize=7.5, fontstyle="italic", color="#888888")
ax1.text(12, 6.9, "Figure 1. Study design and analytical workflow",
         ha="center", fontsize=13, fontweight="bold", color="#222222")
save_fig(fig1, "Figure1_workflow_schematic_v2")
plt.close(fig1)

# ═══════════════════════════════════════════════════════════════
# FIGURE 2E: Graphical overlap/Venn panel
# ═══════════════════════════════════════════════════════════════
print("\n=== Figure 2E: Graphical overlap panel ===")
fig2e, ax2e = plt.subplots(1,1,figsize=(6,5),facecolor="white")
ax2e.set_xlim(0,10); ax2e.set_ylim(0,10); ax2e.axis("off")

# Two circles (Venn-style)
left = Circle((3.65,5.8), 2.2, facecolor="#4472C4", edgecolor="white", lw=2, alpha=0.25)
right = Circle((6.35,5.8), 2.0, facecolor="#ED7D31", edgecolor="white", lw=2, alpha=0.25)
ax2e.add_patch(left); ax2e.add_patch(right)

ax2e.text(2.3,7.5,"TCGA-CHOL",ha="center",fontsize=11,fontweight="bold",color="#4472C4")
ax2e.text(2.3,7.1,"9,291 DEGs",ha="center",fontsize=9.5,color="#333333")
ax2e.text(2.3,6.7,"5,932 up",ha="center",fontsize=8.5,color="#666666")
ax2e.text(2.3,6.35,"3,359 down",ha="center",fontsize=8.5,color="#666666")

ax2e.text(7.7,7.5,"GSE107943",ha="center",fontsize=11,fontweight="bold",color="#ED7D31")
ax2e.text(7.7,7.1,"7,735 DEGs",ha="center",fontsize=9.5,color="#333333")
ax2e.text(7.7,6.7,"3,869 up",ha="center",fontsize=8.5,color="#666666")
ax2e.text(7.7,6.35,"3,866 down",ha="center",fontsize=8.5,color="#666666")

# Overlap center
ax2e.text(5.0,5.85,"4,560",ha="center",fontsize=14,fontweight="bold",color="#2F9B8E")
ax2e.text(5.0,5.45,"common DEGs",ha="center",fontsize=9,color="#2F9B8E")
ax2e.text(5.0,5.1,"4,534 (99.4%)",ha="center",fontsize=9.5,fontweight="bold",color="#C0392B")
ax2e.text(5.0,4.8,"same direction",ha="center",fontsize=8,color="#C0392B")

# Bottom bar
bar = FancyBboxPatch((1.0,2.2), 8.0, 1.6, boxstyle="round,pad=0.1,rounding_size=0.15",
                      facecolor="#27AE60", edgecolor="white", lw=1.5, alpha=0.20)
ax2e.add_patch(bar)
ax2e.text(5.0,3.3,"344 Validated IM/CAF/TAM DEGs",ha="center",fontsize=11,fontweight="bold",color="#27AE60")
ax2e.text(5.0,2.85,"78 upregulated | 266 downregulated",ha="center",fontsize=9,color="#333333")
ax2e.text(5.0,2.5,"ECM organization | Collagen metabolism | Lipid/Steroid/Amino acid metabolism",ha="center",fontsize=7.5,color="#666666")

ax2e.text(5.0,1.2,"Validated in both cohorts with concordant direction",ha="center",fontsize=9,fontstyle="italic",color="#888888")
ax2e.set_title("Cross-cohort DEG overlap",fontsize=12,fontweight="bold")
plt.tight_layout()
img2e = fig_to_img(fig2e)
plt.close(fig2e)
print("  Figure 2E: graphical Venn generated")

# ═══════════════════════════════════════════════════════════════
# FIGURE 2: Assemble A-H with graphical E
# ═══════════════════════════════════════════════════════════════
print("\n=== Figure 2: Re-assemble with graphical E ===")
assemble([
    ("A", SRC/"DEG/Fig2A_TCGA_PCA_tumor_normal.pdf",0,"TCGA PCA"),
    ("B", SRC/"DEG/Fig2B_TCGA_volcano.pdf",0,"TCGA volcano"),
    ("C", SRC/"DEG/Fig2D_GSE107943_PCA_tumor_normal.pdf",0,"GSE PCA"),
    ("D", SRC/"DEG/Fig2E_GSE107943_paired_volcano.pdf",0,"GSE volcano"),
    ("E", None,0,"Overlap Venn", img2e),
    ("F", SRC/"enrichment/Fig3A_GO_BP_up_dotplot.pdf",0,"GO BP up"),
    ("G", SRC/"enrichment/Fig3B_GO_BP_down_dotplot.pdf",0,"GO BP down"),
    ("H", SRC/"enrichment/Fig3E_up_down_pathway_barplot.pdf",0,"Pathway barplot"),
], 4, 2, "Figure 2. Cross-cohort DEG validation and functional enrichment",
    "Figure2_DEG_cross_validation_v2")

# ═══════════════════════════════════════════════════════════════
# FIGURE 3: Fix annotation overlap
# ═══════════════════════════════════════════════════════════════
print("\n=== Figure 3: Fix annot overlap ===")
assemble([
    ("A", SRC/"gsva/Fig4C_TCGA_score_heatmap.pdf",0,"TCGA heatmap"),
    ("B", SRC/"gsva/Fig4D_GSE107943_score_heatmap.pdf",0,"GSE heatmap"),
    ("C", SRC/"gsva/Fig4H_GSE107943_KM_high_low_IM_CAF_TAM_score.pdf",1,"GSE KM"),
    ("D", SRC/"gsva/Fig4I_score_correlation_heatmap.pdf",0,"Score cor"),
    ("E", SRC/"immune/Fig6C_TCGA_aggressive_immune_correlation_heatmap.pdf",0,"Immune cor"),
    ("F", SRC/"immune/Fig6E_TCGA_checkpoint_correlation_heatmap.pdf",0,"Checkpoint cor"),
    ("G", SRC/"immune/Fig6G_aggressive_vs_macrophage_scatter_TCGA.pdf",0,"Macro scatter"),
    ("H", SRC/"immune/Fig6I_aggressive_vs_checkpoint_scatter_TCGA.pdf",0,"Checkpoint scatter"),
], 4, 2, "Figure 3. Aggressive microenvironment score and immune landscape",
    "Figure3_aggressive_score_immune_v2",
    notes={"C":"All survival analyses\nare exploratory"})

# ═══════════════════════════════════════════════════════════════
# FIGURE 4: Better layout
# ═══════════════════════════════════════════════════════════════
print("\n=== Figure 4: Improved layout ===")
fig4 = plt.figure(figsize=(16, 9.5), facecolor="white")
# Top row: 3 UMAPs
gs4 = gridspec.GridSpec(2, 3, figure=fig4, wspace=0.08, hspace=0.18,
                        left=0.03, right=0.97, top=0.935, bottom=0.03)

for idx, (letter, path, pg, desc) in enumerate([
    ("A", SRC/"single_cell/Fig7A_GSE138709_UMAP_clusters.pdf",0,"Clusters"),
    ("B", SRC/"single_cell/Fig7B_GSE138709_UMAP_celltypes.pdf",0,"Cell types"),
    ("C", SRC/"single_cell/Fig7D_GSE138709_aggressive_score_UMAP.pdf",0,"Aggressive score"),
]):
    ax = fig4.add_subplot(gs4[0, idx])
    img = pdf_img(path, pg=pg)
    ax.imshow(img); ax.axis("off")
    add_label(ax, letter)

# Bottom row: D violin (left 1 col), E dotplot (right 2 cols)
ax_d = fig4.add_subplot(gs4[1, 0])
img_d = pdf_img(SRC/"single_cell/Fig7F_GSE138709_score_by_celltype_violin.pdf")
ax_d.imshow(img_d); ax_d.axis("off")
add_label(ax_d, "D")

# E spans columns 1-2
ax_e = fig4.add_subplot(gs4[1, 1:])
img_e = pdf_img(SRC/"single_cell/Fig7H_GSE138709_key_genes_dotplot_by_celltype.pdf")
ax_e.imshow(img_e); ax_e.axis("off")
add_label(ax_e, "E")

fig4.suptitle("Figure 4. Single-cell transcriptomic localization",
              fontsize=12, fontweight="bold", y=0.99)
save_fig(fig4, "Figure4_single_cell_localization_v2")
plt.close(fig4)

# ═══════════════════════════════════════════════════════════════
# FIGURE 5: Copy v2 (already good)
# ═══════════════════════════════════════════════════════════════
print("\n=== Figure 5: Copy existing v2 ===")
for ext in ["pdf","png"]:
    src = BASE / f"submission_CSBJ_UPLOAD_READY/Figures_Main/Figure5_CellChat_communication.{ext}"
    dst = OUT / f"Figure5_CellChat_communication_v2.{ext}"
    import shutil; shutil.copy2(src, dst)
    print(f"  {dst.name} ({dst.stat().st_size//1024} KB)")

# ═══════════════════════════════════════════════════════════════
# FIGURE 6: Fix label and annotation overlap
# ═══════════════════════════════════════════════════════════════
print("\n=== Figure 6: Fix overlap ===")
assemble([
    ("A", SRC/"integrated/Fig9A_integrated_evidence_heatmap.pdf",0,"Evidence heatmap"),
    ("B", SRC/"integrated/Fig9B_top_hub_genes_barplot.pdf",0,"Hub gene barplot"),
    ("C", SRC/"integrated/Fig9C_prioritized_LR_axis_dotplot.pdf",0,"LR axis dotplot"),
    ("D", SRC/"clinical_relevance/Fig10B_axis_scores_survival_forest_GSE107943.pdf",0,"Forest plot"),
], 2, 2, "Figure 6. Integrated hub gene prioritization and exploratory clinical associations",
    "Figure6_hub_genes_clinical_relevance_v2",
    notes={"D":"Exploratory; GSE107943 only; not validated"})

# ═══════════════════════════════════════════════════════════════
# FIGURE 7: Expand to 4 panels
# ═══════════════════════════════════════════════════════════════
print("\n=== Figure 7: Expand to 4 panels ===")

# Panel B: candidate drug summary (use existing drug network)
# Panel D: mechanism schematic (draw clean version)

fig7d, ax7d = plt.subplots(1,1,figsize=(7,6),facecolor="white")
ax7d.set_xlim(0,10); ax7d.set_ylim(0,10); ax7d.axis("off")

# Three cell types as boxes
nodes = [
    ("CAF\nFibroblast", 1.5, 5.5, "#E74C3C"),
    ("Macrophage\nTAM", 5.0, 5.5, "#3498DB"),
    ("Epithelial\nlike", 8.5, 5.5, "#27AE60"),
]
for label, x, y, color in nodes:
    box = FancyBboxPatch((x-1.0,y-0.5), 2.0, 1.0,
                          boxstyle="round,pad=0.08,rounding_size=0.12",
                          facecolor=color, edgecolor="white", lw=1.5, alpha=0.30)
    ax7d.add_patch(box)
    ax7d.text(x, y, label, ha="center", va="center", fontsize=8, fontweight="bold", color=color)

# Arrows with labels
arrows = [
    (2.5,5.5,4.0,5.5,"COL1A1/2-CD44\nMIF-CD74/CXCR4"),
    (4.0,5.0,2.5,5.0,"TGFB1-TGFBR\nPPIA-BSG"),
    (8.5,6.0,6.0,6.0,"MIF-CD74/CXCR4"),
]
for x1,y1,x2,y2,label in arrows:
    ax7d.annotate("", xy=(x2,y2), xytext=(x1,y1),
                  arrowprops=dict(arrowstyle="->",color="#666666",lw=2.0))
    mx, my = (x1+x2)/2, (y1+y2)/2 + 0.25
    ax7d.text(mx, my, label, ha="center", va="bottom", fontsize=6.5, color="#333333",
              bbox=dict(boxstyle="round,pad=0.08",facecolor="white",edgecolor="none",alpha=0.85))

# TAM checkpoint box
ax7d.text(5.0, 3.5, "TAM checkpoint transcripts:\nHAVCR2, IDO1, CD86, PDCD1LG2",
          ha="center", fontsize=7.5, color="#3498DB", fontweight="bold",
          bbox=dict(boxstyle="round,pad=0.2",facecolor="white",edgecolor="#3498DB",lw=1,alpha=0.9))

# Candidate targets at bottom
ax7d.text(5.0, 1.8, "Candidate targets: CXCR4, IDO1, TGFBR1, MIF",
          ha="center", fontsize=7.5, color="#333333",
          bbox=dict(boxstyle="round,pad=0.15",facecolor="#EEEEEE",edgecolor="none",alpha=0.9))
ax7d.text(5.0, 1.2, "Docking = in silico predictions only",
          ha="center", fontsize=7, fontstyle="italic", color="#888888")

ax7d.set_title("Proposed CAF-TAM-Epithelial mechanism", fontsize=10, fontweight="bold")
plt.tight_layout()
img7d = fig_to_img(fig7d)
plt.close(fig7d)

# Assemble Figure 7: A(druggability), B(drug network), C(docking barplot), D(mechanism)
assemble([
    ("A", SRC/"drug_screening/Fig11A_target_druggability_heatmap.pdf",0,"Druggability"),
    ("B", SRC/"drug_screening/Fig11C_candidate_drug_target_network.pdf",0,"Drug-target network"),
    ("C", SRC/"docking/FigS_Docking_affinity_barplot_all_pairs.pdf",0,"Docking affinity"),
    ("D", None,0,"Mechanism model", img7d),
], 2, 2, "Figure 7. Candidate target prioritization and in silico docking",
    "Figure7_therapeutic_implications_v2",
    notes={"C":"In silico predictions only"})

# ═══════════════════════════════════════════════════════════════
print("\n"+"="*60)
print("ALL MAIN FIGURES V2 GENERATED")
print(f"Output: {OUT}")
for f in sorted(OUT.glob("*.pdf")):
    print(f"  {f.name} ({f.stat().st_size//1024} KB)")
print("="*60)
