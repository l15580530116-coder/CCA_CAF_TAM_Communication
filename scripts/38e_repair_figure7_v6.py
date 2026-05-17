#!/usr/bin/env python
"""Step 38e: Repair Figure 7 - fix network, barplot, enlarge mechanism."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path
import fitz, csv, shutil
from PIL import Image
import io

BASE = Path(r"e:\CCA")
SRC = BASE / "figures"
TBL = BASE / "tables"
SRC_V5 = BASE / "figures/main_final_v5"
OUT = BASE / "figures/main_final_v6"
OUT.mkdir(parents=True, exist_ok=True)
DPI = 600
LF = 18  # panel label size
NF = 8

def pdf_img(path, pg=0, dpi=DPI, max_mp=22):
    doc = fitz.open(str(path))
    if len(doc)==0: doc.close(); raise ValueError("0-page")
    if pg>=len(doc): pg=len(doc)-1
    page = doc[pg]; pw,ph = page.rect.width, page.rect.height
    d=dpi
    if (pw*d/72)*(ph*d/72)/1e6>max_mp: d=int(dpi*(max_mp/((pw*dpi/72)*(ph*dpi/72)/1e6))**0.5)
    zoom=d/72.0; mat=fitz.Matrix(zoom,zoom)
    pix=page.get_pixmap(matrix=mat,colorspace=fitz.csRGB)
    img=np.array(Image.open(io.BytesIO(pix.tobytes("png")))); doc.close()
    return img

def add_label(ax, letter, fs=LF):
    ax.text(0.012,0.978,letter,transform=ax.transAxes,fontsize=fs,fontweight="bold",
            color="black",va="top",ha="left",
            bbox=dict(boxstyle="round,pad=0.08",facecolor="white",edgecolor="none",alpha=0.88))

def add_note(ax, text, fs=NF):
    ax.text(0.985,0.015,text,transform=ax.transAxes,fontsize=fs,fontstyle="italic",
            color="#666666",va="bottom",ha="right",
            bbox=dict(boxstyle="round,pad=0.08",facecolor="white",edgecolor="none",alpha=0.88))

def fig_to_img(fig):
    fig.canvas.draw(); img=np.frombuffer(fig.canvas.buffer_rgba(),dtype=np.uint8)
    w,h=fig.canvas.get_width_height(); return img.reshape(h,w,4)[:,:,:3]

# ═══════════════════════════════════════════════════════════════
# LOAD DATA
# ═══════════════════════════════════════════════════════════════

# Drug screening data
with open(TBL/"drug_screening/candidate_drugs_prioritized.csv",newline="") as f:
    drugs = list(csv.DictReader(f))

print(f"Loaded {len(drugs)} drug-target records")

# Per-target drug count (Panel C)
tgt_counts = {}
tgt_drugs = {}
for r in drugs:
    t = r["target_gene"]
    d = r["drug_name"]
    tgt_counts[t] = tgt_counts.get(t, 0) + 1
    if t not in tgt_drugs: tgt_drugs[t] = []
    tgt_drugs[t].append(d)

print("Targets per drug count:")
for t, c in sorted(tgt_counts.items(), key=lambda x: -x[1]):
    print(f"  {t}: {c} ({', '.join(tgt_drugs[t])})")

# Docking residues
with open(TBL/"docking/docking_neighbor_residues_4A.csv",newline="") as f:
    residues = list(csv.DictReader(f))

pair_counts = {}
for r in residues:
    p = r.get("pair","?")
    pair_counts[p] = pair_counts.get(p, 0) + 1

print(f"\nResidues per pair: {pair_counts}, total={sum(pair_counts.values())}")

# ═══════════════════════════════════════════════════════════════
# PANEL B: True bipartite drug-target network
# ═══════════════════════════════════════════════════════════════
print("\n=== Panel B: Bipartite drug-target network ===")

fig_b, ax_b = plt.subplots(1,1,figsize=(8,7),facecolor="white")
ax_b.set_xlim(0,12); ax_b.set_ylim(0,14); ax_b.axis("off")

# Top 14 relationships (by total_priority)
sorted_d = sorted(drugs, key=lambda r: -int(r.get("total_priority","0")))[:14]

# Left column: targets, Right column: drugs
targets_seen = list(dict.fromkeys(r["target_gene"] for r in sorted_d))
drugs_seen = list(dict.fromkeys(r["drug_name"] for r in sorted_d))

# Assign y positions
target_y = {t: 12.5 - i*1.8 for i,t in enumerate(targets_seen)}
drug_y = {d: 12.5 - i*1.8 for i,d in enumerate(drugs_seen)}

# Draw nodes
for t, y in target_y.items():
    box = FancyBboxPatch((0.5, y-0.35), 3.5, 0.7, boxstyle="round,pad=0.05,rounding_size=0.10",
                          facecolor="#E74C3C", edgecolor="white", lw=1.5, alpha=0.25)
    ax_b.add_patch(box)
    ax_b.text(2.25, y, t, ha="center", va="center", fontsize=7.5, fontweight="bold", color="#C0392B")
for d, y in drug_y.items():
    box = FancyBboxPatch((8.0, y-0.35), 3.5, 0.7, boxstyle="round,pad=0.05,rounding_size=0.10",
                          facecolor="#3498DB", edgecolor="white", lw=1.5, alpha=0.25)
    ax_b.add_patch(box)
    ax_b.text(9.75, y, d, ha="center", va="center", fontsize=7, fontweight="bold", color="#2980B9")

# Draw edges
for r in sorted_d:
    t = r["target_gene"]; d = r["drug_name"]
    if t in target_y and d in drug_y:
        yt, yd = target_y[t], drug_y[d]
        ax_b.plot([4.0, 8.0], [yt, yd], color="#AAAAAA", lw=1.5, alpha=0.7)

# Column labels
ax_b.text(2.25, 13.8, "TARGETS", ha="center", fontsize=9, fontweight="bold", color="#E74C3C")
ax_b.text(9.75, 13.8, "DRUGS", ha="center", fontsize=9, fontweight="bold", color="#3498DB")

ax_b.set_title("Candidate drug-target network", fontsize=12, fontweight="bold")
plt.tight_layout()
img_b = fig_to_img(fig_b); plt.close(fig_b)
print("  Panel B: bipartite network rendered")

# ═══════════════════════════════════════════════════════════════
# PANEL C: Drugs per target barplot (FIXED)
# ═══════════════════════════════════════════════════════════════
print("\n=== Panel C: Drugs per target ===")

fig_c, ax_c = plt.subplots(1,1,figsize=(7,6),facecolor="white")
sorted_t = sorted(tgt_counts.items(), key=lambda x: -x[1])[:12]
targets_list = [t for t,c in sorted_t]
counts_list = [c for t,c in sorted_t]
colors_c = ["#E74C3C" if c>=3 else "#3498DB" if c>=2 else "#95A5A6" for c in counts_list]

bars = ax_c.barh(range(len(targets_list)), counts_list, color=colors_c, edgecolor="white", height=0.7)
ax_c.set_yticks(range(len(targets_list))); ax_c.set_yticklabels(targets_list, fontsize=8)
ax_c.set_xlabel("Candidate drugs", fontsize=10)
ax_c.set_title("Drugs per target", fontsize=11, fontweight="bold")
ax_c.invert_yaxis()
# Add count labels
for i, (t, c) in enumerate(zip(targets_list, counts_list)):
    ax_c.text(c+0.05, i, str(c), va="center", fontsize=9, fontweight="bold", color="#333333")
ax_c.set_xlim(0, max(counts_list)+1.5)
plt.tight_layout()
img_c = fig_to_img(fig_c); plt.close(fig_c)
print("  Panel C: per-target barplot rendered")

# ═══════════════════════════════════════════════════════════════
# PANEL D: Docking affinity (use existing source PDF)
# ═══════════════════════════════════════════════════════════════
print("\n=== Panel D: Docking affinity ===")
img_d = pdf_img(SRC/"docking/FigS_Docking_affinity_barplot_all_pairs.pdf")

# ═══════════════════════════════════════════════════════════════
# PANEL E: Docking neighbor residues
# ═══════════════════════════════════════════════════════════════
print("\n=== Panel E: Neighbor residues ===")

fig_e, ax_e = plt.subplots(1,1,figsize=(7,5.5),facecolor="white")
sorted_p = sorted(pair_counts.items(), key=lambda x: -x[1])
pair_names = [p.replace("_","-") for p,c in sorted_p]
res_counts = [c for p,c in sorted_p]
colors_e = ["#9B59B6","#8E44AD","#7D3C98","#6C3483"]
ax_e.barh(range(len(pair_names)), res_counts, color=colors_e, edgecolor="white", height=0.65)
ax_e.set_yticks(range(len(pair_names))); ax_e.set_yticklabels(pair_names, fontsize=8)
ax_e.set_xlabel("Residues within 4 A", fontsize=10)
ax_e.set_title("Docking neighbor residues", fontsize=11, fontweight="bold")
ax_e.invert_yaxis()
for i, c in enumerate(res_counts):
    ax_e.text(c+0.3, i, str(c), va="center", fontsize=9, fontweight="bold", color="#333333")
ax_e.set_xlim(0, max(res_counts)+3)
plt.tight_layout()
img_e = fig_to_img(fig_e); plt.close(fig_e)
print("  Panel E: residue barplot rendered")

# ═══════════════════════════════════════════════════════════════
# PANEL F: LARGE mechanism schematic
# ═══════════════════════════════════════════════════════════════
print("\n=== Panel F: Large mechanism schematic ===")

fig_f, ax_f = plt.subplots(1,1,figsize=(20,9),facecolor="white")
ax_f.set_xlim(0,20); ax_f.set_ylim(0,11); ax_f.axis("off")

# Three large cell boxes
cells = [
    ("CAF\nFibroblast", 3.5, 7.0, "#E74C3C"),
    ("Macrophage\nTAM", 10.0, 7.0, "#3498DB"),
    ("Epithelial\nlike", 16.5, 7.0, "#27AE60"),
]
for label, x, y, color in cells:
    box = FancyBboxPatch((x-2.0, y-1.0), 4.0, 2.0, boxstyle="round,pad=0.1,rounding_size=0.2",
                          facecolor=color, edgecolor="white", lw=2.5, alpha=0.22)
    ax_f.add_patch(box)
    ax_f.text(x, y, label, ha="center", va="center", fontsize=12, fontweight="bold", color=color)

# Arrows with labels
arrows_data = [
    (5.5,7.2,8.0,7.2, "COL1A1/COL1A2 -> CD44"),
    (5.5,6.5,8.0,6.5, "MIF -> CD74/CXCR4"),
    (8.0,5.8,5.5,5.8, "TGFB1 -> TGFBR"),
    (8.0,5.3,5.5,5.3, "PPIA -> BSG"),
    (16.5,8.0,12.0,8.0, "MIF -> CD74/CXCR4"),
]
for x1,y1,x2,y2,label in arrows_data:
    ax_f.annotate("",xy=(x2,y2),xytext=(x1,y1),
                  arrowprops=dict(arrowstyle="->",color="#555555",lw=2.5))
    mx,my = (x1+x2)/2, (y1+y2)/2 + 0.2
    ax_f.text(mx,my,label,ha="center",va="bottom",fontsize=8,color="#333333",
              bbox=dict(boxstyle="round,pad=0.06",facecolor="white",edgecolor="none",alpha=0.85))

# TAM checkpoint box
box_tam = FancyBboxPatch((7.0,3.5),6.0,1.3,boxstyle="round,pad=0.1,rounding_size=0.12",
                           facecolor="#3498DB",edgecolor="white",lw=1.5,alpha=0.15)
ax_f.add_patch(box_tam)
ax_f.text(10.0,4.3,"TAM checkpoint transcripts (mRNA-level)",ha="center",fontsize=9.5,fontweight="bold",color="#2980B9")
ax_f.text(10.0,3.8,"HAVCR2  |  IDO1  |  CD86  |  PDCD1LG2",ha="center",fontsize=8.5,color="#333333")

# Candidate targets box
box_tgt = FancyBboxPatch((7.0,2.0),6.0,1.0,boxstyle="round,pad=0.08,rounding_size=0.10",
                           facecolor="#27AE60",edgecolor="white",lw=1.5,alpha=0.15)
ax_f.add_patch(box_tgt)
ax_f.text(10.0,2.6,"Candidate targets: CXCR4  |  IDO1  |  TGFBR1  |  MIF",ha="center",fontsize=9,fontweight="bold",color="#1E8449")
ax_f.text(10.0,2.25,"Docking = in silico predictions only",ha="center",fontsize=7,fontstyle="italic",color="#888888")

# Title and caveat
ax_f.text(10.0,10.3,"Proposed CAF-TAM-Epithelial communication model",ha="center",fontsize=14,fontweight="bold")
ax_f.text(10.0,0.6,"Computational model; experimental validation required.",ha="center",fontsize=8,fontstyle="italic",color="#888888")

# Legend
leg_elems = [
    mpatches.Patch(color="#E74C3C",alpha=0.22,label="CAF Fibroblast"),
    mpatches.Patch(color="#3498DB",alpha=0.22,label="Macrophage TAM"),
    mpatches.Patch(color="#27AE60",alpha=0.22,label="Epithelial-like"),
]
ax_f.legend(handles=leg_elems,loc="upper right",fontsize=7.5,framealpha=0.9)

plt.tight_layout()
img_f = fig_to_img(fig_f); plt.close(fig_f)
print(f"  Panel F: large mechanism ({img_f.shape[1]}x{img_f.shape[0]})")

# ═══════════════════════════════════════════════════════════════
# ASSEMBLE FIGURE 7 (large canvas: 20 x 18 in)
# ═══════════════════════════════════════════════════════════════
print("\n=== Assembling Figure 7 ===")

fig7 = plt.figure(figsize=(20, 19), facecolor="white")

# Grid: row0=A+B, row1=C+D+E, row2=F spanning full width
gs = gridspec.GridSpec(3, 6, figure=fig7, wspace=0.08, hspace=0.12,
                       left=0.025, right=0.975, top=0.985, bottom=0.015,
                       height_ratios=[1.0, 1.0, 1.6])

# Row 0: A (cols 0-2), B (cols 3-5)
ax_a = fig7.add_subplot(gs[0, 0:3])
img_a = pdf_img(SRC/"drug_screening/Fig11A_target_druggability_heatmap.pdf")
ax_a.imshow(img_a); ax_a.axis("off"); add_label(ax_a, "A")

ax_b = fig7.add_subplot(gs[0, 3:6])
ax_b.imshow(img_b); ax_b.axis("off"); add_label(ax_b, "B")

# Row 1: C (cols 0-1), D (cols 2-3), E (cols 4-5)
ax_c = fig7.add_subplot(gs[1, 0:2])
ax_c.imshow(img_c); ax_c.axis("off"); add_label(ax_c, "C")

ax_d = fig7.add_subplot(gs[1, 2:4])
ax_d.imshow(img_d); ax_d.axis("off"); add_label(ax_d, "D")
add_note(ax_d, "In silico predictions only")

ax_e = fig7.add_subplot(gs[1, 4:6])
ax_e.imshow(img_e); ax_e.axis("off"); add_label(ax_e, "E")

# Row 2: F spans all 6 cols
ax_f = fig7.add_subplot(gs[2, :])
ax_f.imshow(img_f); ax_f.axis("off"); add_label(ax_f, "F", fs=LF)
add_note(ax_f, "Computational model; experimental validation required")

fig7.suptitle("Figure 7. Candidate target prioritization and in silico docking",
              fontsize=14, fontweight="bold", y=0.998)

# Save
for ext in ["pdf", "png"]:
    p = OUT / f"Figure7_therapeutic_implications_v6.{ext}"
    fig7.savefig(p, dpi=DPI, facecolor="white", edgecolor="none", bbox_inches="tight")
    print(f"  {p.name} ({p.stat().st_size//1024} KB)")
plt.close(fig7)

print("\nFigure 7 v6 complete.")
