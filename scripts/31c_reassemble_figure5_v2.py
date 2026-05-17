#!/usr/bin/env python
"""Re-assemble Figure 5 with redrawn Panel A network + existing B/C/D."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import csv
from pathlib import Path
import fitz
from PIL import Image
import io

BASE = Path(r"e:\CCA")
OUT_MAIN = BASE / "figures/main_repaired"
OUT_FINAL = BASE / "figures/main_final"
for d in [OUT_MAIN, OUT_FINAL]:
    d.mkdir(parents=True, exist_ok=True)

DPI = 300
FIGSIZE = 3.5

def pdf_page_to_image(pdf_path, page_idx=0, dpi=DPI):
    doc = fitz.open(str(pdf_path))
    page = doc[page_idx]
    zoom = dpi / 72.0
    mat = fitz.Matrix(zoom, zoom)
    pix = page.get_pixmap(matrix=mat, colorspace=fitz.csRGB)
    img_data = pix.tobytes("png")
    doc.close()
    return np.array(Image.open(io.BytesIO(img_data)))

def add_panel_label(ax, letter, x=0.02, y=0.96, fontsize=14):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=fontsize,
            fontweight="bold", color="black", va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.1", facecolor="white",
                      edgecolor="none", alpha=0.85))

def add_caveat(ax, text, fontsize=7):
    ax.text(0.98, 0.02, text, transform=ax.transAxes, fontsize=fontsize,
            fontstyle="italic", color="#888888", va="bottom", ha="right")

def make_bubble_plot(rows, title, figsize=(7, 5)):
    """Generate a CellChat-style bubble plot (same as Step 31)."""
    rows_sorted = sorted(rows, key=lambda r: float(r["prob"]), reverse=True)[:25]
    rows_sorted = sorted(rows_sorted, key=lambda r: float(r["prob"]))
    pair_labels = [f"{r['ligand']} – {r['receptor']}" for r in rows_sorted]
    probs = [float(r["prob"]) for r in rows_sorted]
    fig, ax = plt.subplots(1, 1, figsize=figsize, facecolor="white")
    sizes = [p * 800 for p in probs]
    scatter = ax.scatter(probs, pair_labels, s=sizes, c=probs, cmap="Reds",
                         edgecolors="grey", linewidth=0.3, alpha=0.85, zorder=3)
    ax.set_xlabel("Communication probability", fontsize=9)
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
    fig.canvas.draw()
    img = np.frombuffer(fig.canvas.buffer_rgba(), dtype=np.uint8)
    w, h = fig.canvas.get_width_height()
    return img.reshape(h, w, 4)[:, :, :3]

# ── Load CellChat data for bubble plots ──
lr_csv = BASE / "tables/cellchat/GSE138709_CAF_TAM_Epithelial_LR_pairs.csv"
with open(lr_csv, newline="") as f:
    all_lr = list(csv.DictReader(f))

caf_to_tam = [r for r in all_lr
              if r["source"] == "CAF_Fibroblast"
              and r["target"] == "Macrophage_TAM"]
tam_to_caf = [r for r in all_lr
              if r["source"] == "Macrophage_TAM"
              and r["target"] == "CAF_Fibroblast"]

# ── Generate panels ──
print("Generating Figure 5 panels...")

# Panel A: new redrawn network
img_a = pdf_page_to_image(
    BASE / "figures/cellchat/Fig8C_CellChat_CAF_TAM_Epithelial_network_redrawn.pdf")
print(f"  Panel A: redrawn network ({img_a.shape[1]}x{img_a.shape[0]})")

# Panel B: existing outgoing/incoming heatmap
img_b = pdf_page_to_image(
    BASE / "figures/cellchat/Fig8D_CellChat_outgoing_incoming_heatmap.pdf")
print(f"  Panel B: heatmap ({img_b.shape[1]}x{img_b.shape[0]})")

# Panel C: CAF→TAM bubble
print(f"  Panel C: CAF→TAM bubble ({len(caf_to_tam)} pairs)")
fig_c = make_bubble_plot(caf_to_tam, "CAF → TAM (inferred communication)")
img_c = fig_to_image(fig_c)
plt.close(fig_c)

# Panel D: TAM→CAF bubble
print(f"  Panel D: TAM→CAF bubble ({len(tam_to_caf)} pairs)")
fig_d = make_bubble_plot(tam_to_caf, "TAM → CAF (inferred communication)")
img_d = fig_to_image(fig_d)
plt.close(fig_d)

# ── Assemble composite ──
print("Assembling composite Figure 5...")
fig, axes = plt.subplots(2, 2, figsize=(8.7, 7.5), facecolor="white")
gs = axes[0, 0].get_gridspec()
for ax_row in axes:
    for ax in ax_row:
        ax.remove()

# Recreate with GridSpec for better control
fig = plt.figure(figsize=(8.7, 7.5), facecolor="white")
gs = gridspec.GridSpec(2, 2, figure=fig, wspace=0.10, hspace=0.15,
                       left=0.03, right=0.97, top=0.94, bottom=0.03)

panels = [
    ("A", img_a, (0, 0)),
    ("B", img_b, (0, 1)),
    ("C", img_c, (1, 0)),
    ("D", img_d, (1, 1)),
]

for letter, img, (r, c) in panels:
    ax = fig.add_subplot(gs[r, c])
    ax.imshow(img)
    ax.axis("off")
    add_panel_label(ax, letter)
    if letter == "C":
        add_caveat(ax, "Inferred communication", fontsize=6.5)
    if letter == "D":
        add_caveat(ax, "Inferred communication", fontsize=6.5)

fig.suptitle("Figure 5. CellChat-inferred CAF–TAM–Epithelial communication network",
             fontsize=11, fontweight="bold", y=0.99)

# Save
stem_main = OUT_MAIN / "Fig5_CellChat_communication_repaired_v2"
stem_final = OUT_FINAL / "Figure5_CellChat_communication"
for stem, label in [(stem_main, "main_repaired"), (stem_final, "main_final")]:
    fig.savefig(f"{stem}.pdf", dpi=DPI, facecolor="white", edgecolor="none",
                bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=DPI, facecolor="white", edgecolor="none",
                bbox_inches="tight")
    print(f"  {label}/: {stem.name}.pdf + .png")

plt.close(fig)
print("Done — Figure 5 re-assembled with redrawn Panel A.")
