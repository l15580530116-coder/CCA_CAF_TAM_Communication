#!/usr/bin/env python
"""Redraw Figure 5A: CellChat-inferred 3-node communication network.
Uses pure matplotlib — no networkx needed.
Data from: tables/cellchat/GSE138709_cellchat_overall_interactions.csv
           tables/cellchat/GSE138709_cellchat_celltype_counts.csv
           tables/cellchat/GSE138709_CAF_TAM_Epithelial_LR_pairs.csv
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle, FancyBboxPatch
import matplotlib.patches as mpatches
import numpy as np
import csv
from pathlib import Path

BASE = Path(r"e:\CCA")
OUT_NET = BASE / "figures/cellchat"
OUT_MAIN = BASE / "figures/main_repaired"
OUT_FINAL = BASE / "figures/main_final"
for d in [OUT_NET, OUT_MAIN, OUT_FINAL]:
    d.mkdir(parents=True, exist_ok=True)

DPI = 300

# ═══════════════════════════════════════════════════════════════
# 1. Load data
# ═══════════════════════════════════════════════════════════════

# Cell type counts (tumor-only where possible)
cell_counts = {}
with open(BASE / "tables/cellchat/GSE138709_cellchat_celltype_counts.csv", newline="") as f:
    for row in csv.DictReader(f):
        cell_counts[row["cell_type"]] = int(row["n_cells"])

print("Cell counts (all samples):")
for ct in ["CAF_Fibroblast", "Macrophage_TAM", "Epithelial_like"]:
    print(f"  {ct}: {cell_counts.get(ct, '?')}")

# Overall interactions (edge weights)
interactions = {}
with open(BASE / "tables/cellchat/GSE138709_cellchat_overall_interactions.csv", newline="") as f:
    for row in csv.DictReader(f):
        src = row["source"]
        tgt = row["target"]
        interactions[(src, tgt)] = {
            "count": int(row["interaction_count"]),
            "weight": float(row["interaction_weight"]),
        }

# Focus: 3 nodes
NODES = ["CAF_Fibroblast", "Macrophage_TAM", "Epithelial_like"]

# Extract relevant edges
edges = {}
for src in NODES:
    for tgt in NODES:
        if src == tgt:
            continue  # skip self-loops for clarity
        key = (src, tgt)
        if key in interactions:
            edges[key] = interactions[key]

print("\nEdges (source → target):")
for (s, t), d in sorted(edges.items()):
    print(f"  {s} → {t}: count={d['count']}, weight={d['weight']:.3f}")

# Top LR pairs per edge (for annotation)
lr_csv = BASE / "tables/cellchat/GSE138709_CAF_TAM_Epithelial_LR_pairs.csv"
edge_top_lr = {}
with open(lr_csv, newline="") as f:
    all_lr = list(csv.DictReader(f))

for (s, t) in edges:
    pairs = [(r, float(r["prob"])) for r in all_lr
             if r["source"] == s and r["target"] == t]
    pairs.sort(key=lambda x: -x[1])
    # Get top 2 unique LR combos
    seen = set()
    top = []
    for r, prob in pairs:
        combo = f"{r['ligand']}–{r['receptor']}"
        if combo not in seen:
            seen.add(combo)
            top.append((combo, prob))
        if len(top) >= 2:
            break
    edge_top_lr[(s, t)] = top

print("\nTop LR pairs per edge:")
for (s, t), tops in edge_top_lr.items():
    names = ", ".join([f"{c} ({p:.3f})" for c, p in tops])
    print(f"  {s} → {t}: {names}")

# ═══════════════════════════════════════════════════════════════
# 2. Draw network
# ═══════════════════════════════════════════════════════════════

# Node positions (triangle layout)
# Epithelial_like at top-center (largest node, tumor cells)
# CAF_Fibroblast at bottom-left
# Macrophage_TAM at bottom-right
pos = {
    "CAF_Fibroblast":    np.array([1.2, 0.3]),
    "Macrophage_TAM":    np.array([3.8, 0.3]),
    "Epithelial_like":   np.array([2.5, 2.0]),
}

# Node metadata
node_meta = {
    "CAF_Fibroblast": {
        "color": "#E74C3C",  # red
        "label": "CAF\nFibroblast",
        "count": cell_counts.get("CAF_Fibroblast", "?"),
    },
    "Macrophage_TAM": {
        "color": "#3498DB",  # blue
        "label": "Macrophage\nTAM",
        "count": cell_counts.get("Macrophage_TAM", "?"),
    },
    "Epithelial_like": {
        "color": "#27AE60",  # green
        "label": "Epithelial\nlike",
        "count": cell_counts.get("Epithelial_like", "?"),
    },
}

# Node sizes based on relative cell count
all_counts = [node_meta[n]["count"] for n in NODES]
min_count, max_count = min(all_counts), max(all_counts)
node_radii = {}
for n in NODES:
    frac = (node_meta[n]["count"] - min_count) / (max_count - min_count + 1)
    node_radii[n] = 0.35 + frac * 0.50  # radius range 0.35–0.85

# Edge thickness based on weight (number of LR pairs)
edge_weights_list = [edges[k]["count"] for k in edges]
min_w, max_w = min(edge_weights_list), max(edge_weights_list)

def edge_lw(count):
    frac = (count - min_w) / max(max_w - min_w, 1)
    return 2.0 + frac * 6.0  # linewidth range 2.0–8.0

# Color map for edges
edge_colors = {
    ("CAF_Fibroblast", "Macrophage_TAM"):   "#C0392B",  # dark red
    ("Macrophage_TAM", "CAF_Fibroblast"):   "#2980B9",  # dark blue
    ("CAF_Fibroblast", "Epithelial_like"):  "#E67E22",  # orange
    ("Epithelial_like", "CAF_Fibroblast"):  "#8E44AD",  # purple
    ("Epithelial_like", "Macrophage_TAM"):  "#16A085",  # teal
    ("Macrophage_TAM", "Epithelial_like"):  "#2C3E50",  # dark grey-blue
}

# ── Create figure ──
fig, ax = plt.subplots(1, 1, figsize=(7.5, 5.5), facecolor="white")
ax.set_xlim(0, 5.0)
ax.set_ylim(-0.2, 2.8)
ax.axis("off")
ax.set_aspect("equal")

# Draw edges
for (src, tgt), d in sorted(edges.items()):
    p1 = pos[src]
    p2 = pos[tgt]
    color = edge_colors.get((src, tgt), "#888888")
    lw = edge_lw(d["count"])
    alpha = 0.75

    # Offset arrow start/end by node radius
    vec = p2 - p1
    dist = np.linalg.norm(vec)
    if dist < 0.01:
        continue
    uvec = vec / dist
    r1 = node_radii[src]
    r2 = node_radii[tgt]
    start = p1 + uvec * (r1 + 0.05)
    end = p2 - uvec * (r2 + 0.10)

    # Draw arrow
    arrow = FancyArrowPatch(
        start, end,
        arrowstyle="->", mutation_scale=18,
        color=color, lw=lw, alpha=alpha,
        connectionstyle="arc3,rad=0.05",
        zorder=1,
    )
    ax.add_patch(arrow)

    # Edge label (LR count + top pair) at midpoint
    mid = (start + end) / 2
    # Perpendicular offset to avoid arrow
    perp = np.array([-uvec[1], uvec[0]]) * 0.18
    label_pos = mid + perp

    top_pairs = edge_top_lr.get((src, tgt), [])
    top_text = ", ".join([c for c, p in top_pairs[:1]]) if top_pairs else ""
    if len(top_text) > 22:
        top_text = top_text[:20] + "…"

    edge_label = f"{d['count']} LR pairs"
    if top_text:
        edge_label += f"\n{top_text}"
    ax.annotate(
        edge_label, xy=label_pos, fontsize=6.0, color="#333333",
        ha="center", va="center", zorder=5,
        bbox=dict(boxstyle="round,pad=0.12", facecolor="white",
                  edgecolor="none", alpha=0.80),
    )

# Draw nodes
for n in NODES:
    x, y = pos[n]
    r = node_radii[n]
    meta = node_meta[n]

    # Node circle
    circle = Circle((x, y), r, facecolor=meta["color"], edgecolor="white",
                    linewidth=1.5, alpha=0.82, zorder=10)
    ax.add_patch(circle)

    # Node label
    ax.text(x, y, meta["label"], ha="center", va="center",
            fontsize=8, fontweight="bold", color="white", zorder=11)

    # Cell count below node
    ax.text(x, y - r - 0.10, f"n = {meta['count']}",
            ha="center", va="top", fontsize=6.5, color="#555555")

# Title and caveat
ax.set_title("CellChat-inferred communication network",
             fontsize=11, fontweight="bold", pad=8)

ax.text(2.5, -0.35,
        "Inferred communication; not experimentally validated",
        ha="center", fontsize=7, fontstyle="italic", color="#888888")

# Legend for node types
legend_elements = [
    mpatches.Patch(color=node_meta["CAF_Fibroblast"]["color"],
                   alpha=0.82, label="CAF Fibroblast"),
    mpatches.Patch(color=node_meta["Macrophage_TAM"]["color"],
                   alpha=0.82, label="Macrophage TAM"),
    mpatches.Patch(color=node_meta["Epithelial_like"]["color"],
                   alpha=0.82, label="Epithelial-like"),
]
leg = ax.legend(handles=legend_elements, loc="lower center",
                ncol=3, fontsize=7, frameon=True, framealpha=0.9,
                bbox_to_anchor=(0.5, -0.32))
leg.set_zorder(20)

plt.tight_layout()

# Save network figure
net_pdf = OUT_NET / "Fig8C_CellChat_CAF_TAM_Epithelial_network_redrawn.pdf"
net_png = OUT_NET / "Fig8C_CellChat_CAF_TAM_Epithelial_network_redrawn.png"
fig.savefig(net_pdf, dpi=DPI, facecolor="white", edgecolor="none",
            bbox_inches="tight")
fig.savefig(net_png, dpi=DPI, facecolor="white", edgecolor="none",
            bbox_inches="tight")
plt.close(fig)

ks_pdf = net_pdf.stat().st_size // 1024
ks_png = net_png.stat().st_size // 1024
print(f"\nFigure 5A network saved:")
print(f"  {net_pdf} ({ks_pdf} KB)")
print(f"  {net_png} ({ks_png} KB)")
print("Done — Figure 5A redrawn.")
