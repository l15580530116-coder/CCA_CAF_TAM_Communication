#!/usr/bin/env python
"""Generate Figure 1: Study Workflow Schematic — 6 horizontal modules."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from pathlib import Path

BASE = Path(r"e:\CCA")
for d in ["figures/main", "figures/editable"]:
    (BASE / d).mkdir(parents=True, exist_ok=True)

# ── Colour palette ──
C = {
    "data":       "#4472C4",   # blue
    "deg":        "#2F9B8E",   # teal
    "score":      "#ED7D31",   # orange
    "single":     "#9B59B6",   # purple
    "cellchat":   "#C0392B",   # red
    "therapeutic":"#27AE60",   # green
    "arrow":      "#555555",
    "bg":         "#FAFAFA",
    "text":       "#222222",
}

plt.rcParams.update({"font.family": "sans-serif", "font.size": 9})

fig, ax = plt.subplots(1, 1, figsize=(22, 8))
ax.set_xlim(0, 22); ax.set_ylim(0, 8)
ax.axis("off")
ax.set_facecolor(C["bg"])

# ── Module positions (x_center, y_center, width, height) ──
modules = [
    ("A", "Data\nCollection",        1.8, 4.0, 3.0, 6.0, C["data"]),
    ("B", "Cross-Cohort\nDEG Discovery",  5.2, 4.0, 3.0, 6.0, C["deg"]),
    ("C", "Microenvironment\nScore & Immune", 8.6, 4.0, 3.0, 6.0, C["score"]),
    ("D", "Single-Cell\nLocalization", 12.0, 4.0, 3.0, 6.0, C["single"]),
    ("E", "CellChat-Inferred\nCommunication", 15.4, 4.0, 3.0, 6.0, C["cellchat"]),
    ("F", "Therapeutic\nImplications", 18.8, 4.0, 3.0, 6.0, C["therapeutic"]),
]

def draw_module(ax, letter, title, cx, cy, w, h, color):
    """Draw a rounded-box module with title and content."""
    x0, y0 = cx - w/2, cy - h/2

    # Module box
    box = FancyBboxPatch((x0, y0), w, h,
                         boxstyle="round,pad=0.15,rounding_size=0.25",
                         facecolor=color, edgecolor="white", linewidth=2, alpha=0.18)
    ax.add_patch(box)

    # Letter circle
    circle = plt.Circle((cx, cy + h/2 - 0.55), 0.35, color=color, ec="white", linewidth=1.5)
    ax.add_patch(circle)
    ax.text(cx, cy + h/2 - 0.55, letter, ha="center", va="center", fontsize=12,
            fontweight="bold", color="white")

    # Title
    ax.text(cx, cy + h/2 - 1.2, title, ha="center", va="top", fontsize=10,
            fontweight="bold", color=color, linespacing=1.3)

    return x0, y0, cx, cy, w, h

# ── Draw all modules ──
for letter, title, cx, cy, w, h, color in modules:
    draw_module(ax, letter, title, cx, cy, w, h, color)

# ── Module content ──
def module_text(ax, cx, cy, y_start, lines, color, fontsize=8.5):
    """Add text lines inside a module."""
    for i, line in enumerate(lines):
        fc = color if "→" not in line and "•" not in line else C["text"]
        fw = "bold" if i == 0 else "normal"
        fs = fontsize if i > 0 else fontsize + 0.5
        ax.text(cx, y_start - i * 0.34, line, ha="center", va="top", fontsize=fs,
                fontweight=fw, color=fc, linespacing=1.15)

# A: Data Collection
module_text(ax, 1.8, 4.0, 6.2, [
    "4 Public Cohorts",
    "TCGA-CHOL (n=44, RNA-seq)",
    "GSE107943 (n=57, RNA-seq)",
    "GSE138709 (32,626 cells, 10X)",
    "GSE26566 (n=169, Microarray)",
    "• 35 + 9 tumor/normal",
    "• 30 + 27 tumor/normal",
    "• 5 iCCA patients",
    "• Expression validation",
], C["data"])

# B: DEG Discovery
module_text(ax, 5.2, 4.0, 6.2, [
    "DEG Cross-Validation",
    "DESeq2 (TCGA)",
    "edgeR paired (GSE107943)",
    "4,534 same-direction DEGs",
    "344 validated IM/CAF/TAM",
    "78 upregulated",
    "266 downregulated",
    "99.4% directional concordance",
], C["deg"])

# C: Score + Immune
module_text(ax, 8.6, 4.0, 6.2, [
    "ssGSEA Scoring + Immune",
    "Aggressive score:",
    "ECM + CAF + TAM − Metab",
    "Score vs CAF: rho = 0.77",
    "Score vs M2 TAM: rho = 0.65",
    "Checkpoint correlations:",
    "HAVCR2 (0.54)  IDO1 (0.48)",
    "CD86 (0.55)  PDCD1LG2 (0.53)",
], C["score"])

# D: Single-Cell
module_text(ax, 12.0, 4.0, 6.2, [
    "scRNA-seq (GSE138709)",
    "8 cell types identified",
    "CAF_Fibroblast:",
    "  Aggressive score = 6.63",
    "  (highest among all types)",
    "Macrophage_TAM expresses:",
    "  HAVCR2, IDO1, CD86",
    "  PDCD1LG2 (transcript level)",
], C["single"])

# E: CellChat
module_text(ax, 15.4, 4.0, 6.2, [
    "CellChat-Inferred Axes",
    "15,824 tumor cells, 6 types",
    "CAF → TAM:",
    "  COL1A1/COL1A2 → CD44",
    "  MIF → CD74/CXCR4",
    "TAM → CAF:",
    "  PPIA → BSG",
    "  TGFB1 → TGFBR",
], C["cellchat"])

# F: Therapeutic
module_text(ax, 18.8, 4.0, 6.2, [
    "Candidate Targets",
    "Hub genes: COL1A1, COL1A2,",
    "  TREM2, POSTN, SPP1",
    "Druggability: 14 highly",
    "  druggable targets",
    "Docking (in silico):",
    "  CXCR4–Plerixafor −8.01",
    "  IDO1–Epacadostat −6.47",
    "  (kcal/mol)",
], C["therapeutic"])

# ── Arrows between modules ──
for i in range(len(modules) - 1):
    x1 = modules[i][2] + modules[i][4]/2 + 0.05
    x2 = modules[i+1][2] - modules[i+1][4]/2 - 0.05
    y_mid = 4.0
    ax.annotate("", xy=(x2, y_mid), xytext=(x1, y_mid),
                arrowprops=dict(arrowstyle="->", color=C["arrow"],
                               lw=2.5, connectionstyle="arc3,rad=0"))

# ── Cohort icons (small boxes above Module A) ──
cohorts = [
    ("TCGA-CHOL", 0.7, 7.0),
    ("GSE107943", 1.8, 7.0),
    ("GSE138709", 2.1, 7.0),
    ("GSE26566", 1.4, 7.0),
]
# Place simplified dataset badges above module A
ax.text(1.8, 7.5, "Input Datasets", ha="center", fontsize=9, fontweight="bold", color=C["data"])

# ── Title ──
ax.text(11, 7.75, "Figure 1. Study Workflow and Analytical Framework",
        ha="center", fontsize=14, fontweight="bold", color=C["text"])

# ── Bottom annotation ──
ax.text(11, 0.35,
        "In silico predictions. All survival analyses are exploratory. Docking scores are computational predictions.",
        ha="center", fontsize=7.5, fontstyle="italic", color="#888888")

# ── Legend bar ──
legend_y = 0.65
categories = [
    ("Data", C["data"]), ("DEG", C["deg"]), ("Score/Immune", C["score"]),
    ("Single-Cell", C["single"]), ("CellChat", C["cellchat"]), ("Therapeutic", C["therapeutic"]),
]
x_start = 6.0
for j, (label, col) in enumerate(categories):
    lx = x_start + j * 2.1
    ax.add_patch(plt.Rectangle((lx-0.12, legend_y-0.06), 0.24, 0.12, color=col, alpha=0.6, ec=None))
    ax.text(lx + 0.18, legend_y, label, fontsize=7, color=C["text"], va="center")

# ── Save ──
for fmt, dpi in [("pdf", 300), ("png", 200), ("svg", None)]:
    path = BASE / f"figures/main/Fig1_workflow_schematic.{fmt}" if fmt != "svg" else BASE / "figures/editable/Fig1_workflow_schematic.svg"
    plt.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white", edgecolor="none")
    print(f"  Saved: {path}")

# Also save high-res PNG to editable
plt.savefig(BASE / "figures/editable/Fig1_workflow_schematic.png", dpi=300, bbox_inches="tight",
            facecolor="white", edgecolor="none")

plt.close()
print("Figure 1 generated successfully")
