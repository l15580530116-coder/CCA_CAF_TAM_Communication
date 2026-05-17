#!/usr/bin/env python
"""Step 38d: Expand Figures 5/6/7 to 6-panel, fix Figure 2E, generate main figures v5."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path
import fitz, csv, io, shutil
from PIL import Image

BASE = Path(r"e:\CCA")
SRC = BASE / "figures"
TBL = BASE / "tables"
OUT = BASE / "figures/main_final_v5"
SRC_V2 = BASE / "figures/main_final_v2"
OUT.mkdir(parents=True, exist_ok=True)
DPI = 600
LF = 17
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
    ax.text(0.015,0.975,letter,transform=ax.transAxes,fontsize=fs,fontweight="bold",
            color="black",va="top",ha="left",
            bbox=dict(boxstyle="round,pad=0.10",facecolor="white",edgecolor="none",alpha=0.88))

def add_note(ax, text, fs=NF):
    ax.text(0.982,0.018,text,transform=ax.transAxes,fontsize=fs,fontstyle="italic",
            color="#555555",va="bottom",ha="right",
            bbox=dict(boxstyle="round,pad=0.10",facecolor="white",edgecolor="none",alpha=0.88))

def save_fig(fig, stem):
    for ext in ["pdf","png"]:
        p=OUT/f"{stem}.{ext}"
        fig.savefig(p,dpi=DPI,facecolor="white",edgecolor="none",bbox_inches="tight")
        print(f"  {stem}.{ext} ({p.stat().st_size//1024} KB)")

def fig_to_img(fig):
    fig.canvas.draw()
    img=np.frombuffer(fig.canvas.buffer_rgba(),dtype=np.uint8)
    w,h=fig.canvas.get_width_height()
    return img.reshape(h,w,4)[:,:,:3]

def assemble(panels, nrows, ncols, title, stem, notes=None, fs=4.0):
    """panels=[(letter, path_or_None, pg, desc, img_override_or_None)]"""
    fig=plt.figure(figsize=(ncols*fs,nrows*fs),facecolor="white")
    gs=gridspec.GridSpec(nrows,ncols,figure=fig,wspace=0.10,hspace=0.20,
                         left=0.03,right=0.97,top=0.935,bottom=0.03)
    for idx,item in enumerate(panels):
        letter,path,pg,desc=item[:4]; override=item[4] if len(item)>4 else None
        row,col=divmod(idx,ncols); ax=fig.add_subplot(gs[row,col])
        try:
            img=override if override is not None else pdf_img(path,pg=pg)
            ax.imshow(img); ax.axis("off"); add_label(ax,letter)
            if notes and letter in notes: add_note(ax,notes[letter])
        except Exception as e:
            ax.text(0.5,0.5,f"Panel {letter}\nNOT AVAILABLE",ha="center",va="center",fontsize=7,color="red")
            ax.axis("off")
    for idx in range(len(panels),nrows*ncols):
        ax=fig.add_subplot(gs[idx//ncols,idx%ncols]); ax.axis("off")
    fig.suptitle(title,fontsize=12,fontweight="bold",y=0.99)
    save_fig(fig,stem); plt.close(fig)

# ═══════════════════════════════════════════════════════════════
# COPY Figures 1, 3, 4 from v2
# ═══════════════════════════════════════════════════════════════
print("=== Copying Figures 1, 3, 4 from v2 ===")
names={1:"workflow_schematic",3:"aggressive_score_immune",4:"single_cell_localization"}
for fid in [1,3,4]:
    n=names[fid]
    for ext in ["pdf","png"]:
        src=SRC_V2/f"Figure{fid}_{n}_v2.{ext}"
        dst=OUT/f"Figure{fid}_{n}_v5.{ext}"
        shutil.copy2(src,dst)
        print(f"  Figure {fid}: copied {ext}")

# ═══════════════════════════════════════════════════════════════
# FIGURE 2E: Clear funnel/flow validation summary
# ═══════════════════════════════════════════════════════════════
print("\n=== Figure 2E: Clear funnel validation summary ===")
fig2e,ax2e=plt.subplots(1,1,figsize=(7,6),facecolor="white")
ax2e.set_xlim(0,10); ax2e.set_ylim(0,12); ax2e.axis("off")

# Funnel-style boxes
boxes=[
    (1.5,10.5,7.0,1.2,"TCGA-CHOL DEGs = 9,291","GSE107943 DEGs = 7,735","#4472C4"),
    (2.0,8.5,6.0,1.2,"Common DEGs = 4,560","Intersection of both cohorts","#2F9B8E"),
    (2.3,6.5,5.4,1.2,"Same-direction = 4,534","99.4% directional concordance","#ED7D31"),
    (2.5,4.5,5.0,1.2,"Validated IM/CAF/TAM DEGs = 344","78 up | 266 down","#27AE60"),
]
for x,y,w,h,line1,line2,color in boxes:
    box=FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.1,rounding_size=0.15",
                        facecolor=color,edgecolor="white",lw=2,alpha=0.18)
    ax2e.add_patch(box)
    ax2e.text(x+w/2,y+h/2+0.15,line1,ha="center",va="center",fontsize=12,fontweight="bold",color=color)
    ax2e.text(x+w/2,y+h/2-0.35,line2,ha="center",va="center",fontsize=9,color="#555555")

# Arrows between boxes
for y1,y2 in [(10.5,9.7),(8.5,7.7),(6.5,5.7)]:
    ax2e.annotate("",xy=(5.0,y2),xytext=(5.0,y1),arrowprops=dict(arrowstyle="->",color="#666666",lw=3))

ax2e.text(5.0,3.2,"Validated in both cohorts with concordant direction\nFDR < 0.05, |log2FC| > 1",ha="center",fontsize=9,fontstyle="italic",color="#888888")
ax2e.set_title("Cross-cohort DEG validation summary",fontsize=13,fontweight="bold")
plt.tight_layout()
img2e=fig_to_img(fig2e); plt.close(fig2e)

# Assemble Figure 2
print("=== Figure 2: Re-assemble with clear E ===")
assemble([
    ("A",SRC/"DEG/Fig2A_TCGA_PCA_tumor_normal.pdf",0,"PCA TCGA"),
    ("B",SRC/"DEG/Fig2B_TCGA_volcano.pdf",0,"Volcano TCGA"),
    ("C",SRC/"DEG/Fig2D_GSE107943_PCA_tumor_normal.pdf",0,"PCA GSE"),
    ("D",SRC/"DEG/Fig2E_GSE107943_paired_volcano.pdf",0,"Volcano GSE"),
    ("E",None,0,"Overlap funnel",img2e),
    ("F",SRC/"enrichment/Fig3A_GO_BP_up_dotplot.pdf",0,"GO BP up"),
    ("G",SRC/"enrichment/Fig3B_GO_BP_down_dotplot.pdf",0,"GO BP down"),
    ("H",SRC/"enrichment/Fig3E_up_down_pathway_barplot.pdf",0,"Pathway barplot"),
],4,2,"Figure 2. Cross-cohort DEG validation and functional enrichment",
    "Figure2_DEG_cross_validation_v5")

# ═══════════════════════════════════════════════════════════════
# FIGURE 5: Expand to A-F (6 panels)
# ═══════════════════════════════════════════════════════════════
print("\n=== Figure 5: Expand to 6 panels ===")

# Load CellChat data for E/F
with open(TBL/"cellchat/GSE138709_CAF_TAM_Epithelial_LR_pairs.csv",newline="") as f:
    lr_all=list(csv.DictReader(f))
with open(TBL/"cellchat/GSE138709_cellchat_overall_interactions.csv",newline="") as f:
    int_all=list(csv.DictReader(f))

# E: Top signaling pathways barplot (by LR pair count summed by pathway)
pathway_counts={}
for r in lr_all:
    pw=r.get("pathway_name","")
    pathway_counts[pw]=pathway_counts.get(pw,0)+1
top_pw=sorted(pathway_counts.items(),key=lambda x:-x[1])[:12]
fig5e,ax5e=plt.subplots(1,1,figsize=(7,5),facecolor="white")
colors=["#C0392B" if i<3 else "#3498DB" if i<6 else "#95A5A6" for i in range(len(top_pw))]
ax5e.barh([p for p,c in reversed(top_pw)],[c for p,c in reversed(top_pw)],
          color=list(reversed(colors)),edgecolor="white",height=0.65)
ax5e.set_xlabel("LR pair count",fontsize=10); ax5e.set_title("Top inferred signaling pathways",fontsize=11,fontweight="bold")
ax5e.tick_params(axis="y",labelsize=7.5); ax5e.tick_params(axis="x",labelsize=8)
plt.tight_layout(); img5e=fig_to_img(fig5e); plt.close(fig5e)

# F: Source-target communication summary (LR count by direction)
st_pairs={}
for r in int_all:
    src,tgt=r["source"],r["target"]
    if src in ["CAF_Fibroblast","Macrophage_TAM","Epithelial_like"] and tgt in ["CAF_Fibroblast","Macrophage_TAM","Epithelial_like"] and src!=tgt:
        key=f"{src} -> {tgt}"; st_pairs[key]=int(r["interaction_count"])
fig5f,ax5f=plt.subplots(1,1,figsize=(7,5),facecolor="white")
st_sorted=sorted(st_pairs.items(),key=lambda x:-x[1])
bar_colors=["#E74C3C" if "CAF" in k.split("->")[0] else "#3498DB" if "TAM" in k.split("->")[0] else "#27AE60" for k,c in st_sorted]
ax5f.barh([k for k,c in reversed(st_sorted)],[c for k,c in reversed(st_sorted)],
          color=list(reversed(bar_colors)),edgecolor="white",height=0.6)
ax5f.set_xlabel("LR pair count",fontsize=10); ax5f.set_title("Source-target communication strength",fontsize=11,fontweight="bold")
ax5f.tick_params(axis="y",labelsize=7.5); ax5f.tick_params(axis="x",labelsize=8)
add_note(ax5f,"CellChat-inferred; not experimentally validated")
plt.tight_layout(); img5f=fig_to_img(fig5f); plt.close(fig5f)

# Regenerate bubbles C/D from CSV
def make_bubble(rows,title,fs=(7.5,6)):
    srt=sorted(rows,key=lambda r:float(r["prob"]),reverse=True)[:25]
    srt=sorted(srt,key=lambda r:float(r["prob"]))
    labels=[f"{r['ligand']} - {r['receptor']}" for r in srt]
    probs=[float(r["prob"]) for r in srt]
    fig,ax=plt.subplots(1,1,figsize=fs,facecolor="white")
    ax.scatter(probs,labels,s=[p*900 for p in probs],c=probs,cmap="Reds",edgecolors="grey",linewidth=0.3,alpha=0.85)
    ax.set_xlabel("Communication probability",fontsize=9); ax.set_title(title,fontsize=11,fontweight="bold")
    ax.tick_params(axis="y",labelsize=7.5); ax.grid(axis="x",alpha=0.3,lw=0.5)
    ax.set_xlim(0,max(probs)*1.15)
    cbar=fig.colorbar(ax.collections[0],ax=ax,shrink=0.6,aspect=20,pad=0.02); cbar.set_label("Probability",fontsize=7)
    plt.tight_layout(); return fig

caf_tam=[r for r in lr_all if r["source"]=="CAF_Fibroblast" and r["target"]=="Macrophage_TAM"]
tam_caf=[r for r in lr_all if r["source"]=="Macrophage_TAM" and r["target"]=="CAF_Fibroblast"]
fc=make_bubble(caf_tam,"CAF -> TAM (inferred)"); img5c=fig_to_img(fc); plt.close(fc)
fd=make_bubble(tam_caf,"TAM -> CAF (inferred)"); img5d=fig_to_img(fd); plt.close(fd)

# Assemble Figure 5: 3 rows x 2 cols
fig5=plt.figure(figsize=(16,18),facecolor="white")
gs5=gridspec.GridSpec(3,2,figure=fig5,wspace=0.10,hspace=0.20,left=0.03,right=0.97,top=0.97,bottom=0.02)
panels5=[
    ("A",pdf_img(SRC/"cellchat/Fig8C_CellChat_CAF_TAM_Epithelial_network_redrawn.pdf")),
    ("B",pdf_img(SRC/"cellchat/Fig8D_CellChat_outgoing_incoming_heatmap.pdf")),
    ("C",img5c), ("D",img5d), ("E",img5e), ("F",img5f),
]
for idx,(letter,img) in enumerate(panels5):
    ax=fig5.add_subplot(gs5[idx//2,idx%2]); ax.imshow(img); ax.axis("off"); add_label(ax,letter)
    if letter in ("C","D","E","F"): add_note(ax,"Inferred communication")
fig5.suptitle("Figure 5. CellChat-inferred CAF-TAM-Epithelial communication network",fontsize=13,fontweight="bold",y=0.995)
save_fig(fig5,"Figure5_CellChat_communication_v5"); plt.close(fig5)

# ═══════════════════════════════════════════════════════════════
# FIGURE 6: Expand to A-F (6 panels)
# ═══════════════════════════════════════════════════════════════
print("\n=== Figure 6: Expand to 6 panels ===")

# E: GSE26566 validation - hub gene expression heatmap (from existing)
# F: Hub gene correlation summary (from existing cor heatmap)

fig6=plt.figure(figsize=(16,18),facecolor="white")
gs6=gridspec.GridSpec(3,2,figure=fig6,wspace=0.10,hspace=0.20,left=0.03,right=0.97,top=0.97,bottom=0.02)
panels6=[
    ("A",SRC/"integrated/Fig9A_integrated_evidence_heatmap.pdf",0,"Evidence heatmap"),
    ("B",SRC/"integrated/Fig9B_top_hub_genes_barplot.pdf",0,"Hub gene ranking"),
    ("C",SRC/"integrated/Fig9C_prioritized_LR_axis_dotplot.pdf",0,"LR axis dotplot"),
    ("D",SRC/"clinical_relevance/Fig10B_axis_scores_survival_forest_GSE107943.pdf",0,"Forest plot"),
    ("E",SRC/"GSE26566_validation/FigS_GSE26566_hub_gene_heatmap.pdf",0,"GSE26566 validation"),
    ("F",SRC/"clinical_relevance/Fig10I_hub_gene_score_correlation_heatmap_TCGA.pdf",0,"Hub gene correlation"),
]
notes6={"D":"Exploratory; GSE107943 only; not validated","E":"Expression-level validation only"}
for idx,(letter,path,pg,desc) in enumerate(panels6):
    ax=fig6.add_subplot(gs6[idx//2,idx%2])
    img=pdf_img(path,pg=pg); ax.imshow(img); ax.axis("off"); add_label(ax,letter)
    if letter in notes6: add_note(ax,notes6[letter])
fig6.suptitle("Figure 6. Integrated hub gene prioritization and exploratory clinical relevance",fontsize=13,fontweight="bold",y=0.995)
save_fig(fig6,"Figure6_hub_genes_clinical_relevance_v5"); plt.close(fig6)

# ═══════════════════════════════════════════════════════════════
# FIGURE 7: Expand to A-F with large mechanism
# ═══════════════════════════════════════════════════════════════
print("\n=== Figure 7: Expand to 6 panels ===")

# B: Drug-target network (from CSV)
drug_data=[]
with open(TBL/"drug_screening/candidate_drugs_prioritized.csv",newline="") as f:
    drug_data=list(csv.DictReader(f))

fig7b,ax7b=plt.subplots(1,1,figsize=(7,6),facecolor="white")
ax7b.set_xlim(0,10); ax7b.set_ylim(0,len(drug_data)*0.9+1)
# Simple bipartite: targets left, drugs right
targets=list(set(r.get("target","?") for r in drug_data[:15]))
drugs=list(set(r.get("drug","?") for r in drug_data[:15]))
# Draw as list
for i,r in enumerate(drug_data[:15]):
    t=r.get("target","?"); d=r.get("drug","?")
    ax7b.text(2,9-i*0.6,t,ha="right",fontsize=7,color="#E74C3C",fontweight="bold")
    ax7b.text(5,9-i*0.6,d,ha="left",fontsize=7,color="#3498DB")
    ax7b.plot([2.2,4.8],[9-i*0.6,9-i*0.6],color="#CCCCCC",lw=1)
ax7b.axis("off"); ax7b.set_title("Candidate drug-target relationships",fontsize=11,fontweight="bold")
plt.tight_layout(); img7b=fig_to_img(fig7b); plt.close(fig7b)

# C: Drugs per target barplot
tgt_counts={}
for r in drug_data:
    t=r.get("target","?"); tgt_counts[t]=tgt_counts.get(t,0)+1
fig7c,ax7c=plt.subplots(1,1,figsize=(7,5),facecolor="white")
sorted_t=sorted(tgt_counts.items(),key=lambda x:-x[1])
ax7c.barh([t for t,c in reversed(sorted_t)],[c for t,c in reversed(sorted_t)],color="#27AE60",edgecolor="white",height=0.6)
ax7c.set_xlabel("Candidate drugs",fontsize=10); ax7c.set_title("Drugs per target",fontsize=11,fontweight="bold")
ax7c.tick_params(axis="y",labelsize=7.5); plt.tight_layout(); img7c=fig_to_img(fig7c); plt.close(fig7c)

# E: Docking neighbor residues summary
res_data=[]
with open(TBL/"docking/docking_neighbor_residues_4A.csv",newline="") as f:
    res_data=list(csv.DictReader(f))
pair_counts={}
for r in res_data:
    pair=r.get("pair","?"); pair_counts[pair]=pair_counts.get(pair,0)+1
fig7e,ax7e=plt.subplots(1,1,figsize=(7,5),facecolor="white")
sorted_p=sorted(pair_counts.items(),key=lambda x:-x[1])
ax7e.barh([p for p,c in reversed(sorted_p)],[c for p,c in reversed(sorted_p)],color="#9B59B6",edgecolor="white",height=0.6)
ax7e.set_xlabel("Residues within 4 A",fontsize=10); ax7e.set_title("Docking neighbor residues",fontsize=11,fontweight="bold")
ax7e.tick_params(axis="y",labelsize=7.5); plt.tight_layout(); img7e=fig_to_img(fig7e); plt.close(fig7e)

# F: Large mechanism schematic
fig7f,ax7f=plt.subplots(1,1,figsize=(12,7),facecolor="white")
ax7f.set_xlim(0,14); ax7f.set_ylim(0,10); ax7f.axis("off")

# Three cell-type boxes (large)
cells=[("CAF\nFibroblast",2.5,5.5,"#E74C3C"),("Macrophage\nTAM",7.0,5.5,"#3498DB"),("Epithelial\nlike",11.5,5.5,"#27AE60")]
for label,x,y,color in cells:
    box=FancyBboxPatch((x-1.3,y-0.7),2.6,1.4,boxstyle="round,pad=0.08,rounding_size=0.15",
                        facecolor=color,edgecolor="white",lw=2,alpha=0.25)
    ax7f.add_patch(box); ax7f.text(x,y,label,ha="center",va="center",fontsize=10,fontweight="bold",color=color)

# Arrows
arrows=[(3.8,5.5,5.7,5.5,"COL1A1/COL1A2 -> CD44\nMIF -> CD74/CXCR4"),
        (5.7,4.8,3.8,4.8,"TGFB1 -> TGFBR\nPPIA -> BSG"),
        (11.5,6.2,8.3,6.2,"MIF -> CD74/CXCR4")]
for x1,y1,x2,y2,label in arrows:
    ax7f.annotate("",xy=(x2,y2),xytext=(x1,y1),arrowprops=dict(arrowstyle="->",color="#666666",lw=2.5))
    mx,my=(x1+x2)/2,(y1+y2)/2+0.3
    ax7f.text(mx,my,label,ha="center",va="bottom",fontsize=7,color="#333333",
              bbox=dict(boxstyle="round,pad=0.08",facecolor="white",edgecolor="none",alpha=0.85))

# TAM checkpoint box
box2=FancyBboxPatch((5.0,2.2),4.0,1.5,boxstyle="round,pad=0.1,rounding_size=0.12",
                     facecolor="#3498DB",edgecolor="white",lw=1.5,alpha=0.15)
ax7f.add_patch(box2)
ax7f.text(7.0,3.2,"TAM checkpoint transcripts",ha="center",fontsize=8.5,fontweight="bold",color="#3498DB")
ax7f.text(7.0,2.7,"HAVCR2 | IDO1 | CD86 | PDCD1LG2",ha="center",fontsize=7.5,color="#333333")
ax7f.text(7.0,2.4,"mRNA-level observations; protein validation required",ha="center",fontsize=6.5,fontstyle="italic",color="#888888")

# Candidate targets
box3=FancyBboxPatch((5.0,0.8),4.0,1.0,boxstyle="round,pad=0.08,rounding_size=0.10",
                     facecolor="#27AE60",edgecolor="white",lw=1.5,alpha=0.15)
ax7f.add_patch(box3)
ax7f.text(7.0,1.4,"Candidate targets: CXCR4 | IDO1 | TGFBR1 | MIF",ha="center",fontsize=8,fontweight="bold",color="#27AE60")
ax7f.text(7.0,1.0,"Docking = in silico predictions only",ha="center",fontsize=6.5,fontstyle="italic",color="#888888")

# Title and caveat
ax7f.text(7.0,9.3,"Proposed CAF-TAM-Epithelial communication model",ha="center",fontsize=12,fontweight="bold")
ax7f.text(7.0,0.2,"Computational model; experimental validation required.",ha="center",fontsize=7.5,fontstyle="italic",color="#888888")
plt.tight_layout(); img7f=fig_to_img(fig7f); plt.close(fig7f)

# Assemble Figure 7: 3 rows = A+B+C, D+E, F-span
fig7=plt.figure(figsize=(20,18),facecolor="white")
gs7=gridspec.GridSpec(3,3,figure=fig7,wspace=0.10,hspace=0.20,left=0.03,right=0.97,top=0.97,bottom=0.02)
row1=[("A",SRC/"drug_screening/Fig11A_target_druggability_heatmap.pdf",0),
      ("B",None,0,img7b),("C",None,0,img7c)]
row2=[("D",SRC/"docking/FigS_Docking_affinity_barplot_all_pairs.pdf",0),
      ("E",None,0,img7e)]

for idx,(letter,path,pg,*rest) in enumerate(row1):
    ax=fig7.add_subplot(gs7[0,idx]); img=rest[0] if rest else pdf_img(path,pg=pg)
    ax.imshow(img); ax.axis("off"); add_label(ax,letter)
    if letter=="D": add_note(ax,"In silico predictions only")
for idx,(letter,path,pg,*rest) in enumerate(row2):
    ax=fig7.add_subplot(gs7[1,idx]); img=rest[0] if rest else pdf_img(path,pg=pg)
    ax.imshow(img); ax.axis("off"); add_label(ax,letter)

# Hide col 2 for row 1 continuation... actually row2 has D+E at col 0,1; col 2 hidden
for idx in range(2,3):
    ax=fig7.add_subplot(gs7[1,idx]); ax.axis("off")

# F spans full row 2
ax_f=fig7.add_subplot(gs7[2,:])
ax_f.imshow(img7f); ax_f.axis("off"); add_label(ax_f,"F")

fig7.suptitle("Figure 7. Candidate target prioritization and in silico docking",fontsize=13,fontweight="bold",y=0.995)
save_fig(fig7,"Figure7_therapeutic_implications_v5"); plt.close(fig7)

# ═══════════════════════════════════════════════════════════════
print("\n"+"="*60)
print("ALL MAIN FIGURES V5 GENERATED")
for f in sorted(OUT.glob("*.pdf")):
    print(f"  {f.name} ({f.stat().st_size//1024} KB)")
