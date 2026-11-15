import json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from collections import Counter
out=json.load(open('../results/full_out.json')); pops=json.load(open('../results/full_pops.json'))
plt.rcParams.update({"font.family":"Liberation Sans","font.size":8,"axes.edgecolor":"#000000","axes.linewidth":0.6,
  "xtick.major.width":0.6,"ytick.major.width":0.6,"xtick.major.size":3,"ytick.major.size":3,"pdf.fonttype":42,"ps.fonttype":42})
BLUE,ORANGE,GRAY="#2a78d6","#eb6834","#8a8984"
yrs=[r['year'] for r in out]; T=[r['T'] for r in out]
fig,(a1,a2)=plt.subplots(2,1,figsize=(5.0,4.6),gridspec_kw=dict(height_ratios=[2.6,1.2],hspace=0.55))
def panel(ax,letter,title):
    ax.text(0,1.06,letter,transform=ax.transAxes,fontsize=8,fontweight="bold",va="bottom",ha="left",
            bbox=dict(boxstyle="square,pad=0.15",fc="white",ec="black",lw=0.6))
    ax.text(0.045,1.06,title,transform=ax.transAxes,fontsize=8,va="bottom",ha="left")
for ax in (a1,a2):
    ax.grid(axis="y",color="#d9d9d9",lw=0.5); ax.set_axisbelow(True)
    for s in ("top","right"): ax.spines[s].set_visible(False)
    ax.set_xlim(2005.5,2026.5); ax.set_xticks(range(2006,2027,2))
a1.axvspan(2020.5,2026.5,color="#efefef",zorder=0,lw=0)
for key,lo,hi,col,lab in (("E","Elo","Ehi",BLUE,"Eligible"),("R","Rlo","Rhi",ORANGE,"Expected to respond")):
    y=[t[key]*100 for t in T]; l=[t[lo]*100 for t in T]; h=[t[hi]*100 for t in T]
    a1.fill_between(yrs,l,h,step="post",color=col,alpha=0.22,lw=0)
    a1.step(yrs,y,where="post",color=col,lw=1.4,label=lab)
a1.plot([2018,2020],[8.82,13.60],ls="none",marker="o",ms=4,mfc="white",mec=BLUE,mew=1.1,label="Published estimate, eligible")
a1.plot([2018,2020],[5.48,7.04],ls="none",marker="s",ms=3.6,mfc="white",mec=ORANGE,mew=1.1,label="Published estimate, respond")
a1.text(2023.5,19.2,"2021-2026",ha="center",fontsize=7,color="#555555")
a1.set_ylim(0,20); a1.set_yticks(range(0,21,5)); a1.set_ylabel("US cancer deaths, %")
a1.legend(loc="upper left",frameon=False,fontsize=7,handlelength=1.8,borderaxespad=0.2)
panel(a1,"A","Eligibility for and expected response to genome-targeted therapy")
cnt=Counter(p["start"] for p in pops if p["cls"]=="T")
a2.bar(yrs,[cnt.get(y,0) for y in yrs],width=0.65,color=BLUE,edgecolor="none")
a2.bar([2026],[-1],width=0.65,color=GRAY,edgecolor="none")
a2.axhline(0,color="black",lw=0.6)
a2.set_ylim(-2,9); a2.set_yticks([0,2,4,6,8]); a2.set_ylabel("Populations, No."); a2.set_xlabel("Year")
panel(a2,"B","Biomarker-defined populations newly eligible each year")
fig.savefig("../figures/Figure1_JAMA.pdf",bbox_inches="tight")
fig.savefig("../figures/Figure1_JAMA.tif",dpi=600,bbox_inches="tight",pil_kwargs={"compression":"tiff_lzw"})
fig.savefig("../figures/Figure1_JAMA.png",dpi=300,bbox_inches="tight")
