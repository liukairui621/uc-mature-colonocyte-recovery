"""Two publication figures with vertically stacked panels and generous gutters."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'runs/005_health_function';OUT=ROOT/'figures/publication'
BLUE='#2F6690';GREEN='#26867C';GRAY='#666666';ORANGE='#B56727'
def style():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.labelsize':10,'xtick.labelsize':9,'ytick.labelsize':9,'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,stem):
    OUT.mkdir(parents=True,exist_ok=True)
    fig.savefig(OUT/(stem+'.pdf'));fig.savefig(OUT/(stem+'.png'),dpi=300);fig.savefig(OUT/(stem+'_preview.png'),dpi=150)
    plt.close(fig)
def health():
    d=pd.read_csv(DATA/'healthy_reference_contrasts.tsv',sep='\t')
    d=d[d.scope.eq('main')&d.outcome.eq('Yes')]
    order=['GSE92415_W6','GSE73661_IFX_W4W6','GSE73661_VDZ_W6','GSE73661_VDZ_W12']
    labels=['GLM / placebo W6 responders (n=32)','IFX W4–6 healed (n=8)','VDZ W6 healed (n=6)','VDZ W12 healed (n=3)']
    fig,axs=plt.subplots(2,1,figsize=(7.6,7.9));fig.subplots_adjust(left=.37,right=.96,top=.93,bottom=.09,hspace=.70)
    for ax,prog,title,color in zip(axs,['CANDIDATE6','COMMON194_INFLAMMATION'],['A  Six-gene colonocyte score','B  Inflammatory-response score'],[BLUE,ORANGE]):
        z=d[d.program.eq(prog)].set_index('context').loc[order]
        y=np.arange(len(z))[::-1];ax.axvline(0,color='0.5',lw=1,ls='--',zorder=0)
        ax.errorbar(z.estimate,y,xerr=[z.estimate-z.lower,z.upper-z.estimate],fmt='o',color=color,ecolor=color,capsize=3,ms=6,lw=1.6)
        ax.set_yticks(y);ax.set_yticklabels(labels);ax.set_ylim(-.6,3.65)
        ax.set_title(title,loc='left',pad=18,fontweight='bold');ax.set_xlabel('Post-treatment minus non-IBD control\n(mean log2-expression score; 95% CI)',labelpad=9)
        ax.tick_params(axis='y',length=0,pad=9);ax.grid(axis='x',alpha=.12)
    save(fig,'Figure5_healthy_reference')
def function():
    d=pd.read_csv(DATA/'ct_program_models.tsv',sep='\t');d=d[d.model.eq('baseline_ANCOVA')].set_index('program')
    c=pd.read_csv(DATA/'ct_candidate_coordination.tsv',sep='\t').set_index('program')
    order=['CT_MARKERS_EXCLUDING_CANDIDATE6','GOBP_INTESTINAL_ABSORPTION','GOBP_BRUSH_BORDER_ASSEMBLY'];labels=['CT markers excluding Candidate6 (16 genes)','Intestinal absorption (44 genes)','Brush-border assembly (7 genes)']
    fig,axs=plt.subplots(2,1,figsize=(7.6,7.7));fig.subplots_adjust(left=.41,right=.94,top=.93,bottom=.10,hspace=.78)
    ax=axs[0];z=d.loc[order];y=np.arange(3)[::-1]
    ax.axvline(0,color='0.5',ls='--',lw=1);ax.errorbar(z.estimate,y,xerr=[z.estimate-z.lower,z.upper-z.estimate],fmt='o',color=BLUE,capsize=3,lw=1.6)
    ax.set_yticks(y);ax.set_yticklabels(labels);ax.set_ylim(-.6,2.65);ax.set_title('A  Baseline-adjusted remission associations',loc='left',fontweight='bold',pad=18)
    ax.set_xlabel('Remission coefficient\n(mean log2-CPM score; HC3 95% CI)',labelpad=10);ax.tick_params(axis='y',length=0,pad=9);ax.grid(axis='x',alpha=.12)
    ax=axs[1];z=c.loc[order];ax.axvline(0,color='0.5',ls='--',lw=1)
    ax.errorbar(z.pearson_r,y,xerr=[z.pearson_r-z.lower,z.upper-z.pearson_r],fmt='o',color=GREEN,capsize=3,lw=1.6)
    ax.set_yticks(y);ax.set_yticklabels(labels);ax.set_ylim(-.6,2.65);ax.set_xlim(-.8,.8)
    ax.set_title('B  Covariation with six-gene score change',loc='left',fontweight='bold',pad=18)
    ax.set_xlabel('Patient-level Pearson correlation\n(Fisher-z 95% CI)',labelpad=10);ax.tick_params(axis='y',length=0,pad=9);ax.grid(axis='x',alpha=.12)
    save(fig,'Figure6_CT_functional_programs')
if __name__=='__main__':style();health();function()
