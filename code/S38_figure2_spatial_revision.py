# -*- coding: utf-8 -*-
"""Redraw Figure 2 from the corrected, versioned spatial outputs.

Inputs are the current post-provenance-correction files:
- S26_spatial_scores.csv (22,915 spots, 16 sections)
- S26b_global_gradient.json
- S29_edge_dist.csv / S29b_gradient_specificity.json
- spatial_donor_sensitivity.json (S36 output)
- S35_selfaudit_redundancy.json
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

ROOT=Path(__file__).resolve().parents[1]
TABLES=ROOT/'processed_tables'
RESULTS=ROOT/'processed_results'
OUTDIR=ROOT/'figures'
OUTDIR.mkdir(parents=True,exist_ok=True)
OUT=OUTDIR/'Figure2_revised.png'

plt.rcParams.update({
    'font.family':'DejaVu Sans','font.size':8.2,
    'axes.linewidth':0.75,'axes.spines.top':False,'axes.spines.right':False,
    'xtick.major.width':0.7,'ytick.major.width':0.7,
    'figure.dpi':300,'savefig.dpi':300,
})

SP=pd.read_csv(TABLES/'S26_spatial_scores.csv')
with open(RESULTS/'S26b_global_gradient.json',encoding='utf-8') as f: GJ=json.load(f)
with open(RESULTS/'S29b_gradient_specificity.json',encoding='utf-8') as f: SJ=json.load(f)
with open(RESULTS/'S36_donor_stratified_spatial_sensitivity.json',encoding='utf-8') as f: DS=json.load(f)
with open(RESULTS/'S35_selfaudit_redundancy.json',encoding='utf-8') as f: SA=json.load(f)

COND=['Skin','Wound1','Wound7','Wound30']
XLAB=['D0 (Skin)','D1 (Wound1)','D7 (Wound7)','D30 (Wound30)']
colors={'F6':'#c51f27','F7':'#1f77b4','F7e':'#169b62','F5':'#68707d'}

# Depth-residualize individual marker z-scores exactly as in S26b.
def depth_res_gene_z(df):
    zcols=[c for c in df.columns if c.endswith('_z')]
    log_n=np.log1p(df['total_counts'].values.astype(float))
    lc=log_n-log_n.mean(); den=(lc**2).sum()+1e-12
    out={}
    for c in zcols:
        v=df[c].values.astype(float)
        beta=(lc*(v-v.mean())).sum()/den
        out[c[:-2]]=v-beta*lc
    return pd.DataFrame(out,index=df.index)

MARK={
 'F6':['CXCL8','MMP3','CHI3L1','CXCL1','CCL2'],
 'F7':['LRRC15','SFRP4','ASPN','POSTN','CTHRC1','ACTA2','TAGLN'],
 'F7e':['COL1A1','COL1A2','COL3A1','FN1','POSTN','LOX'],
 'F1':['COL13A1','COL23A1','SFRP2','CRABP1','WNT5A'],
 'F5':['SCN7A','SOX10','PLP1','NGFR','S100B'],
}

fig=plt.figure(figsize=(14.2,10.0))
gs=fig.add_gridspec(2,2,height_ratios=[1.0,1.05],hspace=0.34,wspace=0.28)

# ---------- Panel A ----------
axA=fig.add_subplot(gs[0,0])
x=np.arange(4)

def one_violin(ax,values,xpos,color,side_width=0.13):
    v=np.asarray(values,float)
    lo,hi=np.quantile(v,[0.01,0.99])
    vc=np.clip(v,lo,hi)
    ys=np.linspace(lo,hi,180)
    sd=float(np.std(vc))
    bw=1.06*sd*(len(vc)**(-1/5)) if sd>0 else 0.05
    bw=max(bw,0.035)
    z=(ys[:,None]-vc[None,:])/bw
    den=np.exp(-0.5*z*z).sum(axis=1)/(len(vc)*bw*np.sqrt(2*np.pi))
    den=den/(den.max()+1e-12)*side_width
    ax.fill_betweenx(ys,xpos-den,xpos+den,color=color,alpha=0.32,lw=0)
    q1,q2,q3=np.quantile(v,[.25,.5,.75])
    ax.plot([xpos,xpos],[q1,q3],color=color,lw=2.3,alpha=.85,zorder=4)
    ax.plot([xpos-.018,xpos+.018],[q2,q2],color='white',lw=1.15,zorder=5)

for i,c in enumerate(COND):
    sub=SP[SP['Condition']==c]
    one_violin(axA,sub['F6_res'],i-.17,colors['F6'])
    one_violin(axA,sub['F7e_res'],i+.17,colors['F7e'])

means={p:[SP.loc[SP['Condition']==c,p+'_res'].mean() for c in COND] for p in ['F6','F7','F7e']}
styles={'F6':dict(color=colors['F6'],marker='o',lw=2.1,ms=5.5,ls='-'),
        'F7':dict(color=colors['F7'],marker='^',lw=1.8,ms=5.2,ls='--'),
        'F7e':dict(color=colors['F7e'],marker='s',lw=2.1,ms=5.2,ls='-')}
for p,st in styles.items():
    axA.plot(x,means[p],zorder=6,**st)
# Label only the key quantitative anchors; full distributions remain visible without numeric clutter.
key_labels=[
    (1,means['F6'][1],1.02,0.49,colors['F6']),
    (1,means['F7e'][1],1.02,-0.82,colors['F7e']),
    (3,means['F7'][3],2.74,0.24,colors['F7']),
    (3,means['F7e'][3],2.72,0.84,colors['F7e']),
    (3,means['F6'][3],2.73,-0.39,colors['F6']),
]
for xx,yy,tx,ty,cc in key_labels:
    axA.text(tx,ty,f'{yy:+.3f}',fontsize=7.0,color=cc,fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.14',fc='white',ec='none',alpha=.82))
axA.axhline(0,color='#555',lw=.65,ls=':',alpha=.75)
counts=[int((SP['Condition']==c).sum()) for c in COND]
axA.set_xticks(x)
axA.set_xticklabels([f'{a}\n(n={b:,})' for a,b in zip(XLAB,counts)],fontsize=7.4)
axA.set_ylabel('Depth-residualized spot score',fontsize=8.7)
axA.set_ylim(-1.23,1.92)
axA.set_xlim(-.48,3.55)
axA.set_title('A. Temporal spot-score distributions across 16 Visium sections',loc='left',fontweight='bold',fontsize=9.6)
# compact legend with inferential statistics
legend_lines=[
    f"F6 inflammatory (peak D1 {means['F6'][1]:+.3f})",
    f"F7e core ECM (trend r={GJ['empirical_p']['F7e']['obs']:+.3f}, p={GJ['empirical_p']['F7e']['p_emp']:.3f})",
    f"F7 broad (trend r={GJ['empirical_p']['F7']['obs']:+.3f}, p={GJ['empirical_p']['F7']['p_emp']:.3f}, ns)",
]
handles=[plt.Line2D([],[],color=styles[p]['color'],marker=styles[p]['marker'],ls=styles[p]['ls'],lw=styles[p]['lw'],markersize=4.8,label=lab)
         for p,lab in zip(['F6','F7e','F7'],legend_lines)]
axA.legend(handles=handles,frameon=True,fontsize=6.6,loc='upper left',edgecolor='#cccccc',framealpha=.92,borderpad=.45,labelspacing=.35)
axA.text(.02,.035,f"Candidate-null tests: F7e p = {GJ['empirical_p']['F7e']['p_emp']:.3f}; "
                   f"F1 p = {GJ['empirical_p']['F1']['p_emp']:.3f}; F5 p = {GJ['empirical_p']['F5']['p_emp']:.3f}",
         transform=axA.transAxes,fontsize=6.9,color='#155e46',fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.28',fc='#effaf4',ec='#7fc6a4',lw=.7))

# ---------- Panel B ----------
axB=fig.add_subplot(gs[0,1])
GZ=depth_res_gene_z(SP)
row_order=['FB1','Smooth_muscle','FB2','Papillary_dermis','Dermis2','Dermis1','Sweatgland_FB','Immune_endo','Immune','Wound_edges','Basal_epi','Suprabasal_epi']
row_lab=['FB1','Smooth muscle','FB2','Papillary dermis','Dermis2','Dermis1','Sweat gland FB','Immune / endothelial','Immune','Wound edge','Basal epidermis','Suprabasal epidermis']
cols=['COL1A1','COL1A2','COL3A1','FN1','LOX','ACTA2','TAGLN','CTHRC1','F7e','F1']
# Compute color values and size values.
prog_scores={
    'F7e':GZ[[g for g in MARK['F7e'] if g in GZ.columns]].mean(axis=1),
    'F1':GZ[[g for g in MARK['F1'] if g in GZ.columns]].mean(axis=1),
}
val=[]; frac=[]
for anno in row_order:
    mask=(SP['AnnoType']==anno).values
    rv=[]; rf=[]
    for col in cols:
        if col in prog_scores:
            s=prog_scores[col].values
            rv.append(float(np.mean(s[mask])))
            rf.append(float(np.mean(s[mask]>0)))
        else:
            rv.append(float(GZ.loc[mask,col].mean()))
            rf.append(float((SP.loc[mask,col].values>0).mean()))
    val.append(rv); frac.append(rf)
val=np.array(val); frac=np.array(frac)
norm=Normalize(vmin=-1.0,vmax=1.0)
cmap=plt.get_cmap('RdBu_r')
for yi in range(len(row_order)):
    axB.axhline(yi,color='#eeeeee',lw=.55,zorder=0)
    for xi in range(len(cols)):
        size=22+250*frac[yi,xi]
        axB.scatter(xi,yi,s=size,c=[cmap(norm(np.clip(val[yi,xi],-1,1)))],edgecolors='#4c4c4c',linewidths=.45,zorder=3)
# group dividers and headers
for xx in [4.5,7.5]:
    axB.axvline(xx,color='#666',lw=.8,ls='--')
axB.text(2,-.85,'Fibrillar ECM\n(FB1-specific)',ha='center',va='bottom',fontsize=7.4,color='#087f5b',fontweight='bold')
axB.text(6,-.85,'Contractile\n(mural spillover)',ha='center',va='bottom',fontsize=7.4,color='#b4232a',fontweight='bold')
axB.text(8.5,-.85,'Score',ha='center',va='bottom',fontsize=7.4,color='#3155a4',fontweight='bold')
axB.set_xticks(range(len(cols))); axB.set_xticklabels(cols,rotation=48,ha='right',fontsize=7.1)
axB.set_yticks(range(len(row_order))); axB.set_yticklabels(row_lab,fontsize=7.2)
axB.invert_yaxis()
axB.set_xlim(-.65,9.72); axB.set_ylim(len(row_order)-.45,-1.35)
axB.set_title('B. Spatial DotPlot: matrix programs versus mural markers',loc='left',fontweight='bold',fontsize=9.6)
# colorbar
cax=axB.inset_axes([1.015,.18,.025,.66])
cb=plt.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=cax)
cb.set_label('Mean depth-res. z',fontsize=6.5,labelpad=2)
cb.ax.tick_params(labelsize=5.8,length=2)

# ---------- Panel C ----------
axC=fig.add_subplot(gs[1,0])
ED=pd.read_csv(TABLES/'S29_edge_dist.csv')
SP['bc']=SP['barcode'].astype(str); ED['bc']=ED['barcode'].astype(str)
M=SP.merge(ED[['bc','h5_stem','edge_dist']],on=['bc','h5_stem'],how='inner')
GZm=depth_res_gene_z(M)
wm=(M['Condition']=='Wound7').values
dist=M.loc[wm,'edge_dist'].values.astype(float)
# Two-unit bins through the complete Day-7 range.
edges=np.arange(0,36,2)
centers=(edges[:-1]+edges[1:])/2
prog_c={}
for p in ['F6','F7','F7e','F5']:
    genes=[g for g in MARK[p] if g in GZm.columns]
    prog_c[p]=GZm.loc[wm,genes].mean(axis=1).values
for p in ['F6','F7e','F7','F5']:
    mu=[]; lo=[]; hi=[]; xs=[]
    for a,b,ce in zip(edges[:-1],edges[1:],centers):
        mm=(dist>=a)&(dist<b)
        if mm.sum()<25: continue
        v=prog_c[p][mm]
        muv=float(np.mean(v)); se=float(np.std(v,ddof=1)/np.sqrt(len(v)))
        xs.append(ce); mu.append(muv); lo.append(muv-1.96*se); hi.append(muv+1.96*se)
    xs=np.array(xs); mu=np.array(mu); lo=np.array(lo); hi=np.array(hi)
    if p in styles:
        st=styles[p]
    else:
        st=dict(color=colors['F5'],marker='D',lw=1.7,ms=4.2,ls=':')
    axC.plot(xs,mu,color=st['color'],lw=st['lw'],ls=st['ls'],marker=st['marker'],ms=st['ms'],label=p,zorder=4)
    axC.fill_between(xs,lo,hi,color=st['color'],alpha=.14,lw=0,zorder=2)
axC.axvspan(0,2.5,color='#c51f27',alpha=.055,lw=0)
axC.axhline(0,color='#555',lw=.65,ls=':',alpha=.75)
axC.set_xlim(0,35); axC.set_ylim(-.48,.36)
axC.set_xlabel('Distance from nearest wound-edge spot (spot-diameter units)',fontsize=8.4)
axC.set_ylabel('Mean depth-residualized score ± 95% CI',fontsize=8.4)
axC.set_title('C. Day-7 continuous spatial gradient from the wound edge',loc='left',fontweight='bold',fontsize=9.6)
# Donor-level statistics from S36.
f6d=DS['wound7_donor_level']['Wound7']['F6']; f7ed=DS['wound7_donor_level']['Wound7']['F7e']
f6r=[v['rho'] for v in f6d.values()]; f7er=[v['rho'] for v in f7ed.values()]
dd=SJ['double_diff']['F6 − F5']['Wound7']
ann=(f"Day-7 donor-stratified F6: 4/4 negative\n"
     f"mean rho = {np.mean(f6r):+.3f} (range {min(f6r):+.3f} to {max(f6r):+.3f})\n"
     f"F6−F5 double difference: rho = {dd['rho']:+.3f}, exploratory p = {dd['p']:.3f}")
axC.text(.035,.965,ann,transform=axC.transAxes,ha='left',va='top',fontsize=6.7,color='#8f171c',fontweight='bold',
         bbox=dict(boxstyle='round,pad=0.38',fc='#fff2f2',ec='#e69a9a',lw=.8))
axC.legend(frameon=False,fontsize=6.7,loc='lower left',ncol=2)
axC.text(.98,.035,'All 6,044 Day-7 wound spots with edge distances',transform=axC.transAxes,ha='right',va='bottom',fontsize=6.5,color='#555')

# ---------- Panel D ----------
axD=fig.add_subplot(gs[1,1])
# Raw program means are used here because the subtype-specificity audit (S35) is defined on category-level program means.
annmeans=SP.groupby('AnnoType').agg(n=('F7','size'),F7=('F7','mean'),F7e=('F7e','mean'))
rows=[r for r in row_order if r in annmeans.index]
ylabel_map=dict(zip(row_order,row_lab))
y=np.arange(len(rows))
for i,r in enumerate(rows):
    a=annmeans.loc[r,'F7']; b=annmeans.loc[r,'F7e']
    axD.plot([a,b],[i,i],color='#aeb9c3',lw=2.0,zorder=2)
    axD.scatter(a,i,s=48,color='#35a6db',edgecolor='#155a78',lw=.6,zorder=4,label='F7 broad' if i==0 else None)
    axD.scatter(b,i,s=48,color='#15965e',edgecolor='#075a38',lw=.6,zorder=4,label='F7e core ECM' if i==0 else None)
axD.axvline(0,color='#444',lw=.8,ls='--')
axD.set_yticks(y); axD.set_yticklabels([f'{ylabel_map[r]} (n = {int(annmeans.loc[r,"n"]):,})' for r in rows],fontsize=7.1)
axD.invert_yaxis()
axD.set_xlim(-.62,1.12)
axD.set_xlabel('Mean spot-level program score (z-scale)',fontsize=8.4)
axD.set_title('D. Subtype specificity: FB1 enrichment versus smooth-muscle dilution',loc='left',fontweight='bold',fontsize=9.6)
q=SA['q1_fb1_vs_otherFIB']; pq=SA['perm_q1']
sm=annmeans.loc['Smooth_muscle']
box=(f"FB1 subtype specificity:\n"
     f"F7e diff. = {q['F7e']['delta']:+.3f} (p = {pq['F7e']['p']:.4f})\n"
     f"F7 diff. = {q['F7']['delta']:+.3f} (p = {pq['F7']['p']:.3f}, ns)")
axD.text(.67,.34,box,transform=axD.transAxes,fontsize=6.8,color='#075a38',fontweight='bold',va='center',
         bbox=dict(boxstyle='round,pad=0.38',fc='#effaf4',ec='#79b794',lw=.8))
smtext=f"Mural smooth muscle:\nF7 = {sm['F7']:+.3f}  →  F7e = {sm['F7e']:+.3f}"
axD.text(.69,.80,smtext,transform=axD.transAxes,fontsize=6.8,color='#9c1f24',fontweight='bold',va='center',
         bbox=dict(boxstyle='round,pad=0.38',fc='#fff2f2',ec='#e69a9a',lw=.8))
sm_i=rows.index('Smooth_muscle')
axD.annotate('',xy=(sm['F7'],sm_i),xytext=(.72,.76),textcoords='axes fraction',
             arrowprops=dict(arrowstyle='->',color='#c51f27',lw=1.0))
axD.legend(frameon=True,fontsize=6.7,loc='lower right',edgecolor='#cccccc')

fig.suptitle('Spatial transcriptomic audit of vascular admixture and matrix compartmentalization (GSE241124)',
             fontsize=12.0,fontweight='bold',y=.985)
fig.text(.5,.006,'Corrected barcode-to-section provenance: 16 sections and 22,915 annotated spots from 4 donors. '
                 'Spatial permutation p-values are exploratory because they do not preserve spatial autocorrelation, section membership, or donor clustering.',
         ha='center',fontsize=6.8,color='#555')
fig.savefig(OUT,bbox_inches='tight',facecolor='white')
fig.savefig(OUT.with_suffix('.pdf'),bbox_inches='tight',facecolor='white')
plt.close(fig)
print('saved',OUT,OUT.stat().st_size)
print('panel A means',means)
print('panel C n wound7',int(wm.sum()),'F6 donor mean',float(np.mean(f6r)),'F7e donor mean',float(np.mean(f7er)))
