"""Complete patient panels; no outcome-dependent case/gene selection."""
from contract import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
COL={'M0':'#666666','M1':'#0072B2','M2':'#009E73','M3':'#D55E00'}
F=B/'figures';F.mkdir(exist_ok=True)
def save(fig,name):fig.savefig(F/(name+'.png'),dpi=170,bbox_inches='tight');plt.close(fig)
def main():
 pt=pd.read_csv(B/'C9B_PATIENT_METRICS.csv');ef=pd.read_csv(B/'C9B_PATIENT_EFFECTS.csv');r=pt[(pt.regime=='V1')&(pt.strategy=='random')];rows=[json.loads(p.read_text()) for p in (B/'logs/runs').glob('*.json')]
 fig,axes=plt.subplots(1,3,figsize=(13,4))
 for ax,s in zip(axes,EVAL):
  rr=next(x for x in rows if x['sample_id']==s and x['regime']=='V1' and x['seed']==0 and x['strategy']=='random' and x['budget']==.2 and not x['control']);sp=RAW/'spatial'/(s+'_spatial');im=Image.open(sp/'tissue_hires_image.png')
  with np.load(CACHE/(s+'_metadata.npz')) as f:xy=f['xy_px']*f['hires_scale']
  with np.load(ROOT/rr['prediction_path']) as f:obs=f['observed_indices']
  ax.imshow(im);ax.scatter(xy[obs,0],xy[obs,1],s=2,c='#00BFC4');ax.set_title(s+' | 20% anchors');ax.axis('off')
 fig.suptitle('Fresh C9 cohort: all held-out patients; random seed0');save(fig,'01_HE_ANCHORS')
 fig,axes=plt.subplots(1,3,figsize=(12,3.5))
 for ax,metric,title in zip(axes,['median_pearson','mean_mse','mean_mae'],['Median gene Pearson','Mean gene MSE','Mean gene MAE']):
  for j,m in enumerate(COL):
   vals=r[(r.budget==.2)&(r.method==m)][metric];ax.scatter(np.repeat(j,len(vals)),vals,color=COL[m]);ax.plot([j-.15,j+.15],[vals.median()]*2,color='black')
  ax.set_xticks(range(4),COL);ax.set_ylabel(title);ax.grid(alpha=.2)
 fig.suptitle('20% V1: points are patients; bars are patient medians');save(fig,'02_METHOD_COMPARISON')
 for metric,title,name in [('median_pearson','Median gene Pearson','03_COVERAGE_PEARSON'),('mean_mse','Mean gene MSE (common units)','04_COVERAGE_MSE')]:
  fig,axes=plt.subplots(1,3,figsize=(12,3.5))
  for ax,s in zip(axes,EVAL):
   for m in COL:
    p=r[(r.patient_id==s)&(r.method==m)].sort_values('budget');ax.plot(p.budget*100,p[metric],marker='o',label=m,color=COL[m],linestyle='--' if m=='M1' else '-')
   ax.set_title(s);ax.set_xlabel('Measured coverage (%)');ax.set_ylabel(title);ax.grid(alpha=.2)
  axes[-1].legend();fig.suptitle('V1 random; M1 and M2 overlap');save(fig,name)
 fig,axes=plt.subplots(1,2,figsize=(10,4))
 for ax,metric in zip(axes,['median_pearson','mean_mse']):
  for j,s in enumerate(EVAL):
   p=pt[(pt.patient_id==s)&(pt.strategy=='random')&(pt.budget==.2)&(pt.method=='M2')].set_index('regime');ax.plot([0,1],p.loc[['V1','V2'],metric],marker='o',label=s)
  ax.set_xticks([0,1],['V1','V2 (3 geometries)']);ax.set_ylabel(metric);ax.grid(alpha=.2)
 axes[-1].legend();fig.suptitle('Matched20% budget; V2 scores compact forbidden40% regions');save(fig,'05_V1_V2')
 d=pd.read_csv(B/'C9B_DISTANCE_TO_ANCHOR.csv');d=d[d.control.isna()&(d.strategy=='random')&(d.budget==.2)&(d.method=='M2')];dp=d.groupby(['patient_id','regime','lower_um','geometry']).mse.mean().groupby(['patient_id','regime','lower_um']).mean().reset_index();dp.to_csv(B/'C9B_DISTANCE_PATIENT_CURVES.csv',index=False)
 fig,axes=plt.subplots(1,3,figsize=(12,3.5))
 for ax,s in zip(axes,EVAL):
  for regime in ['V1','V2']:
   p=dp[(dp.patient_id==s)&(dp.regime==regime)].sort_values('lower_um');ax.plot(p.lower_um,p.mse,marker='o',label=regime)
  ax.set_title(s);ax.set_xlabel('Nearest-anchor bin lower edge (um)');ax.set_ylabel('M2 MSE');ax.grid(alpha=.2)
 axes[-1].legend();fig.suptitle('Distance curves: occupied bins only; errors in common units');save(fig,'06_DISTANCE')
 for comparison,name,title in [('sparse_ST','07_PATIENT_EFFECTS','Sparse ST minus frozen H&E'),('harmonization','08_HARMONIZATION','Harmonization: affine invariance, no measurable gain'),('image','09_IMAGE_GAIN','Incremental image: M3 minus M2')]:
  p=ef[(ef.comparison==comparison)&(ef.strategy=='random')&(ef.regime=='V1')&(ef.budget==.2)];fig,axes=plt.subplots(1,3,figsize=(12,3.7))
  for ax,col in zip(axes,['delta_pearson','relative_mse_reduction','mae_reduction']):
   vals=p[col].to_numpy();vals=np.where(abs(vals)<1e-12,0,vals);ax.bar(p.patient_id,vals,color=['#009E73' if x>=0 else '#D55E00' for x in vals]);ax.axhline(0,color='black',linewidth=.8);ax.tick_params(axis='x',rotation=25);ax.set_ylabel(col);ax.grid(axis='y',alpha=.2)
   if comparison=='harmonization':ax.set_ylim(-.01,.01)
  fig.suptitle(title+' |20% V1; positive is better');save(fig,name)
 tr=pd.read_csv(OUT/'C9_TRANSPORTABILITY_SUMMARY.csv');tr=tr[(tr.coverage==.2)&(tr.mask_regime=='V1')];selection=[('development_C5_primary','T0_E1'),('HER2ST_C7_zero_shot','T0_E1'),('HER2ST_C8A_harmonized','T1_E1'),('fresh_C9_Wu2021','M2')]
 # Retain actual cohort labels from frozen parent table; no invented/missing group.
 pd.DataFrame(tr).to_csv(F/'TRANSPORTABILITY_FIGURE_DATA.csv',index=False)
 fig,axes=plt.subplots(1,2,figsize=(14,4));groups=list(tr.groupby(['cohort','method'],sort=False))
 for ax,metric in zip(axes,['median_pearson','mse']):
  labels=[]
  for j,((cohort,m),g) in enumerate(groups):ax.scatter(np.repeat(j,len(g)),g[metric],s=16);labels.append(cohort.replace('development_','dev_').replace('HER2ST_','HER2_')+'\n'+m)
  ax.set_xticks(range(len(groups)),labels,rotation=65,ha='right',fontsize=6);ax.set_ylabel(metric);ax.grid(alpha=.2)
 fig.suptitle('Same40 genes, common error units; cohorts/roles differ and are not pooled');save(fig,'10_TRANSPORTABILITY')
 print('FIGURES',len(list(F.glob('*.png'))),flush=True)
if __name__=='__main__':main()
