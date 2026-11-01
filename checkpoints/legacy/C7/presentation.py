"""Descriptive tables and figures from locked outcomes; no model or selection changes."""
from support import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
Image.MAX_IMAGE_PIXELS=120_000_000
F=B/'figures';F.mkdir(exist_ok=True)
COLORS={'E0':'#6b7280','E1':'#1674a5','E2':'#d66b17','E3':'#8360a5'}
NAMES={'E0':'E0 H&E','E1':'E1 ST','E2':'E2 fusion','E3':'E3 propagation'}
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':120,'savefig.dpi':160})
def save(fig,name):
 fig.savefig(F/name,bbox_inches='tight',facecolor='white');plt.close(fig)
def main():
 c=pd.read_csv(B/'C7B_COVERAGE_CURVES.csv');e=pd.read_csv(B/'C7B_PATIENT_EFFECTS.csv');d=pd.read_csv(B/'C7B_DISTANCE_TO_ANCHOR.csv');rows=[json.loads(p.read_text()) for p in sorted((B/'logs/runs').glob('*.json'))];clean=[r for r in rows if not r['control']]
 summary=c.groupby(['panel','regime','strategy','nominal_budget','model'])[['median_pearson','mean_mse','mean_mae']].median().reset_index();summary.to_csv(B/'C7B_COHORT_MEDIANS.csv',index=False)
 erows=[]
 for key,g in e.groupby(['panel','regime','strategy','nominal_budget']):
  z=dict(zip(['panel','regime','strategy','nominal_budget'],key))
  for base in ['E0','E1']:
   for metric in ['gain_pearson','gain_mse','gain_mae','relative_mse_reduction']:
    v=g[metric+'_vs_'+base];z[metric+'_vs_'+base+'_median']=v.median();z[metric+'_vs_'+base+'_positive_patients']=int((v>0).sum())
   z['meaningful_vs_'+base+'_patients']=int(g['qualifies_vs_'+base].sum())
  erows.append(z)
 pd.DataFrame(erows).to_csv(B/'C7B_COHORT_EFFECTS.csv',index=False)
 # Equal-patient distance summaries. Within each patient/bin, weight sampled locations.
 d['sse']=d.mse*d.n_hidden;d['sae']=d.mae*d.n_hidden
 dk=['patient_id','panel','regime','strategy','model','lower_um']
 dp=d.groupby(dk,as_index=False)[['sse','sae','n_hidden']].sum();dp['mse']=dp.sse/dp.n_hidden;dp['mae']=dp.sae/dp.n_hidden;dp.to_csv(B/'C7B_DISTANCE_PATIENT_SUMMARY.csv',index=False)
 dist=dp.groupby(['panel','regime','strategy','model','lower_um']).agg(median_patient_mse=('mse','median'),q25=('mse',lambda x:x.quantile(.25)),q75=('mse',lambda x:x.quantile(.75)),patients=('patient_id','nunique'),sampled_locations=('n_hidden','sum')).reset_index();dist.to_csv(B/'C7B_DISTANCE_COHORT_SUMMARY.csv',index=False)
 # Independent distance reconstruction from every saved mask/prediction file.
 n=0
 for r in rows:
  with np.load(ROOT/r['mask_path']) as m,np.load(ROOT/r['prediction_path']) as f:
   actual=cKDTree(m['xy_um'][f['observed_indices']]).query(m['xy_um'][f['hidden_indices']])[0];assert np.array_equal(actual,f['distance_um'])
   for bb in r['distance_bins']:
    selected=(actual>=bb['lower_um'])&(actual<(bb['upper_um'] if bb['upper_um'] is not None else np.inf));assert selected.sum()==bb['n_hidden']
   n+=1
 dump(B/'logs/DISTANCE_AUDIT.json',{'status':'PASS','prediction_files_checked':n,'exact_distance_and_bin_count_reconstruction':True})
 # Descriptive gene-specific gains at the predeclared 20% budget, random policy.
 grows=[];lookup={(r['mask_id'],r['panel'],r['nominal_budget'],r['model']):r for r in clean}
 for r in clean:
  if r['model']!='E2' or r['strategy']!='random' or r['nominal_budget']!=.2:continue
  with np.load(ROOT/r['prediction_path']) as f:
   pred={k:f[k].copy() for k in ['genes','pearson','mse','mae']}
  for base in ['E0','E1']:
   with np.load(ROOT/lookup[(r['mask_id'],r['panel'],.2,base)]['prediction_path']) as f:
    for j,gene in enumerate(pred['genes']):grows.append({'patient_id':r['patient_id'],'panel':r['panel'],'regime':r['regime'],'gene':gene,'comparator':base,'pearson_gain':float(pred['pearson'][j]-f['pearson'][j]),'mse_reduction':float(f['mse'][j]-pred['mse'][j]),'mae_reduction':float(f['mae'][j]-pred['mae'][j])})
 gp=pd.DataFrame(grows).groupby(['patient_id','panel','regime','gene','comparator'],as_index=False)[['pearson_gain','mse_reduction','mae_reduction']].mean();gp.to_csv(B/'C7B_GENE_GAIN_PATIENT.csv',index=False)
 gg=gp.groupby(['panel','regime','gene','comparator']).agg(median_patient_pearson_gain=('pearson_gain','median'),median_patient_mse_reduction=('mse_reduction','median'),median_patient_mae_reduction=('mae_reduction','median'),patients_positive_pearson=('pearson_gain',lambda x:int((x>0).sum())),patients_positive_mse=('mse_reduction',lambda x:int((x>0).sum())),patients_finite_pearson=('pearson_gain','count')).reset_index();gg.to_csv(B/'C7B_GENE_GAIN_SUMMARY.csv',index=False)
 for strategy in ['random','spatial_FPS']:
  for metric,ylabel,short in [('median_pearson','Patient median gene Pearson','pearson'),('mean_mse','Patient mean gene MSE (frozen standardized units)','mse')]:
   fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
   for i,panel in enumerate(panels()):
    for j,regime in enumerate(['V1','V2']):
     ax=axes[i,j];sub=c[(c.panel==panel)&(c.regime==regime)&(c.strategy==strategy)]
     for model,color in COLORS.items():
      tab=sub[sub.model==model].groupby('nominal_budget')[metric];mid=tab.median();ax.plot(mid.index*100,mid.values,'o-',color=color,label=NAMES[model]);ax.fill_between(mid.index*100,tab.quantile(.25),tab.quantile(.75),color=color,alpha=.10)
     ax.set(title=f'{panel} ({len(panels()[panel])} genes) | {regime}',xlabel='Measured locations (%)',ylabel=ylabel,xticks=[10,20,40]);ax.grid(alpha=.15)
     if short=='mse':ax.set_yscale('log')
   axes[0,0].legend(fontsize=9);fig.suptitle(f'{strategy.capitalize()} sampling: median of 8 patients; bands = patient IQR, not confidence intervals\nV1 scored locations change with budget; V2 averages three fixed unavailable regions',fontsize=12);save(fig,f'coverage_{strategy}_{short}.png')
 # Same metric, matched patient comparison at20% (V1 vsV2).
 fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
 for i,panel in enumerate(panels()):
  for j,(metric,label) in enumerate([('median_pearson','Median gene Pearson'),('mean_mse','Mean gene MSE')]):
   ax=axes[i,j]
   for mi,model in enumerate(COLORS):
    for ri,regime in enumerate(['V1','V2']):
     vals=c[(c.panel==panel)&(c.nominal_budget==.2)&(c.strategy=='random')&(c.model==model)&(c.regime==regime)].sort_values('patient_id')[metric].to_numpy();x=mi+(-.16 if ri==0 else .16);ax.scatter(x+np.linspace(-.06,.06,8),vals,c=COLORS[model],marker='o' if ri==0 else '^',s=25,alpha=.75);ax.plot([x-.1,x+.1],[np.median(vals)]*2,color='black',lw=2)
   ax.set(title=f'{panel}: circles V1 / triangles V2',xticks=range(4),xticklabels=list(NAMES.values()),ylabel=label)
   if metric=='mean_mse':ax.set_yscale('log')
 fig.suptitle('20% random anchors | each marker is one patient; bars are medians\nV2: three fixed contiguous regions per patient');save(fig,'models_and_contiguous.png')
 # Paired patient effects.
 for panel in panels():
  fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
  for i,regime in enumerate(['V1','V2']):
   sub=e[(e.panel==panel)&(e.regime==regime)&(e.strategy=='random')&(e.nominal_budget==.2)].sort_values('patient_id')
   for j,(metric,label) in enumerate([('gain_pearson','E2 − comparator: median gene Pearson'),('relative_mse_reduction','Relative MSE reduction (%)')]):
    ax=axes[i,j]
    for offset,base,col in [(-.12,'E0',COLORS['E0']),(.12,'E1',COLORS['E1'])]:ax.scatter(np.arange(8)+offset,sub[metric+'_vs_'+base]*(100 if j else 1),label='vs '+NAMES[base],color=col)
    ax.axhline(0,color='black',lw=.8);ax.set(xticks=range(8),xticklabels=[p[-1] for p in sub.patient_id],xlabel='External patient',ylabel=label,title=regime);ax.legend(fontsize=9)
  fig.suptitle(f'{panel}: paired patient effects at 20% random coverage\nPositive means improvement; no spot-level significance tests');save(fig,'patient_effects_'+panel+'.png')
 fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained');labels=['0–150','150–300','300–600','600–1200','1200–2400','≥2400'];lo=[0,150,300,600,1200,2400]
 for i,panel in enumerate(panels()):
  for j,regime in enumerate(['V1','V2']):
   ax=axes[i,j]
   for model,col in COLORS.items():
    sub=dist[(dist.panel==panel)&(dist.regime==regime)&(dist.strategy=='random')&(dist.model==model)].set_index('lower_um').reindex(lo);ax.plot(range(6),sub.median_patient_mse,'o-',label=NAMES[model],color=col)
   ax.set(title=f'{panel} | {regime}',xticks=range(6),xticklabels=labels,ylabel='Median patient MSE',xlabel='Nearest-anchor distance (µm)');ax.tick_params(axis='x',rotation=25);ax.set_yscale('log');ax.grid(alpha=.15)
 axes[0,0].legend();fig.suptitle('Distance behavior: random policy, all frozen budgets\nWithin-patient location-weighted bins, then equal-patient median; no extrapolation-horizon claim');save(fig,'distance_to_anchor.png')
 # Developmentcomparison identicalgenes, onlyV1samebroadregime.
 comp=pd.read_csv(OUT/'C7_DEVELOPMENT_EXTERNAL_COMPARISON.csv');fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
 for i,panel in enumerate(panels()):
  for j,(metric,label) in enumerate([('median_pearson','Median gene Pearson'),('mse','Mean gene MSE')]):
   ax=axes[i,j]
   for mi,model in enumerate(COLORS):
    sub=comp[(comp.panel==panel)&(comp.model==model)&(comp.budget==.2)&(comp.mask_regime=='V1')&(comp.strategy=='random')];ext=sub[sub.dataset=='external_HER2ST'][metric];dev=sub[sub.dataset=='development_C5_primary'][metric];ax.scatter(mi+np.linspace(-.09,.09,len(ext)),ext,color=COLORS[model],alpha=.7,s=25);ax.scatter([mi+.25]*len(dev),dev,color='black',marker='D',s=40);ax.plot([mi-.12,mi+.12],[ext.median()]*2,color=COLORS[model],lw=3)
   ax.set(title=f'{panel}: identical {len(panels()[panel])} retained genes',xticks=range(4),xticklabels=list(NAMES.values()),ylabel=label)
   if metric=='mse':ax.set_yscale('log')
 fig.suptitle('20% random anchors | development: 1 primary patient (black diamond)\nExternal: 8 independent patients (dots + median); cross-platform errors are not assay-calibrated');save(fig,'development_external.png')
 # Patient H&E overlays and predeclared six genes.
 for regime in ['V1','V2']:
  fig,axes=plt.subplots(2,4,figsize=(15,8),layout='constrained')
  for ax,sample in zip(axes.flat,SAMPLES):
   r=next(r for r in clean if r['sample_id']==sample and r['panel']=='S2-50' and r['regime']==regime and r['geometry']==(0 if regime=='V1' else 1) and r['seed']==0 and r['strategy']=='random' and r['nominal_budget']==.2 and r['model']=='E2')
   with np.load(CACHE/(sample+'_metadata.npz')) as m:xy=m['xy_px'].copy()
   with np.load(ROOT/r['prediction_path']) as p:oi=p['observed_indices'];hi=p['hidden_indices']
   img=Image.open(CACHE/(sample+'.jpg'));size=img.size;img.thumbnail((800,800));ax.imshow(img,extent=(0,size[0],size[1],0));ax.scatter(xy[oi,0],xy[oi,1],s=8,facecolor='#ffdb00',edgecolor='black',lw=.3)
   if regime=='V2':ax.scatter(xy[hi,0],xy[hi,1],s=7,facecolors='none',edgecolors='#00d6ff',lw=.4)
   ax.set_title(f'{sample}: {len(oi)} / {len(xy)} anchors',fontsize=11);ax.axis('off')
  fig.suptitle(f'HER2ST H&E | {regime}, 20% random seed 0\nYellow: measured anchors; cyan outlines (V2): fixed missing region',fontsize=14);save(fig,f'he_anchors_{regime}.png')
  r=next(r for r in clean if r['sample_id']=='A1' and r['panel']=='S2-50' and r['regime']==regime and r['geometry']==(0 if regime=='V1' else 1) and r['seed']==0 and r['strategy']=='random' and r['nominal_budget']==.2 and r['model']=='E2')
  with np.load(ROOT/r['prediction_path']) as p:pr=p['prediction'];h=p['hidden_indices'];genes=list(p['genes'])
  xy=scene('A1','S2-50')['xy'][h];y=truth('A1','S2-50')[h];display=cfg()['figure_genes'];fig,axes=plt.subplots(len(display),3,figsize=(11,17),layout='constrained')
  for i,g in enumerate(display):
   k=genes.index(g);a=y[:,k];b=pr[:,k];err=b-a;vmin=min(a.min(),b.min());vmax=max(a.max(),b.max());lim=max(abs(err).max(),1e-9)
   for j,(v,title,cmap) in enumerate([(a,'Ground truth','viridis'),(b,'Frozen E2 prediction','viridis'),(err,'Prediction − truth','coolwarm')]):
    ax=axes[i,j];im=ax.scatter(xy[:,0],xy[:,1],c=v,s=13,cmap=cmap,vmin=-lim if j==2 else vmin,vmax=lim if j==2 else vmax);ax.set_aspect('equal');ax.invert_yaxis();ax.set_title(f'{g} | {title}',fontsize=10);ax.set(xticks=[],yticks=[]);fig.colorbar(im,ax=ax,shrink=.8)
  fig.suptitle(f'A1 | {regime}, 20% random seed 0 | hidden locations only\nSix genes fixed before evaluation; values in frozen C0 standardized log-count units',fontsize=13);save(fig,f'gene_maps_{regime}.png')
 dump(B/'logs/FIGURE_INDEX.json',{'figures':[str(p.relative_to(ROOT)) for p in sorted(F.glob('*.png'))],'display_gene_selection':'First six retained S2-50 genes, frozen in protocol before external access','statistics':'Patients are units; seeds/geometries averaged within patient; bands are IQR, not confidence intervals.'})
 print('Presentation complete',len(list(F.glob('*.png'))),'figures',flush=True)
if __name__=='__main__':main()
