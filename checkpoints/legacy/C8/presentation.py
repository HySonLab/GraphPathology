"""Patient-level descriptive summaries and figures; no fitting or selection."""
from contract import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':120,'savefig.dpi':150})
COL={'T0':'#767b83','T1':'#157ca8','T2':'#ca8127','F0':'#727881','F1':'#137da6','F2':'#bc8b25','F3':'#8660a5','SELECTED':'#d45828'}
def save(fig,path):fig.savefig(path,bbox_inches='tight',facecolor='white');plt.close(fig)
def main():
 a=pd.read_csv(A/'C8A_PATIENT_METRICS.csv');sel=pd.read_csv(A/'C8A_SELECTED_PATIENT_METRICS.csv');b=pd.read_csv(B/'C8B_PATIENT_METRICS.csv') if (B/'C8B_PATIENT_METRICS.csv').exists() else None
 a.groupby(['transform','model','regime','budget'],as_index=False)[['median_pearson','common_mse','common_mae','mean_mse','mean_target_variance']].median().to_csv(A/'C8A_COHORT_MEDIANS.csv',index=False)
 fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
 for i,regime in enumerate(['V1','V2']):
  for j,(met,label) in enumerate([('median_pearson','Median gene Pearson'),('common_mse','MSE in unchanged C7 units')]):
   ax=axes[i,j]
   for k,kind in enumerate(['T0','T1','T2']):
    vals=[]
    for model in ['E0','E1','E2']:vals.append(a[(a['transform']==kind)&(a.model==model)&(a.regime==regime)&(a.budget==.2)][met].median())
    ax.plot(np.arange(3)+(k-1)*.09,vals,'o-',label=kind,color=COL[kind])
   ax.set(title=regime,xticks=range(3),xticklabels=['E0 H&E','E1 ST','E2 fixed fusion'],ylabel=label)
   if met=='common_mse':ax.set_yscale('log')
 axes[0,0].legend();fig.suptitle('C8A at 20% anchors: median of 8 held-out patient summaries\nT0 reproduces frozen C7; target transforms selected using training patients only');save(fig,A/'figures/target_transforms.png')
 fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
 for i,regime in enumerate(['V1','V2']):
  for j,(met,label) in enumerate([('median_pearson','Median gene Pearson'),('common_mse','Common-unit MSE')]):
   ax=axes[i,j]
   for kind in ['T0','T1','T2']:
    s=a[(a['transform']==kind)&(a.model=='E2')&(a.regime==regime)&(a.budget==.2)].sort_values('sample_id');ax.plot(range(8),s[met],'o-',label=kind,color=COL[kind],lw=1)
   ax.set(title=regime,xticks=range(8),xticklabels=list('ABCDEFGH'),xlabel='Held-out patient',ylabel=label)
   if met=='common_mse':ax.set_yscale('log')
 axes[0,0].legend();fig.suptitle('C8A fixed fusion by patient, 20% anchors\nNo pooled-spot inference; common-unit errors prevent scale-only improvements');save(fig,A/'figures/patient_transforms.png')
 fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
 for j,(met,label) in enumerate([('median_pearson','Median gene Pearson'),('common_mse','Common-unit MSE')]):
  for kind in ['T0','T1','T2']:
   ss=a[(a['transform']==kind)&(a.model=='E2')&(a.regime=='V1')].groupby('budget')[met].median();axes[j].plot(ss.index*100,ss.values,'o-',label=kind,color=COL[kind])
  axes[j].set(xlabel='Observed tissue (%)',ylabel=label,xticks=[10,20,40])
  if j:axes[j].set_yscale('log')
 axes[0].legend();fig.suptitle('C8A V1 coverage response: patient medians');save(fig,A/'figures/target_coverage.png')
 variance=a[(a.model=='E1')&(a.budget==.2)][['fold','sample_id','transform','regime','target_mean','mean_target_variance']];variance.to_csv(A/'C8A_TARGET_SCALE_CONSISTENCY.csv',index=False)
 # Development values re-use saved same-40-gene C5 summaries from C7; no parent rerun.
 old=pd.read_csv(ROOT/'C7/C7_DEVELOPMENT_EXTERNAL_COMPARISON.csv');old=old[(old.dataset=='development_C5_primary')&(old.panel=='S2-50')&(old.model.isin(['E0','E1','E2']))];dev=old.rename(columns={'model':'fusion_method','budget':'coverage'}).copy();dev['target_transform']='T0';dev['target_units']='common_C7';dev['gain_vs_ST_interp']=np.nan;dev['gain_vs_H&E']=np.nan
 for _,inds in dev.groupby(['patient','coverage','mask_regime']).groups.items():
  sub=dev.loc[inds].set_index('fusion_method')
  for ix in inds:dev.loc[ix,'gain_vs_ST_interp']=dev.loc[ix,'median_pearson']-sub.loc['E1','median_pearson'];dev.loc[ix,'gain_vs_H&E']=dev.loc[ix,'median_pearson']-sub.loc['E0','median_pearson']
 cols=['dataset','patient','target_transform','fusion_method','coverage','mask_regime','median_pearson','mse','mae','gain_vs_ST_interp','gain_vs_H&E','target_units'];comps=[dev[cols]]
 for frame,label in [(a,'HER2ST_C8A'),(b,'HER2ST_C8B')]:
  if frame is None:continue
  ext=frame.rename(columns={'sample_id':'patient','transform':'target_transform','model':'fusion_method','budget':'coverage','regime':'mask_regime','common_mse':'mse','common_mae':'mae'}).copy();ext['dataset']=label;ext['target_units']='common_C7'
  for _,inds in ext.groupby(['fold','patient','target_transform','coverage','mask_regime']).groups.items():
   sub=ext.loc[inds].set_index('fusion_method');ref='E1' if label.endswith('A') else 'F0'
   baseim=sub.loc['E0','median_pearson'] if ref=='E1' else float(sel[(sel.sample_id==ext.loc[list(inds)[0],'patient'])&(sel.model=='E0')&(sel.budget==ext.loc[list(inds)[0],'coverage'])&(sel.regime==ext.loc[list(inds)[0],'mask_regime'])].median_pearson.iloc[0])
   for ix in inds:ext.loc[ix,'gain_vs_ST_interp']=ext.loc[ix,'median_pearson']-sub.loc[ref,'median_pearson'];ext.loc[ix,'gain_vs_H&E']=ext.loc[ix,'median_pearson']-baseim
  comps.append(ext[cols])
 comparison=pd.concat(comps,ignore_index=True);comparison.to_csv(OUT/'C8_DEVELOPMENT_TRANSPORTABILITY.csv',index=False)
 if b is None:return
 b.groupby(['model','regime','budget'],as_index=False)[['median_pearson','common_mse','common_mae','mean_mse']].median().to_csv(B/'C8B_COHORT_MEDIANS.csv',index=False)
 for coverage in [False,True]:
  fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
  for i,regime in enumerate(['V1','V2']):
   for j,(met,label) in enumerate([('median_pearson','Median gene Pearson'),('common_mse','Common-unit MSE')]):
    ax=axes[i,j]
    if coverage:
     for model in ['F0','F1','F2','F3','SELECTED']:
      ss=b[(b.model==model)&(b.regime==regime)].groupby('budget')[met];m=ss.median();ax.plot(m.index*100,m,'o-',label=model,color=COL[model]);ax.fill_between(m.index*100,ss.quantile(.25),ss.quantile(.75),color=COL[model],alpha=.07)
     ax.set(xlabel='Observed tissue (%)',xticks=[10,20,40])
    else:
     for k,model in enumerate(['F0','F1','F2','F3','SELECTED']):
      vals=b[(b.model==model)&(b.regime==regime)&(b.budget==.2)][met];ax.scatter(k+np.linspace(-.08,.08,len(vals)),vals,color=COL[model],alpha=.6,s=20);ax.plot([k-.15,k+.15],[vals.median()]*2,color=COL[model],lw=3)
     ax.set(xticks=range(5),xticklabels=['F0 ST','F1 scalar','F2 distance','F3 confidence','Selected'])
    ax.set(title=regime,ylabel=label)
    if met=='common_mse':ax.set_yscale('log')
  if coverage:axes[0,0].legend(ncol=3,fontsize=9)
  fig.suptitle('C8B: '+('coverage response; bands = patient IQR, not confidence intervals\nV2 is evaluated at 20% only unless C8B passes' if coverage else '20% anchors; each dot is a held-out patient\nSELECTED is chosen within each training fold, not from these results'));save(fig,B/'figures'/('coverage_response.png' if coverage else 'fusion_strategies.png'))
 eff=pd.read_csv(B/'C8B_PATIENT_EFFECTS.csv');fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
 for i,regime in enumerate(['V1','V2']):
  ss=eff[(eff.regime==regime)&(eff.budget==.2)].sort_values('sample_id')
  for j,(met,label) in enumerate([('gain_pearson','Paired Pearson gain'),('relative_mse_reduction','Relative common-MSE reduction (%)')]):
   ax=axes[i,j];ax.bar(range(8),ss[met]*(100 if j else 1),color=COL['SELECTED']);ax.axhline(0,color='black',lw=.8);ax.set(title=regime,xticks=range(8),xticklabels=list('ABCDEFGH'),xlabel='Held-out patient',ylabel=label)
 fig.suptitle('Selected calibrated fusion versus ST interpolation at 20%\nPositive means improvement; patients remain the biological units');save(fig,B/'figures/patient_fusion_effects.png')
 # Gene effects, paired within each patient/mask then summarized by patient.
 rows=[json.loads(p.read_text()) for p in (B/'logs/runs').glob('*.json')];clean=[r for r in rows if r['control'] is None];lookup={(r['mask_id'],r['budget'],r['model']):r for r in clean};gg=[];dd=[]
 for r in clean:
  if r['model']!='SELECTED':continue
  base=lookup[(r['mask_id'],r['budget'],'F0')]
  with np.load(ROOT/r['prediction_path']) as f,np.load(ROOT/base['prediction_path']) as z:
   if r['budget']==.2:
    for j,g in enumerate(GENES):gg.append({'patient':r['sample_id'],'regime':r['regime'],'gene':g,'pearson_gain':float(f['pearson'][j]-z['pearson'][j]),'common_mse_reduction':float(z['common_gene_mse'][j]-f['common_gene_mse'][j])})
  bins={x['lower_um']:x for x in base['distance_bins']}
  for x in r['distance_bins']:dd.append({'patient':r['sample_id'],'regime':r['regime'],'budget':r['budget'],'distance_bin':x['lower_um'],'n':x['n'],'fusion_sse':x['common_mse']*x['n'],'ST_sse':bins[x['lower_um']]['common_mse']*x['n']})
 gp=pd.DataFrame(gg).groupby(['patient','regime','gene'],as_index=False)[['pearson_gain','common_mse_reduction']].mean();gp.to_csv(B/'C8B_GENE_PATIENT_EFFECTS.csv',index=False);gs=gp.groupby(['regime','gene']).agg(pearson_gain=('pearson_gain','median'),common_mse_reduction=('common_mse_reduction','median'),patients_positive_pearson=('pearson_gain',lambda x:int((x>0).sum())),patients_positive_MSE=('common_mse_reduction',lambda x:int((x>0).sum()))).reset_index();gs.to_csv(B/'C8B_GENE_GAIN_SUMMARY.csv',index=False)
 fig,axes=plt.subplots(2,1,figsize=(14,8),layout='constrained')
 for ax,regime in zip(axes,['V1','V2']):
  sub=gp[gp.regime==regime].pivot(index='patient',columns='gene',values='pearson_gain').reindex(columns=GENES);lim=max(.001,np.nanmax(abs(sub.values)));im=ax.imshow(sub,cmap='coolwarm',vmin=-lim,vmax=lim,aspect='auto');ax.set(yticks=range(8),yticklabels=list('ABCDEFGH'),xticks=range(40),xticklabels=GENES,title=regime,ylabel='Held-out patient');ax.tick_params(axis='x',labelrotation=90,labelsize=8);fig.colorbar(im,ax=ax,label='Fusion − ST Pearson')
 fig.suptitle('All 40 frozen genes: selected-fusion effects at 20%\nOriginal gene order; no outcome-based gene selection');save(fig,B/'figures/gene_fusion_effects.png')
 ds=pd.DataFrame(dd).groupby(['patient','regime','budget','distance_bin'],as_index=False)[['n','fusion_sse','ST_sse']].sum();ds['fusion_mse']=ds.fusion_sse/ds.n;ds['ST_mse']=ds.ST_sse/ds.n;ds['mse_reduction']=ds.ST_mse-ds.fusion_mse;ds.to_csv(B/'C8B_DISTANCE_PATIENT_EFFECTS.csv',index=False)
 fig,axes=plt.subplots(1,2,figsize=(12,4),layout='constrained')
 for ax,regime in zip(axes,['V1','V2']):
  ss=ds[(ds.regime==regime)&(ds.budget==.2)].groupby('distance_bin').agg(gain=('mse_reduction','median'),patients=('patient','nunique'));x=[0,150,300,600,1200,2400];ss=ss.reindex(x);ax.plot(range(6),ss.gain,'o-',color=COL['SELECTED']);ax.axhline(0,color='black',lw=.7);ax.set(title=regime,xticks=range(6),xticklabels=['0–150','150–300','300–600','600–1200','1200–2400','≥2400'],ylabel='Median paired common-MSE reduction',xlabel='Nearest anchor (µm)');ax.tick_params(axis='x',rotation=25)
 fig.suptitle('Selected fusion versus ST by distance at 20%\nRepeated masks weighted within patient/bin, then patient median; sparse bins are descriptive');save(fig,B/'figures/fusion_distance_gain.png')
 # Negative controls paired to the identical clean run.
 nc=[]
 for r in rows:
  if not r['control']:continue
  base=lookup[(r['mask_id'],r['budget'],r['model'])];nc.append({k:r[k] for k in ['sample_id','model','regime','control','control_seed']}|{'Pearson_change':r['median_pearson']-base['median_pearson'],'common_MSE_change':r['common_mse']-base['common_mse']})
 nc=pd.DataFrame(nc);nc.to_csv(B/'C8B_NEGATIVE_CONTROLS.csv',index=False);nc.groupby(['sample_id','model','regime','control'],as_index=False)[['Pearson_change','common_MSE_change']].mean().to_csv(B/'C8B_NEGATIVE_CONTROL_PATIENTS.csv',index=False)
 fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
 for ax,met,label in zip(axes,['median_pearson','mse'],['Median gene Pearson','MSE in C7 units']):
  cohorts=[('C5 source','development_C5_primary','T0','E2'),('C7 zero-shot','HER2ST_C8A','T0','E2'),('C8 harmonized','HER2ST_C8A','T1','E2'),('C8 calibrated','HER2ST_C8B','T1','SELECTED'),('ST only','HER2ST_C8B','T1','F0')]
  for i,(name,dsname,tr,mod) in enumerate(cohorts):
   vals=comparison[(comparison.dataset==dsname)&(comparison.target_transform==tr)&(comparison.fusion_method==mod)&(comparison.coverage==.2)&(comparison.mask_regime=='V1')][met];ax.scatter(i+np.linspace(-.06,.06,len(vals)),vals,color='#246d88',alpha=.7);ax.plot([i-.12,i+.12],[vals.median()]*2,color='black',lw=2)
  ax.set(xticks=range(5),xticklabels=[r[0] for r in cohorts],ylabel=label);ax.tick_params(axis='x',rotation=15)
  if met=='mse':ax.set_yscale('log')
 fig.suptitle('20% random anchors, identical 40 genes\nC5 has one primary test patient; HER2ST has eight development-CV patients');save(fig,B/'figures/domain_gap.png')
 print('Presentation complete',len(list(OUT.glob('*/figures/*.png'))),'figures')
if __name__=='__main__':main()
