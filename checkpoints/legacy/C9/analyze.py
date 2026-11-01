"""Patient-level, fixed-method summaries and prospective decision-rule evaluation."""
from contract import *

def main():
 check_lock();records=[json.loads(p.read_text()) for p in sorted((B/'logs/runs').glob('*.json'))];assert len(records)==414
 panels=[];gg=[];dd=[];comparables=[];common=parent.panels()['S2-50'];idx=[genes().index(g) for g in common]
 for r in records:
  base={k:r[k] for k in ['run_id','sample_id','patient_id','regime','geometry','strategy','seed','budget','control','control_seed']}
  with np.load(ROOT/r['prediction_path']) as f:
   for m in ['M0','M1','M2','M3']:
    panels.append({**base,'method':m,**r['summary'][m]})
    for j,g in enumerate(genes()):gg.append({**base,'method':m,'gene':g,**{k:float(f[m+'_'+k][j]) for k in ['pearson','spearman','mse','mae']}})
    if not r['control']:comparables.append({**base,'method':m,'median_pearson':core.finite_stat(f[m+'_pearson'][idx],'median'),'mean_mse':float(np.mean(f[m+'_mse'][idx])),'mean_mae':float(np.mean(f[m+'_mae'][idx]))})
  dd.extend({**base,**d} for d in r['distance_bins'])
 df=pd.DataFrame(panels);df.to_csv(B/'C9B_PANEL_METRICS.csv',index=False);pd.DataFrame(gg).to_csv(B/'C9B_GENE_METRICS.csv',index=False);pd.DataFrame(dd).to_csv(B/'C9B_DISTANCE_TO_ANCHOR.csv',index=False)
 clean=df[df.control==''];metrics=list(records[0]['summary']['M0']);keys=['patient_id','regime','strategy','budget','method']
 geo=clean.groupby(keys+['geometry'])[metrics].mean().reset_index();pt=geo.groupby(keys)[metrics].mean().reset_index();pt.to_csv(B/'C9B_PATIENT_METRICS.csv',index=False);geo.to_csv(B/'C9B_GEOMETRY_METRICS.csv',index=False);pt.to_csv(B/'C9B_COVERAGE_CURVES.csv',index=False)
 effects=[]
 for key,g in pt.groupby(['patient_id','regime','strategy','budget']):
  g=g.set_index('method');base=dict(zip(['patient_id','regime','strategy','budget'],key))
  for a,b,label in [('M2','M0','sparse_ST'),('M2','M1','harmonization'),('M3','M2','image')]:
   dr=float(g.loc[a,'median_pearson']-g.loc[b,'median_pearson']);mse=float(g.loc[b,'mean_mse']-g.loc[a,'mean_mse']);mae=float(g.loc[b,'mean_mae']-g.loc[a,'mean_mae']);rel=mse/g.loc[b,'mean_mse']
   effects.append({**base,'comparison':label,'delta_pearson':dr,'mse_reduction':mse,'relative_mse_reduction':rel,'mae_reduction':mae,'meaningful_sparse':dr>=.02 and rel>=.02 and mae>=0})
 ef=pd.DataFrame(effects);ef.to_csv(B/'C9B_PATIENT_EFFECTS.csv',index=False)
 nc=df[df.control!=''].copy();matched=clean[['sample_id','regime','geometry','strategy','seed','budget','method']+metrics].rename(columns={k:'clean_'+k for k in metrics});nc=nc.merge(matched,on=['sample_id','regime','geometry','strategy','seed','budget','method'],validate='many_to_one');nc['delta_pearson']=nc.median_pearson-nc.clean_median_pearson;nc['mse_reduction']=nc.clean_mean_mse-nc.mean_mse;nc['mae_reduction']=nc.clean_mean_mae-nc.mean_mae;nc.to_csv(B/'C9B_NEGATIVE_CONTROLS.csv',index=False)
 primary=ef[(ef.regime=='V1')&(ef.strategy=='random')];st=primary[primary.comparison=='sparse_ST'];support=st.groupby('patient_id').meaningful_sparse.sum();coverage=[]
 for s in EVAL:
  p=pt[(pt.patient_id==s)&(pt.regime=='V1')&(pt.strategy=='random')&(pt.method=='M2')].set_index('budget');coverage.append({'patient_id':s,'delta_pearson_40_minus_10':float(p.loc[.4,'median_pearson']-p.loc[.1,'median_pearson']),'mse_reduction_10_to_40':float(p.loc[.1,'mean_mse']-p.loc[.4,'mean_mse'])})
 cov=pd.DataFrame(coverage);cov.to_csv(B/'C9B_COVERAGE_EFFECTS.csv',index=False);budget_gains=st.groupby('budget')[['delta_pearson','mse_reduction','mae_reduction']].median();core_ok=bool((support>=2).sum()>=2 and (budget_gains>0).all().all() and ((cov.delta_pearson_40_minus_10>0)&(cov.mse_reduction_10_to_40>0)).sum()>=2)
 some_support=int(st[st.meaningful_sparse].patient_id.nunique())>=2
 image=primary[primary.comparison=='image'];imok=(image.delta_pearson>0)&(image.mse_reduction>0)&(image.mae_reduction>0);image_ok=bool(imok.groupby(image.patient_id).sum().ge(2).sum()>=2 and (image.relative_mse_reduction<-.02).mean()<=.25)
 status='MIXED_FRESH_EXTERNAL_VALIDATION' if core_ok or some_support else 'FAIL_FRESH_EXTERNAL_VALIDATION'
 dump(B/'logs/GATE_RESULT.json',{'C9B_status':status,'core_sparse_support':core_ok,'patient_meaningful_budget_counts':support.to_dict(),'coverage_response_patients':int(((cov.delta_pearson_40_minus_10>0)&(cov.mse_reduction_10_to_40>0)).sum()),'harmonization_gain_identifiable':False,'harmonization_pass':False,'reason':'M1/M2 affine invariance; cannot satisfy full-PASS clause3. Core, image and platform limitations reported independently.','fusion_secondary_pass':image_ok,'max_M1_M2_prediction_difference':max(r['affine_max_abs_difference'] for r in records),'runs':len(records),'method_predictions':len(records)*4,'corruption_invariance_checks':sum(r['hidden_corruption_checks'] for r in records),'recommendation':'PROCEED_TO_SECOND_REPLICATION_COHORT' if core_ok else 'NEED_METHOD_REVIEW'})
 # Comparable 40-gene primary subset, matching frozen C7/C8 exactly.
 cp=pd.DataFrame(comparables).groupby(keys+['geometry'])[['median_pearson','mean_mse','mean_mae']].mean().groupby(keys).mean().reset_index();trans=[]
 for _,r in cp[cp.strategy=='random'].iterrows():
  peers=cp[(cp.patient_id==r.patient_id)&(cp.regime==r.regime)&(cp.strategy==r.strategy)&(cp.budget==r.budget)].set_index('method');trans.append({'cohort':'fresh_C9_Wu2021','patient':r.patient_id,'role':'held_out_evaluation','panel':'C7/C8 common S2-50 40 genes','coverage':r.budget,'mask_regime':r.regime,'method':r.method,'median_pearson':r.median_pearson,'mse':r.mean_mse,'mae':r.mean_mae,'gain_vs_H&E':r.median_pearson-peers.loc['M0','median_pearson'],'gain_vs_raw_interpolation':r.median_pearson-peers.loc['M1','median_pearson'],'gain_vs_harmonized_interpolation':r.median_pearson-peers.loc['M2','median_pearson'],'target_units':'common_C0_normalized_log1p_counts'})
 old=pd.read_csv(ROOT/'C8/C8_DEVELOPMENT_TRANSPORTABILITY.csv')
 for _,r in old.iterrows():trans.append({'cohort':r.dataset,'patient':r.patient,'role':'development_or_transportability_development','panel':'C7/C8 common S2-50 40 genes','coverage':r.coverage,'mask_regime':r.mask_regime,'method':r.target_transform+'_'+r.fusion_method,'median_pearson':r.median_pearson,'mse':r.mse,'mae':r.mae,'gain_vs_H&E':r['gain_vs_H&E'],'gain_vs_raw_interpolation':r.gain_vs_ST_interp,'gain_vs_harmonized_interpolation':r.gain_vs_ST_interp,'target_units':'common_C0_normalized_log1p_counts'})
 pd.DataFrame(trans).to_csv(OUT/'C9_TRANSPORTABILITY_SUMMARY.csv',index=False)
 print(json.dumps(json.loads((B/'logs/GATE_RESULT.json').read_text()),indent=2))
if __name__=='__main__':main()
