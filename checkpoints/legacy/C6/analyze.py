"""Frozen fixed-budget comparisons, random distributions, matched-bin contrasts and gates."""
from support import *
from scipy.integrate import trapezoid

def main():
 records=[json.loads(p.read_text()) for p in sorted((A/'logs/runs').glob('*.json'))];records=[r for r in records if r['phase']!='validation'];assert len(records)==4448
 table=pd.DataFrame([{k:v for k,v in r.items() if not isinstance(v,(dict,list))} for r in records]);table.to_csv(A/'C6A_PANEL_METRICS.csv',index=False);table[table.phase=='sensitivity'].to_csv(A/'C6A_PATIENT_SENSITIVITY.csv',index=False);table[table.control.notna()].to_csv(A/'C6A_NEGATIVE_CONTROLS.csv',index=False)
 distance=pd.DataFrame([{**{k:r[k] for k in ['run_id','policy_id','cohort','sample_id','encoder','panel','role','regime','strategy','seed','control','nominal_budget','realized_budget']},**b} for r in records for b in r['distance_bins']]);distance.to_csv(A/'C6A_DISTANCE_METRICS.csv',index=False)
 geometry_cols=['mean_anchor_distance_um','median_anchor_distance_um','p90_anchor_distance_um','max_anchor_distance_um','spatial_coverage_radius_um','unmeasured_components_at_least_50','morphology_mean_nearest_distance','morphology_p90_nearest_distance','observed_morphology_mean_pairwise_squared_distance'];table.drop_duplicates(['policy_id','nominal_budget'])[['policy_id','cohort','sample_id','encoder','regime','strategy','seed','control','nominal_budget','realized_budget','n_observed','n_hidden']+geometry_cols].to_csv(A/'C6A_GEOMETRY_DIAGNOSTICS.csv',index=False)
 ordinary=table[table.control.isna()];keys=['cohort','sample_id','encoder','panel','role','regime','nominal_budget'];summary=ordinary.groupby(keys+['strategy']).agg(realized_budget=('realized_budget','mean'),median_pearson=('median_pearson','mean'),pearson_sd=('median_pearson','std'),mean_pearson=('mean_pearson','mean'),mean_mse=('mean_mse','mean'),mse_sd=('mean_mse','std'),mean_mae=('mean_mae','mean'),n_runs=('run_id','size')).reset_index();summary.to_csv(A/'C6A_BUDGET_CURVES.csv',index=False);ordinary[ordinary.strategy=='A0'].to_csv(A/'C6A_RANDOM_DISTRIBUTION.csv',index=False)
 effects=[]
 for key,d in ordinary.groupby(keys):
  random=d[d.strategy=='A0'];assert len(random)==10
  for strategy in ['A1','A2','A3']:
   z=d[d.strategy==strategy];assert len(z)==1;r=z.iloc[0];rg=float(r.median_pearson-random.median_pearson.mean());mr=float(random.mean_mse.mean()-r.mean_mse);rel=mr/float(random.mean_mse.mean());wins=int(((r.median_pearson>random.median_pearson)&(r.mean_mse<random.mean_mse)).sum());threshold=max(.01,float(random.median_pearson.std()));qual=rg>=threshold and rel>=.02 and wins>=8
   effects.append({**dict(zip(keys,key)),'strategy':strategy,'Pearson_gain':rg,'MSE_reduction':mr,'relative_MSE_reduction':rel,'random_Pearson_SD':float(random.median_pearson.std()),'random_MSE_SD':float(random.mean_mse.std()),'wins_on_both_of_10':wins,'qualifies':bool(qual),'material_change':bool(abs(rg)>=threshold and abs(rel)>=.02)})
 effects=pd.DataFrame(effects);effects.to_csv(A/'C6A_EFFECTS.csv',index=False)
 auc=[]
 for key,d in ordinary.groupby(['cohort','sample_id','encoder','panel','role','regime','strategy','seed']):
  d=d.sort_values('realized_budget');assert len(d)==4;x=d.realized_budget.to_numpy();auc.append({**dict(zip(['cohort','sample_id','encoder','panel','role','regime','strategy','seed'],key)),'Pearson_AUC_coverage':float(trapezoid(d.median_pearson,x)/(x[-1]-x[0])),'MSE_AUC_coverage':float(trapezoid(d.mean_mse,x)/(x[-1]-x[0])),'MAE_AUC_coverage':float(trapezoid(d.mean_mae,x)/(x[-1]-x[0]))})
 pd.DataFrame(auc).to_csv(A/'C6A_AUC_COVERAGE.csv',index=False)
 matched=[];dk=['cohort','sample_id','encoder','panel','role','regime','nominal_budget','lower_um']
 for key,d in distance[distance.control.isna()].groupby(dk):
  r=d[d.strategy=='A0']
  if not len(r):continue
  for strategy in ['A1','A2','A3']:
   z=d[d.strategy==strategy]
   if len(z):q=z.iloc[0];matched.append({**dict(zip(dk,key)),'strategy':strategy,'n_random_masks_with_bin':len(r),'strategy_n_hidden':q.n_hidden,'random_mean_n_hidden':r.n_hidden.mean(),'MSE_reduction_at_distance':r.mse.mean()-q.mse,'MAE_reduction_at_distance':r.mae.mean()-q.mae,'Pearson_gain_at_distance':q.median_pearson-r.median_pearson.mean()})
 pd.DataFrame(matched).to_csv(A/'C6A_MATCHED_DISTANCE.csv',index=False)
 gates=[]
 for (panel,enc,strategy),d in effects[effects.role=='primary'].groupby(['panel','encoder','strategy']):
  p=d[d.cohort=='primary'];fold=d[d.cohort!='primary'];po=p[(p.regime=='open')&p.qualifies];ps=p[(p.regime=='contiguous')&p.qualifies];open_patients=sum(v.qualifies.any() for _,v in fold[fold.regime=='open'].groupby('sample_id'));stress_patients=sum(v.qualifies.any() for _,v in fold[fold.regime=='contiguous'].groupby('sample_id'));no_systematic_mse=bool(d.relative_MSE_reduction.mean()>=0 and (d.relative_MSE_reduction<-.02).mean()<=.25);passed=len(po)>=2 and len(ps)>=1 and open_patients>=2 and stress_patients>=2 and no_systematic_mse
  gates.append({'panel':panel,'encoder':enc,'strategy':strategy,'primary_open_qualifying_budgets':len(po),'primary_stress_qualifying_budgets':len(ps),'conditional_open_supporting_patients':int(open_patients),'conditional_stress_supporting_patients':int(stress_patients),'no_systematic_MSE_harm':no_systematic_mse,'pass':bool(passed)})
 e=effects[effects.role=='primary'];material_patients=int(e[e.cohort!='primary'].groupby('sample_id').material_change.any().sum());status='PASS_SAMPLING_EFFECT' if any(g['pass'] for g in gates) else ('MIXED_SAMPLING_EFFECT' if e.qualifies.any() or material_patients>=2 else 'NO_MEANINGFUL_SAMPLING_EFFECT')
 dump(A/'logs/GATE_RESULT.json',{'C6A_status':status,'rules':gates,'any_qualifying_primary_model_settings':int(e.qualifies.sum()),'patients_with_material_change':material_patients,'technical_audit_required':True,'C6B_authorized_by_gate':status=='PASS_SAMPLING_EFFECT','test_evaluation_runs':len(records)})
 # Full gene metrics remain available without restricting attention to favorable genes.
 frames=[]
 for r in records:
  with np.load(ROOT/r['prediction_path']) as f:frames.append(pd.DataFrame({'gene':f['genes'],**{k:f[k] for k in ['pearson','spearman','mse','mae']}}).assign(**{k:r[k] for k in ['run_id','cohort','sample_id','encoder','panel','role','regime','strategy','seed','control','nominal_budget']}))
 genes=pd.concat(frames,ignore_index=True);genes.to_csv(A/'C6A_GENE_METRICS.csv',index=False);genes.to_csv(A/'C6A_GENE_METRICS.csv.gz',index=False,compression='gzip')
 lock=json.loads((A/'TRAJECTORY_LOCK.json').read_text());policies=[json.loads((A/'policies'/f'{key}.json').read_text()) for key in lock['policies']];pd.DataFrame([{k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in r.items()} for r in policies]).to_csv(A/'C6A_POLICY_MANIFEST.csv',index=False);selections=[]
 for r in policies:
  with np.load(ROOT/r['policy_path']) as f:
   order=f['order'];first=np.searchsorted([count(r['n_locations'],b) for b in BUDGETS],np.arange(1,len(order)+1));selections.append(pd.DataFrame({'policy_id':r['policy_id'],'selection_rank':np.arange(1,len(order)+1),'location_index':order,'barcode':f['selected_barcodes'],'first_budget':np.array(BUDGETS)[first]}))
 selections=pd.concat(selections,ignore_index=True);selections.to_csv(A/'C6A_SELECTIONS.csv',index=False);selections.to_csv(A/'C6A_SELECTIONS.csv.gz',index=False,compression='gzip');print(status,flush=True)
if __name__=='__main__':main()
