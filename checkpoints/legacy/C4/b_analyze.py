"""Apply a locked minimal spatial-specificity gate to C4B predictions."""
from common import *
import math
B=OUT/'C4B'
def main():
 cfg=json.loads((B/'configs/protocol.json').read_text());ev=[json.loads(p.read_text()) for p in (B/'logs/evaluation').glob('*.json')];assert len(ev)==240
 p=pd.DataFrame(ev);g=pd.concat([pd.read_csv(ROOT/r['gene_metrics_path']) for r in ev],ignore_index=True);p.to_csv(B/'C4B_PANEL_METRICS.csv',index=False);g.to_csv(B/'C4B_GENE_METRICS.csv',index=False)
 keys=['phase','sample_id','panel','gene'];means=g.groupby(keys+['control'])[['pearson','mse']].mean().unstack('control');means.columns=[f'{a}_{b}' for a,b in means.columns];means=means.reset_index()
 for name,c in [('spatial','LOCAL'),('matched','MATCHED'),('mean','MEAN')]:means['delta_'+name]=means.pearson_TRUE-means['pearson_'+c];means['delta_MSE_'+name]=means['mse_'+c]-means.mse_TRUE
 seed=g.groupby(keys+['training_seed','control'])[['pearson']].mean().unstack('control');seed.columns=seed.columns.droplevel(0);seed['minimum_seed_matched_gain']=seed.TRUE-seed.MATCHED;seed['minimum_seed_local_gain']=seed.TRUE-seed.LOCAL;seed=seed.groupby(keys)[['minimum_seed_matched_gain','minimum_seed_local_gain']].min().reset_index();means=means.merge(seed,on=keys,validate='one_to_one')
 sd=g[g.control=='MATCHED'].groupby(keys+['graph_seed']).pearson.mean().groupby(keys).std(ddof=1).rename('matched_graph_SD').reset_index();means=means.merge(sd,on=keys,validate='one_to_one');means.to_csv(B/'C4B_GENE_GAIN.csv',index=False);means[means.phase=='lopo'].to_csv(B/'C4B_PATIENT_SENSITIVITY.csv',index=False)
 rule=cfg['gene_rule'];rows=[]
 for r in means[means.phase=='primary'].to_dict('records'):
  folds=means[(means.phase=='lopo')&(means.panel==r['panel'])&(means.gene==r['gene'])];assert len(folds)==4
  good=(folds.delta_spatial>0)&(folds.delta_matched>0)&(folds.mse_TRUE<=1.02*folds.mse_LOCAL)&(folds.mse_TRUE<=1.02*folds.mse_MATCHED)
  r['supporting_patients']=int(good.sum());r['mean_fold_relative_MSE_increase_vs_matched']=float((folds.mse_TRUE/folds.mse_MATCHED-1).mean())
  r['primary_candidate']=bool(r['delta_matched']>=.02 and r['delta_spatial']>0 and r['delta_MSE_spatial']>=0 and r['delta_MSE_matched']>=0)
  r['robust_gene']=bool(min(r['delta_spatial'],r['delta_mean'],r['delta_matched'])>=rule['minimum_primary_mean_Pearson_gain_over_local_mean_and_matched'] and r['minimum_seed_matched_gain']>0 and r['minimum_seed_local_gain']>0 and r['delta_matched']>r['matched_graph_SD'] and r['delta_MSE_spatial']>=0 and r['delta_MSE_matched']>=0 and r['supporting_patients']>=rule['minimum_conditional_patients_positive_local_and_matched_gain_with_MSE_increase_at_most_2percent_against_both'] and r['mean_fold_relative_MSE_increase_vs_matched']<=rule['maximum_mean_fold_relative_MSE_increase_vs_matched']);rows.append(r)
 support=pd.DataFrame(rows);support.to_csv(B/'C4B_GENE_SUPPORT.csv',index=False);panels=[]
 for panel,d in support.groupby('panel'):
  panels.append({'panel':panel,'robust_genes':int(d.robust_gene.sum()),'required_genes':max(3,math.ceil(.02*len(d))),'primary_candidates':int(d.primary_candidate.sum()),'median_spatial_gene_gain':d.delta_spatial.median(),'median_matched_gene_gain':d.delta_matched.median(),'mean_spatial_MSE_reduction':d.delta_MSE_spatial.mean(),'mean_matched_MSE_reduction':d.delta_MSE_matched.mean()})
 status='PASS_SPATIAL_RECHECK' if any(r['robust_genes']>=r['required_genes'] for r in panels) else ('MIXED_SPATIAL_RECHECK' if any(r['primary_candidates']>=3 for r in panels) else 'FAIL_SPATIAL_RECHECK');dump(B/'logs/GATE_RESULT.json',{'status':status,'panels':panels,'descriptive_rules_only':True});pd.DataFrame(panels).to_csv(B/'C4B_GAIN_SUMMARY.csv',index=False)
 # Direct frozen-C3 contrasts; no old fitting, inference or rescoring.
 old=pd.read_csv(ROOT/'C3/C3A/C3A_GENE_METRICS.csv');old=old[old.panel.isin(cfg['panels'])&old.control.isin(['A0','A2','A3','A6'])];oldmean=old.groupby(keys+['control'])[['pearson','mse']].mean().unstack('control');oldmean.columns=[f'{a}_{b}' for a,b in oldmean.columns];oldmean=oldmean.reset_index()
 for name,c in [('spatial','A0'),('matched','A6'),('mean','A2')]:oldmean['C3_delta_'+name]=oldmean.pearson_A3-oldmean['pearson_'+c];oldmean['C3_delta_MSE_'+name]=oldmean['mse_'+c]-oldmean.mse_A3
 comparison=means.merge(oldmean[keys+[c for c in oldmean if c.startswith('C3_')]],on=keys,validate='one_to_one');comparison.to_csv(B/'C4B_VS_FROZEN_C3.csv',index=False)
 pd.DataFrame([json.loads(f.read_text()) for f in (B/'logs/runs').glob('*.json')]).to_csv(B/'C4B_RESOURCE_PROFILE.csv',index=False)
 events=[json.loads(f.read_text()) for f in (B/'logs/registry').glob('*.json')]
 for r in events:
  for k,v in list(r.items()):
   if isinstance(v,(dict,list)):r[k]=json.dumps(v,sort_keys=True)
 pd.DataFrame(events).to_csv(B/'C4B_REGISTRY_EXTENSIONS.csv',index=False);print(status)
if __name__=='__main__':main()
