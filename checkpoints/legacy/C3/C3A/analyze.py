"""Apply the pre-test descriptive gate; retain all genes and randomizations."""
from common import *
import math

def analyze():
 check_lock()
 evals=[json.loads(p.read_text()) for p in sorted((OUT/'logs/evaluation').glob('*.json'))]
 panel=pd.DataFrame(evals);genes=pd.concat([pd.read_csv(ROOT/r['gene_metrics_path']) for r in evals],ignore_index=True)
 panel.to_csv(OUT/'C3A_PANEL_METRICS.csv',index=False);genes.to_csv(OUT/'C3A_GENE_METRICS.csv',index=False)
 randomized=genes[genes.control.isin(['A4','A5','A6'])]
 keys=['phase','panel','sample_id','control','gene','graph_randomization_seed']
 rg=randomized.groupby(keys,as_index=False)[['pearson','spearman','mse','mae']].mean()
 rg.to_csv(OUT/'C3A_RANDOMIZATION_METRICS.csv',index=False)
 uncertainty=rg.groupby(keys[:-1]).agg(pearson_mean=('pearson','mean'),pearson_graph_SD=('pearson','std'),mse_mean=('mse','mean'),mse_graph_SD=('mse','std')).reset_index()
 uncertainty.to_csv(OUT/'C3A_RANDOMIZATION_UNCERTAINTY.csv',index=False)
 key=['phase','panel','sample_id','gene'];means=genes.groupby(key+['control'])[['pearson','mse','mae']].mean()
 wide=means.unstack('control');wide.columns=[f'{metric}_{control}' for metric,control in wide.columns];wide=wide.reset_index()
 for control in ['A0','A1','A2','A4','A5','A6']:
  wide[f'delta_{control}']=wide.pearson_A3-wide[f'pearson_{control}']
  wide[f'mse_change_{control}']=wide.mse_A3-wide[f'mse_{control}']
  wide[f'relative_mse_change_{control}']=wide[f'mse_change_{control}']/wide[f'mse_{control}']
 wide.to_csv(OUT/'C3A_GENE_SPATIAL_GAINS.csv',index=False)
 wide[wide.phase=='lopo'].to_csv(OUT/'C3A_PATIENT_SENSITIVITY.csv',index=False)
 # Training-seed contrasts keep randomization uncertainty distinct from optimizer variability.
 skey=key+['training_seed'];sm=genes.groupby(skey+['control'])[['pearson','mse']].mean().unstack('control');sm.columns=[f'{metric}_{control}' for metric,control in sm.columns];sm=sm.reset_index();sm['delta_A6']=sm.pearson_A3-sm.pearson_A6
 sm.to_csv(OUT/'C3A_SEED_CONTRASTS.csv',index=False)
 rule=CFG['candidate_rule'];robust=CFG['robust_gene_rule'];records=[]
 for row in wide[wide.phase=='primary'].to_dict('records'):
  panelname,gene=row['panel'],row['gene'];folds=wide[(wide.phase=='lopo')&(wide.panel==panelname)&(wide.gene==gene)]
  u=uncertainty[(uncertainty.phase=='primary')&(uncertainty.panel==panelname)&(uncertainty.gene==gene)&(uncertainty.control=='A6')].iloc[0]
  seed=sm[(sm.phase=='primary')&(sm.panel==panelname)&(sm.gene==gene)]
  row.update(lopo_positive_A6_patients=int((folds.delta_A6>0).sum()),lopo_worst_A6_gain=float(folds.delta_A6.min()),lopo_worst_relative_MSE_A0=float(folds.relative_mse_change_A0.max()),lopo_mean_MSE_change_A6=float(folds.mse_change_A6.mean()),A6_graph_randomization_SD=float(u.pearson_graph_SD),minimum_training_seed_A6_gain=float(seed.delta_A6.min()))
  flags={'positive_true_pearson':row['pearson_A3']>rule['true_pearson_gt'], 'beats_local':row['delta_A0']>rule['mean_gain_over_A0_A5_A6_gt'],'beats_degree':row['delta_A5']>rule['mean_gain_over_A0_A5_A6_gt'],'beats_distance':row['delta_A6']>rule['mean_gain_over_A0_A5_A6_gt'],'primary_MSE_acceptable':row['relative_mse_change_A0']<=rule['maximum_relative_primary_MSE_increase_vs_A0'],'patient_direction_consistent':row['lopo_positive_A6_patients']>=rule['minimum_lopo_patients_with_positive_mean_A6_gain'],'no_catastrophic_patient_correlation':row['lopo_worst_A6_gain']>=rule['minimum_worst_patient_mean_A6_gain'],'no_catastrophic_patient_MSE':row['lopo_worst_relative_MSE_A0']<=rule['maximum_worst_patient_relative_MSE_increase_vs_A0']}
  row.update(flags);row['candidate']=all(flags.values())
  row['minimum_gain_over_nonlocal_controls']=min(row[f'delta_{c}'] for c in ['A1','A2','A4','A5','A6'])
  row['robust']=bool(row['candidate'] and row['minimum_gain_over_nonlocal_controls']>=robust['minimum_primary_mean_gain_over_each_A1_A2_A4_A5_A6'] and row['minimum_training_seed_A6_gain']>0 and row['delta_A6']>row['A6_graph_randomization_SD'] and row['mse_change_A0']<=0 and row['mse_change_A6']<=0 and row['lopo_mean_MSE_change_A6']<=0)
  records.append(row)
 candidates=pd.DataFrame(records);candidates.to_csv(OUT/'C3A_CANDIDATE_AUDIT.csv',index=False)
 graphrows=[json.loads(p.read_text()) for p in sorted((OUT/'graphs').glob('*.json'))]
 graphframe=pd.DataFrame(graphrows);graphframe.to_csv(OUT/'C3A_GRAPH_AUDIT.csv',index=False)
 limits=CFG['graph_audit_limits'];issues=[]
 for r in graphrows:
  c=r['graph_control_type']
  if c in ['A4','A5'] and r['changed_edge_fraction']<limits['A4_A5_min_changed_edge_fraction']:issues.append(f"{r['path']}: inadequate edge turnover")
  if c=='A6':
   if r['changed_edge_fraction']<limits['A6_min_changed_edge_fraction']:issues.append(f"{r['path']}: inadequate A6 edge turnover")
   if r['edge_length_ks_rounded_0_001_um']>limits['A6_max_edge_length_KS']:issues.append(f"{r['path']}: A6 length KS")
   if abs(r['relative_mean_edge_length_change'])>limits['A6_max_absolute_relative_mean_length_change']:issues.append(f"{r['path']}: A6 mean length shift")
 passed=candidates[candidates.robust];enough=any(len(passed[passed.panel==p])>=max(3,math.ceil(.02*len(PANELS[p]))) for p in PANELS if p!='S1')
 weak=candidates[(candidates.delta_A6>=.02)&(candidates.mse_change_A6<=0)]
 if issues:status='FAIL_TECHNICAL'
 elif passed.gene.nunique()>=3 and enough:status='PASS_SPATIAL_SPECIFICITY'
 elif candidates.candidate.any() or weak.gene.nunique()>=3:status='MIXED_SPATIAL_SPECIFICITY'
 else:status='FAIL_SPATIAL_SPECIFICITY'
 result={'status':status,'graph_audit_issues':issues,'candidate_gene_panel_count':int(candidates.candidate.sum()),'candidate_unique_genes':sorted(candidates[candidates.candidate].gene.unique()),'robust_gene_panel_count':len(passed),'robust_unique_genes':sorted(passed.gene.unique()),'weak_primary_gene_panel_count':len(weak),'evaluations':len(evals),'new_fits':len(list((OUT/'logs/runs').glob('*.json'))),'C3B_authorized_by_gate':status=='PASS_SPATIAL_SPECIFICITY','spec_sha256':sha(OUT/'C3A_SPEC_LOCK.json')}
 dump(OUT/'logs/GATE_RESULT.json',result)
 registryrows=[json.loads(p.read_text()) for p in sorted((OUT/'logs/registry').glob('*.json'))]
 for r in registryrows:
  for k,v in list(r.items()):
   if isinstance(v,(dict,list)):r[k]=json.dumps(v,sort_keys=True)
 pd.DataFrame(registryrows).to_csv(OUT/'C3A_REGISTRY_EXTENSIONS.csv',index=False)
 print(json.dumps(result,indent=2))
if __name__=='__main__':analyze()
