"""Frozen patient-level analysis and decision rules for the C7 zero-shot test."""
from support import *
import gzip

def analyze():
 rows=[json.loads(p.read_text()) for p in sorted((B/'logs/runs').glob('*.json'))];assert len(rows)==8556
 df=pd.DataFrame([{k:v for k,v in r.items() if k!='distance_bins'} for r in rows]);df.to_csv(B/'C7B_PANEL_METRICS.csv',index=False);clean=df[df.control.isna()].copy();assert len(clean)==8448
 keys=['sample_id','patient_id','panel','regime','strategy','nominal_budget'];metric=['mean_pearson','median_pearson','mean_spearman','median_spearman','mean_mse','median_mse','mean_mae','median_mae'];curves=clean.groupby(keys+['model'],as_index=False)[metric].mean();curves.to_csv(B/'C7B_COVERAGE_CURVES.csv',index=False)
 effects=[]
 for key,g in curves.groupby(keys):
  tab=g.set_index('model');r=dict(zip(keys,key))
  for model in ['E0','E1','E2','E3']:
   for m in metric:r[model+'_'+m]=float(tab.loc[model,m])
  for base in ['E0','E1']:
   r['gain_pearson_vs_'+base]=r['E2_median_pearson']-r[base+'_median_pearson'];r['gain_mse_vs_'+base]=r[base+'_mean_mse']-r['E2_mean_mse'];r['relative_mse_reduction_vs_'+base]=r['gain_mse_vs_'+base]/max(r[base+'_mean_mse'],1e-12);r['gain_mae_vs_'+base]=r[base+'_mean_mae']-r['E2_mean_mae'];r['qualifies_vs_'+base]=r['gain_pearson_vs_'+base]>=(.02 if base=='E0' else .01) and r['relative_mse_reduction_vs_'+base]>=(.02 if base=='E0' else .01)
  effects.append(r)
 effects=pd.DataFrame(effects);effects.to_csv(B/'C7B_PATIENT_EFFECTS.csv',index=False)
 primary=effects[(effects.panel=='S2-50')&(effects.regime=='V1')&(effects.strategy=='random')];assert len(primary)==24
 q0=primary.groupby('patient_id').qualifies_vs_E0.sum();q1=primary.groupby('patient_id').qualifies_vs_E1.sum();support0=int((q0>=2).sum());support1=int((q1>=1).sum());p10=primary[primary.nominal_budget==.1].set_index('patient_id');p40=primary[primary.nominal_budget==.4].set_index('patient_id');coherent=int(((p40.E2_median_pearson>p10.E2_median_pearson)&(p40.E2_mean_mse<p10.E2_mean_mse)).sum());harm=float((primary.relative_mse_reduction_vs_E1<-.02).mean());patientmed=primary.groupby('patient_id')[['gain_pearson_vs_E0','gain_pearson_vs_E1','gain_mse_vs_E0','gain_mse_vs_E1']].mean().median()
 pass_v1=support0>=6 and support1>=4 and (patientmed>0).all() and harm<=.25 and coherent>=6
 stress=effects[(effects.panel=='S2-50')&(effects.regime=='V2')&(effects.strategy=='random')];stress0=int(stress.groupby('patient_id').qualifies_vs_E0.any().sum());stress1=int(stress.groupby('patient_id').qualifies_vs_E1.any().sum());pass_v2=stress0>=4 and stress1>=4
 mixed=int(primary.groupby('patient_id').qualifies_vs_E0.any().sum())>=2 or int(primary.groupby('patient_id').qualifies_vs_E1.any().sum())>=2
 status=('PASS_EXTERNAL_VALIDATION' if pass_v2 else 'PASS_LIMITED_EXTERNAL_VALIDATION') if pass_v1 else ('MIXED_EXTERNAL_VALIDATION' if mixed else 'FAIL_EXTERNAL_VALIDATION')
 gate={'scientific_status_pending_integrity_audit':status,'primary_panel':'S2-50 common 40 genes','external_patients':8,'E2_over_E0_patients_at_least_two_budgets':support0,'E2_over_E1_patients_at_least_one_budget':support1,'coherent_coverage_patients':coherent,'fraction_primary_settings_MSE_harm_vs_E1_gt2pct':harm,'median_patient_gains':patientmed.to_dict(),'V1_pass':bool(pass_v1),'V2_support_over_E0_patients':stress0,'V2_support_over_E1_patients':stress1,'V2_support':bool(pass_v2),'recommendation':'PROCEED_TO_MANUSCRIPT_BENCHMARKING' if pass_v1 and pass_v2 else ('NEED_EXTERNAL_DATASET_EXPANSION' if pass_v1 else 'NEED_METHOD_REVIEW')};dump(B/'logs/GATE_RESULT.json',gate)
 distances=[]
 for r in rows:
  if r['control']:continue
  for b in r['distance_bins']:distances.append({**{k:r[k] for k in ['run_id','sample_id','patient_id','panel','regime','geometry','strategy','seed','nominal_budget','model']},**b})
 pd.DataFrame(distances).to_csv(B/'C7B_DISTANCE_TO_ANCHOR.csv',index=False)
 ncs=[];lookup={(r['mask_id'],r['panel'],r['nominal_budget'],r['model']):r for r in rows if not r['control']}
 for r in rows:
  if not r['control']:continue
  base=lookup[(r['mask_id'],r['panel'],r['nominal_budget'],r['model'])];ncs.append({**{k:r[k] for k in ['run_id','sample_id','panel','regime','model','control','control_seed','nominal_budget']},'pearson_change_control_minus_clean':r['median_pearson']-base['median_pearson'],'mse_change_control_minus_clean':r['mean_mse']-base['mean_mse'],'mae_change_control_minus_clean':r['mean_mae']-base['mean_mae']})
 pd.DataFrame(ncs).to_csv(B/'C7B_NEGATIVE_CONTROLS.csv',index=False)
 with gzip.open(B/'C7B_GENE_METRICS.csv.gz','wt',newline='') as stream:
  fields=['run_id','sample_id','patient_id','panel','regime','geometry','strategy','seed','nominal_budget','model','control','control_seed','gene','pearson','spearman','mse','mae'];writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
  for r in rows:
   with np.load(ROOT/r['prediction_path']) as f:
    for i,g in enumerate(f['genes']):writer.writerow({**{k:r[k] for k in fields if k in r},'gene':g,**{m:float(f[m][i]) for m in ['pearson','spearman','mse','mae']}})
 print(json.dumps(gate,indent=2),flush=True)
if __name__=='__main__':check_lock();analyze()
