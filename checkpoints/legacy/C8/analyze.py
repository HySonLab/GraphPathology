"""Predeclared C8 gates; common-unit errors and patient-level paired effects."""
from contract import *
import argparse,gzip
METRICS=['mean_pearson','median_pearson','mean_spearman','median_spearman','mean_mse','median_mse','mean_mae','median_mae','common_mse','common_mae','target_mean','mean_target_variance']
def load(stage):return pd.DataFrame([{k:v for k,v in json.loads(p.read_text()).items() if k!='distance_bins'} for p in sorted((OUT/stage/'logs/runs').glob('*.json'))])
def patient(df,extra=[]):return df.groupby(['fold','sample_id','patient_id','transform','model','regime','budget']+extra,as_index=False)[METRICS].mean()
def paired(p,a,b):
 keys=['fold','sample_id','regime','budget'];left=p[p.model==a].set_index(keys);right=p[p.model==b].set_index(keys);assert left.index.equals(right.index);out=left[[]].copy()
 for met in ['median_pearson','common_mse','common_mae']:
  out[a+'_'+met]=left[met];out[b+'_'+met]=right[met]
 out['gain_pearson']=left.median_pearson-right.median_pearson;out['mse_reduction']=right.common_mse-left.common_mse;out['relative_mse_reduction']=out.mse_reduction/np.maximum(right.common_mse,1e-12);out['mae_reduction']=right.common_mae-left.common_mae;return out.reset_index()
def analyze_A():
 df=load('C8A');assert len(df)==4320;df.to_csv(A/'C8A_PANEL_METRICS.csv',index=False);p=patient(df);p.to_csv(A/'C8A_PATIENT_METRICS.csv',index=False);lock=json.loads((A/'FIT_LOCK.json').read_text());sels={f:v['selected_transform'] for f,v in lock['folds'].items()};selected=p[p.apply(lambda r:r['transform']==sels[r['fold']],axis=1)].copy();selected.to_csv(A/'C8A_SELECTED_PATIENT_METRICS.csv',index=False)
 sel=selected[selected.model=='E2'].copy();sel['model']='SELECTED';ref=p[(p['transform']=='T0')&(p.model=='E2')].copy();ref['model']='T0';eff=paired(pd.concat([sel,ref]).sort_values(['fold','sample_id','regime','budget']),'SELECTED','T0');eff.to_csv(A/'C8A_TARGET_EFFECTS.csv',index=False);primary=eff[(eff.regime=='V1')&(eff.budget==.2)];stress=eff[eff.regime=='V2'];n_nonzero=sum(v!='T0' for v in sels.values());joint=int(((primary.gain_pearson>=-.01)&(primary.relative_mse_reduction>=.05)).sum());iqr=lambda x:float(x.quantile(.75)-x.quantile(.25));instability=iqr(primary.SELECTED_common_mse)<=iqr(primary.T0_common_mse);budgets=eff[eff.regime=='V1'].groupby('budget')[['gain_pearson','mse_reduction']].median();budget_support=int(((budgets.gain_pearson>=0)&(budgets.mse_reduction>0)).sum());stresssafe=int((stress.relative_mse_reduction>=-.1).sum());passed=n_nonzero>=3 and primary.gain_pearson.median()>=-.005 and primary.relative_mse_reduction.median()>=.1 and joint>=6 and instability and budget_support>=2 and stresssafe>=6 and stress.gain_pearson.median()>=-.02;partial=int(((primary.gain_pearson>=0)&(primary.mse_reduction>0)).sum())>=2
 status='PASS_TARGET_HARMONIZATION' if passed else ('MIXED_TARGET_HARMONIZATION' if partial else 'NO_TARGET_HARMONIZATION_GAIN');r={'status':status,'selected_transforms':sels,'non_T0_folds':n_nonzero,'V1_20_median_paired_Pearson_gain':float(primary.gain_pearson.median()),'V1_20_median_relative_common_MSE_reduction':float(primary.relative_mse_reduction.median()),'V1_20_joint_support_patients':joint,'selected_common_MSE_IQR':iqr(primary.SELECTED_common_mse),'T0_common_MSE_IQR':iqr(primary.T0_common_mse),'error_instability_not_increased':bool(instability),'V1_supported_budgets':budget_support,'V2_safe_MSE_patients':stresssafe,'V2_median_Pearson_gain':float(stress.gain_pearson.median()),'proceed_B':bool(passed)};dump(A/'C8A_STATUS.json',r)
 scales=[]
 for path in sorted((A/'scalers').glob('*.json')):
  ob=json.loads(path.read_text())
  for domain in ['source','destination']:
   for i,g in enumerate(GENES):scales.append({'scaler_path':str(path.relative_to(ROOT)),'transform':ob['transform'],'train_patients':'|'.join(ob['train_patients']),'domain':domain,'gene':g,'center':ob[domain]['center'][i],'scale':ob[domain]['scale'][i],'robust_SD_fallback':ob[domain]['robust_SD_fallback'][i]})
 pd.DataFrame(scales).to_csv(A/'C8A_TARGET_TRANSFORMS.csv',index=False);print(json.dumps(r,indent=2),flush=True)
def analyze_B():
 df=load('C8B');clean=df[df.control.isna()];assert len(clean)==2400 and len(df)-len(clean)==576;df.to_csv(B/'C8B_PANEL_METRICS.csv',index=False);p=patient(clean);p.to_csv(B/'C8B_PATIENT_METRICS.csv',index=False);eff=paired(p,'SELECTED','F0');eff.to_csv(B/'C8B_PATIENT_EFFECTS.csv',index=False);v1=eff[eff.regime=='V1'];v2=eff[eff.regime=='V2'];bud=[]
 for b,g in v1.groupby('budget'):
  count=int(((g.gain_pearson>0)&(g.relative_mse_reduction>=-.01)).sum());supported=g.gain_pearson.median()>=.01 and g.relative_mse_reduction.median()>=0 and count>=5;bud.append({'budget':b,'median_Pearson_gain':float(g.gain_pearson.median()),'median_relative_common_MSE_reduction':float(g.relative_mse_reduction.median()),'joint_improved_patients':count,'supported':bool(supported)})
 geom=patient(clean,['geometry']);gl=[]
 for g,sub in geom[geom.regime=='V2'].groupby('geometry'):
  ee=paired(sub,'SELECTED','F0');gl.append({'geometry':int(g),'median_Pearson_gain':float(ee.gain_pearson.median()),'median_relative_MSE_reduction':float(ee.relative_mse_reduction.median()),'supported':bool(ee.gain_pearson.median()>=0 and ee.relative_mse_reduction.median()>=0)})
 models={f:json.loads((B/'calibration'/f'{f}.json').read_text())['selected_family'] for f in cfg()['folds']};harm=float((v1.relative_mse_reduction<-.02).mean());passed=sum(x['supported'] for x in bud)>=2 and harm<=.25 and v2.gain_pearson.median()>=0 and v2.relative_mse_reduction.median()>=0 and sum(x['supported'] for x in gl)>=2 and sum(v!='F0' for v in models.values())>=2;partial=int(v1.assign(joint=(v1.gain_pearson>=.01)&(v1.mse_reduction>=0)).groupby('sample_id').joint.any().sum())>=2;status='PASS_FUSION_CALIBRATION' if passed else ('MIXED_FUSION_CALIBRATION' if partial else 'NO_FUSION_ADVANTAGE');r={'status':status,'selected_families':models,'budget_support':bud,'geometry_support':gl,'fraction_V1_settings_MSE_harm_gt2pct':harm,'V2_median_Pearson_gain':float(v2.gain_pearson.median()),'V2_median_relative_MSE_reduction':float(v2.relative_mse_reduction.median()),'proceed_C':bool(passed)};dump(B/'C8B_STATUS.json',r);print(json.dumps(r,indent=2),flush=True)
 families=[]
 for fold in cfg()['folds']:
  obj=json.loads((B/'calibration'/f'{fold}.json').read_text())
  for fam,model in obj['family_models'].items():families.append({'fold':fold,'transform':obj['transform'],'family':fam,'selected':fam==obj['selected_family'],'train_patients':'|'.join(obj['train_patients']),'model_json':json.dumps(model,sort_keys=True)})
 pd.DataFrame(families).to_csv(B/'C8B_FUSION_MODELS.csv',index=False)
 with gzip.open(B/'C8B_GENE_METRICS.csv.gz','wt') as f:
  first=True
  for _,rr in df.iterrows():
   with np.load(ROOT/rr.prediction_path) as a:
    frame=pd.DataFrame({k:a[k] for k in ['pearson','spearman','mse','mae','common_gene_mse','common_gene_mae']});frame['gene']=a['genes']
   for k in ['run_id','sample_id','fold','transform','model','regime','geometry','budget','seed','control','control_seed']:frame[k]=rr[k]
   frame.to_csv(f,index=False,header=first);first=False

def analyze_C():
 df=load('C8C');assert len(df)==960;df.to_csv(OUT/'C8C/C8C_PANEL_METRICS.csv',index=False);base=load('C8B');base=base[base.control.isna()&base.model.isin(['SELECTED','F0'])];p=patient(pd.concat([base,df]));p.to_csv(OUT/'C8C/C8C_FOLD_METRICS.csv',index=False);eff=paired(p,'SELECTED','F0');eff.to_csv(OUT/'C8C/C8C_PATIENT_EFFECTS.csv',index=False);v=eff[eff.regime=='V2'].groupby('budget')[['gain_pearson','mse_reduction']].median();support=int(((v.gain_pearson>=0)&(v.mse_reduction>=0)).sum());dump(OUT/'C8C/C8C_STATUS.json',{'status':'PASS_INTERNAL_TRANSPORTABILITY' if support>=2 else 'MIXED_INTERNAL_TRANSPORTABILITY','V2_supported_budgets':support,'overall_pass':support>=2})
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['A','B','C']);arg=p.parse_args();check_lock();globals()['analyze_'+arg.stage]()
