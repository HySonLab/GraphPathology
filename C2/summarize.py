"""Join C2 outcomes to frozen C1 metrics without running C1 models."""
import os,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C2'
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'data/matplotlib'));os.environ.setdefault('XDG_CACHE_HOME',str(ROOT/'data/.cache'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})

def records(folder):return pd.DataFrame([json.loads(p.read_text()) for p in sorted((OUT/folder).glob('*.json'))])
def summarize():
 runs=records('logs/runs');ev=records('logs/evaluation');genes=pd.concat([pd.read_csv(ROOT/p) for p in ev.gene_metrics_path],ignore_index=True)
 frozen=pd.read_csv(ROOT/'C1/C1_GENE_METRICS.csv');frozen=frozen[(frozen.model=='B2')&frozen.phase.isin(['primary','lopo'])]
 keys=['panel','sample_id','seed','gene'];base=frozen[keys+['phase','pearson','spearman','mse','mae']].rename(columns={k:k+'_local' for k in ['pearson','spearman','mse','mae']}).rename(columns={'phase':'baseline_phase'})
 genes['baseline_phase']=np.where(genes.phase=='lopo','lopo','primary')
 gain=genes.merge(base,on=keys+['baseline_phase'],validate='many_to_one')
 assert len(gain)==len(genes)
 for metric in ['pearson','spearman','mse','mae']:gain['delta_'+metric+'_local']=gain[metric]-gain[metric+'_local']
 pool=genes[(genes.model=='MEAN')&genes.phase.isin(['primary','lopo'])][keys+['phase','pearson','mse','mae']].rename(columns={'phase':'baseline_phase',**{k:k+'_mean' for k in ['pearson','mse','mae']}})
 gain=gain.merge(pool,on=keys+['baseline_phase'],validate='many_to_one')
 for metric in ['pearson','mse','mae']:gain['delta_'+metric+'_mean']=gain[metric]-gain[metric+'_mean']
 gain.to_csv(OUT/'C2_GENE_CONTEXTUAL_GAINS.csv',index=False);genes.to_csv(OUT/'C2_GENE_METRICS.csv',index=False)
 ev[ev.phase=='primary'].to_csv(OUT/'C2_METRICS_SUMMARY.csv',index=False);ev.to_csv(OUT/'C2_PATIENT_METRICS.csv',index=False)
 ev[ev.phase=='lopo'].to_csv(OUT/'C2_ROBUSTNESS_LOPO.csv',index=False);ev[ev.phase=='shuffled_graph'].to_csv(OUT/'C2_GRAPH_CONTROLS.csv',index=False)
 resources=runs.merge(ev[['run_id','inference_seconds','inference_peak_rss_mb']],on='run_id',how='left');resources['threads']=4;resources['cpu']='Intel Core i7-14700';resources.to_csv(OUT/'C2_RESOURCE_PROFILE.csv',index=False)
 runs[runs.phase=='selection'][['run_id','panel','model','graph','validation_median_pearson','validation_mse','best_epoch','parameters','local_parameters']].to_csv(OUT/'C2_GRAPH_SENSITIVITY.csv',index=False)
 records('graphs').to_csv(OUT/'C2_GRAPH_STATISTICS.csv',index=False)
 summary=gain.groupby(['phase','panel','model','sample_id','seed']).agg(mean_delta_local=('delta_pearson_local','mean'),median_delta_local=('delta_pearson_local','median'),
     fraction_genes_positive=('delta_pearson_local',lambda v:float((v>0).mean())),fraction_genes_gain_ge_005=('delta_pearson_local',lambda v:float((v>=.05).mean())),
     mean_delta_vs_mean=('delta_pearson_mean','mean'),median_delta_vs_mean=('delta_pearson_mean','median'),mean_mse_change=('delta_mse_local','mean'),mean_mae_change=('delta_mae_local','mean')).reset_index()
 summary.to_csv(OUT/'C2_CONTEXTUAL_GAIN_SUMMARY.csv',index=False)
 shortlist=pd.read_csv(OUT/'C2_VALIDATION_GENE_SHORTLIST.csv')
 primary=gain[(gain.phase=='primary')&gain.model.isin(['GCN','SAGE'])&(gain.panel!='S1')].groupby(['panel','model','gene']).agg(
  primary_delta_local=('delta_pearson_local','mean'),primary_min_seed_delta_local=('delta_pearson_local','min'),primary_delta_mean=('delta_pearson_mean','mean'),primary_mse_change=('delta_mse_local','mean')).reset_index()
 folds=gain[(gain.phase=='lopo')&gain.model.isin(['GCN','SAGE'])&(gain.panel!='S1')].groupby(['panel','model','gene','sample_id'])[['delta_pearson_local','delta_pearson_mean','delta_mse_local']].mean().reset_index()
 folds['positive_both']=(folds.delta_pearson_local>0)&(folds.delta_pearson_mean>0)
 supports=folds.groupby(['panel','model','gene']).agg(lopo_positive_both_patients=('positive_both','sum'),lopo_mean_delta_local=('delta_pearson_local','mean'),lopo_mean_delta_mean=('delta_pearson_mean','mean'),lopo_mean_mse_change=('delta_mse_local','mean')).reset_index()
 candidates=shortlist.merge(primary,on=['panel','model','gene'],validate='one_to_one').merge(supports,on=['panel','model','gene'],validate='one_to_one')
 candidates['supported_candidate']=candidates.validation_candidate&(candidates.primary_delta_local>=.05)&(candidates.primary_min_seed_delta_local>0)&(candidates.primary_delta_mean>0)&(candidates.primary_mse_change<=0)&(candidates.lopo_positive_both_patients>=3)
 candidates.to_csv(OUT/'C2_GENE_CANDIDATE_AUDIT.csv',index=False);folds.to_csv(OUT/'C2_GENE_PATIENT_DEPENDENCE.csv',index=False)
 return runs,ev,gain,summary,candidates

def figures(runs,ev,gain,summary,candidates):
 local=pd.read_csv(ROOT/'C1/C1_METRICS_SUMMARY.csv');local=local[local.model=='B2']
 fig,axes=plt.subplots(1,3,figsize=(13,4))
 for ax,panel in zip(axes,['S2-50','S2-250','S2-280']):
  groups=[local[local.panel==panel].median_pearson]+[ev[(ev.phase=='primary')&(ev.panel==panel)&(ev.model==m)].median_pearson for m in ['MEAN','GCN','SAGE']]
  ax.bar(range(4),[a.mean() for a in groups],yerr=[a.std() for a in groups],color=['#64748b','#059669','#0369a1','#b45309'])
  ax.set_xticks(range(4),['Local','Mean','GCN','SAGE']);ax.set_title(panel);ax.set_ylim(0,1);ax.set_ylabel('Median gene Pearson')
 fig.suptitle('Primary test patient: contextual models versus frozen local MLP\nError bars are seed SD, not patient uncertainty.');fig.tight_layout(rect=(0,0,1,.86));fig.savefig(OUT/'figures/primary_comparison.png',dpi=150,bbox_inches='tight');plt.close(fig)
 for panel in ['S2-50','S2-250','S2-280']:
  fig,axes=plt.subplots(1,3,figsize=(13,4))
  for ax,model in zip(axes,['MEAN','GCN','SAGE']):
   f=gain[(gain.phase=='primary')&(gain.panel==panel)&(gain.model==model)].groupby('gene')[['delta_pearson_local','delta_mse_local']].mean()
   sc=ax.scatter(f.delta_pearson_local,f.delta_mse_local,s=16,alpha=.75,color={'MEAN':'#059669','GCN':'#0369a1','SAGE':'#b45309'}[model]);ax.axhline(0,c='black',lw=.7);ax.axvline(0,c='black',lw=.7)
   ax.set(title=model,xlabel='Gene Pearson change vs local',ylabel='Gene MSE change vs local')
  fig.suptitle(panel+': correlation gains and error changes (seed means)\nLower-right quadrant improves both; all panel genes are shown.');fig.tight_layout(rect=(0,0,1,.84));fig.savefig(OUT/'figures'/f'{panel}_gene_gains.png',dpi=150,bbox_inches='tight');plt.close(fig)
 for panel in ['S2-50','S2-250','S2-280']:
  f=summary[(summary.phase=='lopo')&(summary.panel==panel)]
  fig,ax=plt.subplots(figsize=(10,4));samples=['NCBI783','NCBI785','TENX95','TENX99'];x=np.arange(4)
  for j,model in enumerate(['MEAN','GCN','SAGE']):
   groups=[f[(f.model==model)&(f.sample_id==s)].median_delta_local for s in samples]
   ax.bar(x+(j-1)*.23,[g.mean() for g in groups],.23,yerr=[g.std() for g in groups],label=model)
  ax.axhline(0,c='black',lw=.8);ax.set_xticks(x,samples);ax.set_ylabel('Median per-gene gain vs frozen local');ax.legend()
  ax.set_title(panel+': conditional patient-dependence sensitivity\nNot unbiased validation; seed SD error bars.');fig.tight_layout();fig.savefig(OUT/'figures'/f'{panel}_patient_dependence.png',dpi=150,bbox_inches='tight');plt.close(fig)
 f=runs[runs.phase=='selection'].copy();f['graph_name']=f.graph.map(lambda g:g['name'])
 fig,axes=plt.subplots(1,3,figsize=(14,4))
 for ax,panel in zip(axes,['S2-50','S2-250','S2-280']):
  for model in ['MEAN','GCN','SAGE']:
   g=f[(f.panel==panel)&(f.model==model)].sort_values('graph_name');ax.plot(range(len(g)),g.validation_median_pearson,marker='o',label=model)
  ax.set_xticks(range(len(g)),g.graph_name,rotation=30,ha='right');ax.set(title=panel,ylabel='Validation median gene Pearson');ax.legend(fontsize=8)
 fig.suptitle('Graph sensitivity: validation only, seed 0');fig.tight_layout(rect=(0,0,1,.93));fig.savefig(OUT/'figures/validation_graph_sensitivity.png',dpi=150,bbox_inches='tight');plt.close(fig)
 fig,axes=plt.subplots(1,2,figsize=(11,4))
 for ax,metric in zip(axes,['median_pearson','mean_mse']):
  for j,phase in enumerate(['primary','shuffled_graph']):
   groups=[ev[(ev.phase==phase)&(ev.panel=='S2-50')&(ev.model==m)][metric] for m in ['MEAN','GCN','SAGE']]
   ax.bar(np.arange(3)+(j-.5)*.32,[g.mean() for g in groups],.32,yerr=[g.std() for g in groups],label=phase)
  ax.set_xticks(range(3),['MEAN','GCN','SAGE']);ax.set_ylabel(metric);ax.legend()
 fig.suptitle('S2-50 spatial specificity: real versus vertex-permuted graphs\nSeed SD; graph degree/edge structure preserved under permutation.');fig.tight_layout(rect=(0,0,1,.85));fig.savefig(OUT/'figures/shuffled_graph_controls.png',dpi=150,bbox_inches='tight');plt.close(fig)

if __name__=='__main__':
 results=summarize();figures(*results);print('Assembled',len(results[1]),'evaluations and',len(results[2]),'gene rows')
