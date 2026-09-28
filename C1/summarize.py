"""Assemble descriptive C1 tables and figures from immutable run records."""
import os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'C1'
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'data/matplotlib'))
os.environ.setdefault('XDG_CACHE_HOME',str(ROOT/'data/.cache'))
import argparse,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
from core import correlation_matrix

plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})

def load_records(folder):return [json.loads(p.read_text()) for p in sorted((OUT/folder).glob('*.json'))]

def tables():
 runs=pd.DataFrame(load_records('logs/runs'));evaluation=pd.DataFrame(load_records('logs/evaluation'))
 if evaluation.empty:raise ValueError('No evaluation results')
 primary=evaluation[evaluation.phase=='primary']
 primary.to_csv(OUT/'C1_METRICS_SUMMARY.csv',index=False)
 evaluation.to_csv(OUT/'C1_PATIENT_METRICS.csv',index=False)
 evaluation[evaluation.phase.isin(['NC1','NC2'])].to_csv(OUT/'C1_NEGATIVE_CONTROLS.csv',index=False)
 evaluation[evaluation.phase=='lopo'].to_csv(OUT/'C1_ROBUSTNESS_LOPO.csv',index=False)
 cols=['run_id','phase','panel','model','seed','sample_id','patient_id','coexpression_pearson','coexpression_pair_rmse','n_valid_pairs','n_total_pairs']
 evaluation[cols].to_csv(OUT/'C1_COEXPRESSION_METRICS.csv',index=False)
 genes=pd.concat([pd.read_csv(ROOT/p) for p in evaluation.gene_metrics_path],ignore_index=True)
 genes.to_csv(OUT/'C1_GENE_METRICS.csv',index=False)
 resource=runs.merge(evaluation[['run_id','inference_seconds','inference_peak_rss_mb']],on='run_id',how='left')
 resource['cpu']='Intel Core i7-14700';resource['threads']=4
 resource.to_csv(OUT/'C1_RESOURCE_PROFILE.csv',index=False)
 ranks=[]
 for (panel,model),sub in genes[genes.phase=='primary'].groupby(['panel','model']):
  grouped=sub.groupby('gene',as_index=False).agg(pearson=('pearson','mean'),pearson_seed_sd=('pearson','std'),target_variance=('target_variance','mean'))
  defined=grouped.dropna(subset=['pearson']).sort_values(['pearson','gene'],ascending=[False,True])
  for name,frame in [('top10',defined.head(10)),('bottom10',defined.tail(10).sort_values(['pearson','gene']))]:
   for i,row in enumerate(frame.to_dict('records'),1):ranks.append({'panel':panel,'model':model,'list':name,'rank':i,**row})
 pd.DataFrame(ranks).to_csv(OUT/'C1_GENE_RANKINGS.csv',index=False)
 return runs,evaluation,genes

def plot_gene_distributions(genes):
 for panel in ['S2-50','S2-250','S2-280']:
  fig,axes=plt.subplots(2,2,figsize=(11,8))
  for col,model in enumerate(['B1','B2']):
   data=genes[(genes.phase=='primary')&(genes.panel==panel)&(genes.model==model)]
   data=data.groupby('gene',as_index=False)[['pearson','target_variance']].mean()
   axes[0,col].hist(data.pearson.dropna(),bins=np.linspace(-1,1,31),color=['#0369a1','#b45309'][col],edgecolor='white')
   axes[0,col].set(xlabel='Gene-wise Pearson',ylabel='Genes',title=model+(' (seed mean)' if model=='B2' else ''))
   axes[1,col].scatter(data.target_variance,data.pearson,s=16,alpha=.65,color=['#0369a1','#b45309'][col])
   axes[1,col].set(xlabel='Test target variance (C0 units)',ylabel='Gene-wise Pearson',xscale='symlog')
   axes[0,col].axvline(0,color='black',lw=.8);axes[1,col].axhline(0,color='black',lw=.8)
  fig.suptitle(f'{panel}: descriptive gene predictability on TENX99 (one patient)\nB0 correlations are undefined because its predictions are constant.')
  fig.tight_layout(rect=(0,0,1,.93));fig.savefig(OUT/'figures'/f'{panel}_gene_predictability.png',dpi=150,bbox_inches="tight");plt.close(fig)

def coexpression_figures():
 panel=json.loads((OUT/'SELECTION_LOCK.json').read_text())['panels']['S2-50']
 with np.load(ROOT/'data/cache/TENX99_targets.npz') as f:
  ix=[list(f['genes']).index(g) for g in panel];truth=correlation_matrix(f['normalized_targets'][:,ix])
 for model in ['B0','B1','B2']:
  with np.load(OUT/'predictions'/f'primary_S2-50_{model}_s0_TENX99.npz') as f:pred=correlation_matrix(f['prediction'])
  fig,axes=plt.subplots(1,3,figsize=(14,5));cmap=plt.get_cmap('coolwarm').copy();cmap.set_bad('#dddddd')
  for ax,matrix,title,limit in zip(axes,[truth,pred,pred-truth],['Measured','Predicted','Predicted − measured'],[1,1,2]):
   im=ax.imshow(matrix,vmin=-limit,vmax=limit,cmap=cmap,interpolation='nearest');fig.colorbar(im,ax=ax,shrink=.7)
   ax.set_title(title);ax.set_xlabel('Frozen gene order');ax.set_ylabel('Frozen gene order')
  fig.suptitle(f'S2-50 co-expression: {model}, seed 0, TENX99\nSame frozen gene ordering in all matrices; gray means undefined.')
  fig.tight_layout(rect=(0,0,1,.9));fig.savefig(OUT/'figures'/f'S2-50_coexpression_{model}.png',dpi=150,bbox_inches="tight");plt.close(fig)

def spatial_figures():
 lock=json.loads((OUT/'SELECTION_LOCK.json').read_text())
 with h5py.File(ROOT/'data/raw/IDC/adata/TENX99.h5ad') as f:
  thumb=f['uns/spatial/ST/images/downscaled_fullres'][:]
  sf=float(f['uns/spatial/ST/scalefactors/tissue_downscaled_fullres_scalef'][()])
 manifest=json.loads((ROOT/'reports/C0/C0_DATA_MANIFEST.json').read_text())
 mpp=next(s['microns_per_pixel'] for s in manifest['slides'] if s['sample_id']=='TENX99')
 with np.load(ROOT/'data/cache/TENX99_targets.npz') as f:
  truth=f['normalized_targets'];genes=list(f['genes']);xy=f['xy_um']/1000
 height,width=thumb.shape[:2];extent=[0,width/sf*mpp/1000,height/sf*mpp/1000,0]
 for panel in lock['panels']:
  models=['B0','B1','B2'] if panel=='S1' else ['B1','B2']
  for model in models:
   with np.load(OUT/'predictions'/f'primary_{panel}_{model}_s0_TENX99.npz') as f:pred=f['prediction'];panel_genes=list(f['genes'])
   for gene in lock['map_genes'][panel]:
    y=truth[:,genes.index(gene)];p=pred[:,panel_genes.index(gene)];error=abs(y-p)
    # Full common data range, no independently autoscaled prediction panel.
    lo,hi=float(min(y.min(),p.min())),float(max(y.max(),p.max()))
    fig,axes=plt.subplots(1,4,figsize=(13,7),sharex=True,sharey=True)
    axes[0].imshow(thumb,extent=extent);axes[0].set_title('H&E reference')
    for ax,v,title in zip(axes[1:],[y,p,error],['Measured','Predicted','Absolute error']):
     im=ax.scatter(xy[:,0],xy[:,1],c=v,s=2,cmap='magma' if title=='Absolute error' else 'viridis',
                   vmin=0 if title=='Absolute error' else lo,vmax=float(error.max()) if title=='Absolute error' else hi,rasterized=True)
     ax.set_title(title);fig.colorbar(im,ax=ax,fraction=.04,pad=.02,label='C0 normalized units')
    for ax in axes:ax.set_xlim(extent[:2]);ax.set_ylim(extent[2:]);ax.set_aspect('equal');ax.set_xlabel('x (mm)')
    axes[0].set_ylabel('y (mm)')
    fig.suptitle(f'{panel} · {gene} · {model} seed 0 · TENX99\nGene selected from training variance only; descriptive single-patient map.')
    fig.tight_layout(rect=(0,0,1,.91));fig.savefig(OUT/'figures'/f'{panel}_{gene}_{model}_spatial.png',dpi=130,bbox_inches="tight");plt.close(fig)

def comparison_figures(evaluation):
 data=evaluation[evaluation.panel=='S2-50']
 fig,axes=plt.subplots(1,2,figsize=(11,4))
 for col,model in enumerate(['B1','B2']):
  groups=[data[(data.model==model)&(data.phase==phase)].median_pearson.to_numpy(dtype=float) for phase in ['primary','NC1','NC2']]
  axes[col].bar(range(3),[np.nanmean(a) for a in groups],yerr=[np.nanstd(a,ddof=1) if len(a)>1 else 0 for a in groups],color=['#0369a1','#94a3b8','#cbd5e1'])
  axes[col].set_xticks(range(3),['Correct pairing','Label permutation','Feature permutation']);axes[col].axhline(0,c='black',lw=.8)
  axes[col].set(title=model,ylabel='Median gene-wise Pearson')
 fig.suptitle('S2-50 paired models versus negative controls on TENX99\nError bars: seed SD, not uncertainty across patients.')
 fig.tight_layout(rect=(0,0,1,.86));fig.savefig(OUT/'figures/negative_controls.png',dpi=150,bbox_inches="tight");plt.close(fig)
 lopo=data[data.phase=='lopo']
 if len(lopo):
  fig,ax=plt.subplots(figsize=(10,5));samples=['NCBI783','NCBI785','TENX95','TENX99'];x=np.arange(4)
  for j,model in enumerate(['B1','B2']):
   groups=[lopo[(lopo.model==model)&(lopo.sample_id==s)].median_pearson.to_numpy(dtype=float) for s in samples]
   ax.bar(x+(j-.5)*.32,[np.nanmean(g) for g in groups],.32,yerr=[np.nanstd(g,ddof=1) if len(g)>1 else 0 for g in groups],label=model)
  ax.set_xticks(x,samples);ax.set_ylabel('S2-50 median gene-wise Pearson');ax.axhline(0,color='black',lw=.8);ax.legend()
  ax.set_title('Conditional four-patient robustness sensitivity analysis\nNot unbiased validation; seed SD error bars are not patient uncertainty.')
  fig.tight_layout();fig.savefig(OUT/'figures/conditional_lopo.png',dpi=150,bbox_inches="tight");plt.close(fig)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--figures',action='store_true');args=parser.parse_args()
 _,evaluation,genes=tables()
 if args.figures:
  plot_gene_distributions(genes);coexpression_figures();spatial_figures();comparison_figures(evaluation)
 print(f'Assembled {len(evaluation)} evaluated runs and {len(genes)} gene-metric rows',flush=True)
if __name__=='__main__':main()
