"""Deterministic journal-only graphics from actual, explicitly identified numerical sources."""
import json
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .journal_data import load_patient,original_masks
from .journal_methods import predict_anchors
from .policy import count

GENES=['ERBB2','EPCAM','FOXA1','POSTN']


def initialize_package(root):
    source=root/'manuscript/revised';dest=root/'manuscript/journal_revision'
    if not dest.exists():
        shutil.copytree(source,dest,ignore=shutil.ignore_patterns('main.pdf','VALIDATION_SUMMARY.txt','Codex_Publication_Sync_Instructions.md'))
    (dest/'figures/journal').mkdir(parents=True,exist_ok=True);(dest/'source_data/journal').mkdir(parents=True,exist_ok=True)
    return dest


def molecular_figures(root,data_root,package):
    result=root/'experiments/journal_revision/results/wu_anchors'
    patients=pd.read_csv(result/'PATIENT_METRICS.csv');pairs=pd.read_csv(result/'PAIRED_PATIENT_EFFECTS.csv')
    for name in ['PATIENT_METRICS.csv','PAIRED_PATIENT_EFFECTS.csv','PATIENT_UNCERTAINTY.csv','B2_PATIENT_ANCHOR_DISTANCE.csv','GEOMETRY_METRICS.csv','VALID_CORRELATION_COUNTS.csv']:
        shutil.copyfile(result/name,package/'source_data/journal'/('WU_'+name))
    old=root/'experiments/journal_revision/results/original_v1_reanalysis_r2'
    for name in ['ORIGINAL_V1_PAIRED_PATIENT_EFFECTS.csv','ORIGINAL_V1_PATIENT_BOOTSTRAP.csv']:shutil.copyfile(old/name,package/'source_data/journal'/name)
    order=['1160920F','CID4290','CID44971'];colors=['#0072B2','#D55E00','#009E73']
    fig,axes=plt.subplots(2,1,figsize=(7,6.5),layout='constrained')
    a=patients[(patients.budget==.2)&patients.regime.isin(['V1','V2'])&patients.method.isin(['idw_k16_p2','uniform_k16'])]
    for patient,color in zip(order,colors):
        values=[]
        for regime,method in [('V1','idw_k16_p2'),('V1','uniform_k16'),('V2','idw_k16_p2'),('V2','uniform_k16')]:values.append(float(a[(a.patient==patient)&(a.regime==regime)&(a.method==method)].median_pearson.iloc[0]))
        axes[0].plot(range(4),values,'o-',color=color,label=patient)
    axes[0].set_xticks(range(4),['V1 IDW','V1 uniform','V2 IDW','V2 uniform'],rotation=15);axes[0].set_ylabel('Patient median gene Pearson r');axes[0].set_title('A  Operator x historical protocol (20%)');axes[0].legend(fontsize=8)
    b=pairs[(pairs.budget==.2)&(pairs.origin=='new_B2_identical_targets')]
    for j,(method,color) in enumerate([('idw_k16_p2','#0072B2'),('uniform_k16','#D55E00')]):
        values=[float(b[(b.patient==p)&b.contrast.str.startswith(method)].delta_pearson.iloc[0]) for p in order]
        axes[1].plot(np.arange(3)+(j-.5)*.12,values,'o',color=color,label=method)
    axes[1].axhline(0,color='gray',lw=.8);axes[1].set_xticks(range(3),order);axes[1].set_ylabel('Region-exclusion minus distributed r');axes[1].set_title('B  Same H and anchor count (20%)');axes[1].legend(fontsize=8)
    for ax in axes:ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.2)
    fig.savefig(package/'figures/journal/JournalFigure7.pdf');fig.savefig(package/'figures/journal/JournalFigure7.png',dpi=150);plt.close(fig)
    data=load_patient(root,data_root,'wu','1160920F');xy=data['xy'];y=data['targets'];n=len(xy)
    for regime in ['V1','V2']:
        record=next(r for r in original_masks(root,'wu','1160920F') if r['regime']==regime and r['seed']==0 and r['geometry']==(0 if regime=='V1' else 1))
        with np.load(record['absolute_mask_path']) as f:obs=f['order'][:count(n,.2)];hidden=np.setdiff1d(np.arange(n),obs) if regime=='V1' else np.flatnonzero(f['forbidden'])
        config={'kind':'idw','k':16,'p':2} if regime=='V1' else {'kind':'uniform','k':16}
        pred,_=predict_anchors(xy,y,obs,hidden,config)
        fig,axes=plt.subplots(4,4,figsize=(7,8.5),layout='constrained');source=[]
        for row,gene in enumerate(GENES):
            g=list(data['genes']).index(gene);truth=y[hidden,g];yp=pred[:,g];low=min(truth.min(),yp.min());high=max(truth.max(),yp.max());err=abs(truth-yp)
            constant_span = high-low <= 1e-12
            if constant_span:
                center=(low+high)/2;low=center-.01;high=center+.01
            for col,ax in enumerate(axes[row]):
                ax.scatter(xy[:,0],xy[:,1],s=2,color='#e3e3e3',rasterized=True)
                if col==1:artist=ax.scatter(xy[obs,0],xy[obs,1],c=y[obs,g],s=10,marker='x',cmap='viridis',vmin=low,vmax=high,rasterized=True)
                else:
                    values=truth if col==0 else yp if col==2 else err
                    artist=ax.scatter(xy[hidden,0],xy[hidden,1],c=values,s=4,cmap='magma' if col==3 else 'viridis',vmin=0 if col==3 else low,vmax=max(err.max(),1e-8) if col==3 else high,rasterized=True)
                clip_low=col==1 and y[obs,g].min()<low;clip_high=col==1 and y[obs,g].max()>high
                extension='both' if clip_low and clip_high else 'min' if clip_low else 'max' if clip_high else 'neither'
                bar=fig.colorbar(artist,ax=ax,shrink=.65,pad=.015,extend=extension);bar.ax.tick_params(labelsize=7.5);bar.ax.yaxis.get_offset_text().set_fontsize(7);ax.set_aspect('equal');ax.invert_yaxis();ax.set_xticks([]);ax.set_yticks([])
                if row==0:ax.set_title(['Truth\nscored sites','Observed\nanchors (x)','Prediction\nscored sites','Absolute\nerror'][col],fontsize=8)
                if col==0:ax.set_ylabel(gene+('\nconstant field' if constant_span else ''),fontsize=9)
            for index,i in enumerate(hidden):source.append({'regime':regime,'mask_sha256':record['mask_sha256'],'patient':'1160920F','barcode':data['barcodes'][i],'gene':gene,'x_um':xy[i,0],'y_um':xy[i,1],'role':'scored','truth':truth[index],'prediction':yp[index],'absolute_error':err[index]})
            context=np.setdiff1d(np.arange(n),np.union1d(obs,hidden))
            for i in context:source.append({'regime':regime,'mask_sha256':record['mask_sha256'],'patient':'1160920F','barcode':data['barcodes'][i],'gene':gene,'x_um':xy[i,0],'y_um':xy[i,1],'role':'context','truth':np.nan,'prediction':np.nan,'absolute_error':np.nan})
            for i in obs:source.append({'regime':regime,'mask_sha256':record['mask_sha256'],'patient':'1160920F','barcode':data['barcodes'][i],'gene':gene,'x_um':xy[i,0],'y_um':xy[i,1],'role':'observed','truth':y[i,g],'prediction':np.nan,'absolute_error':np.nan})
        fig.suptitle(f'Wu 1160920F | {regime} | 20% seed0 | standardized log1p counts\nScored truth/prediction share scales per gene; gray = other registered sites; no tissue annotations',fontsize=9)
        fig.savefig(package/f'figures/journal/JournalMap_{regime}.pdf');fig.savefig(package/f'figures/journal/JournalMap_{regime}.png',dpi=120);plt.close(fig)
        pd.DataFrame(source).to_csv(package/f'source_data/journal/WU_REAL_MAP_{regime}.csv.gz',index=False)
    return {'figures':['JournalFigure7','JournalMap_V1','JournalMap_V2'],'map_patient':'1160920F','unavailable_map_patient':'BASS_P01; source acquisition blocked; no substituted map'}


def learned_figure(root,package):
    result=root/'experiments/journal_revision/results/hest_learned'
    for name in ['SEED_METRICS.csv','TRAINING_SEED_SUMMARY.csv','EVALUATION_RUN_METRICS.csv']:
        shutil.copyfile(result/name,package/'source_data/journal'/('HEST_'+name))
    x=pd.read_csv(result/'SEED_METRICS.csv');x=x[x.budget==.2]
    configs=[('local','local')]+[(m,c) for m in ['MEAN','GCN','SAGE'] for c in ['spatial','self','shuffled']]
    labels=['Local']+[m+'\n'+c for m,c in configs[1:]]
    fig,axes=plt.subplots(2,1,figsize=(7,6.5))
    for row,route in enumerate(['image','exploratory_residual_anchor']):
        ax=axes[row]
        for j,(model,control) in enumerate(configs):
            a=x[(x.model==model)&(x.control==control)&(x.route==route)].sort_values('training_seed')
            if len(a)!=3:raise ValueError('Missing matched training seeds in figure')
            ax.scatter(j+np.array([-.12,0,.12]),a['median_pearson'],color=['#0072B2','#D55E00','#009E73'],s=28,zorder=3)
        anchor=x[(x.model=='local')&(x.route=='anchor_idw_k16_p2')]['median_pearson'].iloc[0]
        ax.axhline(anchor,color='#555555',ls='--',lw=1)
        ax.set_xticks(range(len(configs)),labels,fontsize=8);ax.set_ylabel('Median gene Pearson r',fontsize=10)
        ax.set_title('Image-only' if row==0 else 'Exploratory residual + anchors',fontsize=10)
        ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.2)
    handles=[plt.Line2D([],[],marker='o',linestyle='',color=c,label=f'Training seed {s}') for s,c in zip(range(3),['#0072B2','#D55E00','#009E73'])]
    handles.append(plt.Line2D([],[],ls='--',color='#555555',label='IDW reference'))
    fig.legend(handles=handles,loc='lower center',ncol=4,fontsize=8)
    fig.suptitle('New journal_v2 fits | One-patient TENX99 development evaluation\n20% coverage; points average three fixed masks',fontsize=11)
    fig.tight_layout(rect=(0,.04,1,.93))
    fig.savefig(package/'figures/journal/JournalFigure8.pdf');fig.savefig(package/'figures/journal/JournalFigure8.png',dpi=140);plt.close(fig)
    return {'configurations':10,'training_seeds':3,'evaluation_patients':1,'ci_claim':False}
