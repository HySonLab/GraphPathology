"""Reader-facing labels only, using immutable validated patient-summary tables."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

P=Path(__file__).resolve().parents[1]
NAMES={'gps_full':'GPS full','gps_local':'GPS local','gps_no_relative':'GPS no distances','gps_self':'GPS self','gps_shuffled':'GPS shuffled','gps_degree_length':'GPS rewiring attempt','distance_fixed':'Distance fixed','egnn_latent':'EGNN latent','distance_disabled':'Distance disabled','egnn_self':'EGNN self','egnn_shuffled':'EGNN shuffled','egnn_degree_length':'EGNN rewiring attempt','matched_local':'Local MLP','matched_mean':'Mean','matched_gcn':'GCN','matched_sage':'GraphSAGE','idw_k16_p2':'IDW (16, power2)','uniform_k16':'Uniform (16)'}
ROUTES={'image':'Image only','fixed_20pct_image':'Fixed fusion','exploratory_residual_anchor':'Observed-residual interpolation'}
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
patient=pd.read_csv(P/'source_data/extension/PATIENT_METRICS.csv');full=patient[(patient.design=='historical_V1')&(patient.budget==.2)]
order=[('new',v) for v in list(NAMES)[:16]]+[('frozen',v) for v in sorted(full[full.origin=='frozen'].variant.unique())]+[('molecular','idw_k16_p2'),('molecular','uniform_k16')]
def name(origin,v):
 if origin!='frozen':return NAMES[v]
 family,control=v.split('_',1);family={'local':'Local MLP','MEAN':'Mean','GCN':'GCN','SAGE':'GraphSAGE'}[family]
 return family+' / '+control+' [minibatch]'
def save(fig,stem):
 for suffix in ['.pdf','.png']:fig.savefig(P/'figures/extension'/(stem+suffix),bbox_inches='tight',**({'dpi':160} if suffix=='.png' else {}))
 plt.close(fig)
for cohort in ['Development','Wu','Bassiouni']:
 for route,title in ROUTES.items():
  fig,axs=plt.subplots(1,3,figsize=(9,10),sharey=True);subset=full[(full.cohort==cohort)&((full.route==route)|(full.origin=='molecular'))]
  for y,(origin,v) in enumerate(order):
   q=subset[(subset.origin==origin)&(subset.variant==v)]
   for ax,metric in zip(axs,['median_pearson','mean_mse','mean_mae']):
    for k,row in enumerate(q.itertuples()):ax.plot(getattr(row,metric),y+(.1*(k-1) if cohort=='Wu' else 0),'o',ms=3,color=plt.cm.tab10(k) if cohort=='Wu' else '.7',alpha=.8)
    ax.plot(q[metric].median(),y,'|',color='black',ms=10)
  axs[0].set_yticks(range(len(order)),[name(o,v) for o,v in order]);axs[0].invert_yaxis()
  for ax,label in zip(axs,['Median gene Pearson','Mean MSE','Mean MAE']):ax.set_xlabel(label);ax.grid(axis='x',alpha=.2)
  axs[0].set_title(('HEST development' if cohort=='Development' else cohort)+' / '+title,loc='left')
  if cohort=='Wu':
   for k,p in enumerate(['1160920F','CID4290','CID44971']):axs[0].plot([],[],'o',color=plt.cm.tab10(k),label=p)
   fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',ncol=3,fontsize=11)
  fig.tight_layout(rect=(0,.04,1,1) if cohort=='Wu' else (0,0,1,1));save(fig,'ExtensionMatrix_'+cohort+'_'+route)
# Preserve the original distance-curve aggregation and all18 series; change titles/legend only.
bm=pd.read_csv(P/'source_data/extension/BASSIOUNI_DISTANCE_PATIENT_METRICS.csv');fig,axs=plt.subplots(3,1,figsize=(9,10),sharex=True)
for ax,(route,title) in zip(axs,ROUTES.items()):
 for j,v in enumerate(NAMES):
  q=bm[(bm.variant==v)&((bm.route==route)|(bm.origin=='molecular'))].groupby('low_um')[['mean_mse']].median()
  ax.plot(q.index,q.mean_mse,'o-',ms=3,color=plt.cm.tab20(j),label=NAMES[v])
 ax.set_title(title,fontsize=12);ax.set_ylabel('Median patient mean MSE');ax.grid(alpha=.2)
axs[-1].set_xlabel('Nearest-anchor distance-bin lower bound (um)');fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',ncol=4,fontsize=9);fig.tight_layout(rect=(0,.12,1,1));save(fig,'ExtensionDistanceBins')

# Original minibatch points and IDW reference; editorial labels only.
def minibatch_labels():
    x=pd.read_csv(P/'source_data/journal/HEST_SEED_METRICS.csv');x=x[x.budget==.2]
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
        ax.set_title('Image-only' if row==0 else 'Observed-residual interpolation',fontsize=10)
        ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.2)
    handles=[plt.Line2D([],[],marker='o',linestyle='',color=c,label=f'Training seed {s}') for s,c in zip(range(3),['#0072B2','#D55E00','#009E73'])]
    handles.append(plt.Line2D([],[],ls='--',color='#555555',label='IDW reference'))
    fig.legend(handles=handles,loc='lower center',ncol=4,fontsize=8)
    fig.suptitle('Minibatch neural comparisons | TENX99 development (n=1)\n20% coverage; points average three fixed masks',fontsize=11)
    fig.tight_layout(rect=(0,.04,1,.93))
    fig.savefig(P/'figures/journal/JournalFigure8.pdf');fig.savefig(P/'figures/journal/JournalFigure8.png',dpi=140);plt.close(fig)

minibatch_labels()
