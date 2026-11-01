"""Static descriptive figures; no biological annotation or test-driven tuning."""
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import sparse
from graph_controls import edge_array,lengths
LABELS={'A0':'Local MLP','A1':'Self-only SAGE','A2':'Mean + MLP','A3':'True SAGE','A4':'Shuffled','A5':'Degree-matched','A6':'Distance-matched'}
COLORS={'A0':'#777777','A1':'#8888bb','A2':'#d69a32','A3':'#197b58','A4':'#ba6c87','A5':'#528cba','A6':'#8f65a3'}

def save(fig,name):
 fig.savefig(OUT/'figures'/name,dpi=170,bbox_inches='tight');plt.close(fig)

def main():
 p=pd.read_csv(OUT/'C3A_PANEL_METRICS.csv');g=pd.read_csv(OUT/'C3A_GENE_SPATIAL_GAINS.csv');c=pd.read_csv(OUT/'C3A_CANDIDATE_AUDIT.csv');primary=p[p.phase=='primary'];pg=g[g.phase=='primary'];panels=list(PANELS)
 fig,axes=plt.subplots(2,4,figsize=(18,8))
 for j,panel in enumerate(panels):
  for metric,ax in zip(['median_pearson','mean_mse'],axes[:,j]):
   for i,control in enumerate(LABELS):
    d=primary[(primary.panel==panel)&(primary.control==control)];vals=d.groupby('graph_randomization_seed')[metric].mean() if control in ['A4','A5','A6'] else d[metric]
    ax.bar(i,vals.mean(),color=COLORS[control],alpha=.65);ax.scatter(i+np.linspace(-.16,.16,len(vals)),vals,s=19,color='black',zorder=3)
   ax.set_xticks(range(7),LABELS,rotation=45);ax.set_title(panel);ax.set_ylabel('Median gene Pearson' if metric=='median_pearson' else 'Mean gene MSE')
 fig.suptitle('Primary patient: all controls; randomized dots = five graph seeds (three training seeds averaged)\nNonrandomized dots = three training seeds; uncertainty is developmental, not patient confidence')
 fig.tight_layout();save(fig,'panel_controls.png')
 for control,tag in [('A4','shuffle'),('A5','degree'),('A6','distance')]:
  fig,axes=plt.subplots(1,4,figsize=(15,3.5))
  for ax,panel in zip(axes,panels):
   d=pg[pg.panel==panel][f'delta_{control}'];ax.hist(d,bins=20 if len(d)>1 else 1,color=COLORS[control]);ax.axvline(0,color='black',lw=1);ax.set_title(f'{panel}; median {d.median():+.3f}');ax.set_xlabel(f'True − {LABELS[control]} Pearson');ax.set_ylabel('Genes')
  fig.suptitle('Per-gene differences; training/randomization means, primary patient');fig.tight_layout();save(fig,f'gain_distribution_{tag}.png')
 fig,axes=plt.subplots(1,3,figsize=(16,4.5))
 for ax,panel in zip(axes,panels[1:]):
  d=pg[pg.panel==panel].sort_values('delta_A6');ax.scatter(np.arange(len(d)),d.delta_A6,s=12,c=np.where(d.gene.isin(c[c.candidate].gene),'#d69a32','#197b58'));ax.axhline(0,color='black',lw=.7);ax.set_title(panel);ax.set_xlabel('All genes ranked by distance-control gain');ax.set_ylabel('True − distance-matched Pearson')
 fig.suptitle('All genes retained; gold marks candidate identities (criteria remain panel-specific)');fig.tight_layout();save(fig,'gene_spatial_gains.png')
 fig,axes=plt.subplots(1,2,figsize=(12,4))
 lopo=g[g.phase=='lopo']
 for ax,metric,title in zip(axes,['delta_A6','mse_change_A6'],['Median per-gene true − distance Pearson','Mean true − distance MSE (negative is better)']):
  table=lopo.pivot_table(index='sample_id',columns='panel',values=metric,aggfunc='median' if metric=='delta_A6' else 'mean').reindex(index=IDS,columns=panels)
  bound=max(abs(table.to_numpy()).max(),1e-6);im=ax.imshow(table,vmin=-bound,vmax=bound,cmap='RdBu' if metric=='delta_A6' else 'RdBu_r',aspect='auto')
  ax.set_xticks(range(4),panels);ax.set_yticks(range(4),IDS);ax.set_title(title)
  for i in range(4):
   for j in range(4):ax.text(j,i,f'{table.iloc[i,j]:+.3f}',ha='center',va='center',fontsize=9,color='black')
  fig.colorbar(im,ax=ax)
 fig.suptitle('Conditional patient sensitivity; overlapping folds, not unbiased cross-validation');fig.tight_layout();save(fig,'patient_comparison.png')
 fig,axes=plt.subplots(1,3,figsize=(15,4))
 for ax,panel in zip(axes,panels[1:]):
  d=pg[pg.panel==panel];ax.scatter(d.delta_A6,d.mse_change_A6,s=14,alpha=.65,color='#197b58');ax.axhline(0,color='black',lw=.8);ax.axvline(0,color='black',lw=.8);ax.set_title(panel);ax.set_xlabel('True − distance Pearson');ax.set_ylabel('True − distance MSE')
 fig.suptitle('Correlation and error can disagree; lower-right quadrant improves both');fig.tight_layout();save(fig,'correlation_vs_error.png')
 fig,axes=plt.subplots(2,4,figsize=(17,7))
 for j,sample in enumerate(IDS):
  xy=load(sample)['xy']
  for control in ['A3','A4','A5','A6']:
   matches=[json.loads(f.read_text()) for f in (OUT/'graphs').glob(f'{sample}_delaunay_cap250_{control}_*.json')]
   a=sparse.load_npz(ROOT/matches[0]['path']);le=lengths(edge_array(a),xy);deg=np.diff(a.indptr)
   axes[0,j].hist(deg,bins=np.arange(0,12)-.5,density=True,histtype='step',label=LABELS[control],color=COLORS[control])
   values=np.sort(le);axes[1,j].plot(values,np.arange(1,len(values)+1)/len(values),label=LABELS[control],color=COLORS[control])
  axes[0,j].set_title(sample);axes[0,j].set_xlabel('Node degree');axes[1,j].set_xlabel('Physical edge length (µm)');axes[1,j].set_xscale('log');axes[1,j].set_ylabel('Empirical CDF')
 axes[0,0].legend(fontsize=8);fig.suptitle('Delaunay graph diagnostics: first randomization shown; all graphs audited in CSV');fig.tight_layout();save(fig,'graph_diagnostics.png')
 # Descriptive maps for every candidate gene-panel row (mean over all relevant seeds).
 raw=pd.read_csv(OUT/'C3A_GENE_METRICS.csv');sample=SPLITS['splits']['test'][0];data=load(sample)
 for row in c[c.candidate].itertuples():
  gene=row.gene;panel=row.panel;idx=PANELS[panel].index(gene);truth=data['y'][:,ix(panel)[idx]];arrays=[truth]
  for control in ['A0','A3','A6']:
   recs=primary[(primary.panel==panel)&(primary.control==control)];vectors=[]
   for pp in recs.prediction_path:
    with np.load(ROOT/pp) as f:vectors.append(f['prediction'][:,idx])
   arrays.append(np.mean(vectors,axis=0))
  vmin=min(np.min(v) for v in arrays);vmax=max(np.max(v) for v in arrays);fig,axes=plt.subplots(1,4,figsize=(16,4))
  for ax,values,label in zip(axes,arrays,['Measured','Frozen local','True SAGE','Distance control']):
   im=ax.scatter(data['xy'][:,0],data['xy'][:,1],c=values,s=2,cmap='viridis',vmin=vmin,vmax=vmax,rasterized=True);ax.set_aspect('equal');ax.invert_yaxis();ax.set_title(label);ax.set_axis_off()
  fig.suptitle(f'{panel} / {gene}: primary-patient descriptive map, common color scale');fig.colorbar(im,ax=axes.tolist(),shrink=.7,label='Frozen normalized expression');save(fig,f'map_{panel}_{gene}.png')
if __name__=='__main__':main()
