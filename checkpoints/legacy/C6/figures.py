"""Scientific figures with fixed examples, raw random outcomes and separate geometries."""
from support import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy.spatial import ConvexHull
import h5py
STRATEGIES=['A0','A1','A2','A3'];NAMES={'A0':'Random','A1':'Spatial FPS','A2':'Morphology FPS','A3':'Joint FPS'};COLORS={'A0':'#777777','A1':'#176f9d','A2':'#bd6a25','A3':'#568540'}
def save(fig,name):fig.savefig(A/'figures'/name,dpi=145,bbox_inches='tight');plt.close(fig)
def main():
 curves=pd.read_csv(A/'C6A_BUDGET_CURVES.csv');runs=pd.read_csv(A/'C6A_PANEL_METRICS.csv');effects=pd.read_csv(A/'C6A_EFFECTS.csv');distance=pd.read_csv(A/'C6A_DISTANCE_METRICS.csv');geometry=pd.read_csv(A/'C6A_GEOMETRY_DIAGNOSTICS.csv');primary=curves[curves.cohort=='primary']
 for enc in ENCODERS:
  short='ResNet18' if enc==ENCODERS[0] else 'Phikon'
  for regime in REGIMES:
   for role in ['primary','baseline']:
    fig,axes=plt.subplots(2,2,figsize=(11,8))
    for col,panel in enumerate(PANELS):
     d=primary[(primary.encoder==enc)&(primary.regime==regime)&(primary.role==role)&(primary.panel==panel)]
     for ax,metric,sd,label in zip(axes[:,col],['median_pearson','mean_mse'],['pearson_sd','mse_sd'],['Median gene Pearson','Mean gene MSE']):
      for strategy in STRATEGIES:
       q=d[d.strategy==strategy].sort_values('realized_budget');x=100*q.realized_budget;ax.plot(x,q[metric],'o-',label=NAMES[strategy],color=COLORS[strategy],ms=4)
       if strategy=='A0':ax.fill_between(x,q[metric]-q[sd],q[metric]+q[sd],alpha=.15,color=COLORS[strategy])
      ax.set_title(panel);ax.set_xlabel('Realized measured locations (%)');ax.set_ylabel(label);ax.grid(alpha=.2)
    fig.legend(*axes[0,0].get_legend_handles_labels(),loc='lower center',ncol=4);fig.suptitle(f'{short} / {regime} / {role} reconstruction\nPrimary patient; random mean ± SD over 10 seeds; deterministic policies have one trace');fig.tight_layout(rect=[0,.05,1,.92]);save(fig,f'curves_{short}_{regime}_{role}.png')
  for role in ['primary','baseline']:
   fig,axes=plt.subplots(2,2,figsize=(12,8))
   for row,regime in enumerate(REGIMES):
    for col,panel in enumerate(PANELS):
     ax=axes[row,col];d=distance[(distance.cohort=='primary')&(distance.encoder==enc)&(distance.role==role)&(distance.panel==panel)&(distance.regime==regime)&(distance.nominal_budget==.2)&distance.control.isna()]
     for strategy in STRATEGIES:
      q=d[d.strategy==strategy].groupby('lower_um').mse.mean().reindex([0,150,300,600,1200,2400]);ax.plot(range(6),q,'o-',label=NAMES[strategy],color=COLORS[strategy])
     ax.set_xticks(range(6),['0–150','150–300','300–600','600–1200','1200–2400','≥2400'],rotation=35,ha='right',fontsize=8);ax.set_title(f'{panel} / {regime}');ax.set_ylabel('Mean gene MSE');ax.set_xlabel('Nearest-anchor distance (estimated µm)');ax.grid(alpha=.2)
   fig.legend(*axes[0,0].get_legend_handles_labels(),loc='lower center',ncol=4);fig.suptitle(f'{short} / {role} / nominal 20% / primary patient\nRandom bin means over available seeds; gaps are empty bins');fig.tight_layout(rect=[0,.05,1,.92]);save(fig,f'distance_{short}_{role}.png')
  fig,axes=plt.subplots(2,2,figsize=(11,8))
  for col,regime in enumerate(REGIMES):
   d=geometry[(geometry.cohort=='primary')&(geometry.encoder==enc)&(geometry.regime==regime)&geometry.control.isna()]
   for ax,metric,label in zip(axes[:,col],['p90_anchor_distance_um','morphology_mean_nearest_distance'],['90th percentile nearest-anchor distance (µm)','Mean nearest-anchor morphology distance']):
    for strategy in STRATEGIES:
     q=d[d.strategy==strategy].groupby('realized_budget')[metric].mean();ax.plot(100*q.index,q,'o-',label=NAMES[strategy],color=COLORS[strategy])
    ax.set_xlabel('Realized measured locations (%)');ax.set_title(regime);ax.set_ylabel(label);ax.grid(alpha=.2)
  fig.legend(*axes[0,0].get_legend_handles_labels(),loc='lower center',ncol=4);fig.suptitle(f'{short}: spatial and morphology coverage diagnostics\nDistances over scored hidden sites; diagnostics alone do not establish reconstruction gain');fig.tight_layout(rect=[0,.05,1,.92]);save(fig,f'geometry_{short}.png')
 # All random observations; no seed selected by outcome.
 for regime in REGIMES:
  fig,axes=plt.subplots(2,4,figsize=(16,8))
  for col,(panel,budget) in enumerate([(p,b) for p in PANELS for b in [.1,.2]]):
   d=runs[(runs.cohort=='primary')&(runs.encoder==ENCODERS[0])&(runs.role=='primary')&(runs.regime==regime)&(runs.panel==panel)&(runs.nominal_budget==budget)&runs.control.isna()]
   for ax,metric in zip(axes[:,col],['median_pearson','mean_mse']):
    for j,strategy in enumerate(STRATEGIES):
     q=d[d.strategy==strategy].sort_values('seed');jitter=np.linspace(-.13,.13,len(q)) if len(q)>1 else np.zeros(1);ax.scatter(j+jitter,q[metric],color=COLORS[strategy],s=20,alpha=.8);ax.plot([j-.2,j+.2],[q[metric].mean()]*2,color='black',lw=1)
    ax.set_xticks(range(4),STRATEGIES);ax.set_title(f'{panel} / {budget:.0%}');ax.set_ylabel(metric);ax.grid(axis='y',alpha=.2)
  fig.suptitle(f'{regime}: every random seed versus deterministic strategies\nPrimary patient / ResNet18 / primary reconstruction; bars mark means');fig.tight_layout(rect=[0,0,1,.92]);save(fig,f'random_distribution_{regime}.png')
 for panel in PANELS:
  fig,axes=plt.subplots(4,2,figsize=(16,15));colkeys=[(s,b) for s in ['A1','A2','A3'] for b in BUDGETS]
  for row,(enc,regime) in enumerate([(e,g) for e in ENCODERS for g in REGIMES]):
   d=effects[(effects.cohort!='primary')&(effects.encoder==enc)&(effects.panel==panel)&(effects.role=='primary')&(effects.regime==regime)]
   for ax,metric in zip(axes[row],['Pearson_gain','relative_MSE_reduction']):
    v=d.pivot(index='sample_id',columns=['strategy','nominal_budget'],values=metric).reindex(index=IDS,columns=pd.MultiIndex.from_tuples(colkeys));limit=max(.02,float(np.nanmax(abs(v.to_numpy()))));im=ax.imshow(v,aspect='auto',cmap='RdBu',vmin=-limit,vmax=limit);ax.set_yticks(range(4),IDS);ax.set_xticks(range(12),[f'{s}\n{b:.0%}' for s,b in colkeys],fontsize=8);ax.set_title(f'{"ResNet18" if enc==ENCODERS[0] else "Phikon"} / {regime} / {metric}');fig.colorbar(im,ax=ax,fraction=.025)
  fig.suptitle(f'{panel}: conditional patient sensitivity / positive = improvement over random mean\nNot unbiased external validation; each row uses its indicated symmetric color scale');fig.tight_layout(rect=[0,0,1,.95]);save(fig,f'patients_{panel}.png')
 # Fixed H&E overlay example from the pretest protocol, with barcode-joined pixel coordinates.
 lock=json.loads((A/'TRAJECTORY_LOCK.json').read_text());policies=[json.loads((A/'policies'/f'{pid}.json').read_text()) for pid in lock['policies']];sample='TENX99'
 with h5py.File(ROOT/'data/raw/IDC/adata'/f'{sample}.h5ad') as f:thumb=f['uns/spatial/ST/images/downscaled_fullres'][:];scale=float(f['uns/spatial/ST/scalefactors/tissue_downscaled_fullres_scalef'][()])
 with h5py.File(ROOT/'data/raw/IDC/patches'/f'{sample}.h5') as f:xy=f['coords'][:][:,::-1]*scale;barcodes=np.asarray(f['barcode'][:]).reshape(-1).astype(str)
 assert np.array_equal(barcodes,c5.image_data(ENCODERS[0],sample)['barcodes']);fig,axes=plt.subplots(2,4,figsize=(15,13))
 for row,regime in enumerate(REGIMES):
  for col,strategy in enumerate(STRATEGIES):
   ax=axes[row,col];r=next(r for r in policies if r['cohort']=='primary' and r['encoder']==ENCODERS[0] and r['regime']==regime and r['strategy']==strategy and r['seed']==0 and r['control'] is None)
   with np.load(ROOT/r['policy_path']) as f:selected=f['order'][:count(len(xy),.1)];forbidden=f['forbidden']
   ax.imshow(thumb);ax.scatter(xy[selected,0],xy[selected,1],s=1,c='#ffff00',alpha=.8)
   if forbidden.any():points=xy[forbidden];hull=ConvexHull(points);ax.add_patch(Polygon(points[hull.vertices],facecolor='#00ffff',alpha=.15,edgecolor='#00ffff',linewidth=2))
   ax.set_xlim(0,thumb.shape[1]);ax.set_ylim(thumb.shape[0],0);ax.set_title(f'{regime} / {NAMES[strategy]}\n{len(selected)} measured');ax.axis('off')
 fig.suptitle('Predeclared TENX99 H&E overlay / ResNet18 / 10% budget / A0 seed 0\nYellow: measured locations. Cyan: unavailable stress region; no anchors enter it.');fig.tight_layout(rect=[0,0,1,.94]);save(fig,'selected_locations_HE.png')
 fig,axes=plt.subplots(1,2,figsize=(12,5))
 for ax,regime in zip(axes,REGIMES):
  for strategy in STRATEGIES:
   r=runs[(runs.cohort=='primary')&(runs.encoder==ENCODERS[0])&(runs.panel=='S2-50')&(runs.role=='primary')&(runs.regime==regime)&(runs.nominal_budget==.1)&(runs.strategy==strategy)&(runs.seed==0)&runs.control.isna()].iloc[0]
   with np.load(ROOT/r.prediction_path) as f:d=np.sort(f['distance_um'])
   ax.plot(d,np.arange(1,len(d)+1)/len(d),label=NAMES[strategy],color=COLORS[strategy])
  ax.set_title(regime);ax.set_xlabel('Nearest-anchor distance (estimated µm)');ax.set_ylabel('Fraction of scored hidden locations');ax.grid(alpha=.2)
 fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=4);fig.suptitle('Fixed TENX99 / ResNet18 / 10% budget / A0 seed 0 distance distributions');fig.tight_layout(rect=[0,.08,1,.93]);save(fig,'distance_distributions.png')
if __name__=='__main__':main()
