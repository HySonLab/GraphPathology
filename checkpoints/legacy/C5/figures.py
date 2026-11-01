"""Predeclared geometry-specific coverage, distance, patient and spatial diagnostics."""
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
LABELS=['B0','B1/R','B1/P','B2a','B2b','B3/R','B3/P','B4/R','B4/P'];COLORS=dict(zip(LABELS,['#999999','#537aac','#95b4d7','#ab6c28','#d49e59','#257c58','#69b595','#775a9c','#b395cd']))
TITLES={'random':'Random observation: interpolation','block':'Contiguous hidden region: extrapolation','capture_window':'Contiguous observed region: local-to-global extrapolation'}
def annotate(d):
 d=d.copy();d['label']=d['model']+d['encoder'].map({'resnet18_reference':'/R','pathology_encoder_1':'/P','none':''}).fillna('');return d

def save(fig,name):fig.savefig(OUT/'figures'/name,dpi=155,bbox_inches='tight');plt.close(fig)
def main():
 curves=annotate(pd.read_csv(OUT/'C5_COVERAGE_CURVES.csv'));primary=curves[curves.cohort=='primary'];distance=annotate(pd.read_csv(OUT/'C5_DISTANCE_TO_ANCHOR.csv'));distance=distance[(distance.cohort=='primary')&distance.control.isna()&(distance.nominal_coverage==.2)];panels=list(PANELS)
 for regime in REGIMES:
  fig,axes=plt.subplots(2,3,figsize=(16,9))
  for j,panel in enumerate(panels):
   d=primary[(primary.regime==regime)&(primary.panel==panel)]
   for ax,metric,ylabel in zip(axes[:,j],['median_pearson','mean_mse'],['Median gene Pearson','Mean gene MSE']):
    for label in LABELS:
     q=d[d.label==label].sort_values('realized_coverage');ax.plot(100*q.realized_coverage,q[metric],'o',ms=3,lw=1.6,color=COLORS[label],label=label,ls='--' if label.endswith('/P') else '-');
    ax.set_title(panel);ax.set_xlabel('Realized observed locations (%)');ax.set_ylabel(ylabel);ax.grid(alpha=.2)
  handles,legend_labels=axes[1,0].get_legend_handles_labels();fig.legend(handles,legend_labels,loc='lower center',ncol=9,fontsize=9);fig.suptitle(TITLES[regime]+'\nPrimary patient; mean over three frozen masks; R = ResNet18, P = Phikon');fig.tight_layout(rect=[0,.05,1,.93]);save(fig,'coverage_'+regime+'.png')
  fig,axes=plt.subplots(2,3,figsize=(16,9));ticks=[0,150,300,600,1200,2400];ticklabels=['0–150','150–300','300–600','600–1200','1200–2400','≥2400']
  for j,panel in enumerate(panels):
   d=distance[(distance.regime==regime)&(distance.panel==panel)].groupby(['label','lower_um'])[['mse','mae','median_pearson','n_hidden']].mean().reset_index()
   for ax,metric,ylabel in zip(axes[:,j],['mse','median_pearson'],['Mean gene MSE','Median gene Pearson (≥20 locations)']):
    for label in LABELS:
     q=d[d.label==label].set_index('lower_um').reindex(ticks);ax.plot(range(6),q[metric],'o',ms=3,color=COLORS[label],label=label,ls='--' if label.endswith('/P') else '-')
    ax.set_xticks(range(6),ticklabels,rotation=35,ha='right',fontsize=8);ax.set_xlabel('Distance to nearest ST anchor (estimated µm)');ax.set_title(panel);ax.set_ylabel(ylabel);ax.grid(alpha=.2)
  handles,legend_labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,legend_labels,loc='lower center',ncol=9,fontsize=9);fig.suptitle(TITLES[regime]+'\nNominal 20% coverage; equal-mask bin means; gaps = no available locations');fig.tight_layout(rect=[0,.05,1,.93]);save(fig,'distance_'+regime+'.png')
 # Matched primary 20% endpoints make B1/B2/B3/B4 and encoder comparisons explicit.
 q=primary[primary.nominal_coverage==.2];fig,axes=plt.subplots(2,3,figsize=(17,9))
 for j,regime in enumerate(REGIMES):
  for i,(panel,metric) in enumerate([('S2-50','median_pearson'),('S2-280','mean_mse')]):
   ax=axes[i,j];d=q[(q.regime==regime)&(q.panel==panel)].set_index('label').reindex(LABELS);ax.bar(range(9),d[metric],color=[COLORS[l] for l in LABELS]);ax.set_xticks(range(9),LABELS,rotation=45,ha='right');ax.set_title(f'{regime} / {panel}');ax.set_ylabel(metric)
 fig.suptitle('Primary comparison at nominal 20% coverage; three-mask means; B0 correlation undefined');fig.tight_layout(rect=[0,0,1,.94]);save(fig,'model_encoder_comparison.png')
 # Full frozen gene order; no outcome-ranked example selection.
 genes=pd.read_csv(OUT/'C5_GENE_COVERAGE.csv.gz');selection=json.loads((OUT/'SELECTION_LOCK.json').read_text())
 for regime in REGIMES:
  fig,axes=plt.subplots(1,4,figsize=(15,12),sharey=True)
  for ax,(model,enc,title) in zip(axes,[('B1',ENCODERS[0],'H&E only'),(selection['simple']['S2-50'][regime]['best_family'],'none','Selected interpolation'),('B3',ENCODERS[0],'Simple fusion'),('B4',ENCODERS[0],'Learned propagation')]):
   d=genes[(genes.cohort=='primary')&(genes.panel=='S2-50')&(genes.regime==regime)&(genes.model==model)&(genes.encoder==enc)];v=d.pivot(index='gene',columns='nominal_coverage',values='pearson').reindex(index=PANELS['S2-50'],columns=[.05,.1,.2,.4,.6,.8]);im=ax.imshow(v,vmin=-.2,vmax=1,cmap='viridis',aspect='auto');ax.set_xticks(range(6),['5','10','20','40','60','80']);ax.set_xlabel('Nominal coverage (%)');ax.set_title(title);ax.set_yticks(range(50),PANELS['S2-50'],fontsize=6)
  fig.colorbar(im,ax=axes.ravel().tolist(),fraction=.02,label='Mean per-gene Pearson');fig.suptitle(TITLES[regime]+'\nS2-50, ResNet18; frozen gene order, no inferred minimum-coverage threshold');save(fig,'gene_coverage_'+regime+'.png')
 for panel in ['S2-50','S2-280']:
  fig,axes=plt.subplots(2,3,figsize=(17,8));mse_limit=float(curves[(curves.cohort!='primary')&(curves.panel==panel)&(curves.nominal_coverage==.2)].mean_mse.max())
  for j,regime in enumerate(REGIMES):
   d=curves[(curves.cohort!='primary')&(curves.panel==panel)&(curves.regime==regime)&(curves.nominal_coverage==.2)]
   for ax,metric in zip(axes[:,j],['median_pearson','mean_mse']):
    v=d.pivot(index='sample_id',columns='label',values=metric).reindex(index=IDS,columns=LABELS);im=ax.imshow(v,aspect='auto',cmap='viridis' if metric=='median_pearson' else 'magma_r',vmin=-.2 if metric=='median_pearson' else 0,vmax=1 if metric=='median_pearson' else mse_limit);ax.set_xticks(range(9),LABELS,rotation=45,ha='right');ax.set_yticks(range(4),IDS);ax.set_title(regime+' / '+metric);fig.colorbar(im,ax=ax,fraction=.035)
  fig.suptitle(panel+' conditional patient sensitivity at nominal 20% coverage; not unbiased cross-validation; shared scales across geometries');fig.tight_layout(rect=[0,0,1,.94]);save(fig,'patients_'+panel+'.png')
 # Fixed test sample, gene, coverage, replicate and model chosen in the spec.
 rows=[json.loads(f.read_text()) for f in (OUT/'logs/runs').glob('*.json')];sample='TENX99';panel='S2-50';gene=PANELS[panel][0]
 with np.load(ROOT/'C4/embeddings/resnet18_reference'/f'{sample}.npz') as f:xy=f['xy_um']
 truth=labels(sample,panel)[:,0];fig,axes=plt.subplots(3,4,figsize=(16,12));maskidx={r['mask_id']:r for r in json.loads((OUT/'configs/masks.json').read_text())}
 for row,regime in enumerate(REGIMES):
  r=next(r for r in rows if r['cohort']=='primary' and r['panel']==panel and r['model']=='B3' and r['encoder']==ENCODERS[0] and r['regime']==regime and r['nominal_coverage']==.2 and r['replicate']==0 and r['control'] is None);obs=mask_array(maskidx[r['mask_id']]);hidden=~obs
  with np.load(ROOT/r['prediction_path']) as f:pred=f['prediction'][:,0]
  lo,hi=np.quantile(truth,[.01,.99]);datasets=[(np.ones(len(xy),bool),obs.astype(float),'Observed mask',None,None),(np.ones(len(xy),bool),truth,'Ground truth',lo,hi),(hidden,pred,'Hidden prediction',lo,hi),(hidden,abs(truth[hidden]-pred),'Hidden absolute error',0,None)]
  for ax,(take,values,title,vmin,vmax) in zip(axes[row],datasets):
   im=ax.scatter(xy[take,0],xy[take,1],c=values,s=2,cmap='coolwarm' if title=='Observed mask' else 'viridis',vmin=vmin,vmax=vmax,rasterized=True);ax.set_title(regime+' / '+title,fontsize=9);ax.set_aspect('equal');ax.set_xlim(xy[:,0].min()-500,xy[:,0].max()+500);ax.set_ylim(xy[:,1].max()+500,xy[:,1].min()-500);ax.set_xlabel('Estimated µm');cb=fig.colorbar(im,ax=ax,fraction=.04)
   if title=='Observed mask':cb.set_label('0 hidden; 1 observed',fontsize=7)
 fig.suptitle(f'Predeclared maps: TENX99 / {gene} / nominal 20% / replicate 0 / ResNet18 B3\nAll maps share spatial bounds; truth/prediction share color limits; observed sites are not credited as reconstructions');fig.tight_layout(rect=[0,0,1,.94]);save(fig,'spatial_examples.png')
if __name__=='__main__':main()
