"""C4A correlation/error, negative-control and patient-sensitivity diagnostics."""
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def save(fig,name):fig.savefig(A/'figures'/name,dpi=170,bbox_inches='tight');plt.close(fig)
def main():
 p=pd.read_csv(A/'C4A_PANEL_METRICS.csv');s=pd.read_csv(A/'C4A_GENE_SUPPORT.csv');primary=p[p.phase=='primary'];panels=list(PANELS);encoders=['resnet18_reference','pathology_encoder_1'];colors=['#777777','#267b60']
 fig,axes=plt.subplots(2,4,figsize=(17,7))
 for j,panel in enumerate(panels):
  for ax,metric in zip(axes[:,j],['median_pearson','mean_mse']):
   for k,model in enumerate(['D0','D1']):
    for i,encoder in enumerate(encoders):
     vals=primary[(primary.panel==panel)&(primary.model==model)&(primary.encoder==encoder)][metric];pos=3*k+i;ax.bar(pos,vals.mean(),color=colors[i],alpha=.7);ax.scatter(pos+np.linspace(-.1,.1,len(vals)),vals,color='black',s=18,zorder=3)
   ax.set_xticks([0,1,3,4],['ResNet\nRidge','Phikon\nRidge','ResNet\nMLP','Phikon\nMLP'],fontsize=8);ax.set_title(panel);ax.set_ylabel('Median gene Pearson' if metric=='median_pearson' else 'Mean gene MSE')
 fig.suptitle('Primary test patient; fixed targets, matched tuning budget; dots = training seeds (Ridge deterministic)');fig.tight_layout();save(fig,'primary_comparison.png')
 fig,axes=plt.subplots(2,3,figsize=(15,8))
 for i,model in enumerate(['D0','D1']):
  for j,panel in enumerate(panels[1:]):
   ax=axes[i,j];d=s[(s.phase=='primary')&(s.panel==panel)&(s.model==model)];ax.scatter(d.delta_encoder,d.delta_mse,c=np.where(d.strong_gene,'#d89b30','#267b60'),s=18,alpha=.7);ax.axhline(0,color='black',lw=.8);ax.axvline(0,color='black',lw=.8);ax.set_title(f'{panel} / {model}');ax.set_xlabel('Phikon − ResNet gene Pearson');ax.set_ylabel('ResNet MSE − Phikon MSE')
 fig.suptitle('Encoder gains: upper-right improves both; gold = strong descriptive gene criteria');fig.tight_layout();save(fig,'gene_gain_vs_error.png')
 fig,axes=plt.subplots(1,2,figsize=(12,4));fold=s[(s.phase=='lopo')&(s.model=='D1')]
 for ax,metric,title in zip(axes,['delta_encoder','delta_mse'],['Median paired gene Pearson gain','Mean per-gene MSE reduction']):
  t=fold.pivot_table(index='sample_id',columns='panel',values=metric,aggfunc='median' if metric=='delta_encoder' else 'mean').reindex(index=IDS,columns=panels);v=t.to_numpy();bound=max(abs(v).max(),1e-6);im=ax.imshow(v,vmin=-bound,vmax=bound,cmap='RdBu',aspect='auto');ax.set_xticks(range(4),panels);ax.set_yticks(range(4),IDS);ax.set_title(title)
  for i in range(4):
   for j in range(4):ax.text(j,i,f'{v[i,j]:+.3f}',ha='center',va='center',fontsize=9,color='white' if abs(v[i,j])/bound>.5 else 'black')
  fig.colorbar(im,ax=ax)
 fig.suptitle('MLP conditional patient sensitivity; positive means Phikon improves; not unbiased cross-validation');fig.tight_layout();save(fig,'patient_sensitivity.png')
 fig,axes=plt.subplots(1,2,figsize=(10,4))
 for ax,model in zip(axes,['D0','D1']):
  groups=[primary[(primary.encoder=='pathology_encoder_1')&(primary.model==model)&(primary.panel=='S2-50')].median_pearson]
  groups += [p[(p.phase=='negative')&(p.model==model)&(p.control==c)].median_pearson for c in ['NC1','NC2']]
  for i,vals in enumerate(groups):ax.bar(i,vals.mean(),alpha=.6,color=['#267b60','#c56c71','#7297bb'][i]);ax.scatter(i+np.linspace(-.12,.12,len(vals)),vals,color='black')
  ax.axhline(.2,color='#aa4444',ls='--',lw=.8);ax.set_xticks(range(3),['Paired','Labels shuffled','Locations shuffled'],rotation=15);ax.set_title(model);ax.set_ylabel('Median gene Pearson')
 fig.suptitle('Phikon S2-50 negative controls; held-out targets/features unchanged');fig.tight_layout();save(fig,'negative_controls.png')
 fig,axes=plt.subplots(1,3,figsize=(15,4))
 for ax,panel in zip(axes,panels[1:]):
  for model,color in zip(['D0','D1'],['#7297bb','#267b60']):
   d=s[(s.phase=='primary')&(s.panel==panel)&(s.model==model)];ax.hist(d.delta_encoder,bins=25,histtype='step',lw=1.5,label=model,color=color)
  ax.axvline(0,color='black',lw=.8);ax.set_title(panel);ax.set_xlabel('Per-gene encoder Pearson gain');ax.set_ylabel('Genes')
 axes[0].legend();fig.suptitle('All genes retained, including negative encoder gains');fig.tight_layout();save(fig,'gain_distributions.png')
if __name__=='__main__':main()
