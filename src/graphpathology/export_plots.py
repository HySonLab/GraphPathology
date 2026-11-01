"""Regenerate publication figures from exact exported CSVs, without inference."""
from pathlib import Path
import argparse,json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=None;M=None;D=None
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','savefig.dpi':300})
COL={'M0':'#7b7b7b','M1':'#087e8b','M3':'#b55628'};LAB={'M0':'H&E only','M1':'Sparse ST','M3':'H&E + ST'}
NAMES={1:'Fig1_study_overview',2:'Fig2_development_reconstruction',3:'Fig3_geometry_distance',4:'Fig4_HER2ST_transportability',5:'Fig5_fresh_external_Wu',6:'Fig6_second_replication',7:'Fig7_cross_cohort_synthesis'}
def load(name):return pd.read_csv(D/(name+'.csv'))
def finish(fig,axs,name,supp=False):
 trace=[]
 for ax in fig.axes:
  record={'title':ax.get_title(),'lines':[],'collections':[],'yscale':ax.get_yscale()}
  for line in ax.lines:record['lines'].append({'x':np.asarray(line.get_xdata()).tolist(),'y':np.asarray(line.get_ydata()).tolist()})
  for col in ax.collections:
   record['collections'].append({'offsets':np.asarray(col.get_offsets()).tolist(),'values':None if col.get_array() is None else np.ma.filled(col.get_array(),np.nan).tolist(),'clim':list(col.get_clim())})
  trace.append(record)
 (M/(name+'_plotted_values.json')).write_text(json.dumps(trace,indent=2)+'\n')
 for i,ax in enumerate(np.array(axs).ravel()):
  if name!='Fig3_geometry_distance' or i<2:ax.text(-.13,1.10,chr(65+i),transform=ax.transAxes,fontweight='bold',fontsize=12,va='top')
 fig.tight_layout(pad=2.4)
 if name=='Fig3_geometry_distance':
  for i,ax in enumerate(np.array(axs).ravel()):
   if i>=2:
    pos=ax.get_subplotspec().get_position(fig);fig.text(pos.x0-.025,pos.y1+.012,chr(65+i),fontweight='bold',fontsize=12,va='top')
 folder=M/('supp_figures' if supp else 'figures');fig.savefig(folder/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(folder/(name+'.svg'),bbox_inches='tight');svg=folder/(name+'.svg');svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n');fig.savefig(folder/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
def curves(ax,d,metric,title):
 for m,g in d.groupby('method'):
  if m not in COL:continue
  x=g.groupby('coverage')[metric].median();ax.plot(x.index*100,x.values,'o-',color=COL[m],label=LAB[m],ms=4)
 ax.set(xlabel='Measured locations (%)',ylabel={'median_pearson':'Median gene Pearson r','mean_mse':'Mean gene MSE','mean_mae':'Mean gene MAE'}[metric],title=title);ax.legend(frameon=False,fontsize=8)
def distance(ax,cohort):
 d=pd.read_csv(M/'supp_tables/SuppTable7_distance_bins.csv',low_memory=False);d=d[(d.cohort==cohort)&(d.coverage==.2)&(d.method=='M1')&d.regime.isin(['V1','V2'])]
 for reg,grp in d.groupby('regime'):
  p=grp.groupby(['patient','geometry','lower_um']).mse.mean().groupby(['patient','lower_um']).mean().reset_index();g=p.groupby('lower_um').mse.median();ax.plot(g.index,g.values,'o-',label=reg,color='#087e8b' if reg=='V1' else '#8d4ca6')
 ax.set(xlabel='Distance-bin lower bound (µm)',ylabel='Mean squared error',title='Nearest-anchor distance · ST only');ax.legend(frameon=False)
def fig1():
 fig,axs=plt.subplots(1,3,figsize=(12,3.6),gridspec_kw={'width_ratios':[1.35,1,1]});ax=axs[0];ax.axis('off')
 for y,txt,fc,ec in [(.90,'Measured ST + tissue coordinates\nOptional frozen H&E baseline','#e4f2f1','#087e8b'),(.48,'Frozen interpolation / fusion\nPredict hidden expression','#f4eee9','#b55628')]:ax.text(.5,y,txt,ha='center',va='center',bbox={'boxstyle':'round,pad=.5','fc':fc,'ec':ec},transform=ax.transAxes)
 for lo,hi in [(.57,.76),(.20,.34)]:ax.annotate('',xy=(.5,lo),xytext=(.5,hi),xycoords='axes fraction',arrowprops={'arrowstyle':'->'})
 ax.text(.5,.08,'Coverage · distance · geometry\nPatient-level correlation and error',ha='center',va='center',transform=ax.transAxes);d=load('Fig1_source')
 for ax,reg in zip(axs[1:],['V1','V2']):
  p=d[d.regime==reg];ax.scatter(p.x,p.y,c=np.where(p.observed,'#087e8b',np.where(p.forbidden,'#8d4ca6','#c9ced3')),s=35);ax.set_aspect('equal');ax.set_title(reg+(': distributed anchors' if reg=='V1' else ': contiguous missing region'));ax.set_xticks([]);ax.set_yticks([]);ax.set_xlabel('Illustrative geometry; teal = measured')
 finish(fig,axs,NAMES[1])
def fig2():
 d=load('development_patient');fig,axs=plt.subplots(2,2,figsize=(9,6.5))
 for j,reg in enumerate(['V1','V2']):
  for i,met in enumerate(['median_pearson','mean_mse']):curves(axs[i,j],d[d.regime==reg],met,reg+' · development primary patient')
 finish(fig,axs,NAMES[2])
def fig3():
 fig,axs=plt.subplots(2,2,figsize=(10,7));distance(axs[0,0],'Development');d=load('SFig6_source');d=d[(d.phase=='primary')&(d.panel=='S2-50')&(d.encoder=='resnet18_reference')&(d.role=='baseline')&d.control.isna()&(d.regime=='open')&d.strategy.isin(['A0','A1'])]
 for label,g in d.groupby('strategy'):
  q=g.groupby('nominal_budget').mean_mse.mean();axs[0,1].plot(q.index*100,q,'o-',label={'A0':'Random · 10-seed mean','A1':'Spatial FPS · deterministic'}[label])
 axs[0,1].set(xlabel='Measured locations (%)',ylabel='Mean gene MSE',title='Secondary C6 interpolation-only sampling');axs[0,1].legend(frameon=False,fontsize=8)
 d=load('Fig3_map_source');lim=d.absolute_error.quantile(.98)
 for ax,reg in zip(axs[1],['random','block']):
  p=d[d.regime==reg];ob=p[p.observed];h=p[~p.observed];s=ax.scatter(h.x_um,h.y_um,c=h.absolute_error,s=1,cmap='magma',vmin=0,vmax=lim,rasterized=True);ax.scatter(ob.x_um,ob.y_um,s=.4,c='#53d5d0',alpha=.4,rasterized=True);ax.set_aspect('equal');ax.invert_yaxis();ax.set(xlabel='x (µm)',ylabel='y (µm)',title=f'{reg} · {h.gene.iloc[0]} · fixed fusion');ax.title.set_position((.5,1.09));fig.colorbar(s,ax=ax,label='Absolute error (98th percentile cap)')
 finish(fig,axs,NAMES[3])
def fig4():
 d=load('Fig4_source');d=d[(d.coverage==.2)&(d.cohort=='HER2ST_C8A')];fig,axs=plt.subplots(2,2,figsize=(10,7))
 for j,reg in enumerate(['V1','V2']):
  for i,(metric,title) in enumerate([('median_pearson','Median gene Pearson r'),('mean_mse','Mean gene MSE · common C0 units')]):
   ax=axs[i,j]
   for k,trans in enumerate(['T0','T1']):
    x=[]
    for z,model in enumerate(['E0','E1','E2']):
     q=d[(d.regime==reg)&(d.method==trans+'_'+model)][metric];x.append(q.median());ax.scatter(np.repeat(z+(k-.5)*.24,len(q)),q,s=12,alpha=.4,color=['#7b7b7b','#b55628'][k])
    ax.plot(np.arange(3)+(k-.5)*.24,x,'o-',label='Zero-shot' if k==0 else 'T1 calibrated',color=['#7b7b7b','#b55628'][k])
   ax.set_xticks(range(3),['H&E','ST','H&E + ST']);ax.set(ylabel=title,title=reg+' · 8 patients · 40 genes');ax.legend(frameon=False)
   if metric=='mean_mse':ax.set_yscale('log')
 finish(fig,axs,NAMES[4])
def external(n,cohort):
 d=load('Fig5_source' if n==5 else 'Fig6_source');fig,axs=plt.subplots(2,3,figsize=(13,7));p=d[(d.regime=='V1')&(d.coverage==.2)];ax=axs[0,0]
 for j,m in enumerate(['M0','M1','M3']):
  q=p[p.method==m];ax.scatter(np.repeat(j,len(q)),q.median_pearson,s=20,c=COL[m],alpha=.65);ax.plot([j-.2,j+.2],[q.median_pearson.median()]*2,lw=3,color=COL[m])
 ax.set_xticks(range(3),['H&E','ST','H&E + ST']);ax.set(ylabel='Median gene Pearson r',title=f'20% V1 · {p.patient.nunique()} patients · 50 genes')
 curves(axs[0,1],d[d.regime=='V1'],'median_pearson','V1 coverage response');curves(axs[0,2],d[d.regime=='V1'],'mean_mse','V1 error response');axs[0,2].set_yscale('log')
 w=p.pivot(index='patient',columns='method',values='median_pearson');ax=axs[1,0]
 for j,(a,b) in enumerate([('M1','M0'),('M3','M1')]):
  val=w[a]-w[b];ax.scatter(np.repeat(j,len(val)),val,s=20,alpha=.6,c=COL[a]);ax.plot([j-.2,j+.2],[val.median()]*2,color=COL[a],lw=3)
 ax.axhline(0,color='k',lw=.7);ax.set_xticks([0,1],['ST − H&E','Fusion − ST']);ax.set(ylabel='Paired Δ median Pearson',title='Patient effects · 20% V1')
 ax=axs[1,1];p=d[(d.method=='M1')&(d.coverage==.2)].pivot(index='patient',columns='regime',values='median_pearson')
 for _,r in p.iterrows():ax.plot([0,1],r[['V1','V2']],'o-',alpha=.4,color='#087e8b',ms=3)
 ax.set_xticks([0,1],['V1','V2']);ax.set(ylabel='Median gene Pearson r',title='Paired patient protocol contrast');distance(axs[1,2],cohort);finish(fig,axs,NAMES[n])
def fig7():
 d=pd.read_csv(ROOT/'C10/CROSS_COHORT_SUMMARY.csv');cohorts=['Development','HER2ST zero-shot','HER2ST calibrated','Wu fresh validation','Bassiouni second replication'];fig,axs=plt.subplots(5,2,figsize=(10,18))
 for j,c in enumerate(cohorts):
  p=d[(d.cohort==c)&(d.coverage==.2)&(d.regime=='V1')].pivot(index='patient',columns='method',values='median_pearson');ax=axs[j,0]
  for k,(a,b) in enumerate([('M1','M0'),('M3','M1')]):
   val=p[a]-p[b];ax.scatter(np.repeat(k,len(val)),val,s=15,c=COL[a],alpha=.65);ax.plot([k-.2,k+.2],[val.median()]*2,c=COL[a],lw=2)
  ax.axhline(0,lw=.6,c='k');ax.set_xticks([0,1],['ST−H&E','Fusion−ST'],rotation=20);ax.set_title(c.replace(' ','\n',1)+f'\n{40 if c.startswith("HER2") else 50} genes · n={len(p)}',fontsize=9);ax.set_ylabel('Paired Δ median Pearson')
  p=d[(d.cohort==c)&(d.method=='M1')&(d.regime=='V1')&d.coverage.isin([.1,.2,.4])];ax=axs[j,1]
  for _,g in p.groupby('patient'):ax.plot(g.coverage*100,g.median_pearson,'o-',c='#087e8b',alpha=.35,ms=3)
  ax.set(xlabel='Measured locations (%)',ylabel='Median gene Pearson r');ax.set_xticks([10,20,40])
 finish(fig,axs,NAMES[7])
def grouped(ax,d,group,metric,ylabel):
 vals=[];labs=[]
 for lab,g in d.groupby(group,dropna=False):
  q=g[metric].dropna()
  if len(q):vals.append(q);labs.append(str(lab))
 for i,q in enumerate(vals):ax.scatter(np.repeat(i,len(q)),q,s=15,alpha=.5,color='#087e8b');ax.plot([i-.2,i+.2],[q.mean()]*2,c='#b55628',lw=2)
 ax.set_xticks(range(len(labs)),labs,rotation=25,ha='right');ax.set_ylabel(ylabel)
def supplement(n):
 fig,axs=plt.subplots(1,3 if n==10 else 2,figsize=(15 if n==10 else 11,4));d=load(f'SFig{n}_source');title=''
 if n in [1,2]:
  d=d[d.panel=='S2-50'];a=d.copy();a['label']=a.model
  if n==2:
   base=load('SFig1_source');base=base[(base.panel=='S2-50')&(base.model=='B2')].copy();base['label']='Local MLP';a=pd.concat([base,a])
  grouped(axs[0],a,'label','median_pearson','Median gene Pearson r');c=load(f'SFig{n}_controls_source');c=c[c.panel=='S2-50'].copy();c['label']=c.phase.astype(str)+' / '+c.model.astype(str);grouped(axs[1],c,'label','mean_mse','Mean gene MSE');title=f'C{n} primary patient · 50 genes · model seeds / graph controls'
 elif n==3:
  d=load('SFig3_panel_source');d=d[(d.phase=='primary')&(d.panel=='S2-50')];p=d.groupby(['control','graph_randomization_seed'])[['median_pearson','mean_mse']].mean().reset_index();p['label']=p.control.map({'A0':'Local','A1':'Self-only','A2':'Mean pool','A3':'True SAGE','A4':'Vertex perm.','A5':'Degree rewire','A6':'Degree+length'});grouped(axs[0],p,'label','median_pearson','Median gene Pearson r');grouped(axs[1],p,'label','mean_mse','Mean gene MSE');title='C3 primary patient · 50 genes · training-seed means; randomized controls: five graph seeds'
 elif n==4:
  d=d[(d.panel=='S2-50')&(d.phase=='primary')&(d.model=='D1')].copy();d['label']=d.encoder.map({'pathology_encoder_1':'Phikon','resnet18_reference':'ResNet18'});grouped(axs[0],d,'label','median_pearson','Median gene Pearson r');grouped(axs[1],d,'label','mean_mse','Mean gene MSE');title='C4 primary patient · local MLP · 50 genes · three training seeds'
 elif n==5:
  d=d[(d.panel=='S2-50')&(d.encoder=='resnet18_reference')&(d.nominal_coverage==.2)].copy();d['label']=d.control.astype(str)+' / '+d.model.astype(str);grouped(axs[0],d,'label','median_pearson','Median gene Pearson r');grouped(axs[1],d,'label','mean_mse','Mean gene MSE');title='Sparse-reconstruction controls · 20% · 50 genes · one random mask, three permutation seeds'
 elif n==6:
  d=d[(d.panel=='S2-50')&(d.phase=='primary')&(d.encoder=='resnet18_reference')&(d.nominal_budget==.2)&(d.role=='baseline')&d.control.isna()]
  for ax,reg in zip(axs,['open','contiguous']):grouped(ax,d[d.regime==reg],'strategy','median_pearson','Median gene Pearson r');ax.set_title(reg)
  title='C6 interpolation baseline · A0 random, A1 spatial, A2 morphology, A3 joint · 20%'
 elif n==7:
  d=d[(d.panel=='S2-50')&(d.nominal_budget==.2)&(d.strategy=='random')&d.model.isin(['E0','E1','E2'])]
  for ax,reg in zip(axs,['V1','V2']):grouped(ax,d[d.regime==reg],'model','median_pearson','Median gene Pearson r');ax.set_title(reg)
  title='C7 zero-shot · 8 patients · common 40 genes · E0 H&E, E1 ST, E2 fusion'
 elif n==8:
  d=d[d.budget==.2];grouped(axs[0],d,'regime','gain_pearson','Selected fusion − ST: Δ median Pearson');grouped(axs[1],d,'regime','relative_mse_reduction','Selected fusion: relative MSE reduction');title='C8 selected fusion · 8 transportability-development patients · 40 genes · 20%'
 elif n in [9,10]:
  d=d[d.method.isin(['M1','M3'])].copy();d['label']=d.control+' / '+d.method;grouped(axs[0],d,'label','delta_pearson','Control − matched clean: Δ Pearson');grouped(axs[1],d,'label','mse_reduction','Clean − control: mean MSE');title=f'C{n} negative controls · one patient · 20% · three seeds × two regimes'
  if n==10:
   q=load('SFig10_QC_source');axs[2].scatter(range(len(q)),q.native_patch_pixels,s=20,color='#087e8b');axs[2].set(xlabel='Patient (source order)',ylabel='Native pixels per 112-µm patch',title='Public image resolution');axs[2].set_xticks(range(len(q)),q.patient_id.str.replace('BASS_',''),rotation=90);title+='; H&E QC: all 22 patients'
 fig.suptitle(title,fontsize=10);finish(fig,axs,f'SFig{n}_checkpoint_diagnostics',True)
