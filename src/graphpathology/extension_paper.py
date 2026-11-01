"""Separate journal-only extension package assembled strictly from completed exports."""
import json,re,shutil,subprocess
from pathlib import Path
import numpy as np
import pandas as pd
from .assets import sha
from .extension_models import VARIANTS
from .journal_paper import table,esc

NAMES={'gps_full':'GPS full','gps_local':'GPS local','gps_no_relative':'GPS no distances','gps_self':'GPS self','gps_shuffled':'GPS shuffled','gps_degree_length':'GPS rewired','distance_fixed':'Distance fixed','egnn_latent':'EGNN latent','distance_disabled':'Distance disabled','egnn_self':'EGNN self','egnn_shuffled':'EGNN shuffled','egnn_degree_length':'EGNN rewired','matched_local':'Local MLP','matched_mean':'Mean','matched_gcn':'GCN','matched_sage':'GraphSAGE'}
COHORTS=['Development','Wu','Bassiouni'];ROUTES=['image','fixed_20pct_image','exploratory_residual_anchor']


def clean_prose(text):
 protected={}
 pattern=r'\\cite\{[^}]*\}|TENX95|TENX99|NCBI783|NCBI785|ResNet18'
 def protect(match):
  key='PROTECTED'+chr(65+len(protected));protected[key]=match.group();return key
 text=re.sub(pattern,protect,text)
 text=re.sub(r'([A-Za-z]{2,})([0-9])',r'\1 \2',text)
 text=re.sub(r'([0-9])([A-Za-z]{2,})',r'\1 \2',text)
 for key,value in protected.items():text=text.replace(key,value)
 return text

def figure(path,caption,label,width=r'\textwidth'):
 return '\n'.join([r'\begin{figure}[!htbp]',r'\centering\includegraphics[width='+width+']{'+path+'}',r'\caption{'+caption+'}',r'\label{'+label+'}',r'\end{figure}',r'\FloatBarrier',''])

def plot_all(root,package):
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 from matplotlib.patches import FancyBboxPatch
 plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
 base=root/'experiments/journal_revision/transformer_egnn';src=package/'source_data/extension';out=package/'figures/extension';out.mkdir(parents=True,exist_ok=True);figures=[]
 def save(fig,name):
  fig.savefig(out/(name+'.pdf'),bbox_inches='tight');fig.savefig(out/(name+'.png'),dpi=160,bbox_inches='tight');plt.close(fig);figures.append(name)
 # Architecture diagram describes real code, not a learned biological interpretation.
 fig,ax=plt.subplots(figsize=(10,6));ax.set(xlim=(0,10),ylim=(0,6));ax.axis('off')
 def box(x,y,w,h,text,color):
  ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.08',facecolor=color,edgecolor='#283747'));ax.text(x+w/2,y+h/2,text,ha='center',va='center',fontsize=13)
 def arrow(a,b):ax.annotate('',xy=b,xytext=a,arrowprops={'arrowstyle':'->','color':'#283747','lw':1.5})
 box(.2,4.6,2.4,1,'Frozen 512-D images\nOriginal physical xy','#eaf2f8');box(3.3,4.6,2.8,1,'Two 64-D blocks\nLocal RBF-gated mean\n+32-landmark attention','#e8f8f5');box(6.8,4.6,2.9,1,'50 scalar outputs g\nNo molecular labels\nin image network','#fef9e7');arrow((2.7,5.1),(3.2,5.1));arrow((6.2,5.1),(6.7,5.1))
 box(3.3,2.6,2.8,1.2,'Alternative:2D EGNN\nRelative squared distances\nPrivate latent p updates\nFixed-edge, two64-D blocks','#f5eef8');arrow((1.4,4.5),(3.2,3.2));arrow((6.2,3.2),(8.2,4.5))
 box(.2,.5,2.5,1.4,'Observed labels y[O]\nOriginal xy + fixed O/H\nIDW k16, p2','#fdebd0');box(3.5,.5,5.9,1.4,'Three separately evaluated routes\nImage: g[H]\nFixed:0.2 g[H]+0.8 IDW(y[O])\nResidual:g[H]+IDW(y[O]-g[O])','#f2f3f4');arrow((2.8,1.2),(3.4,1.2));arrow((8.2,4.5),(8.2,2.0))
 ax.text(.2,0,'Latent coordinates never change masks, anchor weights or physical-distance diagnostics.',fontsize=10);save(fig,'ExtensionArchitecture')
 # Entire20-estimator Bassiouni grid, no selected winner.
 a=pd.read_csv(src/'BASSIOUNI_PATIENT_METRICS.csv');a=a[a.experiment=='A'];methods=list(a.method.unique());ref=a[a.method=='idw_k16_p2'].set_index(['patient','budget']);fig,axs=plt.subplots(1,2,figsize=(9,8),sharey=True)
 for j,budget in enumerate([.1,.2,.4]):
  corr=[];mse=[]
  for method in methods:
   g=a[(a.method==method)&(a.budget==budget)].set_index(['patient','budget']);r=ref.loc[g.index];corr.append(np.median(g.median_pearson-r.median_pearson));mse.append(100*np.median((g.mean_mse-r.mean_mse)/r.mean_mse))
  axs[0].plot(corr,np.arange(len(methods))+(j-1)*.18,'o',label=f'{int(budget*100)}%');axs[1].plot(mse,np.arange(len(methods))+(j-1)*.18,'o')
 axs[0].set_yticks(range(len(methods)),methods);axs[0].invert_yaxis();axs[0].legend(title='Anchor budget');axs[0].set_xlabel('Median paired Pearson change vs IDW');axs[1].set_xlabel('Median paired relative MSE (%)')
 for ax in axs:ax.axvline(0,color='.5',ls='--');ax.grid(axis='x',alpha=.2)
 fig.tight_layout();save(fig,'ExtensionBassiouniEstimators')
 pairs=pd.read_csv(src/'BASSIOUNI_PAIRED_PATIENT_EFFECTS.csv');b2=pairs[pairs.origin=='new_B2_identical_targets'];fig,axs=plt.subplots(1,3,figsize=(10,3.8))
 for ax,budget in zip(axs,[.1,.2,.4]):
  for method,marker in [('idw_k16_p2','o'),('uniform_k16','^')]:
   q=b2[(b2.budget==budget)&b2.contrast.str.startswith(method)];ax.scatter(q.delta_pearson,100*q.relative_mse,s=25,marker=marker,label=method)
  ax.set(title=f'{int(100*budget)}% anchors',xlabel='Region-exclusion minus distributed r',ylabel='Paired relative MSE (%)');ax.axhline(0,color='.7');ax.axvline(0,color='.7')
 axs[0].legend(fontsize=9);fig.tight_layout();save(fig,'ExtensionBassiouniPlacement')
 uncertainty=pd.read_csv(src/'PATIENT_UNCERTAINTY.csv');fig,axs=plt.subplots(1,2,figsize=(9,7.8),sharey=True);variants=list(VARIANTS)
 for ax,effect,scale,label in zip(axs,['delta_pearson','relative_mse'],[1,100],['Paired Pearson change vs IDW','Paired relative MSE (%)']):
  for j,(design,condition,color) in enumerate([('historical_V1','historical','#176a87'),('B2_identical_H','region_exclusion','#bd6426')]):
   q=uncertainty[(uncertainty.cohort=='Bassiouni')&(uncertainty.design==design)&(uncertainty.condition==condition)&(uncertainty.budget==.2)&(uncertainty.reference=='idw_k16_p2')&(uncertainty.route=='exploratory_residual_anchor')&(uncertainty.effect==effect)&(uncertainty.origin=='extension_new')].set_index('candidate').loc[variants]
   y=np.arange(len(variants))+(j-.5)*.2;ax.errorbar(q.median_effect*scale,y,xerr=np.vstack([(q.median_effect-q.ci_low)*scale,(q.ci_high-q.median_effect)*scale]),fmt='o',ms=4,color=color,label=design+': '+condition)
  ax.axvline(0,color='.5',ls='--');ax.set_xlabel(label);ax.grid(axis='x',alpha=.2)
 axs[0].set_yticks(range(len(variants)),[NAMES[v] for v in variants]);axs[0].invert_yaxis();axs[0].legend(fontsize=9,loc='lower left');fig.tight_layout();save(fig,'ExtensionNeuralPaired')
 # Full baseline/control matrix on each route: no architecture omitted from these plots.
 patient=pd.read_csv(src/'PATIENT_METRICS.csv');full=patient[(patient.design=='historical_V1')&(patient.budget==.2)]
 order=[('new',v) for v in variants]+[('frozen',v) for v in sorted(full[full.origin=='frozen'].variant.unique())]+[('molecular','idw_k16_p2'),('molecular','uniform_k16')]
 for cohort in COHORTS:
  for route in ROUTES:
   fig,axs=plt.subplots(1,3,figsize=(9,10),sharey=True);subset=full[(full.cohort==cohort)&((full.route==route)|(full.origin=='molecular'))]
   for y,(origin,variant) in enumerate(order):
    q=subset[(subset.origin==origin)&(subset.variant==variant)]
    for ax,metric in zip(axs,METRICS):
     for k,row in enumerate(q.itertuples()):ax.plot(getattr(row,metric),y+(.1*(k-1) if cohort=='Wu' else 0),'o',ms=3,color=plt.cm.tab10(k) if cohort=='Wu' else '.7',alpha=.8)
     ax.plot(q[metric].median(),y,'|',color='black',ms=10)
   labels=[NAMES.get(v,v)+( ' [old]' if o=='frozen' else '') for o,v in order];axs[0].set_yticks(range(len(order)),labels);axs[0].invert_yaxis()
   for ax,label in zip(axs,['Median gene Pearson','Mean MSE','Mean MAE']):ax.set_xlabel(label);ax.grid(axis='x',alpha=.2)
   axs[0].set_title(cohort+' / '+route,loc='left')
   if cohort=='Wu':
    for k,p in enumerate(['1160920F','CID4290','CID44971']):axs[0].plot([],[],'o',color=plt.cm.tab10(k),label=p)
    fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',ncol=3,fontsize=11)
   fig.tight_layout(rect=(0,.04,1,1) if cohort=='Wu' else (0,0,1,1));save(fig,'ExtensionMatrix_'+cohort+'_'+route)
 # Reach is a software diagnostic, not evidence of molecular biology.
 reach=pd.read_csv(src/'ATTENTION_REACH.csv' if (src/'ATTENTION_REACH.csv').exists() else base/'diagnostics/ATTENTION_REACH.csv');reach.to_csv(src/'ATTENTION_REACH.csv',index=False) if not (src/'ATTENTION_REACH.csv').exists() else None
 return figures


def prepare(root):
 package=root/'manuscript/journal_revision/transformer_egnn'
 if package.exists():raise ValueError('Refuse to overwrite an existing versioned manuscript; build into a clean checkout/output instead')
 shutil.copytree(root/'manuscript/journal_revision',package,ignore=lambda p,n:[x for x in n if x=='transformer_egnn' or (x.startswith('GraphPathology_') and x.endswith('.pdf'))])
 for name in ['JOURNAL_BUILD_REPORT.json']:
  if (package/name).exists():(package/name).unlink()
 src=package/'source_data/extension';src.mkdir(parents=True)
 base=root/'experiments/journal_revision/transformer_egnn'
 for p in (base/'summary').glob('*.csv'):shutil.copy2(p,src/p.name)
 for p in (base/'bassiouni_anchors').glob('*.csv'):shutil.copy2(p,src/('BASSIOUNI_'+p.name))
 for name in ['ATTENTION_REACH.csv','GEOMETRIC_INVARIANCE.csv','DISTANT_FEATURE_PERTURBATION.csv','CPU_BENCHMARK.csv']:shutil.copy2(base/'diagnostics'/name,src/name)
 shutil.copytree(base/'fits',src/'fits')
 shutil.copytree(root/'src/graphpathology',package/'repository_documentation/extension_code/graphpathology',ignore=shutil.ignore_patterns('__pycache__'))
 for name in ['PROTOCOL.md','DATA_RECOVERY.md','FROZEN_HEST_REGRESSION.json']:shutil.copy2(base/name,package/'repository_documentation'/name)
 (package/'scripts/regenerate_extension_figures.py').write_text('from pathlib import Path\nimport sys\np=Path(__file__).resolve().parents[1]\nsys.path.insert(0,str(p/\"repository_documentation/extension_code\"))\nfrom graphpathology.extension_paper import plot_all,plot_more\nplot_all(p,p);plot_more(p,p)\n')
 return package

METRICS=['median_pearson','mean_mse','mean_mae']


def sections(root,package,figures):
 src=package/'source_data/extension';u=pd.read_csv(src/'PATIENT_UNCERTAINTY.csv');b=pd.read_csv(src/'BASSIOUNI_PATIENT_UNCERTAINTY.csv');patient=pd.read_csv(src/'PATIENT_METRICS.csv');fits=pd.read_csv(src/'FIT_INVENTORY.csv')
 def effect(variant,reference,metric,design='historical_V1',condition='historical',route='exploratory_residual_anchor'):
  q=u[(u.cohort=='Bassiouni')&(u.origin=='extension_new')&(u.candidate==variant)&(u.reference==reference)&(u.effect==metric)&(u.budget==.2)&(u.design==design)&(u.condition==condition)&(u.route==route)]
  if len(q)!=1:raise ValueError('Ambiguous/missing scientific effect '+str((variant,reference,metric,len(q))))
  return q.iloc[0]
 def fmt(x,scale=1):return f'{x.median_effect*scale:+.4f} [{x.ci_low*scale:+.4f}, {x.ci_high*scale:+.4f}]'
 rows=[]
 for method,name in [('idw_k16_p2','IDW'),('uniform_k16','Uniform')]:
  vals=[]
  for metric,scale in [('delta_pearson',1),('relative_mse',100),('delta_mae',1)]:
   x=b[(b.origin=='new_B2_identical_targets')&(b.budget==.2)&b.contrast.str.startswith(method)&(b.effect==metric)].iloc[0];vals.append(fmt(x,scale))
  rows.extend([[name,label,value] for label,value in zip([r'$\Delta r$',r'MSE\%',r'$\Delta$ MAE'],vals)])
 (package/'tables/ExtensionPlacement.tex').write_text(table(r'\textbf{New exact-input Bassiouni placement effects.} Region-exclusion minus distributed anchors at20\% coverage, identical H and anchor counts. Median paired patient effects with descriptive95\% percentile intervals from10,000 patient resamples. Positive error changes deteriorate. All22 patients are retained; masks are averaged within patient. The full20-estimator and three-budget exports accompany the paper.','tab:extension_placement',['Operator','Effect',r'Median [95\% interval]'],rows))
 controls=[('gps_full','gps_local'),('gps_full','gps_no_relative'),('egnn_latent','distance_fixed'),('egnn_latent','distance_disabled')];rows=[]
 for candidate,reference in controls:
  vals=[fmt(effect(candidate,reference,m),s) for m,s in [('delta_pearson',1),('relative_mse',100),('delta_mae',1)]]
  rows.extend([[NAMES[candidate]+' vs '+NAMES[reference],label,value] for label,value in zip([r'$\Delta r$',r'MSE\%',r'$\Delta$ MAE'],vals)])
 (package/'tables/ExtensionControls.tex').write_text(table(r'\textbf{Predefined architecture contrasts in recovered Bassiouni patients.} At20\% historical V1 coverage, the exploratory residual route uses identical measured anchors for both models. Effects average all three training seeds and three mask seeds within each patient before pairing; relative MSE is normalized by the reference patient mean MSE after these averages. Intervals resample22 whole patients, not genes, masks or training seeds. They are descriptive, unadjusted intervals, not simultaneous significance claims. Other routes, controls, budgets and individual Wu effects are fully exported.','tab:extension_controls',['Candidate vs control','Effect',r'Median [95\% interval]'],rows))
 primary=[]
 for variant in ['gps_full','egnn_latent']:
  x=effect(variant,'idw_k16_p2','delta_pearson');e=effect(variant,'idw_k16_p2','relative_mse');ma=effect(variant,'idw_k16_p2','delta_mae')
  primary.append(NAMES[variant]+r' changed patient correlation by '+fmt(x)+r', relative MSE by '+fmt(e,100)+r'\% and MAE by '+fmt(ma)+'.')
 results=r'''
\subsection{Recovered exact inputs completed the22-patient matched baseline study}
The previously failed individual GEO transfers did not establish that the data were unavailable. Selective byte ranges from the official TAR recovered all88 required original files with exact historical SHA256 matches, transferring291,569,893 bytes rather than the37,690,060,800-byte whole bundle. All22 patients passed the original primary numerical replay before new scoring. We then completed19,800 new anchor-only runs and39,600 exact NaN/random substitutions under the existing A/B1/B2 protocol, in a separate output directory. The20-estimator grid remains fully reported rather than used to select a winner (Fig.~\ref{fig:extension_baselines}).

The controlled placement comparison holds the hidden set H and total anchor count identical, unlike the historical V1/V2 comparison. Each patient contributes one seed-averaged effect. Table~\ref{tab:extension_placement} and Fig.~\ref{fig:extension_placement} report both fixed estimators, all22 patients and all budgets. The uncertainty is descriptive patient resampling within this previously inspected study; it is not fresh validation or a universal biological distance scale.
'''+figure('figures/extension/ExtensionBassiouniEstimators.pdf',r'\textbf{Complete matched-information baseline grid in22 previously inspected Bassiouni patients.} Every declared estimator and budget is shown. Points are medians of paired patient effects relative to fixed IDW k16/p2 after averaging ten original mask trajectories. Correlation increases and error decreases improve. No estimator was selected from these outcomes.','fig:extension_baselines')+figure('figures/extension/ExtensionBassiouniPlacement.pdf',r'\textbf{Identical-target placement effects in all22 recovered patients.} Each point is one patient, averaged over ten mask seeds; shapes distinguish fixed IDW and uniform k16. Both anchor placements score the same H with equal anchor counts. Panels retain all three budgets. These are not historical V2-minus-V1 effects.','fig:extension_placement')+r'\input{tables/ExtensionPlacement}'+r'''

\subsection{Bounded nonlocal and latent-coordinate models were tested without external selection}
We committed an outcome-blind local protocol before fitting16 specified configurations with seeds0/1/2, retaining all48 new fits alongside the30 frozen journal-v2 heads. This is a post-v1 exploratory protocol, not formal preregistration. Every new model used the same two source training patients, validation-only stopping, full-section update budget and frozen image embeddings. The GraphGPS-inspired implementation combines local gated message passing with at most32 spatial landmarks; the latent-coordinate EGNN is compared directly with a fixed-coordinate distance model (Fig.~\ref{fig:extension_architecture}). These are new implementations inspired by published architectures, not official GraphGPS benchmark reproductions\cite{rampasek2022,satorras2021}.

The complete78-model matrix crossed identical masks and information routes in TENX99 (one previously inspected development patient), three previously inspected Wu patients and22 previously inspected Bassiouni patients. External image outputs used the same frozen C9 T1 affine transformation with no destination refitting. Image-only comparisons against IDW are unequal-information comparisons. Fixed and residual routes provide the same measured anchors as IDW and test incremental image/context information. No evaluation patient, architecture outcome, epoch or seed selected a new fit, coefficient or display subset.

'''+ ' '.join(primary)+r'''
These predefined contrasts retain their error and correlation directions even when they disagree. Figure~\ref{fig:extension_neural} shows every new architecture/control against IDW for distributed and region-exclusion settings. Table~\ref{tab:extension_controls} tests the nonlocal and coordinate-update components against their specified controls. Per-patient effects, all metrics, all masks, every training seed and all three routes are supplied in source data. The three Wu patients remain individually visible; no primary three-patient bootstrap claim is made. Software success and geometric tests do not themselves establish an accuracy benefit, spatial specificity, superiority over published systems or learned long-range molecular correlation.
'''+figure('figures/extension/ExtensionArchitecture.pdf',r'\textbf{Implemented extension and information boundaries.} The image network receives frozen image embeddings and coordinate-defined structure. The GraphGPS-inspired branch combines local gated aggregation with bounded attention to32 coordinate-only landmarks. The alternative EGNN changes private latent coordinates through relative displacements; the distance-only model fixes them. Original physical coordinates always define masks, anchor interpolation and physical-distance diagnostics. No hidden expression is an image-network input. The optional learned anchor-conditioned architecture was not executed.','fig:extension_architecture')+figure('figures/extension/ExtensionNeuralPaired.pdf',r'\textbf{All16 new architectures/controls versus molecular IDW in22 Bassiouni patients.} Exploratory residual route at20\% coverage; colors distinguish historical distributed V1 and identical-H region exclusion. Bars are descriptive95\% patient percentile intervals for median paired effects after mask/training-seed averages. Positive correlation and negative relative MSE improve. The two designs have different scored support and are not subtracted here. All models are shown regardless of direction; nominally equal parameter counts do not imply equal effective function classes.','fig:extension_neural')+r'\input{tables/ExtensionControls}'+r'''

\subsection{Geometric and reach checks delimit architectural claims}
All384 saved-checkpoint geometric checks passed. They establish scalar invariance and latent-coordinate equivariance conditional on fixed embeddings and the documented graph/landmark construction or transformed fixed graph. Float64 primary-model tolerances were10^{-8}; float32 and preaggregated baseline tolerances were10^{-4}. Physical latent tolerances include the250\,\um{} scale. Recomputed-graph unit tests use nondegenerate synthetic coordinates without distance ties. No claim of rotation/reflection invariance of the frozen CNN or complete image pipeline follows.

Connectivity fractions and actual attention mass beyond250/600/1200\,\um{}, together with a distant-feature perturbation diagnostic, distinguish functional nonlocal access from a verbal long-range claim. Neither establishes a biological correlation length. Degree/length controls report actual successful rewiring; zero or sparse swaps are limited controls rather than independent positive evidence. Full route/control matrices, occupied distance bins with support, parameter counts, validation-selected epochs and resource records accompany the Supplementary Information.
'''
 methods=r'''
\subsection{Post-v1 bounded Transformer and latent-coordinate extension}
The timestamped local protocol commit preceded all48 fits. AdamW used learning rate10^{-3}, weight decay10^{-4}, maximum80 epochs, patience10 and minimum validation-MSE improvement10^{-4}; there was one full-section update per source training patient per epoch, in seeded random order. TENX95 alone selected epochs. Local/mean/GCN widths were128; concatenated GraphSAGE used67. The older30 heads used their original pooled512-node minibatch regime and remain frozen references, not optimization-matched controls. New full-section baselines supply the matched-update-budget comparison. Models differ in parameter count and compute; fixed epoch/update budgets are not matched floating-point-operation budgets.

GraphGPS-inspired blocks use a512-to64 projection, two64-wide blocks, four attention heads, dropout0.1, pre-LayerNorm, local linear messages gated by eight relative-distance radial features, parallel bounded attention, residual additions and a two-layer GELU feed-forward path. The local graph is Delaunay with a250\,\um{} cap. Global keys/values are at most32 deterministic coordinate-only farthest-point landmarks, giving O(NL) attention storage with L\leq32; there is no dense N-by-N allocation. The nonlocal-disabled model attends self plus up to seven nearest graph neighbors, padding unused slots with self. Thus full/local changes attention range and key count as well as connectivity, despite identical121,018 nominal parameters. Local radial features use log(1+d/250) with eight centers0--1 and width0.15; attention relative bias uses eight centers0--4 and width0.5. The distance-disabled version removes explicit gates/biases but retains coordinate-defined graph/landmark selection and is not coordinate-free.

The invariant fixed-coordinate and EGNN alternatives use two64-wide message blocks. With centered normalized p=(x-mean(x))/250, m_{ij}=phi_e(h_i,h_j,RBF(||p_i-p_j||^2)), node updates aggregate the mean message through phi_h, using SiLU, LayerNorm and dropout0.1. Eight radial centers span0--4 with width0.5. EGNN additionally updates private coordinates as
\[
p_i^{l+1}=p_i^l+\frac{1}{\max(1,|\mathcal N_i|)}\sum_{j\in\mathcal N_i}0.1\tanh\!\big(\phi_x(m_{ij})\big)(p_i^l-p_j^l).
\]
Each phi_x is a bias-free64-to1 map initialized at standard deviation0.001. Only the first coordinate update affects the second layer's scalar readout; the last update has no downstream scalar loss. EGNN has87,282 nominal parameters versus87,154 for fixed distance (128 additional coordinate coefficients, of which the last64 do not affect prediction). Both are E(2)-compatible conditional on frozen features, so their accuracy difference cannot identify equivariance itself. Latent positions never change original masks, anchors or IDW distances.

Self-only and coordinate/graph-shuffled variants retain nominal parameter tensors but have different effective function classes; inactive attention or edge parameters do not make them capacity-equivalent. Conjugating the graph and assigned coordinates together preserves edge-length geometry while breaking image-to-location assignment. Degree/length rewiring preserves each node's degree and global50\,\um{} edge-length bins below the same250\,\um{} cap, targeting E/4 accepted double swaps with at most20E attempts. Actual successes are reported, including infeasible or weak rewiring. No hidden expression enters any graph, landmark or distance encoding.

Development uses nine original random masks (three seeds at each budget); external cohorts use original V1 and central V2 geometry1 with seeds0/1/2. All78 models also use B2 seeds0/1/2 at10/20/40\% budgets. This bounded extension does not replace the original ten-seed/three-geometry analyses. B2 uses identical H and equal anchor counts in both placements. For every fixed/residual route, all unobserved labels, including unscored labels, were replaced with NaNs and independent random numbers and predictions required exact equality. Image-only B2 metrics must be identical across placements, budgets and mask seeds because H is fixed. Full per-gene metrics and defined-correlation counts are retained. Seed averages precede patient effects, cohorts remain separate, and the22-patient bootstrap uses10,000 resamples and seed20261101; intervals are descriptive and unadjusted.

Exact Bassiouni raw recovery activated the protocol's conditional cohort analysis without model or calibration changes. All88 original SHA256 checks, source barcode joins and independently derived coordinates were verified; source coordinates agreed with preserved label-free coordinates within10^{-8}\,\um{}, which then governed exact neighbor ties. External image patches used the identical112\,\um{} field, nearest-integer side/top-left, LANCZOS224 resize, frozen ResNet18 ImageNet weights and their published transform; no labels entered feature extraction. The existing C9 T1 map was applied to all image predictions before both image and observed-residual routes, with no new calibration fit. Both external cohorts were previously inspected. Low native image resolution and cross-platform scale/domain shift constrain morphology conclusions.

Conditional geometric checks apply rotations, reflections, translations and node permutations at float64/float32, recording scalar and physical latent-coordinate errors separately. Attention diagnostics report physical and assigned-coordinate connectivity fractions and actual mass beyond250/600/1200\,\um{} per head/layer/seed/patient. The distant perturbation zeros one coordinate-selected remote landmark feature; it tests functional dependency only. Optional learned anchor-conditioned models and complete-image equivariance experiments were skipped by protocol. CPU timings, memory, exact environments, optimizer states, complete histories and data/checkpoint hashes are retained; model weights and raw datasets are excluded from Git.
'''
 return clean_prose(results),clean_prose(methods)


def supplementary(root,package,figures):
 src=package/'source_data/extension';base=root/'experiments/journal_revision/transformer_egnn';chunks=[r'\clearpage\section{Post-v1 bounded Transformer and EGNN extension}',r'These are exploratory analyses in previously inspected cohorts. Full per-gene/mask/seed exports, hashes, loss histories and paired patient effects accompany the source package. The patient is the biological unit. Black ticks in the complete matrices are cohort medians; dots are individual patient summaries after seed averages, with separate Wu patient colors. No spread of dots is a confidence interval. Molecular rows receive same-section measured labels; image-only rows do not. Older frozen models are marked [old] and were trained under a different update regime.']
 for cohort in COHORTS:
  for route in ROUTES:
   name='ExtensionMatrix_'+cohort+'_'+route
   chunks.append(figure('figures/extension/'+name+'.pdf',r'\textbf{Complete '+esc(cohort)+' '+esc(route)+r' matrix at20\% historical V1 coverage.} All16 new architectures/controls, all10 older configuration families (three fits each), and both fixed molecular comparators are shown. Model summaries average training/mask seeds within patient. Pearson, MSE and MAE are retained together; no external winner is selected. Fixed and residual routes use identical anchors as their molecular references. The HEST development plot has one patient and no population interval.','fig:'+name))
 diag=pd.read_csv(src/'GEOMETRIC_INVARIANCE.csv');rows=[]
 def sci(x):
  if x==0:return '0'
  a,b=f'{x:.2e}'.split('e');return '$'+a+r'\times10^{'+str(int(b))+'}$'
 for (variant,precision),g in diag.groupby(['variant','precision'],sort=False):
  rows.append([NAMES[variant],precision,sci(g.max_scalar_abs_delta.max()),sci(g.max_latent_um_abs_delta.max()) if g.max_latent_um_abs_delta.notna().any() else '--',str(len(g))])
 chunks.append(r'\begin{longtable}{@{}llrrr@{}}\caption{\textbf{All384 saved-checkpoint conditional geometric checks.} Maxima across all three seeds and recorded transforms. Physical latent errors are in micrometers. Tests are conditional on fixed image features, not the complete CNN pipeline. Baseline preaggregation remains float32 even when the head is tested in float64.}\label{tab:extension_geometric}\\\toprule Model & Precision & Max scalar error & Max latent error & Checks \\\midrule\endfirsthead\toprule Model & Precision & Max scalar error & Max latent error & Checks \\\midrule\endhead'+'\n'+'\n'.join(' & '.join(r)+r' \\' for r in rows)+r'\bottomrule\end{longtable}')
 fits=pd.read_csv(src/'FIT_INVENTORY.csv');bench=pd.read_csv(src/'CPU_BENCHMARK.csv');rows=[]
 for variant in VARIANTS:
  f=fits[fits.variant==variant];b=bench[bench.variant==variant];rows.append([NAMES[variant],str(int(f.parameters.iloc[0])),f'{f.selected_epoch.min():.0f}--{f.selected_epoch.max():.0f}',f'{f.elapsed_seconds.median():.2f}',f'{1000*b.median_forward_seconds.median():.1f}',f'{b.peak_process_rss_mb.max():.0f}'])
 chunks.append(table(r'\textbf{Resource accounting for every new configuration.} All three seeds are included. Training wall times were measured under the documented concurrent workload and are not a controlled throughput comparison. Standardized inference uses the same7,760-node TENX95 section, four CPU threads, two warmups and median of five forwards, including model-specific feature aggregation/Tensor conversion but excluding frozen CNN and graph construction. Peak process RSS includes interpreter and data, not only model activations. No GPU or paid external service was used. Full per-seed values, graph construction times and environments are exported.','tab:extension_resources',['Model','Parameters','Best epoch range','Fit s','Forward ms','RSS MiB'],rows))
 rows=[]
 for r in fits.sort_values(['variant','seed']).itertuples():rows.append([NAMES[r.variant],str(r.seed),str(r.selected_epoch),str(r.epochs_ran),f'{r.elapsed_seconds:.2f}',f'{r.peak_rss_mb:.0f}'])
 chunks.append(r'\begin{longtable}{@{}lrrrrr@{}}\caption{\textbf{Every new fit, without seed selection.} Histories, optimizer states, input hashes and complete checkpoint SHA256 values are available in the fit JSON files and inventory CSV; weights are distributed separately from Git.}\label{tab:extension_fits}\\\toprule Model & Seed & Best epoch & Ran & Seconds & RSS MiB \\\midrule\endfirsthead\toprule Model & Seed & Best epoch & Ran & Seconds & RSS MiB \\\midrule\endhead'+'\n'+'\n'.join(' & '.join(r)+r' \\' for r in rows)+r'\bottomrule\end{longtable}')
 chunks.append(r'\paragraph{Source interpretation.} Undefined correlations remain undefined; valid-gene counts accompany every panel row. Raw Spearman is secondary because near-tied ranks are floating-point sensitive. Complete distance-bin errors retain occupied-bin support and omit empty bins. Native external images are heavily downsampled, the model is trained in only two source patients, and calibration is frozen. The original selected primary heads are still missing; these48 new states and30 older journal-v2 states do not replace them. Optional anchor-conditioned learning and complete-image equivariance tests have status SKIPPED\_PROTOCOL. No null or missing analysis is converted into a negative biological finding.')
 return clean_prose('\n\n'.join(chunks))


def assemble(root,package,figures):
 result,methods=sections(root,package,figures)
 for exponent in [-3,-4,-8]:
  for name,value in [('result',result),('methods',methods)]:
   value=value.replace('10^{'+str(exponent)+'}',r'$10^{'+str(exponent)+'}$')
   if name=='result':result=value
   else:methods=value
 methods=methods.replace(r'L\leq32',r'$L\leq32$').replace('p=(x-mean(x))/250',r'$p=(x-\mathrm{mean}(x))/250$').replace('m_{ij}=phi_e(h_i,h_j,RBF(||p_i-p_j||^2))',r'$m_{ij}=\phi_e(h_i,h_j,\mathrm{RBF}(\|p_i-p_j\|^2))$').replace('through phi_h',r'through $\phi_h$').replace('Each phi_x',r'Each $\phi_x$')
 article=(package/'article_content.tex').read_text().replace('separately preregistered','protocol-specified exploratory').replace('preregistered exploratory','protocol-specified exploratory').replace('preregistered','protocol-specified')
 article=article.replace(r'\section{Discussion}',result+'\n'+r'\section{Discussion}',1).replace(r'\section{Ethics}',methods+'\n'+r'\section{Ethics}',1)
 article=article.replace('New Bassiouni matched-estimator, identical-target and spatial-map analyses are therefore unavailable, rather than inferred from old summaries.','At that earlier revision, the new Bassiouni analyses were unavailable. The separate later extension below recovered exact original inputs and completed matched-estimator, identical-target and model comparisons; no historical summary substituted for a new prediction.')
 article=article.replace('In the present analysis, the defensible conclusion is that the entire contiguous-completion protocol is harder, not that geometry alone explains the measured decline.','The historical protocol comparison remains compound. The new identical-H experiments instead isolate anchor placement for the specified central geometry and equal counts; they do not establish a universal spatial effect.')
 article=article.replace('That test did not show a reproducible benefit.','That original fixed-fusion test did not show a reproducible benefit. The separately reported exploratory residual-anchor extension evaluates new models and does not replace that historical result.')
 article=article.replace('New Bassiouni raw-data comparisons remained blocked by failed source transfers.','The later extension recovered exact Bassiouni inputs and completed matched baselines and geometry-aware neural controls.')
 article=article.replace('but analogous new Bassiouni analyses remain blocked by source acquisition','and the later exact-input Bassiouni extension completes the analogous baseline study')
 article=article.replace('Bassiouni intervals reanalyze original exports and cannot replace unavailable new predictions.','Those earlier intervals reanalyzed original exports; the separate later extension also reports intervals from new exact-input Bassiouni predictions.')
 article=article.replace('During this revision, the selected GEO transfers returned errors or empty bodies despite advertised nonzero sizes; new Bassiouni raw-data analyses were not performed.','Initial individual GEO transfers failed; the later extension recovered all88 exact original files through bounded HTTP byte ranges from the official TAR and performed the separately reported new analyses.')
 # Abstract is generated only after every required result is complete.
 u=pd.read_csv(package/'source_data/extension/PATIENT_UNCERTAINTY.csv')
 def delta(v):
  q=u[(u.cohort=='Bassiouni')&(u.candidate==v)&(u.reference=='idw_k16_p2')&(u.origin=='extension_new')&(u.effect=='delta_pearson')&(u.design=='historical_V1')&(u.condition=='historical')&(u.budget==.2)&(u.route=='exploratory_residual_anchor')].iloc[0];return f'{q.median_effect:+.3f}'
 abstract=r'The cost of spatial transcriptomics motivates reconstruction from sparse molecular observations. GraphPathology separates image context from interpolation of measured anchors. Original50-gene results in25 evaluation patients are preserved. A separate post-v1 exploratory study completed matched20-estimator and identical-target geometry comparisons in22 Bassiouni patients after exact input recovery, supplementing the prior three-patient Wu analysis. We trained48 new GraphGPS-inspired, fixed-distance, latent-coordinate EGNN and baseline/control models with a common source-only update budget, and evaluated them alongside30 frozen journal-v2 heads in26 previously inspected patients. At20\% distributed coverage in Bassiouni, residual-anchor GPS and EGNN changed median paired patient correlation relative to fixed IDW by '+delta('gps_full')+' and '+delta('egnn_latent')+r', respectively. Correlation, MSE and MAE, component ablations and individual patient effects are retained together. Conditional geometric tests and bounded attention reach validate software properties rather than establish biological long-range correlation. No new fit replaces the missing historical primary states; no external patient selected architecture, calibration or seed. Limited image resolution and cross-platform transport constrain morphology conclusions. Fivefold fewer measured locations is conditional variable-cost arithmetic, not demonstrated total assay savings.'
 article=re.sub(r'(\\section\*\{Abstract\}\n).*?(\n\n\\noindent)',lambda m:m.group(1)+clean_prose(abstract)+m.group(2),article,count=1,flags=re.S)
 (package/'article_content.tex').write_text(article)
 supplement=(package/'supplementary_content.tex').read_text().replace('preregistered','protocol-specified').replace('New Bassiouni prediction/map analyses remain blocked by failed canonical-source transfers; no old summary substitutes for a new prediction.','At that earlier revision, Bassiouni prediction/map analyses were blocked by failed individual transfers. The separate later extension recovered exact inputs through bounded TAR ranges and reports new predictions; no old summary substitutes for a new prediction.')
 (package/'supplementary_content.tex').write_text(supplement+'\n'+supplementary(root,package,figures))
 refs=r'''
@inproceedings{rampasek2022,
 title={Recipe for a General, Powerful, Scalable Graph Transformer},
 author={Rampasek, Ladislav and Galkin, Mikhail and Dwivedi, Vijay Prakash and Luu, Anh Tuan and Wolf, Guy and Beaini, Dominique},
 booktitle={Advances in Neural Information Processing Systems},volume={35},pages={14501--14515},year={2022}}
@inproceedings{satorras2021,
 title={E(n) Equivariant Graph Neural Networks},
 author={Satorras, Victor Garcia and Hoogeboom, Emiel and Welling, Max},
 booktitle={Proceedings of the 38th International Conference on Machine Learning},
 series={Proceedings of Machine Learning Research},volume={139},pages={9323--9332},year={2021}}
'''
 with (package/'references.bib').open('a') as f:f.write(refs)
 (package/'README.md').write_text('# Separate post-v1 GraphTransformer/EGNN journal version\n\nThis version extends, and does not overwrite, manuscript/journal_revision or manuscript/revised. Every result comes from the completed extension exports. See the outcome-blind local protocol, recovery record, complete fit inventory and validation report under experiments/journal_revision/transformer_egnn. This is exploratory work in previously inspected cohorts, not formal preregistration or fresh blind validation. Run `python -m graphpathology.extension paper` to rebuild from committed exports without downloads/training.\n')


def plot_more(root,package):
 import matplotlib.pyplot as plt
 base=root/'experiments/journal_revision/transformer_egnn';out=package/'figures/extension';src=package/'source_data/extension';figures=[]
 def save(fig,name):
  fig.savefig(out/(name+'.pdf'),bbox_inches='tight');fig.savefig(out/(name+'.png'),dpi=160,bbox_inches='tight');plt.close(fig);figures.append(name)
 reaches=[pd.read_csv(base/'diagnostics/ATTENTION_REACH.csv')]+[pd.read_csv(p) for p in sorted((base/'diagnostics/recovered').glob('*/ATTENTION_REACH.csv'))] if not (src/'ALL_ATTENTION_REACH.csv').exists() else [pd.read_csv(src/'ALL_ATTENTION_REACH.csv')]
 reach=pd.concat(reaches,ignore_index=True);reach['cohort']=np.where(reach.patient=='TENX99','Development',np.where(reach.patient.str.startswith('BASS_'),'Bassiouni','Wu'));reach.to_csv(src/'ALL_ATTENTION_REACH.csv',index=False) if not (src/'ALL_ATTENTION_REACH.csv').exists() else None
 metrics=['physical_fraction_beyond','assigned_fraction_beyond','attention_mass_beyond'];variants=[v for v in VARIANTS if v.startswith('gps_')]
 r=reach.groupby(['cohort','patient','variant','threshold_um','seed'])[metrics].mean().groupby(['cohort','patient','variant','threshold_um']).mean().reset_index()
 fig,axs=plt.subplots(3,3,figsize=(10,8),sharex=True,sharey=True)
 for row,cohort in enumerate(COHORTS):
  for col,metric in enumerate(metrics):
   ax=axs[row,col]
   for j,threshold in enumerate([250,600,1200]):
    for k,variant in enumerate(variants):
     q=r[(r.cohort==cohort)&(r.variant==variant)&(r.threshold_um==threshold)][metric];ax.plot(q,np.full(len(q),k+(j-1)*.18),'.',color=plt.cm.tab10(j),alpha=.25);ax.plot(q.median(),k+(j-1)*.18,'o',ms=4,color=plt.cm.tab10(j))
   ax.set_yticks(range(len(variants)),[NAMES[v] for v in variants]);ax.set_xlim(-.02,1.02);ax.grid(axis='x',alpha=.2)
   if row==0:ax.set_title(['Physical far fraction','Assigned-coordinate far fraction','Attention mass physically far'][col],fontsize=11)
   if col==0:ax.set_ylabel(cohort)
 axs[0,0].invert_yaxis()
 for j,threshold in enumerate([250,600,1200]):axs[-1,0].plot([],[],'o',color=plt.cm.tab10(j),label=f'>{threshold} um')
 axs[-1,0].legend(fontsize=9);fig.tight_layout();save(fig,'ExtensionAttentionReach')
 # Loss curves include all source-only fits and all seeds; no test results enter selection.
 fitrecords=[json.loads(p.read_text()) for p in ((src/'fits') if (src/'fits').exists() else (base/'fits')).glob('*.json') if not p.name.endswith('_FAILED.json')];fig,axs=plt.subplots(4,4,figsize=(10,8))
 for ax,variant in zip(axs.flat,VARIANTS):
  for record in fitrecords:
   if record['variant']!=variant:continue
   h=pd.DataFrame(record['history']);color=plt.cm.tab10(record['seed']);ax.plot(h.epoch,h.training_mse,color=color,lw=1);ax.plot(h.epoch,h.validation_mse,color=color,lw=1,ls='--')
  ax.set_title(NAMES[variant],fontsize=11);ax.set_xlabel('Epoch',fontsize=10);ax.set_ylabel('MSE',fontsize=10);ax.tick_params(labelsize=9)
 fig.tight_layout();save(fig,'ExtensionLearningCurves')
 # Distance bins: equal mask means, then equal training-seed means, then patient medians.
 frames=[]
 for p in sorted((base/'evaluation').glob('BASS_*/*/BINS.csv.gz')):
  d=pd.read_csv(p)
  if not d.empty and d.origin.iloc[0]=='new':frames.append(d)
 for p in sorted((base/'evaluation').glob('BASS_*/MOLECULAR_BINS.csv.gz')):
  d=pd.read_csv(p);d['variant']=d.model_id;d['origin']='molecular';d['training_seed']=-1;d['route']='anchor_only';frames.append(d)
 if (src/'BASSIOUNI_DISTANCE_PATIENT_METRICS.csv').exists():bm=pd.read_csv(src/'BASSIOUNI_DISTANCE_PATIENT_METRICS.csv')
 else:
  bins=pd.concat(frames,ignore_index=True);bins=bins[(bins.design=='B2_identical_H')&(bins.condition=='region_exclusion')&(bins.budget==.2)]
  keys=['patient','variant','origin','route','training_seed','low_um'];bm=bins.groupby(keys)[['mean_mse','mean_mae','n_locations']].mean().groupby([k for k in keys if k!='training_seed']).mean().reset_index();bm.to_csv(src/'BASSIOUNI_DISTANCE_PATIENT_METRICS.csv',index=False)
 support=bm[(bm.variant=='idw_k16_p2')].groupby('low_um').agg(n_patients=('patient','nunique'),median_locations_per_mask=('n_locations','median'),min_locations_per_mask=('n_locations','min'),max_locations_per_mask=('n_locations','max'));support.to_csv(src/'BASSIOUNI_DISTANCE_SUPPORT.csv') if not (src/'BASSIOUNI_DISTANCE_SUPPORT.csv').exists() else None
 fig,axs=plt.subplots(3,1,figsize=(9,10),sharex=True)
 for ax,route in zip(axs,ROUTES):
  for j,variant in enumerate(list(VARIANTS)+['idw_k16_p2','uniform_k16']):
   q=bm[(bm.variant==variant)&((bm.route==route)|(bm.origin=='molecular'))].groupby('low_um').mean(numeric_only=True)
   # Plot cohort median patient errors, not a pooled spot error.
   q=bm[(bm.variant==variant)&((bm.route==route)|(bm.origin=='molecular'))].groupby('low_um')[['mean_mse']].median()
   ax.plot(q.index,q.mean_mse,'o-',ms=3,color=plt.cm.tab20(j),label=NAMES.get(variant,variant))
  ax.set_title(route,fontsize=12);ax.set_ylabel('Median patient mean MSE');ax.grid(alpha=.2)
 axs[-1].set_xlabel('Nearest-anchor distance-bin lower bound (um)');fig.legend(*axs[0].get_legend_handles_labels(),loc='lower center',ncol=4,fontsize=9);fig.tight_layout(rect=(0,.12,1,1));save(fig,'ExtensionDistanceBins')
 return figures


def build(root,output=None,bibliography=True):
 from .publication import compile_documents
 base=root/'experiments/journal_revision/transformer_egnn';status=json.loads((base/'summary/STATUS.json').read_text());assert status['status']=='COMPLETE' and status['patients']==26
 package=prepare(root) if output is None else Path(output)
 if output is not None:
  if package.exists():raise ValueError('Build destination must not exist')
  shutil.copytree(root/'manuscript/journal_revision/transformer_egnn',package)
 else:
  figures=plot_all(root,package)+plot_more(root,package);assemble(root,package,figures)
  text=(package/'supplementary_content.tex').read_text()
  for name,title,caption in [
   ('ExtensionAttentionReach','Conditional attention reach','Dots are individual patient summaries after equal seed/head/layer means; large dots are cohort medians. Colors identify250/600/1200micrometer thresholds. This is a software/connectivity diagnostic, not evidence of molecular correlation. All GPS variants, seeds and patients are exported.'),
   ('ExtensionLearningCurves','All48 source-only training histories','Every one of16 configurations and all three seeds are retained. Solid curves show training MSE, dashed curves validation MSE. Color identifies seeds0/1/2. Only TENX95 validation MSE selected epochs; evaluation outcomes did not select architecture or seeds.'),
   ('ExtensionDistanceBins','Distance-resolved recovered Bassiouni reconstruction','All16 new models and both molecular baselines,20percent region-exclusion B2. Curves are medians of per-patient conditional-bin mean errors after equal mask and training-seed averages. Occupied bins only; empty bins are omitted. Source tables retain number of patients and min/median/max mean locations per mask in each bin. Different bins need not contain the same patients or locations and do not identify a biological correlation length.')]:
   text+='\n'+figure('figures/extension/'+name+'.pdf',r'\textbf{'+title+'.} '+caption,'fig:'+name)
  (package/'supplementary_content.tex').write_text(text)
 reports=compile_documents(package,bibliography=bibliography)
 names={'main':'GraphPathology_Combined_Manuscript_Transformer_EGNN.pdf','article_only':'GraphPathology_Article_Only_Transformer_EGNN.pdf','supplementary':'GraphPathology_Supplementary_Transformer_EGNN.pdf'}
 for stem,name in names.items():shutil.copy2(package/(stem+'.pdf'),package/name)
 text=subprocess.check_output(['pdftotext','-layout',str(package/'main.pdf'),'-'],text=True)
 figs=set(re.findall(r'\bFigure\s+(S?\d+)\s*:',text));tables=set(re.findall(r'\bTable\s+(S?\d+)\s*:',text));refs=len(re.findall(r'\\bibitem',(package/'main.bbl').read_text()))
 if refs!=31:raise ValueError('Expected all29 prior references plus GraphGPS and EGNN')
 if figs!={str(i) for i in range(1,13)}|{'S'+str(i) for i in range(1,27)}:raise ValueError('Incomplete/incorrect figure numbering: '+str(figs))
 if tables!={str(i) for i in range(1,6)}|{'S'+str(i) for i in range(1,17)}:raise ValueError('Incomplete/incorrect table numbering: '+str(tables))
 pages=text.split('\f')
 for i,page in enumerate(pages):
  if page.strip() and page.strip().splitlines()[-1].strip()!=str(i+1):raise ValueError('Discontinuous combined page number')
 report={'status':'COMPLETE','label':'post_v1_exploratory','documents':reports,'figures':len(figs),'tables':len(tables),'references':refs,'no_final_tex_diagnostics':True,'continuous_pages':True,'outputs':{name:sha(package/name) for name in names.values()},'original_package':'manuscript/journal_revision; unchanged','source_status_sha256':sha(base/'summary/STATUS.json')}
 (package/'EXTENSION_BUILD_REPORT.json').write_text(json.dumps(report,indent=2)+'\n');return report
