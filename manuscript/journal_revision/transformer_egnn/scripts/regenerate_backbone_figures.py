# Regenerate matched displays from recorded exports only; no training or research downloads.
import json,tempfile,shutil
from pathlib import Path
import numpy as np
import pandas as pd
STUDY='backbone_graph_comparison'
LABEL={'resnet18': 'ResNet18', 'phikon': 'Phikon'}
HEAD={'GCN_shuffled': 'GCN shuffled', 'SAGE_shuffled': 'GraphSAGE shuffled', 'local_local': 'local', 'GCN_spatial': 'GCN', 'SAGE_spatial': 'GraphSAGE', 'MEAN_spatial': 'mean', 'idw_k16_p2_none': 'IDW'}

def select_primary(frame):
    return frame[(frame.design=='historical_V1')&(frame.budget==.2)]

def effects(summary,backbone,cohort,route,candidate,reference):
    q=summary[(summary.backbone==backbone)&(summary.cohort==cohort)&(summary.route==route)&(summary.candidate==candidate)&(summary.reference==reference)].set_index('effect')
    assert len(q)==3,(backbone,cohort,route,candidate,reference)
    return q

def value(x,precision):
    return '--' if not np.isfinite(x) else f'{x:.{precision}f}'

def make_table(patient,summary):
    rows=[r'\begin{table}[p]',r'\centering',r'\caption{\textbf{Matched graph-versus-local and anchor-reference comparisons.} Primary 20\% distributed coverage uses the same three mask seeds and hidden sites, with three training seeds averaged within patient. Both backbones received the same two-rate validation-only selection opportunity. Entries separate image-only prediction (A) from observed-residual reconstruction (B). Each metric pair gives candidate/reference cohort medians; effects are medians of paired patient changes, so they need not equal differences of these marginal medians. Brackets give improved/valid patients for that endpoint. Relative MSE is computed within patient after seed averaging. External image outputs use the same frozen C9 T1 map once; no same-section molecular values enter A, and all B comparators use identical anchors. Bassiouni intervals appear in Fig.~\ref{fig:matched_backbone_effects}; Wu effects are reported individually in the source tables. Local denotes the non-graph MLP on the same backbone.}',r'\label{tab:matched_backbone}',r'\scriptsize\setlength{\tabcolsep}{2pt}\renewcommand{\arraystretch}{1.08}',r'\begin{tabular}{lrrrrrr}\toprule',r'Candidate/reference & $r_C/r_R$ & MSE$_C$/MSE$_R$ & MAE$_C$/MAE$_R$ & $\Delta r$ [$n$] & $\Delta$MSE\% [$n$] & $\Delta$MAE [$n$]\\\midrule']
    for panel,route in [('A. Image-only prediction','image'),('B. Observed-residual reconstruction','exploratory_residual_anchor')]:
        rows.append(r'\multicolumn{7}{l}{\textbf{'+panel+r'}}\\')
        for backbone in ['resnet18','phikon']:
            for co,n in [('Wu',3),('Bassiouni',22)]:
                rows.append(r'\multicolumn{7}{l}{'+f'{LABEL[backbone]}; {co}, $n={n}$'+r'}\\')
                pairs=[('GCN_spatial','local_local'),('SAGE_spatial','local_local')]
                if route!='image':pairs += [('local_local','idw_k16_p2_none'),('GCN_spatial','idw_k16_p2_none'),('SAGE_spatial','idw_k16_p2_none')]
                for candidate,reference in pairs:
                    cm,cc=candidate.split('_',1)
                    c=patient[(patient.backbone==backbone)&(patient.cohort==co)&(patient.route==route)&(patient.model==cm)&(patient.control==cc)]
                    if reference=='idw_k16_p2_none':r=patient[(patient.backbone=='molecular')&(patient.cohort==co)&(patient.route=='anchor_only')&(patient.model=='idw_k16_p2')]
                    else:r=patient[(patient.backbone==backbone)&(patient.cohort==co)&(patient.route==route)&(patient.model=='local')&(patient.control=='local')]
                    assert len(c)==len(r)==n
                    q=effects(summary,backbone,co,route,candidate,reference)
                    absolute=[value(c[k].median(),digits)+'/'+value(r[k].median(),digits) for k,digits in [('median_pearson',3),('mean_mse',4),('mean_mae',4)]]
                    paired=[]
                    for endpoint,multiplier,digits in [('delta_pearson',1,4),('relative_mse',100,2),('delta_mae',1,4)]:
                        v=q.loc[endpoint];number=value(v.median_effect*multiplier,digits);paired.append(number+f" [{int(v.n_improved)}/{int(v.n_patients)}]")
                    rows.append(f"{HEAD[candidate]}/{HEAD[reference]} & "+' & '.join(absolute+paired)+r'\\')
        rows.append(r'\midrule')
    rows[-1]=r'\bottomrule';rows += [r'\end{tabular}',r'\end{table}']
    return '\n'.join(rows)+'\n'

def plot(patient_effects,summary,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    rows=[('image','GCN_spatial','local_local','Image: GCN / local'),('image','SAGE_spatial','local_local','Image: GraphSAGE / local'),('exploratory_residual_anchor','GCN_spatial','local_local','Residual: GCN / local'),('exploratory_residual_anchor','SAGE_spatial','local_local','Residual: GraphSAGE / local'),('exploratory_residual_anchor','local_local','idw_k16_p2_none','Residual: local / IDW'),('exploratory_residual_anchor','GCN_spatial','idw_k16_p2_none','Residual: GCN / IDW'),('exploratory_residual_anchor','SAGE_spatial','idw_k16_p2_none','Residual: GraphSAGE / IDW')]
    fig,axes=plt.subplots(2,3,figsize=(11,7.8),sharey=True)
    for row_index,backbone in enumerate(['resnet18','phikon']):
        for ax,(metric,multiplier,xlabel) in zip(axes[row_index],[('delta_pearson',1,'Paired change in Pearson r'),('relative_mse',100,'Paired relative MSE change (%)'),('delta_mae',1,'Paired MAE change')]):
            for index,(route,candidate,reference,label) in enumerate(rows):
                q=patient_effects[(patient_effects.backbone==backbone)&(patient_effects.cohort=='Bassiouni')&(patient_effects.route==route)&(patient_effects.candidate==candidate)&(patient_effects.reference==reference)]
                assert len(q)==22
                vals=q.sort_values('patient')[metric].to_numpy(float)*multiplier;ax.plot(vals,np.repeat(index,len(vals)),'.',color='.65',ms=4,alpha=.8)
                s=effects(summary,backbone,'Bassiouni',route,candidate,reference).loc[metric];m=s.median_effect*multiplier;lo=s.ci_low*multiplier;hi=s.ci_high*multiplier
                ax.errorbar(m,index,xerr=[[m-lo],[hi-m]],fmt='o',color='#0072b2' if row_index==0 else '#d55e00',ms=4,capsize=2,lw=1.2)
            ax.axvline(0,color='black',lw=.7);ax.axhline(1.5,color='.8',lw=.5);ax.axhline(3.5,color='.8',lw=.5);ax.set_xlabel(xlabel);ax.grid(axis='x',alpha=.15)
            ax.set_yticks(range(7),[r[3] for r in rows]);ax.set_ylim(6.6,-.6)
        axes[row_index,0].set_title(LABEL[backbone]+' — Bassiouni, 22 patients',loc='left',fontsize=11)
    fig.tight_layout()
    for extension in ['pdf','png']:fig.savefig(out/f'BackbonePatientEffects.{extension}',bbox_inches='tight',**({'dpi':180} if extension=='png' else {}))
    plt.close(fig)

def resource_table(root,package):
    base=root/f'experiments/journal_revision/{STUDY}'
    chosen=json.loads((base/'SELECTED_FIT_MANIFEST.json').read_text());allfits=json.loads((base/'ALL_CANDIDATE_FIT_MANIFEST.json').read_text())
    lines=[r'\begin{table}[p]',r'\centering',r'\caption{\textbf{Matched-backbone fitting and validation selection.} Each configuration used the same two candidate learning rates (0.001 and 0.0001) with seeds 0/1/2. One rate was selected per backbone/configuration using mean best qualifying TENX95 MSE across all three seeds; no evaluation outcome selected it. Epoch triplets retain seed order 0/1/2. Updates count all minibatches, including the final partial batch. Times sum the three selected fits; peak RSS is their maximum and includes loaded data. These observed execution costs are not a controlled speed comparison across encoders or historical execution sessions. All 72 candidate histories, validation losses, state/configuration/input hashes and selected flags are machine-readable.}',r'\label{tab:backbone_resources}',r'\scriptsize\setlength{\tabcolsep}{3pt}',r'\begin{tabular}{llrrrrrr}\toprule',r'Backbone & Head/control & Parameters & Rate & Best epochs & Run epochs & Updates & Seconds / MiB\\\midrule']
    for backbone in ['resnet18','phikon']:
        for model,control in [('local','local'),('MEAN','spatial'),('GCN','spatial'),('SAGE','spatial'),('GCN','shuffled'),('SAGE','shuffled')]:
            rows=sorted([r for r in chosen if r['backbone']==backbone and r['model']==model and r['control']==control],key=lambda r:r['seed']);assert len(rows)==3
            name={'local':'local','MEAN':'mean','GCN':'GCN','SAGE':'GraphSAGE'}[model]+(' shuffled' if control=='shuffled' else '')
            lines.append(f"{LABEL[backbone]} & {name} & {rows[0]['parameters']:,} & {rows[0]['learning_rate']:.4f} & "+'/'.join(str(r['selected_epoch']) for r in rows)+' & '+'/'.join(str(r['epochs_ran']) for r in rows)+f" & {sum(r.get('optimizer_updates',14*r['epochs_ran']) for r in rows)} & {sum(r['fit_seconds'] for r in rows):.1f} / {max(r['peak_rss_mb'] for r in rows):.0f}"+r'\\')
    lines += [r'\bottomrule\end{tabular}',r'\end{table}']
    (package/'tables/BackboneResources.tex').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':
    package=Path(__file__).resolve().parents[1];source=package/'source_data/backbone'
    patient=select_primary(pd.read_csv(source/'PATIENT_METRICS.csv'));summary=select_primary(pd.read_csv(source/'PAIRED_SUMMARY.csv'));paired=select_primary(pd.read_csv(source/'PAIRED_PATIENT_EFFECTS.csv'))
    (package/'tables/BackboneComparison.tex').write_text(make_table(patient,summary))
    plot(paired,summary,package/'figures/backbone')
    with tempfile.TemporaryDirectory() as directory:
        root=Path(directory);base=root/f'experiments/journal_revision/{STUDY}';base.mkdir(parents=True)
        for name in ['ALL_CANDIDATE_FIT_MANIFEST.json','SELECTED_FIT_MANIFEST.json']:shutil.copy2(source/name,base/name)
        resource_table(root,package)
    print('Regenerated paired table, patient-effect figure and fitting-resource table from recorded exports.')
