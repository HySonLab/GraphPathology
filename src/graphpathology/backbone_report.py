"""Publication tables and a paired-effect plot from completed, locked comparisons."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .assets import sha
from .backbone_data import STUDY,record,write

LABEL={'resnet18':'ResNet18','phikon':'Phikon'}
HEAD={'GCN_shuffled':'GCN shuffled','SAGE_shuffled':'GraphSAGE shuffled','local_local':'local','GCN_spatial':'GCN','SAGE_spatial':'GraphSAGE','MEAN_spatial':'mean','idw_k16_p2_none':'IDW'}


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
                    for endpoint,multiplier,digits in [('delta_pearson',1,4),('relative_mse',100,1),('delta_mae',1,4)]:
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


def report(root):
    base=root/f'experiments/journal_revision/{STUDY}';status=json.loads((base/'ANALYSIS_STATUS.json').read_text())
    assert status['status']=='COMPLETE','No completed Phikon study: do not generate a final comparison table from expected outcomes.'
    summary=select_primary(pd.read_csv(base/'summary/PAIRED_SUMMARY.csv'));patient=select_primary(pd.read_csv(base/'summary/PATIENT_METRICS.csv'));paired=select_primary(pd.read_csv(base/'summary/PAIRED_PATIENT_EFFECTS.csv'))
    package=root/'manuscript/journal_revision/transformer_egnn';figures=package/'figures/backbone';figures.mkdir(parents=True,exist_ok=True)
    (package/'tables/BackboneComparison.tex').write_text(make_table(patient,summary));plot(paired,summary,figures)
    source=package/'source_data/backbone';source.mkdir(parents=True,exist_ok=True)
    import shutil
    for name in ['PATIENT_METRICS.csv','PAIRED_PATIENT_EFFECTS.csv','PAIRED_SUMMARY.csv','DEFINED_GENE_COUNTS.csv','COMMON_VALID_PATIENT_EFFECTS.csv']:shutil.copy2(base/'summary'/name,source/name)
    selected=json.loads((base/'SELECTED_FIT_MANIFEST.json').read_text());candidates=json.loads((base/'ALL_CANDIDATE_FIT_MANIFEST.json').read_text())
    rows=[]
    selected_keys={(r['backbone'],r['model'],r['control'],r['seed'],r['learning_rate']) for r in selected}
    for r in candidates:
        rows.append({'backbone':r['backbone'],'model':r['model'],'control':r['control'],'seed':r['seed'],'learning_rate':r['learning_rate'],'selected_for_evaluation':(r['backbone'],r['model'],r['control'],r['seed'],r['learning_rate']) in selected_keys,'status':r['status'],'parameters':r['parameters'],'selected_epoch':r['selected_epoch'],'epochs_ran':r['epochs_ran'],'optimizer_updates':r.get('optimizer_updates',14*r['epochs_ran']),'validation_mse':r['validation_mse'],'fit_seconds':r['fit_seconds'],'peak_rss_mb':r['peak_rss_mb'],'checkpoint':r['checkpoint'],'checkpoint_sha256':r['checkpoint_sha256']})
    fits=pd.DataFrame(rows);fits.to_csv(base/'MODEL_MANIFEST.csv',index=False);fits.to_csv(source/'MODEL_MANIFEST.csv',index=False)
    history_dir=base/'histories';history_dir.mkdir(exist_ok=True)
    for r in candidates:
        history=pd.DataFrame(r['history']) if 'history' in r else pd.read_csv(root/Path(r['source_fit_record']).with_name(Path(r['source_fit_record']).stem+'_history.csv'))
        assert len(history)==r['epochs_ran'] and np.isfinite(history[['training_mse','validation_mse']].to_numpy()).all()
        assert abs(history.iloc[r['selected_epoch']-1].validation_mse-r['validation_mse'])<1e-12
        history.to_csv(source/(r['run_id']+'_history.csv'),index=False);history.to_csv(history_dir/(r['run_id']+'_history.csv'),index=False)
    resource_table(root,package)
    title=['# Matched backbone comparison report','','Status: COMPLETE. ResNet18 selected heads/results: REUSED_VERIFIED; Phikon selected heads/results: COMPLETE. All 72 source-validation candidates are preserved; 36 validation-selected predictors are evaluated. This exploratory follow-up uses previously inspected patients, not fresh blind validation.','','Primary condition: historical distributed 20%, exact masks0/1/2; three training seeds averaged within patient. External predictions use the fixed C9 T1 contract once. Pearson, MSE and MAE are interpreted jointly.','','## Question A — graph context versus a matched local image head','']
    for backbone in ['resnet18','phikon']:
        title += [f'### {LABEL[backbone]}','', '| Cohort | Comparison | Median Δr | Median relative MSE (%) | Median ΔMAE | Improved patients (r/MSE/MAE) |','|---|---|---:|---:|---:|---|']
        for co in ['Development','Wu','Bassiouni']:
            for candidate in ['GCN_spatial','SAGE_spatial']:
                q=effects(summary,backbone,co,'image',candidate,'local_local');r=q.loc['delta_pearson'];m=q.loc['relative_mse'];a=q.loc['delta_mae'];title += [f"| {co} | {HEAD[candidate]} − local | {r.median_effect:+.4f} | {100*m.median_effect:+.1f} | {a.median_effect:+.4f} | {int(r.n_improved)}/{int(r.n_patients)}; {int(m.n_improved)}/{int(m.n_patients)}; {int(a.n_improved)}/{int(a.n_patients)} |"]
        title+=['','Opposite Pearson/error directions represent a tradeoff. Mean and shuffled controls are reported separately below; these primary contrasts alone do not establish topology specificity.','']
    title+=['## Question B — observed-residual reconstruction, then molecular support','','All contrasts below use the same anchors and IDW weights within a run. Gains versus local and versus IDW answer separate questions.','']
    for backbone in ['resnet18','phikon']:
        title += [f'### {LABEL[backbone]}','', '| Cohort | Candidate / reference | Median Δr | Median relative MSE (%) | Median ΔMAE | Improved patients (r/MSE/MAE) |','|---|---|---:|---:|---:|---|']
        for co in ['Development','Wu','Bassiouni']:
            for candidate,ref in [('GCN_spatial','local_local'),('SAGE_spatial','local_local'),('local_local','idw_k16_p2_none'),('GCN_spatial','idw_k16_p2_none'),('SAGE_spatial','idw_k16_p2_none')]:
                q=effects(summary,backbone,co,'exploratory_residual_anchor',candidate,ref);r=q.loc['delta_pearson'];m=q.loc['relative_mse'];a=q.loc['delta_mae'];title += [f"| {co} | {HEAD[candidate]} / {HEAD[ref]} | {r.median_effect:+.4f} | {100*m.median_effect:+.1f} | {a.median_effect:+.4f} | {int(r.n_improved)}/{int(r.n_patients)}; {int(m.n_improved)}/{int(m.n_patients)}; {int(a.n_improved)}/{int(a.n_patients)} |"]
        title+=['']
    wu=paired[(paired.cohort=='Wu')&paired.candidate.isin(['GCN_spatial','SAGE_spatial'])&paired.reference.isin(['local_local','idw_k16_p2_none'])&paired.route.isin(['image','exploratory_residual_anchor'])]
    wu.to_csv(base/'summary/WU_INDIVIDUAL_PRIMARY_EFFECTS.csv',index=False);wu.to_csv(source/'WU_INDIVIDUAL_PRIMARY_EFFECTS.csv',index=False)
    title += ['## Aggregation and topology controls','','These Bassiouni primary-condition contrasts assess specificity separately from the graph-versus-local question. Positive Pearson and negative errors favor the candidate. Mixed directions are tradeoffs.','','| Backbone | Route | Candidate / reference | Median Δr | Relative MSE (%) | ΔMAE | Improved patients (r/MSE/MAE) |','|---|---|---|---:|---:|---:|---|']
    for backbone in ['resnet18','phikon']:
        for route in ['image','exploratory_residual_anchor']:
            for candidate in ['GCN_spatial','SAGE_spatial']:
                for reference in ['MEAN_spatial',candidate.replace('spatial','shuffled')]:
                    q=effects(summary,backbone,'Bassiouni',route,candidate,reference);r=q.loc['delta_pearson'];m=q.loc['relative_mse'];a=q.loc['delta_mae']
                    title += [f"| {LABEL[backbone]} | {'Image' if route=='image' else 'Residual'} | {HEAD[candidate]} / {HEAD[reference]} | {r.median_effect:+.4f} | {100*m.median_effect:+.1f} | {a.median_effect:+.4f} | {int(r.n_improved)}/22; {int(m.n_improved)}/22; {int(a.n_improved)}/22 |"]
    common=select_primary(pd.read_csv(base/'summary/COMMON_VALID_PATIENT_EFFECTS.csv'))
    title += ['','## Common-valid-gene sensitivity','','Primary definitions remain unchanged. The following separately labelled sensitivity intersects the finite-correlation genes of each paired run before patient averaging. Gene counts are patient means over mask/training seeds, not extra biological replicates.','','| Backbone | Route | Candidate / reference | Primary Δr | Common-valid Δr | Mean valid genes, patient range |','|---|---|---|---:|---:|---|']
    for backbone in ['resnet18','phikon']:
        for route in ['image','exploratory_residual_anchor']:
            for candidate in ['GCN_spatial','SAGE_spatial']:
                for reference in (['local_local'] if route=='image' else ['local_local','idw_k16_p2_none']):
                    c=common[(common.backbone==backbone)&(common.cohort=='Bassiouni')&(common.route==route)&(common.candidate==candidate)&(common.reference==reference)];assert len(c)==22
                    primary=effects(summary,backbone,'Bassiouni',route,candidate,reference).loc['delta_pearson'].median_effect
                    title += [f"| {LABEL[backbone]} | {'Image' if route=='image' else 'Residual'} | {HEAD[candidate]} / {HEAD[reference]} | {primary:+.4f} | {c.common_delta_pearson.median():+.4f} | {c.common_valid_genes.min():.2f}–{c.common_valid_genes.max():.2f} |"]
    title += ['','## Controls, uncertainty and complete supporting results','','`PAIRED_SUMMARY.csv` includes all prespecified spatial/mean/shuffled, molecular, route and budget contrasts with descriptive 10,000-patient-resampling intervals in Bassiouni. `PAIRED_PATIENT_EFFECTS.csv` retains all three individual Wu effects and the single development patient; cohorts are never pooled. `DEFINED_GENE_COUNTS.csv` and the separately labelled common-valid-gene sensitivity preserve undefined-correlation information. Full per-run/per-gene data or checksum-verified prior pointers remain available. Fixed fusion, 10/40% coverage, identical-H placement and historical central-V2 are secondary conditions, not outcome-selected alternatives.','','## Provenance and limits','','The model manifest records every seed, state hash, selected/run epoch, update count and actual observed cost. Each configuration received the same two learning-rate trials, selected on mean TENX95 MSE across all three seeds before new selected outcomes. The validation-guided refinement and initial early-overfitting diagnostic remain documented. Head capacities are approximately matched within each backbone, not between encoders. Graph operators are fixed and encoders are frozen. Cross-backbone differences also include pretraining, feature dimension and effective field of view. The prior 48 advanced fits and original 30 minibatch fits remain intact, with their distinct update schedules. Original primary image-head binaries remain missing. All new full exports/states are author-local; no public release or journal submission occurred.']
    (base/'COMPARISON_REPORT.md').write_text('\n'.join(title)+'\n')
    write(base/'REPORT_STATUS.json',record(root,status='COMPLETE',model_configurations=36,new_candidate_fits=54,reused_candidate_fits=18,all_candidate_fits=72,all_planned_seeds=True,source_hashes={f.name:sha(f) for f in source.iterdir() if f.is_file()}))
    return {'status':'COMPLETE','table':str(package/'tables/BackboneComparison.tex')}

if __name__=='__main__':
    print(json.dumps(report(Path(__file__).resolve().parents[2]),indent=2))
