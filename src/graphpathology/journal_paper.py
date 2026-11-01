"""Assemble the separate journal manuscript using actual frozen journal result tables."""
from pathlib import Path
import json
import re
import shutil
import subprocess
import pandas as pd
from .assets import sha


def esc(x):return str(x).replace('_',r'\_').replace('%',r'\%')

def table(caption,label,headers,rows,widths=None):
    cols=widths or ('l'+'r'*(len(headers)-1))
    return '\n'.join([r'\begin{table}[!htbp]',r'\caption{'+caption+'}',r'\label{'+label+'}',r'\centering\footnotesize',r'\begin{tabular}{@{}'+cols+'@{}}',r'\toprule',' & '.join(headers)+r' \\',r'\midrule',*[' & '.join(row)+r' \\' for row in rows],r'\bottomrule',r'\end{tabular}',r'\end{table}',r'\FloatBarrier',''])


def scientific_sections(root,package):
    src=package/'source_data/journal';stats=pd.read_csv(src/'ORIGINAL_V1_PATIENT_BOOTSTRAP.csv')
    paired=pd.read_csv(src/'WU_PAIRED_PATIENT_EFFECTS.csv');wu=pd.read_csv(src/'WU_PATIENT_METRICS.csv')
    # All numerical fields come from computed exports, without selecting an external winner.
    r=[]
    for comparison,name in [('anchors_vs_image','Anchors minus image'),('fusion_vs_anchors','Fusion minus anchors')]:
        for metric,label,scale in [('delta_pearson',r'$\Delta r$',1),('relative_mse',r'Relative MSE (\%)',100),('delta_mae',r'$\Delta$ MAE',1)]:
            x=stats[(stats.cohort=='Bassiouni')&(stats.coverage==.2)&(stats.regime=='V1')&(stats.contrast==comparison)&(stats.effect==metric)].iloc[0]
            r.append([name,label,f'{x.median_effect*scale:.4f}',f'[{x.ci_low*scale:.4f}, {x.ci_high*scale:.4f}]',f'{x.fraction_improved*22:.0f}/22'])
    t3=table(r'\textbf{Post-v1 patient uncertainty reanalysis of original Bassiouni exports.} Median paired effects at 20\% V1 coverage, with 10,000 patient-resampling percentile 95\% intervals (seed 20261101). Negative error changes improve; positive correlation changes improve. These are descriptive post-v1 intervals, not simultaneous bands or p-values. Anchors versus image is an unequal-information contrast; fusion versus anchors tests incremental image contribution. All original point estimates and source bytes are unchanged.', 'tab:journal_uncertainty',['Contrast','Effect','Median','95\% interval','Improved'],r)
    (package/'tables/JournalTable3.tex').write_text(t3)
    patients=['1160920F','CID4290','CID44971'];rows=[]
    from .journal_methods import candidates
    for config in candidates():
        method=config['method'];a=wu[(wu.experiment=='A')&(wu.budget==.2)&(wu.method==method)].set_index('patient')
        ref=wu[(wu.experiment=='A')&(wu.budget==.2)&(wu.method=='idw_k16_p2')].set_index('patient')
        rows.append([esc(method)]+[f'{a.loc[p,"median_pearson"]:.3f}' for p in patients]+[f'{100*(a.loc[p,"mean_mse"]-ref.loc[p,"mean_mse"])/ref.loc[p,"mean_mse"]:+.1f}' for p in patients])
    s10=table(r'\textbf{All 20 matched-information estimators in Wu at 20\% V1 coverage.} Each patient value averages ten original random trajectories. P1=1160920F, P2=CID4290, P3=CID44971. Relative MSE is candidate minus fixed IDW k16/p2, divided by that reference; negative values improve. No external winner was selected. All 10/20/40\% results, MAE and individual paired effects accompany the CSV exports. Gaussian bandwidth uses observed-anchor geometry only.', 'tab:journal_baselines',['Estimator',r'$r$ P1',r'$r$ P2',r'$r$ P3',r'MSE\% P1',r'MSE\% P2',r'MSE\% P3'],rows)
    (package/'tables/JournalTableS10.tex').write_text(s10)
    rows=[]
    for patient in patients:
        for regime in ['V1','V2']:
            vals=[]
            for method in ['idw_k16_p2','uniform_k16']:
                x=wu[(wu.patient==patient)&(wu.regime==regime)&(wu.budget==.2)&(wu.method==method)].iloc[0]
                vals += [f'{x.median_pearson:.3f}',f'{x.mean_mse:.4f}']
            rows.append([patient,regime]+vals)
    (package/'tables/JournalTableS11.tex').write_text(table(r'\textbf{Operator by historical protocol at 20\% coverage.} IDW and uniform operators share each original mask and scored set. V2 averages geometries equally after averaging seeds. The V1/V2 contrast still changes spatial support and is not a causal geometry estimate. All budgets are supplied in source data.', 'tab:journal_factorial',['Patient','Protocol','IDW $r$','IDW MSE','Uniform $r$','Uniform MSE'],rows))
    rows=[]
    for patient in patients:
        for method in ['idw_k16_p2','uniform_k16']:
            x=paired[(paired.patient==patient)&(paired.budget==.2)&(paired.origin=='new_B2_identical_targets')&paired.contrast.str.startswith(method)].iloc[0]
            rows.append([patient,'IDW' if method.startswith('idw') else 'Uniform',f'{x.delta_pearson:+.4f}',f'{100*x.relative_mse:+.1f}',f'{x.delta_mae:+.4f}'])
    (package/'tables/JournalTableS12.tex').write_text(table(r'\textbf{Identical-target anchor-placement effects at 20\% coverage in Wu.} Region-exclusion minus distributed anchors, with the same H and total anchor count. Each row averages ten seeds. The holdout H is fixed across conditions and budgets. All three patients deteriorate in Pearson, MSE and MAE with both estimators. This central-region result is post-v1 robustness evidence in previously inspected patients, not new blind validation.', 'tab:journal_identical',['Patient','Operator',r'$\Delta r$',r'Relative MSE\%',r'$\Delta$ MAE'],rows))
    # New learned results must exist; never print placeholders or infer missing training outcomes.
    fit=root/'experiments/journal_revision/results/hest_learned';status=json.loads((fit/'EVALUATION_STATUS.json').read_text())
    if status['status']!='COMPLETE_DEVELOPMENT_ONLY' or status['fit_models']!=30:raise ValueError('All preregistered learned fits must be complete')
    learned=pd.read_csv(src/'HEST_SEED_METRICS.csv');learned=learned[learned.budget==.2];rows=[]
    configs=[('local','local')]+[(m,c) for m in ['MEAN','GCN','SAGE'] for c in ['spatial','self','shuffled']]
    for model,control in configs:
        for seed in range(3):
            a=learned[(learned.model==model)&(learned.control==control)&(learned.training_seed==seed)]
            im=a[a.route=='image'].iloc[0];res=a[a.route=='exploratory_residual_anchor'].iloc[0]
            rows.append([esc(model+' '+control),str(seed),f'{im.median_pearson:.3f}',f'{im.mean_mse:.4f}',f'{res.median_pearson:.3f}',f'{res.mean_mse:.4f}'])
    (package/'tables/JournalTableS13.tex').write_text(table(r'\textbf{All new controlled neural fits at 20\% TENX99 coverage.} Each value averages three identical original masks. Models fit NCBI783/NCBI785 and use TENX95 alone for stopping. Seeds are training variability, not patients. Residual is the exploratory image plus interpolated observed-residual predictor. All budgets, fixed 20\% fusion and anchor-only values are in CSVs. No architecture or seed was selected using TENX99. Nominal parameter counts are 72,114 (local/mean/GCN) and 72,075 (SAGE); SAGE self-only duplicates self inputs and has a smaller effective function class.', 'tab:journal_neural',['Configuration','Seed','Image $r$','Image MSE','Residual $r$','Residual MSE'],rows))
    pure=[]
    for model in ['GCN','SAGE']:
        a=learned[(learned.model==model)&(learned.control=='spatial')&(learned.route=='image')]
        pure.append(model+f' mean correlation {a.median_pearson.mean():.3f} and MSE {a.mean_mse.mean():.4f}')
    local=learned[(learned.model=='local')&(learned.route=='image')]
    neural_text=('All 30 preregistered new heads completed training. At 20\\% coverage, the local image head had training-seed mean correlation '
                 f'{local.median_pearson.mean():.3f} and MSE {local.mean_mse.mean():.4f}; '+('; '.join(pure))+'. These values describe one previously used development patient, not population generalization. All spatial, self-only and shuffled configurations and seeds are shown, rather than a selected best model (Fig.~\\ref{fig:journal_learned}; Supplementary Table~\\ref{tab:journal_neural}). Local, mean-self and GCN-self predictions agreed bitwise for every matched seed, as required by their identical inputs and parameterization. The fixed-fusion and exploratory residual-anchor results are reported separately from image-only prediction and from the missing historical fitted heads.')
    result_text=r'''
\subsection{Post-v1 matched-information and identical-target checks}
The following analyses were specified in a committed protocol before new comparator outcomes and remain separate from the immutable arXiv v1 record. The Wu and Bassiouni studies had already been inspected; no new blind external cohort is claimed. Exact public Wu archives and pinned HEST inputs were recovered with checksum verification. GEO HTTPS/API requests failed and direct anonymous FTP returned empty bodies for the selected Bassiouni files; checksum checks rejected those transfers. New Bassiouni matched-estimator, identical-target and spatial-map analyses are therefore unavailable, rather than inferred from old summaries.

In the three Wu evaluation patients, we evaluated all 20 prespecified anchor-only estimators, spanning nearest-anchor, uniform, inverse-distance and Gaussian reconstruction, on identical original V1 masks at 10/20/40\% coverage (Supplementary Table~\ref{tab:journal_baselines}). The grid is a descriptive same-information comparison; it was not used to select an external winner. Separately, IDW k16/p2 and uniform k16 were crossed with every original random V1/V2 mask and scoring support (Fig.~\ref{fig:journal_controls}a; Supplementary Table~\ref{tab:journal_factorial}). Changing the estimator is distinguishable from changing the historical protocol, but the original V1/V2 contrast remains compound.

The identical-target experiment fixed a coordinate-defined central region R and a single hidden set H, then compared equal numbers of distributed anchors and anchors restricted outside R. At 20\% coverage, region exclusion reduced IDW correlation by 0.106--0.245 and increased relative MSE by 20.1--45.4\% across the three patients. Uniform reconstruction also deteriorated in each patient: correlation fell by 0.093--0.188 and relative MSE increased by 21.7--38.9\% (Fig.~\ref{fig:journal_controls}b; Supplementary Table~\ref{tab:journal_identical}). Both conditions scored the same H, unlike the historical V1/V2 contrast. These controlled placement effects concern a single predeclared central geometry in three previously inspected sections; they do not establish a universal missing-region effect. Nearest-anchor distances and occupied-bin errors/support are retained for every budget. Predetermined real-data gene maps include constant and poorly reconstructed fields rather than selected favorable examples (Supplementary Figs.~\ref{fig:journal_map_v1} and \ref{fig:journal_map_v2}).

Across 2,700 evaluated Wu runs, all 5,400 NaN/random substitutions of every unobserved label were invariant. Initial replay identified physically negligible coordinate recalibration differences that changed neighbor ties. Preserved label-free mask coordinates, independently checked against source-derived coordinates within $10^{-8}$~\um{}, restored all 360 original random-reference runs (54,000 per-gene Pearson/MSE/MAE values) within the declared $10^{-10}$ tolerances. Raw Spearman ranks remained sensitive to floating-point near ties and are explicitly excluded from numerical-exactness and improvement claims; the diagnostic and raw exports are retained.

\begin{figure}[!htbp]
\centering\includegraphics[width=\textwidth]{figures/journal/JournalFigure7.pdf}
\caption{\textbf{Post-v1 estimator and placement controls in previously inspected Wu patients.} (a) Both operators share each historical mask and scored set; V2 values average the three geometries. (b) Region-exclusion minus distributed correlation on the identical H, with equal total anchor count and ten random seeds. The three individual patients are shown; masks and genes are not biological replicates. All points use 20\% coverage and the frozen 50-gene panel.}
\label{fig:journal_controls}
\end{figure}
\FloatBarrier

\subsection{Post-v1 patient uncertainty preserved the original effects}
The new uncertainty analysis resamples patients from the original frozen patient-level exports, not individual genes, spots, masks or model seeds. At 20\% V1 coverage in Bassiouni, the median paired anchor-versus-image correlation gain was 0.329 (95\% patient percentile interval 0.304--0.438). Fusion minus anchors instead had a median correlation change of $-0.0108$ (interval $-0.0187$ to $-0.0081$), with median relative MSE increase 10.7\% (6.6--14.5\%). These post-v1 descriptive intervals reinforce the original direction without converting an unequal-information image comparison into an algorithmic benchmark (Table~\ref{tab:journal_uncertainty}). All three Wu patient effects are tabulated individually; no primary three-patient bootstrap claim is made. The 10,000 resamples use patient units and fixed seed 20261101, and cohorts are never pooled.
\input{tables/JournalTable3}

\subsection{New controlled neural development fits were distinct from historical models}
'''+neural_text+r'''
\begin{figure}[!htbp]
\centering\includegraphics[width=\textwidth]{figures/journal/JournalFigure8.pdf}
\caption{\textbf{New journal-v2 neural development analyses.} Three training seeds for each local, mean, GCN and GraphSAGE configuration, including self-only and node-shuffled controls; all evaluated on TENX99, a single previously used development patient. Points average three fixed masks at 20\% coverage. Dashed lines give identical-mask IDW reference values. Top: image-only prediction; bottom: the preregistered exploratory residual-anchor baseline. Color identifies training seed, not an independent patient. These newly trained weights and caches are not the missing original states and were never tuned on Wu/Bassiouni.}
\label{fig:journal_learned}
\end{figure}
\FloatBarrier
'''
    return result_text


def assemble(root,package):
    results=scientific_sections(root,package)
    # Standalone article uses stable explicit SI numbering, as did the original author sources.
    for label,value in [('tab:journal_baselines','S10'),('tab:journal_factorial','S11'),('tab:journal_identical','S12'),('tab:journal_neural','S13'),('fig:journal_map_v1','S13'),('fig:journal_map_v2','S14')]:results=results.replace(r'\ref{'+label+'}',value)
    article=(root/'manuscript/revised/article_content.tex').read_text()
    abstract='''The cost of spatial transcriptomics motivates reconstruction from sparse molecular observations. GraphPathology separates coordinate-defined histological context from interpolation of measured expression anchors. The original evaluation covered a fixed 50-gene panel in 25 patients from two breast cancer studies; at 20\\% distributed coverage, molecular reconstruction achieved median patient correlations of 0.359 and 0.416. We preserved those results and conducted separately preregistered post-v1 analyses. Twenty anchor-only estimators were compared under matched information in three previously inspected Wu patients. With identical hidden targets and equal anchor counts, excluding a central region reduced correlation and increased error in every patient for both fixed estimators. Patient-resampling intervals from original 22-patient exports supported the original negative fixed-fusion result. New local, mean, GCN and GraphSAGE heads were trained with matched spatial, self-only and shuffled controls, and evaluated as one-patient development analyses rather than recovered historical models or fresh validation. New Bassiouni raw-data comparisons remained blocked by failed source transfers. The evidence supports useful distributed reconstruction while limiting claims about morphology, geometry and generalization. Fivefold fewer measured locations represents conditional variable-cost arithmetic, not demonstrated total assay savings.'''
    article=re.sub(r'(\\section\*\{Abstract\}\n).*?(\n\n\\noindent)',lambda m:m.group(1)+abstract+m.group(2),article,count=1,flags=re.S)
    article=article.replace(r'\section{Results}',r'\noindent\textit{Journal revision: original arXiv v1 results are retained below. Added analyses are explicitly post-v1 and use previously inspected cohorts.}'+'\n\n'+r'\section{Results}',1)
    article=article.replace(r'\section{Discussion}',results+'\n'+r'\section{Discussion}',1)
    article=article.replace('A future factorial study holding the interpolation rule and evaluation support fixed would be needed to estimate the isolated contribution of missing-region geometry.', 'The post-v1 operator and identical-target analyses now address these distinctions in the three Wu patients, but analogous new Bassiouni analyses remain blocked by source acquisition. They do not retroactively make the original V1/V2 contrast a pure geometry intervention.')
    article=article.replace('We report descriptive effects and patient-support counts without new hypothesis tests, spot-level significance claims or post-hoc population confidence intervals.', 'The original v1 record reported descriptive effects and patient-support counts without hypothesis tests or spot-level significance claims. The added journal analysis reports separately labelled descriptive patient bootstrap intervals; it does not reinterpret seeds or spots as biological replicates.')
    article=article.replace('The original public datasets remain available from their source repositories.', 'Canonical source links are listed above. During this revision, the selected GEO transfers returned errors or empty bodies despite advertised nonzero sizes; new Bassiouni raw-data analyses were not performed.')
    article=article.replace('the schematic spatial examples are synthetic.', 'the original schematic spatial examples are synthetic. The added journal maps instead show actual Wu counts, anchors, predictions and errors, with no inferred tissue annotations.')
    methods=r'''
\subsection{Post-v1 journal protocols, estimators and neural controls}
A committed protocol preceded new comparator outcomes. Original arXiv sources, numerical targets and archived checkpoint records were immutable. Public data and encoder weights were downloaded into a new ignored cache, with source terms, full SHA256 records and published Zenodo MD5 checks. New learned states were labelled journal-v2 and stored separately. No paid computing service was used.

Experiment A used all ten original random V1 trajectories at each of 10/20/40\% coverage for nearest anchor; uniform k=4,8,16,32; IDW at those k with p=1,2,3; and Gaussian k=8,16,32. We retained cKDTree ordering, six-decimal distances and a $10^{-6}$~\um{} inverse-weight floor. Gaussian bandwidth was the median observed-anchor distance to its kth other observed neighbor, with k limited by available anchors; it never used hidden expression. Every declared estimator is reported. Experiment B1 held each original mask and scored set fixed while changing IDW k16/p2 versus uniform k16, for both V1 and all three V2 geometries.

Experiment B2 used only the predeclared central PCA-coordinate region R (40\% of locations). H was a fixed random 20\% subset of R (about 8\% of all locations), selected using SHA256(``H|patient'') first eight bytes, little-endian. Distributed anchors were drawn from all non-H sites; region-exclusion anchors from outside R. Seeds 0--9 gave nested prefixes with equal total counts at every budget. R, H, barcodes and anchor trajectories were persisted before scoring. Both estimators scored exactly H; no expression-based geometry or exclusion was used. Seed means preceded geometry means and patient pairing. Both primary error and correlation effects and nearest-anchor distances were exported, with occupied-bin errors and support.

Uncertainty resampled whole paired patients 10,000 times (seed 20261101) and took percentile 95\% intervals for the median paired correlation change, MAE change and relative MSE. Fractions improved used the same finite paired patients, whose number is explicit. No cohorts were pooled and no p-values or simultaneous-band claims were made. Wu's three individual effects are primary; Bassiouni intervals reanalyze original exports and cannot replace unavailable new predictions. Raw Spearman was retained descriptively with explicit near-tie sensitivity, not used for numerical-exactness or improvement claims.

New neural fits used the original HEST 50-gene panel and C0 training normalization, frozen ImageNet ResNet18 embeddings with the original PIL resize/crop transform, and a feature scaler fit only on NCBI783/NCBI785. Coordinate-only Delaunay edges were capped at 250~\um{} within each section. A single learned linear/GELU/dropout/readout layer followed mean or symmetric GCN aggregation; GraphSAGE concatenated self and neighbor mean. Widths were 128 for local/mean/GCN and 67 for SAGE (72,114 versus 72,075 parameters). Spatial, self-only and degree-preserving node-shuffled controls shared each family's nominal parameters, labels and training settings. SAGE self-only concatenated X with X and has a smaller effective function class. AdamW used lr $10^{-3}$, decay $10^{-4}$, batch size 512, dropout0.1, maximum 80 epochs, patience 10 and minimum validation improvement $10^{-4}$; TENX95 alone selected stopping. Seeds 0/1/2 were all retained. Deterministic float32 CPU execution used four torch/BLAS threads; caches and checkpoints have independent new hashes.

TENX99 evaluation used three exact original random masks per budget. Alongside image and IDW predictions, we evaluated fixed $0.2\,g_\theta+0.8\,\mathrm{IDW}(y_{\mathcal O})$ and the separately preregistered exploratory residual baseline
\[
\widehat y_i=g_\theta(i)+\sum_{j\in\mathcal O}w_{ij}\big[y_j-g_\theta(j)\big].
\]
Only observed residual vectors entered interpolation; replacing every unobserved label with NaNs or random values had to leave predictions unchanged. This is residual interpolation, not a novelty claim. No external cohort selected a new model, coefficient or seed, and no biological confidence interval was inferred from one patient's masks or training seeds.
'''
    article=article.replace(r'\section{Ethics}',methods+'\n'+r'\section{Ethics}',1)
    (package/'article_content.tex').write_text(article)
    supplement=(root/'manuscript/revised/supplementary_content.tex').read_text()
    supplement=supplement.replace('No confidence intervals or significance tests have been added.', 'The original v1 displays contain no confidence intervals or significance tests. Separate post-v1 patient bootstrap intervals are now reported in main Table 3 and the journal CSVs; no significance tests were added.')
    supplement=supplement.replace('twelve supplementary figures and nine supplementary tables', 'the original twelve supplementary figures and nine supplementary tables, followed by separately labelled post-v1 analyses, two real-data maps and four new tables')
    supplement+=r'''
\clearpage
\section{Post-v1 journal analyses: source-data and limitations}
New source files are under source\_data/journal and results under experiments/journal\_revision. These analyses reuse previously inspected Wu patients and original Bassiouni exports. New Bassiouni prediction/map analyses remain blocked by failed canonical-source transfers; no old summary substitutes for a new prediction. Primary real-data replay passes per-gene Pearson/MSE/MAE tolerances; raw Spearman near-tie sensitivity is explicit. Every new Wu run was invariant to corruption of every unobserved label. The B2 H and anchor trajectories are label-free, persisted and hashed. All estimators, budgets, patient effects, bin supports and new training seeds are retained in machine-readable exports.
\input{tables/JournalTableS10}
\input{tables/JournalTableS11}
\input{tables/JournalTableS12}
\input{tables/JournalTableS13}
\clearpage
\begin{figure}[!htbp]
\centering\includegraphics[width=\textwidth]{figures/journal/JournalMap_V1.pdf}
\caption{\textbf{Predetermined real Wu V1 gene maps.} Patient 1160920F, 20\% random seed0, original mask, fixed IDW k16/p2. Genes ERBB2, EPCAM, FOXA1 and POSTN were designated before comparisons. Truth/prediction use shared per-gene scales over scored sites; errors have separate scales. Anchors are crosses and gray sites are other registered locations. Standardized log1p counts are not absolute transcript abundance. Constant scored fields are labelled rather than omitted or given misleading contrast. Source: Wu et al., Zenodo 4739739, CC-BY4.0; no tumor/cell annotations are inferred.}
\label{fig:journal_map_v1}
\end{figure}
\FloatBarrier
\clearpage
\begin{figure}[!htbp]
\centering\includegraphics[width=\textwidth]{figures/journal/JournalMap_V2.pdf}
\caption{\textbf{Predetermined real Wu V2 gene maps.} Same patient, genes, 20\% budget and seed0, central original V2 geometry 1, historical uniform k16 operator. Colored truth/prediction/error points are the original forbidden scoring region; anchor crosses lie outside. Scales and source terms match the preceding figure. FOXA1 has a constant scored field; numerical residuals do not establish biological recovery. These maps retain unsuccessful/constant examples and are not independent validation. The requested BASS\_P01 map was not generated because its raw inputs could not be acquired.}
\label{fig:journal_map_v2}
\end{figure}
\FloatBarrier
'''
    (package/'supplementary_content.tex').write_text(supplement)
    return {'abstract_words':len(abstract.split()),'new_main_figures':2,'new_main_tables':1,'new_supplement_figures':2,'new_supplement_tables':4}


def build(root,output):
    """Compile in a new staging directory; historical publisher sources are never written."""
    from .publication import compile_documents
    source=root/'manuscript/journal_revision';output=Path(output).resolve()
    if output.exists():raise FileExistsError('Use a new journal build staging directory')
    shutil.copytree(source,output,ignore=shutil.ignore_patterns('*.aux','*.log','*.out','*.bbl','*.blg','*.fls','*.pass*.txt','GraphPathology_*Journal_Revision.pdf'))
    report=compile_documents(output)
    text=subprocess.check_output(['pdftotext','-layout',str(output/'main.pdf'),'-'],text=True)
    if '??' in text or 'PLACEHOLDER' in text:raise ValueError('Unresolved or placeholder manuscript content')
    figures=set(re.findall(r'\bFigure\s+(S?\d+)\s*:',text));tables=set(re.findall(r'\bTable\s+(S?\d+)\s*:',text))
    if figures!={str(n) for n in range(1,9)}|{'S'+str(n) for n in range(1,15)}:raise ValueError('Journal figure numbering mismatch: '+str(figures))
    if tables!={'1','2','3'}|{'S'+str(n) for n in range(1,14)}:raise ValueError('Journal table numbering mismatch: '+str(tables))
    for short,name in [('main','Combined_Manuscript'),('article_only','Article_Only'),('supplementary','Supplementary')]:
        target=source/('GraphPathology_'+name+'_Journal_Revision.pdf');shutil.copyfile(output/(short+'.pdf'),target)
    result={'status':'BUILT_AND_TEXT_CHECKED','documents':report,'figures':len(figures),'tables':len(tables),'new_results':'post-v1, previously inspected cohorts; new Bassiouni raw experiments blocked',
            'source_hashes':{str(p.relative_to(source)):sha(p) for p in source.rglob('*') if p.is_file() and p.suffix in ['.tex','.csv','.gz','.pdf']},'build_staging':str(output)}
    (source/'JOURNAL_BUILD_REPORT.json').write_text(json.dumps(result,indent=2)+'\n');return result
