"""Final readable report and legends; scientific values come from exported tables."""
from synthesis import *

def main():
 t=pd.read_csv(M/'tables/Table2_reconstruction.csv');p=t[(t.cohort=='Bassiouni second replication')&(t.coverage==.2)&(t.regime=='V1')]
 report='''# C10 — Second replication and manuscript package

**PASS_LIMITED_SECOND_REPLICATION**
**Recommendation: READY_WITH_LIMITED_SECOND_REPLICATION**

The fresh GSE210616 cohort contributes 22 TNBC patients and 30,003 retained tissue locations, using the first numbered section for every published patient. All 50 original primary genes are present. No C10 patient was used for training, calibration or selection. Public H&E inputs contain only 24–29 native pixels per 112-µm patch; this limits conclusions about the image baseline and fusion.

## Primary result: 20% distributed coverage

'''+p[['method','median_pearson','mean_mse','mean_mae','patient_support_vs_M0_all_metrics']].to_markdown(index=False)+'''

M0 is H&E only, M1 is sparse-ST interpolation, and M3 is fixed fusion. Each value is a median across patient summaries; random seeds are averaged within patients. All 22 patients met the prespecified meaningful ST criterion at all three budgets: Pearson gain at least 0.02, relative MSE reduction at least 2%, and non-worsening MAE. All 22 improved correlation and MSE between 10% and 40% V1 coverage.

Fusion improved 20% V1 correlation in only 3/22 patients and improved MSE or MAE in none. The median paired correlation change was −0.010830 and the median relative MSE increase was 10.7078%. This is negative evidence for the tested fixed fusion rule, not proof that morphology can never help. M3 deliberately reuses the frozen Wu calibration object; it was not refit on C10.

## Coverage, geometry and distance

At 10%, 20% and 40% V1 coverage, ST median patient correlations were 0.380374, 0.415862 and 0.465266; MSE values were 0.100076, 0.092256 and 0.085249. MAE was 0.210553, 0.203274 and 0.196547.

At 20% V2, ST correlation/MSE/MAE were 0.170747/0.149598/0.262427. Every patient worsened in all three metrics relative to V1. The median paired correlation change was −0.266923 and relative MSE increase was 52.8206%. This contrast includes the scoring population and the frozen interpolation family as well as missing-region geometry. C5 block masks are a different development protocol.

Distance replication is partial. Farthest occupied-bin MSE exceeded nearest-bin MSE in 15/22 V1 and 20/22 V2 patients. Strict increases across every occupied bin occurred in only 6/22 and 13/22 patients. Sparse remote bins, tissue composition and changing patient support constrain pooled curves. No universal distance cutoff is claimed. Complete counts and eligible bin correlations are retained in C10B_DISTANCE_TO_ANCHOR.csv and Supplementary Table 7.

Secondary spatial FPS effects are reported in MANUSCRIPT/data/C10_FPS_effects.csv; they do not reopen sampling optimization. Anchor-expression and anchor-position permutations degraded reconstruction on the prespecified control patient, whereas image permutation produced much smaller changes. Matched control values are in C10B_NEGATIVE_CONTROLS.csv. Descriptive per-gene gains are in C10B_GENE_PATIENT_GAINS.csv and C10B_GENE_GAIN_SUMMARY.csv; they are not biological interpretations or outcome-selected gene panels.

## Integrity

The audit reproduced 8,766 method predictions and 35,064 per-gene metric arrays exactly, replayed 968 label-free mask trajectories, and verified 17,532 recorded hidden-label corruption comparisons. Raw inputs, image caches, parent models, scalers, barcodes and the pre-C10 registry prefix matched their hashes. The protocol and mask commits precede target access. Affine M1/M2 equivalence was checked at every inference, and duplicate M2 metric rows were omitted.

C0–C9 were verified read-only through the full chained manifests; their experiments were neither rerun nor retuned. Before/after verification logs are under C10A/logs and C10B/logs. Original registry bytes are preserved. Qualification, downloads, physical-scale checks, resource records and exact run identities are retained with the scientific locks.

## Manuscript handoff

The package contains seven main and ten supplementary figures, each as a 300-dpi PNG and SVG, with source CSVs and figure-specific plotting entrypoints. Four main and eight supplementary tables are available as CSV and Markdown; the large per-gene CSV is compressed. The draft, legends, citation worklist, source index, numerical ledger and claim-evidence matrix retain the negative findings and separate 40-gene HER2ST transportability development from 50-gene fresh validation.

The manuscript is ready for human writing under the stated limitations. It does not establish general H&E benefit, essential graph topology, optimal sampling, whole-transcriptome recovery, cost savings, clinical utility or causal biology. Native-resolution external H&E remains a gap for stronger morphology claims, but does not block the limited sparse-ST manuscript. Stop after C10.
'''
 (OUT/'C10_REPORT.md').write_text(report);(B/'C10B_REPORT.md').write_text(report.replace('# C10 — Second replication and manuscript package','# C10B — Locked second replication',1))
 legends='''# Figure legends

All figures are descriptive. Errors use common C0 training-normalized log1p count units. Patients are biological units; mask or model seeds quantify algorithmic variation. No confidence intervals or significance tests are shown. Distance bins are [0,150), [150,300), [300,600), [600,1200), [1200,2400), and [2400,infinity) µm. Only occupied bins contribute; patient support may change between bins. Bin Pearson requires at least 20 locations per run, while error summaries retain smaller bins. Exact counts are in Supplementary Table 7.

**Figure 1. Study overview.** (A) Measured ST and coordinates support frozen reconstruction of hidden expression, with an optional H&E baseline. (B) Distributed anchors (teal) and unmeasured locations (gray). (C) A contiguous forbidden region (purple), with anchors outside it. These are synthetic illustrative coordinates with 20% anchors and random seed 17, not experimental results. No biological n or uncertainty estimate applies.

**Figure 2. Development reconstruction.** One primary development patient, TENX99, and 50 genes. Panels A/C show V1 distributed sampling; B/D show the earlier C5 block protocol, labeled V2. A/B: median gene Pearson; C/D: mean gene MSE. Each point averages three frozen masks at 5%, 10%, 20%, 40%, 60%, or 80% coverage. Gray: ResNet18 H&E; teal: validation-selected interpolation; brown: fixed fusion. The C5 block protocol differs from the fixed 40% forbidden-region protocol in C9/C10. Lines connect tested settings; no error bars. Fusion's correlation increase need not imply lower error than ST.

**Figure 3. Geometry, distance and development reliability.** One development patient, TENX99, 50 genes. (A) ST-only MSE by nearest-anchor distance at 20% coverage, averaged over masks, for V1 and C5 block V2. (B) Secondary C6 open interpolation-only comparison: random A0, mean over ten seeds, versus deterministic spatial FPS A1 at 5/10/20/40%. This does not establish general sampling optimality. (C/D) Frozen C5 ResNet18 fusion at 20%, replicate 0, for the first original panel gene, FOXA1, under random/block masks. Hidden-location absolute errors are colored; observed anchors are teal. Both maps use the same maximum, the 98th percentile across all three prespecified map regimes, including capture-window source data. Larger errors saturate; no locations are deleted. Exact values, barcodes and prediction hashes are retained. Spatial axes have equal physical aspect ratio. No outcome-based gene selection or biological annotations were used.

**Figure 4. HER2ST transportability development.** Eight patients and the exact common 40-gene intersection. A/C: V1; B/D: V2; 20% random coverage. A/B show median gene Pearson; C/D show mean gene MSE in common C0 units on a log axis. Dots are patient summaries after averaging seeds/geometries; connected markers are medians across patients. Gray is T0 zero-shot; brown is T1 calibrated. H&E, ST and fixed fusion are separate methods. T0 exactly replays C7. T1 belongs to C8 transportability development, not fresh validation. Affine interpolation maps back to the same prediction, so calibration is not a distinct interpolation advantage.

**Figure 5. Fresh Wu validation.** Three held-out local Wu patients, 50 genes, with calibration confined to three other patients. (A) Patient correlations at 20% V1; horizontal bars are patient medians. (B/C) Median patient correlation/MSE at 10/20/40% V1; the MSE axis is logarithmic. (D) Paired patient ST−H&E and fusion−ST correlation effects at 20% V1; bars are median paired effects, which need not equal differences of cohort medians. (E) Paired ST-only V1/V2 patient correlations at 20%. (F) ST-only distance-bin MSE at 20%, averaging seeds within geometry, geometries within patient, then taking the median of available patients. Gray M0: H&E; teal M1: ST; brown M3: fixed 0.2 image fusion. Equivalent M2 is omitted from duplicate curves; no harmonization gain is claimed. Public images contain 15 or 25 native pixels per 112-µm patch.

**Figure 6. Second fresh replication.** Twenty-two Bassiouni patients, one first section each, 50 genes. All patients are evaluation patients; there is no C10 fitting. Panels match Figure 5: (A) 20% V1 correlations; (B/C) coverage correlation/error; (D) paired ST and fusion effects; (E) paired V1/V2; (F) distance-bin error. Ten random seeds are averaged within geometry, and three V2 geometries are averaged equally within patient. Bars and cohort curves show patient medians. M3 transports the frozen Wu scaler and image weight. Public images contain 24–29 native pixels per 112-µm patch. Distance curves are descriptive and can be nonmonotonic; variable support and location composition preclude a universal cutoff.

**Figure 7. Cross-cohort synthesis without pooling panels.** Columns show development primary n=1, HER2ST zero-shot n=8, HER2ST calibrated n=8, fresh Wu n=3, and fresh Bassiouni n=22. Both HER2ST columns reuse the same patients. Development/Wu/Bassiouni use 50 genes; HER2ST uses 40. (A–E) Patient effects at 20% V1: ST−H&E and fusion−ST; dots are patients and bars are medians. (F–J) Individual patient ST correlation curves at 10/20/40% V1. Cohort-specific axes and labels preserve panel and analytical-role differences; no pooled scalar is calculated.

## Supplementary figures

**Supplementary Figure 1. Morphology-only baselines and controls.** C1, one primary TENX99 patient, 50 genes. (A) Frozen B1 ridge and B2 MLP median Pearson; B0 constant predictions have undefined correlations and are absent. (B) NC1/NC2 MSE by model; definitions remain in the C1 report. Dots are saved training-seed runs, bars are means. These are not independent patients.

**Supplementary Figure 2. Graph context.** C2, one TENX99 patient, 50 genes. (A) Local C1 MLP, mean-pooling MLP, GCN and GraphSAGE correlation; (B) shuffled-graph MSE. Dots are training-seed runs and bars are means. No graph selection or C1 rerun was performed for C10. Contextual gain alone does not establish topology-specific information.

**Supplementary Figure 3. Topology controls.** C3, one TENX99 patient, 50 genes. (A) Median Pearson; (B) mean MSE. Local, self-only, mean-pooling and true SAGE controls accompany vertex permutation, degree rewiring and degree-plus-length rewiring. Training seeds are averaged within graph seed; randomized controls have five graph seeds. Dots are those averages; bars are means across displayed dots. No new models or graph controls were run.

**Supplementary Figure 4. Encoders.** C4, one primary TENX99 patient, 50 genes, local D1 MLP. ResNet18 and Phikon are historical frozen runs. Dots represent three training seeds, bars their means; (A) Pearson, (B) MSE. Phikon is not introduced into C10.

**Supplementary Figure 5. Development controls.** C5, one primary patient, 50 genes, 20% coverage, ResNet18 branch, recorded B3/B4 NC1/NC2/NC3 settings. (A) Pearson; (B) MSE. Dots are recorded control runs across geometries; bars are means. These are diagnostic distributions, not paired estimates of geometry effects. The frozen C5 audits separately document hidden-label invariance.

**Supplementary Figure 6. Sampling strategies.** C6, one primary patient, 50 genes, 20% coverage, ResNet18 interpolation-only baseline. (A) Open; (B) contiguous. A0 is ten-seed random sampling; A1 spatial FPS, A2 morphology FPS and A3 joint FPS are deterministic locked policies. Dots are policy runs; bars are means. This limited baseline diagnostic does not summarize the selected learned/fusion families or establish robust optimization.

**Supplementary Figure 7. External patient distributions.** C7, eight HER2ST patients, common 40 genes, 20% random coverage. (A) V1; (B) V2. E0: H&E; E1: ST; E2: fusion. Dots are patient means across seeds/geometries; bars are arithmetic patient means, unlike Figure 4's medians. These are retrospectively displayed zero-shot results. HER2ST is not fresh after C7 inspection.

**Supplementary Figure 8. Calibration and selected fusion.** C8, eight HER2ST patients, 40 genes, 20% coverage. (A) Selected fusion−ST Pearson; (B) relative MSE reduction, positive meaning improvement. Dots are patient/fold results grouped by V1/V2; bars are arithmetic means. Frozen fold selection is unchanged. This is transportability development, with no new weight selection.

**Supplementary Figure 9. Fresh-validation controls.** C9, one prespecified held-out patient, 1160920F, 50 genes, 20% coverage, random mask seed 0, V1 geometry 0 and V2 geometry 1, three control seeds. (A) Control−matched-clean Pearson; (B) clean−control MSE. NC1 permutes anchor vectors, NC2 anchor positions, NC3 image predictions. Dots are matched control settings; bars are means; M1 and M3 are shown. Negative MSE reduction indicates increased control error.

**Supplementary Figure 10. Replication controls and image QC.** C10 uses the same control structure on BASS_P01, 50 genes, 20% coverage. (A/B) Matched Pearson/MSE effects for M1/M3, three control seeds and two fixed regimes; dots are control settings and bars are means. (C) Native pixels per 112-µm patch for all 22 patients in source order. This is image-input QC, not expression performance. Coordinate overlays for all patients are in C10A/logs/COORDINATE_QA.png.

## Reproduction and source data

Run the figure-specific Python entrypoint beside an output, or `python C10/MANUSCRIPT/figures/plot_figures.py --figure all`. Each figure has a `_sources.csv` listing exact input CSVs and SHA-256 hashes. Source row numbers remain in the exported data and SOURCE_INDEX.json. PNGs are 300 dpi; SVG text remains editable. The original source CSVs include unrounded values; manuscript rounding does not alter plotted values.
'''
 (M/'FIGURE_LEGENDS.md').write_text(legends)
 print('READABLE REPORT AND LEGENDS COMPLETE')
if __name__=='__main__':main()
