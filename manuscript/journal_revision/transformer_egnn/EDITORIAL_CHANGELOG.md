# Editorial change log

Scope: editorial integration only, starting from PR #3 head `d1a0d1bb4e9d802b92186c010b1f86f3efdff2d6`. No acquisition, fitting, patient prediction or numerical-export regeneration was performed. PR #3 remains draft and unmerged. The [complete pre-edit source/PDF snapshot](https://github.com/HySonLab/GraphPathology/tree/d1a0d1bb4e9d802b92186c010b1f86f3efdff2d6/manuscript/journal_revision/transformer_egnn) remains immutable in Git history; the earlier journal package and arXiv v1 are unchanged.

## Before / after outline

| Before | After |
| --- | --- |
| Abstract organized around successive fits and named advanced architectures | 172-word findings-centered abstract: framework, 25-patient reconstruction, controlled 22-patient placement effect, high-level GNN result and conditional cost interpretation; requested code URL retained |
| Results organized by original/follow-up/recovered-data execution stages | Six scientific questions: study/information access; distributed reconstruction; controlled placement; matched estimators; graph representations; failure regimes/calibration |
| Controlled 22-patient experiment near the end of Results | Central geometry result immediately follows distributed reconstruction, with identical-H/equal-count explanation, both estimators, paired effects and patient intervals |
| Two overlapping coordinate/graph-construction subsections | One graph-definition section with physical coordinates and a single set of mean/GCN/SAGE equations |
| Training repeated across multiple sections | One section distinguishes historical primary ensemble, 30 additional minibatch heads and 48 full-section fits; update schedules and target scaling remain distinct |
| Molecular/residual descriptions interleaved | One interpolation/calibration/residual section; development fitted ridge correction is explicitly distinct and fully specified in Supplementary Methods S6 |
| Detailed head searches, Phikon, conditional holdouts and HER2ST gates repeated or scattered | Unique details consolidated in Supplementary Methods S6, with main-text pointers |
| Discussion starts with model/history language and repeats Results | Findings-led synthesis: placement, molecular support, neural negative findings, transport limitations, cost bounds and prospective validation |
| Nine generic matrix captions mention HEST even on other cohorts | Nine route/cohort-specific captions with patient counts, information access, seed averaging and replication unit; common matrix interpretation stated once |
| Code identifiers appear in plot titles | Image only / Fixed fusion / Observed-residual interpolation; frozen models labeled as minibatch references; weak rewiring labeled as an attempt |
| Availability mixes failed transfers and recovered data | Current data recovery stated once; author-local exports/states and missing original heads accurately distinguished |
| No dedicated AI contribution paragraph | Distinct factual Methods disclosure attributes scientific ideas/direction/oversight to Truong-Son Hy and computational/editorial assistance to Codex/ChatGPT |

Approximate prose counts (same script, excluding figure environments, equation blocks and citation/input commands) decreased from 4,891 to 2,431 words for Introduction/Results/Discussion and from 6,009 to 3,975 for main Methods. These are reproducible editorial estimates, not the publisher's final word count. Unique technical details were moved or retained, not deleted to meet a page target. Main Methods subsections decreased from 18 to 9.

## Scientific and display preservation

- All 40 figures, 22 tables and 31 references remain. No display-item cap was imposed.
- Scientific CSV/gzip exports, protocols, fit records and prior results are hash-protected in `experiments/editorial/BASELINE.json`; every numerical table body is checked against the pre-edit commit.
- The 11 relabeled figure assets (nine matrices, one distance figure and the minibatch development comparison) were redrawn only from unchanged validated summary tables, with the same selection/aggregation, point values, line series and color mapping. No patient-level predictions were run.
- Full seed, distance/support, calibration, identity-control, geometry and runtime details remain in the main Methods, SI and source tables.
- Main figures and tables were reordered without removing evidence; stable LaTeX labels preserve cross-references. Controlled placement is now main Figure 6/Table 3; negative graph comparisons are Figures 9--11/Table 4; uncertainty is Table 5.

## Journal requirements and departures

Checked [Scientific Reports submission guidance](https://www.nature.com/srep/author-instructions/submission-guidelines) and [AI policy](https://www.nature.com/srep/journal-policies/ai) on 2026-11-01. The abstract is unstructured, at most 200 words, and has no bibliographic citation or figure reference. The requested code URL is retained: the guidance does not explicitly prohibit a code URL. AI assistance is disclosed in Methods; AI tools are not authors. Sole-human authorship and the Editorial Board competing-interest disclosure are preserved.

The rich display set is intentionally retained for author review, as instructed. Final display placement, submission-portal limits and wording of tool disclosure remain author submission decisions; no journal submission or arXiv change occurred. Author-local assets were not uploaded and no new release/account was created. See `AUTHOR_RELEASE_CHECKLIST.md` and `CLAIM_EVIDENCE_AUDIT.md`.
