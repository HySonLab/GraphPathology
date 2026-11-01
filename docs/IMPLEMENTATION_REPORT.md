# Implementation and handoff

A local revision is implemented on codex/publication-reproducibility. It is **partial against the full revision brief**:software migration, numerical package,available-export views and bounded tests are complete;original model replay,full raw-data/research orchestration,the revised source baseline and both document builds remain unavailable/unimplemented as detailed below.

- Archived all 74,253 starting tracked files under checkpoints/legacy using git mv,including reports/C0 and C1–C10. Preserved shared source/configuration/environment/registry and the obsolete publication package as evidence. Byte verification reports zero changed/missing files.
- Inventoried 31,824 missing referenced paths,including 1,372 learned-state and 10 fitted PCA paths;no original fitted states were installed. Feature caches,predictions,masks,preprocessing and fitted states are distinguished. Every missing original has an explicit manifest entry;primary head paths/checksums are in DATA_AND_MODELS.md.
- Added package-qualified graph/core/external/residual/gating/selection/sampling/vision modules;safe asset resolution,archive verification,doctor,synthetic smoke,bounded embedding inference and publication entrypoints. Copied numerical formulas remain separately versioned from archived hashed source. Bounded inference checks original identities/hashes and never falls back to new weights.
- Replaced C0-only README with manuscript scope/cohort roles,portable setup,tested quick start and truthful limitations. Added methods,target,migration,data/model,environment and test documentation;CPU CI configuration added.
- Tested clean installation and 84 bounded CPU checks. Generated 26 of 29 display export views:Figures 3–5,S1–S12;Tables 1–2,S1–S9. Figures 1, 2 and 6 are blocked. 66 publication checks passed,including printed Table2 values and main external curves/paired effects. These are older-export views;exact revised source-array/layout equivalence is unverified.
- No revised TeX/source_data/methods_data/renderers/templates were supplied or found in GitHub. Both document builds are blocked. No original-state prediction replay or full research job was executed. No claims of full reproduction are made.

## Commands

From this repository after README installation:

```bash
python -m graphpathology smoke --output outputs/quickstart
python -m pytest -q
python -m graphpathology verify-archive --historical --report outputs/new-audit.json
python -m graphpathology publication --available-only --output outputs/new-publication
python scripts/validate_export_views.py --output outputs/new-publication
python -m graphpathology doctor --require-models
python -m graphpathology infer --config configs/replay.example.json --dry-run
python -m graphpathology build-documents --dry-run
```

The last three commands currently report missing originals/source prerequisites. Do not interpret their expected failure as a scientific test pass. See TEST_REPORT.md for measured results and historical test disposition;PUBLICATION_STATUS.json lists each of 29 targets.

## GitHub handoff

No remote branch/PR exists yet:public cloning succeeded,but a push dry-run failed because no GitHub write credentials are configured. After authenticating Git for HySonLab/GraphPathology:

```bash
git push -u origin codex/publication-reproducibility
```

Open a PR from that branch to the repository's default branch. Use a title such as "Archive checkpoint evidence and add bounded manuscript-aligned reproduction tools". Describe the missing assets and partial publication status above;do not mark the full reproducibility task complete. Do not force-push or rewrite history. The local commit identifier is provided in the chat handoff.
