# Matched visual-backbone follow-up

This exploratory follow-up separates two questions: whether spatial context improves a local image head within one frozen encoder, and whether the same improvement remains useful with molecular anchors and against anchor-only interpolation. Previously inspected patients are not a new blind validation. The original arXiv v1, prior numerical exports and all earlier experiments remain preserved.

The [committed protocol](../experiments/journal_revision/backbone_graph_comparison/PROTOCOL.md), [existing-evidence audit](../experiments/journal_revision/backbone_graph_comparison/EXISTING_BACKBONE_COMPARISON_AUDIT.md), [source/validation-guided refinement](../experiments/journal_revision/backbone_graph_comparison/PROTOCOL_V2_TRAINING_REFINEMENT.md) and [resource record](../experiments/journal_revision/backbone_graph_comparison/RESOURCE_ADJUSTMENT.md) describe the sequence and limits. The initial protocol was committed before new fitting. A documented two-rate refinement followed source/validation learning-curve diagnostics before any Phikon development or external outcome was evaluated; it was not the original protocol and was not selected on external outcomes. No further tuning is planned.

## Design and availability

ResNet18 and Phikon each have local, neighborhood-mean, GCN, mean-GraphSAGE, GCN-shuffled and GraphSAGE-shuffled heads, with seeds 0/1/2. The operators aggregate fixed image features; the encoder is frozen. Head capacities are approximately matched within each encoder. Two learning rates receive identical validation-only opportunities. Across 72 candidates, 18 ResNet18 fits are reused and 54 fits are new. Mean TENX95 MSE across all three seeds chooses one rate per configuration, yielding 36 selected heads. All candidate histories and states are retained.

Fitting uses NCBI783/NCBI785; validation uses TENX95. TENX99 is one development patient. Wu's three evaluation patients and the locked 22 Bassiouni patients are analyzed separately. Source-only feature scaling, the ordered 50-gene targets and the exact historical neural masks are shared within each pair. ResNet18's effective crop is approximately 98 μm, while Phikon uses the full registered approximately 112 μm patch. Cross-encoder differences include preprocessing, pretraining and feature dimension; their differences are not isolated graph effects.

Image-only prediction receives no same-section molecular labels. Anchor-assisted reconstruction applies the same observed-residual rule and interpolation weights to every head. The IDW molecular reference remains present. External image predictions receive the frozen C9 T1 transport map once. Primary coverage is 20% distributed with exact mask seeds 0/1/2, distinct from the original ten-mask molecular headline summaries. Supporting analyses retain 10/40% coverage, identical-hidden-target placement and historical central V2. Historical V1/V2 changes are not isolated placement interventions.

The manifests distinguish `COMPLETE` and `REUSED_VERIFIED`; interim files explicitly identify pending extraction or scoring. A final report is generated only after all 936 selected patient/model bundles are complete. The comparison report, selected and candidate manifests, learning histories, paired patient summaries, finite-correlation gene counts and common-valid-gene sensitivity live in [the study directory](../experiments/journal_revision/backbone_graph_comparison/). Patient effects average mask seeds, geometries where applicable and training seeds before pairing. Descriptive Bassiouni intervals resample patients 10,000 times. Wu effects remain individual; TENX99 receives no population interval.

## Reproduction

The ordinary installation and software tests do not download research data or models. Phikon extraction additionally requires the pinned `transformers==4.57.3` environment recorded in the study environment manifest. Downloaded encoders retain their own terms. The public Phikon checkpoint is pinned by revision and SHA-256; neither raw tissue inputs nor pretrained encoder weights are redistributed.

Research commands are explicit and can be costly. The acquisition/preparation command resumes completed patients and flushed partial feature arrays. Training commands reuse completed candidates, preserving their state hashes and histories. The refinement is bounded to the documented two rates; do not add trials based on evaluation results.

```bash
python -m graphpathology.backbone_data prepare
python -m graphpathology.backbone_train train
python -m graphpathology.backbone_refinement train
python -m graphpathology.backbone_refinement select
python -m graphpathology.backbone_selected_evaluate
python -m graphpathology.backbone_analysis
python -m graphpathology.backbone_report
python -m graphpathology.backbone_selected_verify --output outputs/backbone-integrity.json
```

Check each module's `--help` before execution. Public-checkout checks deliberately skip absent author-local states and per-gene arrays, identify those skips and verify their committed provenance instead. They do not claim a new neural replay:

```bash
python -m graphpathology.backbone_selected_verify --committed-only --output outputs/backbone-integrity.json
python -m pytest -q
python -m graphpathology.backbone_paper build --output outputs/matched-paper
python -m graphpathology.backbone_paper package --output outputs/matched-source.zip --check-directory outputs/matched-source-check
```

Build output and clean-extraction directories must be new. The maintained journal sources drive document builds; historical generation cannot restore the removed conceptual schematic or overwrite the expanded opening. Full author-local numerical and head-state bundles include checksums and restoration paths. A manifest is not a public release or DOI. The original three primary image-head binaries remain missing; newly trained heads do not recover them. Draft PR #3 remains unmerged for author review.
