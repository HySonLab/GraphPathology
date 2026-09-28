# C0 development checkpoint

Status: **PASS_DEVELOPMENT_ONLY**

Source commit: `2b54db3379500fec4603a61751145c929db9a72a`. Exact source hashes are in C0_SPEC_LOCK.json.

Prepared 35,536 paired locations across four HEST IDC patients, 280 genes, and 216 frozen observation masks.

## Verified

- Unique barcode joins and patch/expression coordinate agreement within 1.01 full-resolution pixels.
- Patient-disjoint train/validation/test assignments before target preprocessing.
- Gene identifiers intersected and reordered by name; held-out expression values never used for panel ranking.
- Training-only log-count normalization and gene ranking; raw counts retained separately.
- Every frozen mask passed hidden-label NaN corruption invariance.
- Raw assets and feature/target caches SHA-256 checked; feature barcode order verified.

## Limitations

- Four patients only; no publication-level statistical claim
- 280 shared biological targets after control exclusion and identifier harmonization of 541 assay features; 1000/3000-gene panels unavailable
- Generic ImageNet ResNet18 reference, not pathology foundation embeddings
- No validated cell/nucleus labels; RGB and texture summaries only
- Simulated square/circular masks; tissue-area coverage unavailable
- Thumbnail alignment audit cannot establish cell-level registration accuracy

## Next checkpoint

C1 is not executed. Review C0 cohort/split suitability and encoder choice before model development, as specified by the supplied execution plan.
The one-gene study is the top training log-variance gene. The 50/250 panels are available; the 1000-gene request is explicitly unsupported by this assay.
No expression prediction, performance, biological discovery, or publication claims are made.

## Reproduction

See README.md. All masks are in C0_MASKS/ and their geometry/checksums in C0_MASK_INDEX.json.
C0_SPEC_LOCK.json refuses silent replacement. Later revisions must use a separately versioned checkpoint.
