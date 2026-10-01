# C9 fresh external validation

C0-C8 remain immutable. C9 uses the six-patient **local Visium cohort** of Wu et al. (2021), not the HER2ST cohort reanalyzed in that publication. Source: https://zenodo.org/records/4739739 .

Run from the repository root with `.venv/bin/python`:

1. `C9/acquire.py`, `C9/acquire_images.py`: checksum-verified public archives.
2. `C9/qualify.py`: barcode/gene/coordinate metadata only.
3. Inspect `pdfimages -list` and rendered PDF; `C9/prepare_images.py`: frozen image encoder and heads, no expression matrix access.
4. `C9/prepare_protocol.py`, then commit prospective source; `C9/lock.py`, then commit the lock.
5. `C9/run.py masks`, `C9/run.py fit`; commit the mask/calibration locks before evaluation.
6. `C9/run.py evaluate`, `C9/analyze.py`, `C9/verify.py`, figures/reports, parent verification and final freeze.

Stages with freeze files are deliberately non-overwriting. Do not rerun a frozen checkpoint in place. Raw assets, image feature caches and predictions are local ignored files, retained under C9 and SHA-256 manifest-covered; all source, configs, masks, metrics, calibration objects, logs and reports are committed.

The primary panel is the original 50 genes, unchanged. The cross-checkpoint comparison uses the exact common C7/C8 40-gene subset, explicitly labeled. Per-gene z-score harmonization commutes with the frozen linear ST interpolation: M1/M2 common-unit predictions are mathematically identical. This is an identifiability limitation, not a new positive harmonization result. No numerical roundoff or native-unit error reduction can satisfy that claim.
