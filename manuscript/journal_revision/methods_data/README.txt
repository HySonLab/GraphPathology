VERIFIED MODEL AND CALIBRATION DETAILS

These files are publication-oriented transcriptions of the public implementation,
actual fitted-run records and saved preprocessing objects at source revision
03506defbe205cabd13865780d6a1224b348a781. No models were fitted or predictions
recomputed to create these files. See METHODS_SOURCE_AUDIT.txt for source links.

image_feature_scaler.csv: all512 feature means and scales, in encoder output order.
target_and_calibration_constants.csv: original log-count normalization, additional
head target scaler and source/Wu calibration moments for the50 primary genes.
primary_model_details.json: exact layer, optimization, selection and ensemble rules,
with the recorded epochs and SHA256 identifiers of the three fitted heads.
software_versions.json: recorded research environments and editorial drawing version.
implementation_sources.csv: exact public source links, locations and SHA256 values.

Numerical constants are written using full Python float precision without
presentation rounding. These are settings/constants, not additional experiment
results. The fitted neural parameter binaries are not included.
