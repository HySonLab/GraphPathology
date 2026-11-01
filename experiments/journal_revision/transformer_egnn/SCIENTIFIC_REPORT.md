# Post-v1 GraphTransformer/EGNN scientific report

All16 configurations and seeds0/1/2 were fixed in local protocol commit2e03fae636d6f6527d25db58f4cbab72f7c1a286 before fitting. This is exploratory work in previously inspected cohorts, not formal preregistration or fresh blind validation. No external label selected an architecture, epoch, seed or calibration coefficient.

Exact GEO byte ranges recovered all88 Bassiouni input files with original SHA256 matches (291,569,893 bytes transferred; no37.7GB whole-series transfer). The matched A/B1/B2 study completed19,800 runs and39,600 exact all-unobserved-label corruption checks. At20% B2 coverage, all22 patients deteriorated under region exclusion for both fixed estimators: IDW median paired Δr=-0.154080 and relative MSE=+44.9864%; uniform Δr=-0.150195 and relative MSE=+39.9128%. H and anchor counts are identical between placements. See [full baseline source tables](bassiouni_anchors/PATIENT_UNCERTAINTY.csv).

The neural matrix has48 new fits plus30 frozen journal-v2 fits evaluated in26 patients (development1, Wu3, Bassiouni22), using image, fixed0.2/0.8 and exploratory residual routes. All routes use the same targets/masks; only the image-only versus molecular comparison is unequal-information. Mask means precede geometry/patient means, then training seeds are averaged within patient. Relative MSE is normalized after these averages. Wu patients are shown individually below; only Bassiouni uses10,000 patient bootstrap resamples (seed20261101). No cohort is pooled.

The primary neural results are negative. At20% Bassiouni V1 coverage, GPS full residual reconstruction changed median paired correlation by−0.050421 and MSE by+13.4264%; EGNN latent changed correlation by−0.047994 and MSE by+13.0083%. Neither primary model improved Pearson/MSE/MAE in any of the22 patients. Nonlocal GPS versus its local ablation and latent EGNN versus fixed distance also worsened the primary point estimates. This does not support a stronger neural novelty claim under the tested conditions.

## Validation and resources

| Configuration | Parameters | Selected val MSE mean ± SD | Selected epoch range | Fit wall s median | Peak process MiB max |
| --- | --- | --- | --- | --- | --- |
| GPS full | 121018 | 0.76267 ± 0.06881 | 5–60 | 5.67 | 629.1 |
| GPS local | 121018 | 0.77027 ± 0.04794 | 5–56 | 7.04 | 691.5 |
| GPS no distances | 121018 | 0.75976 ± 0.05823 | 5–42 | 5.45 | 600.4 |
| GPS self | 121018 | 0.74145 ± 0.01672 | 47–53 | 15.76 | 546.4 |
| GPS shuffled | 121018 | 0.73777 ± 0.07506 | 6–76 | 22.15 | 634.2 |
| GPS rewired | 121018 | 0.76267 ± 0.06881 | 5–60 | 3.78 | 633.0 |
| Distance fixed | 87154 | 0.76781 ± 0.03505 | 7–53 | 12.47 | 638.2 |
| EGNN latent | 87282 | 0.76041 ± 0.03299 | 7–57 | 13.95 | 647.5 |
| Distance disabled | 87154 | 0.76690 ± 0.03599 | 7–74 | 12.10 | 650.2 |
| EGNN self | 87282 | 0.80009 ± 0.04158 | 8–35 | 1.91 | 506.1 |
| EGNN shuffled | 87282 | 0.77522 ± 0.02525 | 7–40 | 14.76 | 643.2 |
| EGNN rewired | 87282 | 0.76041 ± 0.03299 | 7–57 | 13.95 | 651.1 |
| Local MLP | 72114 | 0.70188 ± 0.05036 | 8–77 | 5.30 | 497.4 |
| Mean | 72114 | 0.77901 ± 0.03597 | 8–9 | 1.24 | 522.4 |
| GCN | 72114 | 0.78088 ± 0.03561 | 8–9 | 1.26 | 524.3 |
| GraphSAGE | 72075 | 0.67721 ± 0.00617 | 39–60 | 3.77 | 569.7 |

Validation uses TENX95 only. All new models use two full-section source updates per epoch; older30 reference fits retain their different original minibatch regime. Parameter counts are nominal, not equal effective capacity. Training wall time includes concurrent workload and is not a controlled speed ranking. [Standardized CPU inference](diagnostics/CPU_BENCHMARK.csv) reports all seeds, four threads, two warmups and five timed forwards on the same7,760-node validation section, excluding frozen CNN/graph setup. No GPU or paid service was used.

## All architecture test metrics at20% historical V1 coverage

Each reported value is a patient summary after all three mask/training seeds; Bassiouni rows are medians across22 such summaries. Wu remains individual. Full10/20/40%/central-V2/B2/per-gene exports are retained, including undefined-correlation counts and secondary near-tie-sensitive Spearman. No model is selected from these tables.

### Development / TENX99 / image

| Configuration | Pearson | MSE | MAE | Mean defined genes range |
| --- | --- | --- | --- | --- |
| GPS full | 0.60413 | 1.121255 | 0.845918 | 50.0–50.0 |
| GPS local | 0.59643 | 1.086989 | 0.831929 | 50.0–50.0 |
| GPS no distances | 0.61462 | 1.148482 | 0.853632 | 50.0–50.0 |
| GPS self | 0.66827 | 1.045583 | 0.813314 | 50.0–50.0 |
| GPS shuffled | 0.65973 | 1.142689 | 0.856138 | 50.0–50.0 |
| GPS rewired | 0.60413 | 1.121255 | 0.845918 | 50.0–50.0 |
| Distance fixed | 0.59624 | 1.016555 | 0.805821 | 50.0–50.0 |
| EGNN latent | 0.62487 | 1.017355 | 0.804067 | 50.0–50.0 |
| Distance disabled | 0.60605 | 1.042494 | 0.814132 | 50.0–50.0 |
| EGNN self | 0.56196 | 1.049511 | 0.824017 | 50.0–50.0 |
| EGNN shuffled | 0.60800 | 1.048835 | 0.821126 | 50.0–50.0 |
| EGNN rewired | 0.62487 | 1.017355 | 0.804067 | 50.0–50.0 |
| Local MLP | 0.64674 | 0.938263 | 0.782851 | 50.0–50.0 |
| Mean | 0.61480 | 1.046753 | 0.840870 | 50.0–50.0 |
| GCN | 0.61172 | 1.047919 | 0.841626 | 50.0–50.0 |
| GraphSAGE | 0.67496 | 0.892144 | 0.760090 | 50.0–50.0 |
| idw_k16_p2 | 0.79225 | 0.268379 | 0.390382 | 50.0–50.0 |
| uniform_k16 | 0.72706 | 0.332362 | 0.440470 | 50.0–50.0 |

### Development / TENX99 / fixed_20pct_image

| Configuration | Pearson | MSE | MAE | Mean defined genes range |
| --- | --- | --- | --- | --- |
| GPS full | 0.81454 | 0.281903 | 0.416981 | 50.0–50.0 |
| GPS local | 0.81347 | 0.280745 | 0.415815 | 50.0–50.0 |
| GPS no distances | 0.81424 | 0.284259 | 0.418610 | 50.0–50.0 |
| GPS self | 0.82145 | 0.267462 | 0.403949 | 50.0–50.0 |
| GPS shuffled | 0.81955 | 0.277708 | 0.412576 | 50.0–50.0 |
| GPS rewired | 0.81454 | 0.281903 | 0.416981 | 50.0–50.0 |
| Distance fixed | 0.81585 | 0.272321 | 0.408867 | 50.0–50.0 |
| EGNN latent | 0.81812 | 0.271096 | 0.406997 | 50.0–50.0 |
| Distance disabled | 0.81596 | 0.273464 | 0.409066 | 50.0–50.0 |
| EGNN self | 0.81525 | 0.276161 | 0.413063 | 50.0–50.0 |
| EGNN shuffled | 0.81896 | 0.272431 | 0.408983 | 50.0–50.0 |
| EGNN rewired | 0.81812 | 0.271096 | 0.406997 | 50.0–50.0 |
| Local MLP | 0.81778 | 0.272980 | 0.411566 | 50.0–50.0 |
| Mean | 0.80581 | 0.290565 | 0.426500 | 50.0–50.0 |
| GCN | 0.80575 | 0.290717 | 0.426639 | 50.0–50.0 |
| GraphSAGE | 0.81467 | 0.272387 | 0.410136 | 50.0–50.0 |
| idw_k16_p2 | 0.79225 | 0.268379 | 0.390382 | 50.0–50.0 |
| uniform_k16 | 0.72706 | 0.332362 | 0.440470 | 50.0–50.0 |

### Development / TENX99 / exploratory_residual_anchor

| Configuration | Pearson | MSE | MAE | Mean defined genes range |
| --- | --- | --- | --- | --- |
| GPS full | 0.81586 | 0.255848 | 0.385681 | 50.0–50.0 |
| GPS local | 0.81370 | 0.265751 | 0.392153 | 50.0–50.0 |
| GPS no distances | 0.82087 | 0.246104 | 0.376436 | 50.0–50.0 |
| GPS self | 0.80211 | 0.326182 | 0.438063 | 50.0–50.0 |
| GPS shuffled | 0.82107 | 0.268477 | 0.396821 | 50.0–50.0 |
| GPS rewired | 0.81586 | 0.255848 | 0.385681 | 50.0–50.0 |
| Distance fixed | 0.80198 | 0.298254 | 0.420282 | 50.0–50.0 |
| EGNN latent | 0.80538 | 0.290332 | 0.413771 | 50.0–50.0 |
| Distance disabled | 0.80600 | 0.286757 | 0.410835 | 50.0–50.0 |
| EGNN self | 0.78702 | 0.312026 | 0.431884 | 50.0–50.0 |
| EGNN shuffled | 0.79004 | 0.318677 | 0.436586 | 50.0–50.0 |
| EGNN rewired | 0.80538 | 0.290332 | 0.413771 | 50.0–50.0 |
| Local MLP | 0.82004 | 0.253188 | 0.386051 | 50.0–50.0 |
| Mean | 0.82338 | 0.237146 | 0.367640 | 50.0–50.0 |
| GCN | 0.82123 | 0.239227 | 0.369745 | 50.0–50.0 |
| GraphSAGE | 0.82576 | 0.232583 | 0.367041 | 50.0–50.0 |
| idw_k16_p2 | 0.79225 | 0.268379 | 0.390382 | 50.0–50.0 |
| uniform_k16 | 0.72706 | 0.332362 | 0.440470 | 50.0–50.0 |

### Wu / 1160920F / image

| Configuration | Pearson | MSE | MAE | Mean defined genes range |
| --- | --- | --- | --- | --- |
| GPS full | -0.01118 | 0.303322 | 0.388559 | 50.0–50.0 |
| GPS local | -0.00935 | 0.311209 | 0.391638 | 50.0–50.0 |
| GPS no distances | -0.01641 | 0.309291 | 0.391411 | 50.0–50.0 |
| GPS self | -0.02039 | 0.295444 | 0.370622 | 50.0–50.0 |
| GPS shuffled | -0.01439 | 0.310505 | 0.389658 | 50.0–50.0 |
| GPS rewired | -0.01118 | 0.303322 | 0.388559 | 50.0–50.0 |
| Distance fixed | -0.02556 | 0.301003 | 0.373698 | 50.0–50.0 |
| EGNN latent | -0.01714 | 0.320603 | 0.389553 | 50.0–50.0 |
| Distance disabled | -0.01892 | 0.300440 | 0.371380 | 50.0–50.0 |
| EGNN self | -0.01618 | 0.266237 | 0.362373 | 50.0–50.0 |
| EGNN shuffled | -0.00529 | 0.311462 | 0.382754 | 50.0–50.0 |
| EGNN rewired | -0.01714 | 0.320603 | 0.389553 | 50.0–50.0 |
| Local MLP | -0.04819 | 0.424797 | 0.452876 | 50.0–50.0 |
| Mean | -0.04253 | 0.316934 | 0.399982 | 50.0–50.0 |
| GCN | -0.04136 | 0.317968 | 0.400521 | 50.0–50.0 |
| GraphSAGE | -0.06854 | 0.374522 | 0.425233 | 50.0–50.0 |
| idw_k16_p2 | 0.35716 | 0.081742 | 0.195417 | 49.7–49.7 |
| uniform_k16 | 0.34457 | 0.081339 | 0.198511 | 49.7–49.7 |

### Wu / 1160920F / fixed_20pct_image

| Configuration | Pearson | MSE | MAE | Mean defined genes range |
| --- | --- | --- | --- | --- |
| GPS full | 0.35182 | 0.088080 | 0.213966 | 50.0–50.0 |
| GPS local | 0.35120 | 0.088388 | 0.213328 | 50.0–50.0 |
| GPS no distances | 0.35224 | 0.088310 | 0.214458 | 50.0–50.0 |
| GPS self | 0.34254 | 0.087762 | 0.207554 | 50.0–50.0 |
| GPS shuffled | 0.35054 | 0.088381 | 0.211384 | 50.0–50.0 |
| GPS rewired | 0.35182 | 0.088080 | 0.213966 | 50.0–50.0 |
| Distance fixed | 0.34921 | 0.087946 | 0.210906 | 50.0–50.0 |
| EGNN latent | 0.34972 | 0.088763 | 0.211509 | 50.0–50.0 |
| Distance disabled | 0.34971 | 0.087913 | 0.210216 | 50.0–50.0 |
| EGNN self | 0.34839 | 0.086654 | 0.213013 | 50.0–50.0 |
| EGNN shuffled | 0.35109 | 0.088394 | 0.211392 | 50.0–50.0 |
| EGNN rewired | 0.34972 | 0.088763 | 0.211509 | 50.0–50.0 |
| Local MLP | 0.33832 | 0.093067 | 0.212888 | 50.0–50.0 |
| Mean | 0.34958 | 0.088709 | 0.214627 | 50.0–50.0 |
| GCN | 0.34950 | 0.088751 | 0.214633 | 50.0–50.0 |
| GraphSAGE | 0.33825 | 0.091026 | 0.210801 | 50.0–50.0 |
| idw_k16_p2 | 0.35716 | 0.081742 | 0.195417 | 49.7–49.7 |
| uniform_k16 | 0.34457 | 0.081339 | 0.198511 | 49.7–49.7 |

### Wu / 1160920F / exploratory_residual_anchor

| Configuration | Pearson | MSE | MAE | Mean defined genes range |
| --- | --- | --- | --- | --- |
| GPS full | 0.33716 | 0.083410 | 0.200977 | 50.0–50.0 |
| GPS local | 0.33721 | 0.083153 | 0.200444 | 50.0–50.0 |
| GPS no distances | 0.34101 | 0.082933 | 0.199901 | 50.0–50.0 |
| GPS self | 0.29508 | 0.089939 | 0.213263 | 50.0–50.0 |
| GPS shuffled | 0.32584 | 0.085039 | 0.204315 | 50.0–50.0 |
| GPS rewired | 0.33716 | 0.083410 | 0.200977 | 50.0–50.0 |
| Distance fixed | 0.32140 | 0.084143 | 0.202754 | 50.0–50.0 |
| EGNN latent | 0.32353 | 0.084182 | 0.202874 | 50.0–50.0 |
| Distance disabled | 0.32326 | 0.083957 | 0.202331 | 50.0–50.0 |
| EGNN self | 0.31417 | 0.084969 | 0.204678 | 50.0–50.0 |
| EGNN shuffled | 0.32600 | 0.084137 | 0.202880 | 50.0–50.0 |
| EGNN rewired | 0.32353 | 0.084182 | 0.202874 | 50.0–50.0 |
| Local MLP | 0.28719 | 0.093378 | 0.218361 | 50.0–50.0 |
| Mean | 0.33710 | 0.082447 | 0.198454 | 50.0–50.0 |
| GCN | 0.33638 | 0.082506 | 0.198611 | 50.0–50.0 |
| GraphSAGE | 0.31006 | 0.086413 | 0.206907 | 50.0–50.0 |
| idw_k16_p2 | 0.35716 | 0.081742 | 0.195417 | 49.7–49.7 |
| uniform_k16 | 0.34457 | 0.081339 | 0.198511 | 49.7–49.7 |

### Wu / CID4290 / image

| Configuration | Pearson | MSE | MAE | Mean defined genes range |
| --- | --- | --- | --- | --- |
| GPS full | 0.13290 | 0.262233 | 0.384991 | 50.0–50.0 |
| GPS local | 0.13294 | 0.280385 | 0.390614 | 50.0–50.0 |
| GPS no distances | 0.12143 | 0.274279 | 0.388259 | 50.0–50.0 |
| GPS self | 0.15908 | 0.283411 | 0.404700 | 50.0–50.0 |
| GPS shuffled | 0.10183 | 0.295704 | 0.414170 | 50.0–50.0 |
| GPS rewired | 0.13290 | 0.262233 | 0.384991 | 50.0–50.0 |
| Distance fixed | 0.14812 | 0.276719 | 0.388303 | 50.0–50.0 |
| EGNN latent | 0.13943 | 0.265691 | 0.384306 | 50.0–50.0 |
| Distance disabled | 0.15079 | 0.272500 | 0.382509 | 50.0–50.0 |
| EGNN self | 0.12799 | 0.233543 | 0.366013 | 50.0–50.0 |
| EGNN shuffled | 0.11737 | 0.261608 | 0.380685 | 50.0–50.0 |
| EGNN rewired | 0.13943 | 0.265691 | 0.384306 | 50.0–50.0 |
| Local MLP | 0.15034 | 0.230900 | 0.358067 | 50.0–50.0 |
| Mean | 0.15685 | 0.222737 | 0.348965 | 50.0–50.0 |
| GCN | 0.15669 | 0.222821 | 0.348921 | 50.0–50.0 |
| GraphSAGE | 0.20504 | 0.229585 | 0.368000 | 50.0–50.0 |
| idw_k16_p2 | 0.31464 | 0.098518 | 0.217633 | 50.0–50.0 |
| uniform_k16 | 0.26140 | 0.105694 | 0.229913 | 50.0–50.0 |

### Wu / CID4290 / fixed_20pct_image

| Configuration | Pearson | MSE | MAE | Mean defined genes range |
| --- | --- | --- | --- | --- |
| GPS full | 0.32458 | 0.101222 | 0.234055 | 50.0–50.0 |
| GPS local | 0.32303 | 0.102082 | 0.232726 | 50.0–50.0 |
| GPS no distances | 0.32237 | 0.101993 | 0.233413 | 50.0–50.0 |
| GPS self | 0.32926 | 0.101162 | 0.236335 | 50.0–50.0 |
| GPS shuffled | 0.32058 | 0.103204 | 0.238841 | 50.0–50.0 |
| GPS rewired | 0.32458 | 0.101222 | 0.234055 | 50.0–50.0 |
| Distance fixed | 0.32604 | 0.101730 | 0.233904 | 50.0–50.0 |
| EGNN latent | 0.32479 | 0.101079 | 0.232527 | 50.0–50.0 |
| Distance disabled | 0.32628 | 0.101377 | 0.231864 | 50.0–50.0 |
| EGNN self | 0.32353 | 0.100261 | 0.233438 | 50.0–50.0 |
| EGNN shuffled | 0.32195 | 0.101281 | 0.232578 | 50.0–50.0 |
| EGNN rewired | 0.32479 | 0.101079 | 0.232527 | 50.0–50.0 |
| Local MLP | 0.32594 | 0.099763 | 0.230682 | 50.0–50.0 |
| Mean | 0.32120 | 0.100354 | 0.230739 | 50.0–50.0 |
| GCN | 0.32134 | 0.100345 | 0.230695 | 50.0–50.0 |
| GraphSAGE | 0.32899 | 0.099782 | 0.233273 | 50.0–50.0 |
| idw_k16_p2 | 0.31464 | 0.098518 | 0.217633 | 50.0–50.0 |
| uniform_k16 | 0.26140 | 0.105694 | 0.229913 | 50.0–50.0 |

### Wu / CID4290 / exploratory_residual_anchor

| Configuration | Pearson | MSE | MAE | Mean defined genes range |
| --- | --- | --- | --- | --- |
| GPS full | 0.30251 | 0.115389 | 0.246671 | 50.0–50.0 |
| GPS local | 0.30624 | 0.108904 | 0.238940 | 50.0–50.0 |
| GPS no distances | 0.30755 | 0.107679 | 0.237373 | 50.0–50.0 |
| GPS self | 0.28224 | 0.132921 | 0.270150 | 50.0–50.0 |
| GPS shuffled | 0.30179 | 0.106249 | 0.236628 | 50.0–50.0 |
| GPS rewired | 0.30251 | 0.115389 | 0.246671 | 50.0–50.0 |
| Distance fixed | 0.31294 | 0.106065 | 0.237255 | 50.0–50.0 |
| EGNN latent | 0.30507 | 0.109820 | 0.241708 | 50.0–50.0 |
| Distance disabled | 0.31442 | 0.106048 | 0.237260 | 50.0–50.0 |
| EGNN self | 0.30327 | 0.107664 | 0.238494 | 50.0–50.0 |
| EGNN shuffled | 0.29640 | 0.111763 | 0.243358 | 50.0–50.0 |
| EGNN rewired | 0.30507 | 0.109820 | 0.241708 | 50.0–50.0 |
| Local MLP | 0.29775 | 0.112944 | 0.246245 | 50.0–50.0 |
| Mean | 0.32177 | 0.098133 | 0.220704 | 50.0–50.0 |
| GCN | 0.32171 | 0.098133 | 0.220838 | 50.0–50.0 |
| GraphSAGE | 0.32569 | 0.101229 | 0.230317 | 50.0–50.0 |
| idw_k16_p2 | 0.31464 | 0.098518 | 0.217633 | 50.0–50.0 |
| uniform_k16 | 0.26140 | 0.105694 | 0.229913 | 50.0–50.0 |

### Wu / CID44971 / image

| Configuration | Pearson | MSE | MAE | Mean defined genes range |
| --- | --- | --- | --- | --- |
| GPS full | 0.04364 | 0.269663 | 0.389781 | 50.0–50.0 |
| GPS local | 0.04446 | 0.283433 | 0.393738 | 50.0–50.0 |
| GPS no distances | 0.04382 | 0.268818 | 0.384897 | 50.0–50.0 |
| GPS self | 0.07540 | 0.310459 | 0.423154 | 50.0–50.0 |
| GPS shuffled | 0.04742 | 0.296236 | 0.411315 | 50.0–50.0 |
| GPS rewired | 0.04364 | 0.269663 | 0.389781 | 50.0–50.0 |
| Distance fixed | 0.06763 | 0.276050 | 0.382730 | 50.0–50.0 |
| EGNN latent | 0.05793 | 0.268120 | 0.383714 | 50.0–50.0 |
| Distance disabled | 0.07049 | 0.271587 | 0.377141 | 50.0–50.0 |
| EGNN self | 0.04916 | 0.261641 | 0.385136 | 50.0–50.0 |
| EGNN shuffled | 0.05217 | 0.275320 | 0.387284 | 50.0–50.0 |
| EGNN rewired | 0.05793 | 0.268120 | 0.383714 | 50.0–50.0 |
| Local MLP | 0.04964 | 0.270975 | 0.384140 | 50.0–50.0 |
| Mean | 0.05740 | 0.249649 | 0.370853 | 50.0–50.0 |
| GCN | 0.05716 | 0.249865 | 0.370863 | 50.0–50.0 |
| GraphSAGE | 0.08688 | 0.259498 | 0.383293 | 50.0–50.0 |
| idw_k16_p2 | 0.50332 | 0.113984 | 0.243434 | 50.0–50.0 |
| uniform_k16 | 0.43519 | 0.123362 | 0.257315 | 50.0–50.0 |

### Wu / CID44971 / fixed_20pct_image

| Configuration | Pearson | MSE | MAE | Mean defined genes range |
| --- | --- | --- | --- | --- |
| GPS full | 0.49018 | 0.117447 | 0.256012 | 50.0–50.0 |
| GPS local | 0.49148 | 0.117908 | 0.254743 | 50.0–50.0 |
| GPS no distances | 0.49432 | 0.117259 | 0.254639 | 50.0–50.0 |
| GPS self | 0.48213 | 0.119247 | 0.258335 | 50.0–50.0 |
| GPS shuffled | 0.49840 | 0.118578 | 0.258152 | 50.0–50.0 |
| GPS rewired | 0.49018 | 0.117447 | 0.256012 | 50.0–50.0 |
| Distance fixed | 0.49636 | 0.117520 | 0.253405 | 50.0–50.0 |
| EGNN latent | 0.49330 | 0.117143 | 0.253135 | 50.0–50.0 |
| Distance disabled | 0.49541 | 0.117243 | 0.251978 | 50.0–50.0 |
| EGNN self | 0.49502 | 0.117261 | 0.256372 | 50.0–50.0 |
| EGNN shuffled | 0.49663 | 0.117531 | 0.253712 | 50.0–50.0 |
| EGNN rewired | 0.49330 | 0.117143 | 0.253135 | 50.0–50.0 |
| Local MLP | 0.49281 | 0.117485 | 0.254861 | 50.0–50.0 |
| Mean | 0.50049 | 0.116892 | 0.254783 | 50.0–50.0 |
| GCN | 0.50040 | 0.116897 | 0.254751 | 50.0–50.0 |
| GraphSAGE | 0.49818 | 0.116989 | 0.256206 | 50.0–50.0 |
| idw_k16_p2 | 0.50332 | 0.113984 | 0.243434 | 50.0–50.0 |
| uniform_k16 | 0.43519 | 0.123362 | 0.257315 | 50.0–50.0 |

### Wu / CID44971 / exploratory_residual_anchor

| Configuration | Pearson | MSE | MAE | Mean defined genes range |
| --- | --- | --- | --- | --- |
| GPS full | 0.42301 | 0.138394 | 0.270458 | 50.0–50.0 |
| GPS local | 0.43638 | 0.130247 | 0.262871 | 50.0–50.0 |
| GPS no distances | 0.44276 | 0.128571 | 0.260952 | 50.0–50.0 |
| GPS self | 0.37238 | 0.159599 | 0.294264 | 50.0–50.0 |
| GPS shuffled | 0.45203 | 0.124919 | 0.258287 | 50.0–50.0 |
| GPS rewired | 0.42301 | 0.138394 | 0.270458 | 50.0–50.0 |
| Distance fixed | 0.44387 | 0.126620 | 0.260374 | 50.0–50.0 |
| EGNN latent | 0.43045 | 0.130075 | 0.263801 | 50.0–50.0 |
| Distance disabled | 0.44385 | 0.126612 | 0.260129 | 50.0–50.0 |
| EGNN self | 0.43218 | 0.129282 | 0.263629 | 50.0–50.0 |
| EGNN shuffled | 0.42999 | 0.130229 | 0.264717 | 50.0–50.0 |
| EGNN rewired | 0.43045 | 0.130075 | 0.263801 | 50.0–50.0 |
| Local MLP | 0.41691 | 0.136554 | 0.270941 | 50.0–50.0 |
| Mean | 0.48135 | 0.117064 | 0.247879 | 50.0–50.0 |
| GCN | 0.48115 | 0.117098 | 0.248078 | 50.0–50.0 |
| GraphSAGE | 0.45165 | 0.123443 | 0.255977 | 50.0–50.0 |
| idw_k16_p2 | 0.50332 | 0.113984 | 0.243434 | 50.0–50.0 |
| uniform_k16 | 0.43519 | 0.123362 | 0.257315 | 50.0–50.0 |

### Bassiouni / 22 patient medians / image

| Configuration | Pearson | MSE | MAE | Mean defined genes range |
| --- | --- | --- | --- | --- |
| GPS full | 0.04649 | 0.353163 | 0.421550 | 46.7–50.0 |
| GPS local | 0.04202 | 0.376169 | 0.434357 | 46.7–50.0 |
| GPS no distances | 0.04606 | 0.360647 | 0.425807 | 46.7–50.0 |
| GPS self | 0.06406 | 0.389232 | 0.460542 | 46.7–50.0 |
| GPS shuffled | 0.05420 | 0.395130 | 0.451556 | 46.7–50.0 |
| GPS rewired | 0.04649 | 0.353163 | 0.421550 | 46.7–50.0 |
| Distance fixed | 0.04657 | 0.368355 | 0.427388 | 46.7–50.0 |
| EGNN latent | 0.05246 | 0.356763 | 0.426782 | 46.7–50.0 |
| Distance disabled | 0.04547 | 0.359811 | 0.419420 | 46.7–50.0 |
| EGNN self | 0.04419 | 0.349633 | 0.422734 | 46.7–50.0 |
| EGNN shuffled | 0.03844 | 0.359976 | 0.424057 | 46.7–50.0 |
| EGNN rewired | 0.05246 | 0.356763 | 0.426782 | 46.7–50.0 |
| Local MLP | 0.05754 | 0.342391 | 0.401985 | 46.7–50.0 |
| Mean | 0.04337 | 0.336414 | 0.404488 | 46.7–50.0 |
| GCN | 0.04204 | 0.337036 | 0.404531 | 46.7–50.0 |
| GraphSAGE | 0.06260 | 0.322166 | 0.398320 | 46.7–50.0 |
| idw_k16_p2 | 0.41786 | 0.092172 | 0.202715 | 44.7–50.0 |
| uniform_k16 | 0.39844 | 0.094582 | 0.209673 | 44.7–50.0 |

### Bassiouni / 22 patient medians / fixed_20pct_image

| Configuration | Pearson | MSE | MAE | Mean defined genes range |
| --- | --- | --- | --- | --- |
| GPS full | 0.40504 | 0.099674 | 0.227799 | 46.7–50.0 |
| GPS local | 0.40532 | 0.100098 | 0.228110 | 46.7–50.0 |
| GPS no distances | 0.40550 | 0.100158 | 0.227796 | 46.7–50.0 |
| GPS self | 0.40296 | 0.100710 | 0.230021 | 46.7–50.0 |
| GPS shuffled | 0.40734 | 0.099704 | 0.230743 | 46.7–50.0 |
| GPS rewired | 0.40503 | 0.099674 | 0.227799 | 46.7–50.0 |
| Distance fixed | 0.40483 | 0.099840 | 0.227357 | 46.7–50.0 |
| EGNN latent | 0.40318 | 0.099953 | 0.227500 | 46.7–50.0 |
| Distance disabled | 0.40441 | 0.099916 | 0.227104 | 46.7–50.0 |
| EGNN self | 0.40511 | 0.098067 | 0.227033 | 46.7–50.0 |
| EGNN shuffled | 0.40449 | 0.100108 | 0.226586 | 46.7–50.0 |
| EGNN rewired | 0.40318 | 0.099953 | 0.227500 | 46.7–50.0 |
| Local MLP | 0.40955 | 0.101023 | 0.227017 | 46.7–50.0 |
| Mean | 0.40913 | 0.099750 | 0.225692 | 46.7–50.0 |
| GCN | 0.40861 | 0.099827 | 0.225723 | 46.7–50.0 |
| GraphSAGE | 0.41176 | 0.100721 | 0.228851 | 46.7–50.0 |
| idw_k16_p2 | 0.41786 | 0.092172 | 0.202715 | 44.7–50.0 |
| uniform_k16 | 0.39844 | 0.094582 | 0.209673 | 44.7–50.0 |

### Bassiouni / 22 patient medians / exploratory_residual_anchor

| Configuration | Pearson | MSE | MAE | Mean defined genes range |
| --- | --- | --- | --- | --- |
| GPS full | 0.36274 | 0.104425 | 0.226324 | 46.7–50.0 |
| GPS local | 0.36720 | 0.101909 | 0.223967 | 46.7–50.0 |
| GPS no distances | 0.37071 | 0.100160 | 0.222204 | 46.7–50.0 |
| GPS self | 0.31315 | 0.133329 | 0.266050 | 46.7–50.0 |
| GPS shuffled | 0.36342 | 0.105295 | 0.233173 | 46.7–50.0 |
| GPS rewired | 0.36274 | 0.104425 | 0.226324 | 46.7–50.0 |
| Distance fixed | 0.36323 | 0.103006 | 0.226427 | 46.7–50.0 |
| EGNN latent | 0.36142 | 0.104720 | 0.228689 | 46.7–50.0 |
| Distance disabled | 0.36642 | 0.102399 | 0.225718 | 46.7–50.0 |
| EGNN self | 0.36551 | 0.104725 | 0.230463 | 46.7–50.0 |
| EGNN shuffled | 0.35853 | 0.106645 | 0.229545 | 46.7–50.0 |
| EGNN rewired | 0.36142 | 0.104720 | 0.228689 | 46.7–50.0 |
| Local MLP | 0.35329 | 0.116472 | 0.246472 | 46.7–50.0 |
| Mean | 0.39515 | 0.094289 | 0.214188 | 46.7–50.0 |
| GCN | 0.39366 | 0.094480 | 0.214674 | 46.7–50.0 |
| GraphSAGE | 0.37959 | 0.100870 | 0.227195 | 46.7–50.0 |
| idw_k16_p2 | 0.41786 | 0.092172 | 0.202715 | 44.7–50.0 |
| uniform_k16 | 0.39844 | 0.094582 | 0.209673 | 44.7–50.0 |

## Matched patient effects and component controls

The following Bassiouni residual-route contrasts use the same20% V1 anchors and scored support. Effects are median paired patient differences; brackets are descriptive95% percentile intervals, unadjusted for multiple contrasts. They are not p-values or simultaneous bands. Every model/control/route/budget and each Wu patient effect is in the [paired exports](summary/PAIRED_PATIENT_EFFECTS.csv) and [uncertainty table](summary/PATIENT_UNCERTAINTY.csv).

| Candidate | Reference | Δr [interval] | Relative MSE% [interval] | ΔMAE [interval] |
| --- | --- | --- | --- | --- |
| GPS full | idw_k16_p2 | -0.0504 [-0.0629, -0.0389] | +13.4264 [+11.2710, +18.3865] | +0.0244 [+0.0204, +0.0295] |
| GPS full | Local MLP | +0.0217 [+0.0132, +0.0353] | -5.4759 [-9.0280, -4.6610] | -0.0106 [-0.0163, -0.0085] |
| GPS full | GCN | -0.0298 [-0.0352, -0.0238] | +10.8518 [+9.0284, +13.4617] | +0.0159 [+0.0143, +0.0196] |
| GPS full | GraphSAGE | -0.0134 [-0.0186, -0.0060] | +5.3568 [+3.4670, +7.5308] | +0.0068 [+0.0044, +0.0102] |
| GPS full | GPS local | -0.0050 [-0.0074, -0.0044] | +2.7228 [+2.0494, +3.6169] | +0.0039 [+0.0030, +0.0054] |
| GPS full | GPS no distances | -0.0121 [-0.0135, -0.0088] | +4.9138 [+3.6647, +5.7686] | +0.0064 [+0.0056, +0.0078] |
| GPS full | GPS self | +0.0561 [+0.0511, +0.0689] | -19.2594 [-22.4784, -16.7362] | -0.0352 [-0.0402, -0.0303] |
| GPS full | GPS shuffled | -0.0006 [-0.0057, +0.0036] | +1.8206 [-0.6598, +3.3526] | +0.0019 [-0.0019, +0.0030] |
| GPS full | GPS rewired | +0.0000 [-0.0000, +0.0000] | -0.0000 [-0.0000, +0.0000] | +0.0000 [-0.0000, +0.0000] |
| GPS local | idw_k16_p2 | -0.0418 [-0.0542, -0.0311] | +10.5093 [+8.4339, +14.3485] | +0.0204 [+0.0158, +0.0262] |
| GPS local | Local MLP | +0.0281 [+0.0219, +0.0390] | -8.2650 [-13.9130, -7.5101] | -0.0147 [-0.0223, -0.0128] |
| GPS local | GCN | -0.0233 [-0.0282, -0.0183] | +8.3845 [+5.7592, +9.9819] | +0.0128 [+0.0098, +0.0155] |
| GPS local | GraphSAGE | -0.0052 [-0.0127, -0.0017] | +2.0764 [+0.6364, +4.6618] | +0.0026 [+0.0008, +0.0061] |
| GPS no distances | idw_k16_p2 | -0.0359 [-0.0513, -0.0272] | +8.4534 [+6.0964, +12.1110] | +0.0175 [+0.0129, +0.0238] |
| GPS no distances | Local MLP | +0.0322 [+0.0260, +0.0461] | -10.3746 [-14.2648, -8.7761] | -0.0186 [-0.0236, -0.0151] |
| GPS no distances | GCN | -0.0196 [-0.0241, -0.0124] | +6.2582 [+4.2489, +7.2460] | +0.0095 [+0.0074, +0.0121] |
| GPS no distances | GraphSAGE | -0.0045 [-0.0085, +0.0031] | +0.5525 [-1.3328, +2.3452] | +0.0002 [-0.0025, +0.0030] |
| GPS self | idw_k16_p2 | -0.1043 [-0.1358, -0.0884] | +41.0551 [+36.1532, +51.0692] | +0.0603 [+0.0535, +0.0665] |
| GPS self | Local MLP | -0.0313 [-0.0487, -0.0208] | +16.4714 [+13.6775, +18.6281] | +0.0241 [+0.0196, +0.0271] |
| GPS self | GCN | -0.0845 [-0.1047, -0.0725] | +38.3237 [+34.2734, +43.7396] | +0.0515 [+0.0476, +0.0563] |
| GPS self | GraphSAGE | -0.0689 [-0.0813, -0.0544] | +30.3671 [+28.4758, +36.1643] | +0.0418 [+0.0372, +0.0466] |
| GPS shuffled | idw_k16_p2 | -0.0515 [-0.0620, -0.0362] | +13.5205 [+9.1462, +17.2020] | +0.0237 [+0.0191, +0.0321] |
| GPS shuffled | Local MLP | +0.0218 [+0.0191, +0.0282] | -7.8771 [-10.0197, -6.5667] | -0.0123 [-0.0161, -0.0112] |
| GPS shuffled | GCN | -0.0299 [-0.0359, -0.0216] | +10.6893 [+8.2578, +12.1268] | +0.0158 [+0.0133, +0.0213] |
| GPS shuffled | GraphSAGE | -0.0120 [-0.0200, -0.0062] | +3.3983 [+2.3276, +7.1799] | +0.0055 [+0.0036, +0.0101] |
| GPS rewired | idw_k16_p2 | -0.0504 [-0.0629, -0.0389] | +13.4264 [+11.2710, +18.3862] | +0.0244 [+0.0204, +0.0295] |
| GPS rewired | Local MLP | +0.0217 [+0.0132, +0.0353] | -5.4759 [-9.0280, -4.6610] | -0.0106 [-0.0163, -0.0085] |
| GPS rewired | GCN | -0.0298 [-0.0352, -0.0238] | +10.8518 [+9.0284, +13.4617] | +0.0159 [+0.0143, +0.0196] |
| GPS rewired | GraphSAGE | -0.0134 [-0.0186, -0.0060] | +5.3568 [+3.4670, +7.5308] | +0.0068 [+0.0044, +0.0102] |
| Distance fixed | idw_k16_p2 | -0.0424 [-0.0579, -0.0341] | +10.3942 [+7.2434, +14.7061] | +0.0222 [+0.0165, +0.0274] |
| Distance fixed | Local MLP | +0.0287 [+0.0198, +0.0372] | -9.3348 [-12.8840, -7.8668] | -0.0143 [-0.0193, -0.0128] |
| Distance fixed | GCN | -0.0223 [-0.0357, -0.0187] | +7.6922 [+6.5873, +9.6929] | +0.0137 [+0.0113, +0.0153] |
| Distance fixed | GraphSAGE | -0.0083 [-0.0164, -0.0050] | +1.8643 [+1.1658, +3.3732] | +0.0030 [+0.0020, +0.0051] |
| EGNN latent | idw_k16_p2 | -0.0480 [-0.0698, -0.0368] | +13.0083 [+9.0839, +18.3071] | +0.0252 [+0.0188, +0.0310] |
| EGNN latent | Local MLP | +0.0200 [+0.0162, +0.0260] | -7.3900 [-10.7409, -5.5995] | -0.0116 [-0.0160, -0.0088] |
| EGNN latent | GCN | -0.0283 [-0.0403, -0.0228] | +10.3474 [+7.2166, +13.1317] | +0.0177 [+0.0128, +0.0203] |
| EGNN latent | GraphSAGE | -0.0138 [-0.0214, -0.0072] | +4.3146 [+2.4545, +4.9327] | +0.0069 [+0.0044, +0.0086] |
| EGNN latent | Distance fixed | -0.0048 [-0.0074, -0.0008] | +1.6673 [+0.7107, +2.8664] | +0.0026 [+0.0014, +0.0046] |
| EGNN latent | Distance disabled | -0.0065 [-0.0096, -0.0024] | +1.9736 [+1.0100, +3.3448] | +0.0028 [+0.0020, +0.0053] |
| EGNN latent | EGNN self | +0.0010 [+0.0006, +0.0037] | -0.1804 [-0.7275, +0.5186] | +0.0000 [-0.0012, +0.0016] |
| EGNN latent | EGNN shuffled | +0.0039 [-0.0001, +0.0066] | -0.2937 [-0.9599, -0.1112] | -0.0004 [-0.0014, -0.0000] |
| EGNN latent | EGNN rewired | +0.0000 [+0.0000, +0.0000] | +0.0000 [-0.0000, +0.0000] | +0.0000 [-0.0000, +0.0000] |
| Distance disabled | idw_k16_p2 | -0.0405 [-0.0597, -0.0335] | +9.8403 [+7.1189, +14.1857] | +0.0213 [+0.0160, +0.0271] |
| Distance disabled | Local MLP | +0.0311 [+0.0208, +0.0400] | -10.0686 [-13.6626, -8.1110] | -0.0159 [-0.0202, -0.0134] |
| Distance disabled | GCN | -0.0216 [-0.0344, -0.0189] | +7.4531 [+5.9922, +9.1361] | +0.0134 [+0.0106, +0.0145] |
| Distance disabled | GraphSAGE | -0.0079 [-0.0150, -0.0027] | +1.5127 [+0.8827, +2.9378] | +0.0025 [+0.0017, +0.0045] |
| EGNN self | idw_k16_p2 | -0.0511 [-0.0677, -0.0441] | +12.3835 [+10.1729, +17.5820] | +0.0241 [+0.0209, +0.0306] |
| EGNN self | Local MLP | +0.0206 [+0.0144, +0.0254] | -7.7348 [-10.4805, -5.7935] | -0.0122 [-0.0162, -0.0097] |
| EGNN self | GCN | -0.0336 [-0.0402, -0.0265] | +10.1021 [+8.6704, +12.3428] | +0.0170 [+0.0146, +0.0195] |
| EGNN self | GraphSAGE | -0.0155 [-0.0241, -0.0080] | +3.7574 [+3.1709, +5.6761] | +0.0060 [+0.0054, +0.0090] |
| EGNN shuffled | idw_k16_p2 | -0.0493 [-0.0718, -0.0404] | +13.6558 [+9.9737, +17.2445] | +0.0256 [+0.0199, +0.0323] |
| EGNN shuffled | Local MLP | +0.0167 [+0.0116, +0.0242] | -6.0580 [-10.5926, -5.4036] | -0.0096 [-0.0157, -0.0080] |
| EGNN shuffled | GCN | -0.0316 [-0.0402, -0.0246] | +10.7360 [+7.9536, +13.7217] | +0.0166 [+0.0143, +0.0209] |
| EGNN shuffled | GraphSAGE | -0.0195 [-0.0267, -0.0071] | +4.6112 [+2.3100, +6.9704] | +0.0073 [+0.0043, +0.0099] |
| EGNN rewired | idw_k16_p2 | -0.0480 [-0.0698, -0.0368] | +13.0083 [+9.0839, +18.3071] | +0.0252 [+0.0188, +0.0310] |
| EGNN rewired | Local MLP | +0.0200 [+0.0162, +0.0260] | -7.3900 [-10.7409, -5.5995] | -0.0116 [-0.0160, -0.0088] |
| EGNN rewired | GCN | -0.0283 [-0.0403, -0.0228] | +10.3474 [+7.2166, +13.1317] | +0.0177 [+0.0128, +0.0203] |
| EGNN rewired | GraphSAGE | -0.0138 [-0.0214, -0.0072] | +4.3146 [+2.4545, +4.9328] | +0.0069 [+0.0044, +0.0086] |
| Local MLP | idw_k16_p2 | -0.0741 [-0.0873, -0.0549] | +23.6273 [+18.6730, +30.3872] | +0.0379 [+0.0318, +0.0462] |
| Local MLP | Local MLP | +0.0000 [+0.0000, +0.0000] | +0.0000 [+0.0000, +0.0000] | +0.0000 [+0.0000, +0.0000] |
| Local MLP | GCN | -0.0490 [-0.0640, -0.0415] | +20.9459 [+16.3577, +23.6678] | +0.0304 [+0.0257, +0.0349] |
| Local MLP | GraphSAGE | -0.0333 [-0.0420, -0.0273] | +13.0721 [+10.6217, +15.4627] | +0.0191 [+0.0159, +0.0241] |
| Mean | idw_k16_p2 | -0.0183 [-0.0277, -0.0134] | +2.4808 [+1.2463, +4.6721] | +0.0070 [+0.0057, +0.0114] |
| Mean | Local MLP | +0.0498 [+0.0423, +0.0658] | -17.4837 [-19.4211, -14.1923] | -0.0306 [-0.0355, -0.0261] |
| Mean | GCN | +0.0014 [+0.0007, +0.0021] | -0.2209 [-0.2954, -0.1637] | -0.0004 [-0.0005, -0.0003] |
| Mean | GraphSAGE | +0.0179 [+0.0113, +0.0238] | -5.8468 [-7.6043, -4.5686] | -0.0110 [-0.0135, -0.0095] |
| GCN | idw_k16_p2 | -0.0198 [-0.0286, -0.0150] | +2.5437 [+1.4677, +4.9562] | +0.0074 [+0.0060, +0.0120] |
| GCN | Local MLP | +0.0490 [+0.0415, +0.0640] | -17.3094 [-19.1382, -14.0564] | -0.0304 [-0.0349, -0.0257] |
| GCN | GCN | +0.0000 [+0.0000, +0.0000] | +0.0000 [+0.0000, +0.0000] | +0.0000 [+0.0000, +0.0000] |
| GCN | GraphSAGE | +0.0148 [+0.0103, +0.0221] | -5.6644 [-7.3673, -4.4134] | -0.0107 [-0.0129, -0.0091] |
| GraphSAGE | idw_k16_p2 | -0.0414 [-0.0538, -0.0234] | +8.2149 [+6.0927, +12.7503] | +0.0179 [+0.0156, +0.0239] |
| GraphSAGE | Local MLP | +0.0333 [+0.0273, +0.0420] | -11.5603 [-13.3920, -9.6000] | -0.0191 [-0.0241, -0.0159] |
| GraphSAGE | GCN | -0.0148 [-0.0221, -0.0103] | +6.0046 [+4.6178, +7.9532] | +0.0107 [+0.0091, +0.0129] |
| GraphSAGE | GraphSAGE | +0.0000 [+0.0000, +0.0000] | +0.0000 [+0.0000, +0.0000] | +0.0000 [+0.0000, +0.0000] |

## Geometric/long-range controls and limitations

All384 checkpoint-level rotation/reflection/translation/permutation checks passed with explicit precision-dependent tolerances. They are conditional on fixed embeddings and documented graph/landmark handling; no end-to-end CNN equivariance is claimed. The latent EGNN really updates private coordinates; its final coordinate update does not affect scalar prediction, and both the fixed-distance and latent models are E(2)-compatible. The latent model has128 extra nominal coefficients (last64 inactive for scalar loss). Full/local GPS has the same121,018 nominal parameters but changes key count and connection structure as well as range. Self/no-distance controls change effective function class; no-distance retains implicit graph/landmark geometry.

[Attention/reach and distant-feature diagnostics](diagnostics/) report physical versus assigned connectivity and attention mass beyond250/600/1200um. A remote perturbation tests dependency, not learned molecular correlation. Degree/length controls retain actual accepted-swap counts; zero/low swaps are limited evidence. Low-resolution external images (nominal112um native crop plus published encoder center crop), generic frozen features, source training in only two patients and transported calibration constrain morphology claims.

## Main text, supplement and novelty

The22-patient identical-H/equal-count placement experiment is the clearest new main-text evidence: it separates anchor placement from estimator choice and changed scored support. Keep all20 estimators, all neural routes/seeds, individual Wu effects, distance supports, resources and geometric checks fully available in the supplement/source tables. Keep a compact negative neural paired-effect/control summary in the main text to avoid hiding the failed hypotheses; move the full matrices and diagnostics to the supplement; it must retain errors alongside correlation and cannot establish biological long-range learning from attention reach.

The implemented architectures follow existing GraphGPS/EGNN ideas. This study supports a rigorous information-controlled reconstruction benchmark and conditional software properties, not a new equivariance theorem, general state-of-the-art claim, clinical benefit or demonstrated total assay savings. Any stronger geometry-specific neural novelty claim must be supported by the full ablation/error/patient matrix and future independent native-resolution cohorts.

## Preservation, statuses and author bundles

All407 protected prior files and74,253 archived files remain unchanged. The original29 publication targets/71 checks and810 frozen-model development primary-metric rows pass; new raw Bassiouni references are separately audited. The original selected primary heads remain missing (BLOCKED_INPUTS), and the78 post-v1 states are not replacements. Optional anchor-conditioned learning and full-image equivariance have SKIPPED_PROTOCOL status. All completed analyses remain visible regardless of direction. The [experiment matrix](summary/EXPERIMENT_MATRIX.csv) separates these statuses from negative findings.

Raw data/model binaries are excluded from Git. Full per-gene exports and48 new best/resume/optimizer states are preserved in author-local ZIP bundles with a [release-ready manifest](RELEASE_READY_MANIFEST.json). No public release upload or new credential exposure was performed; the manifest provides exact author paths/checksums for an approved release destination. The existing30-state bundle remains separate. PR3 stays draft and unmerged for author review.


## Numerical diagnostic summary

In the22 Bassiouni patients, median patient GPS-full physical connection fractions beyond250/600/1200um were0.9901/0.9424/0.7867; corresponding attention mass fractions were0.9905/0.9470/0.8002 (head/layer/seed means within patient). The coordinate-selected remote perturbation was23,349.7um away in TENX99: maximum output change was0.01283 for GPS full and exactly0 for GPS local and EGNN latent. This proves functional access only, not biological correlation.

The recovered-cohort degree/length control achieved zero accepted swaps in 64/66 patient–seed graphs (maximum1). It is therefore a limited control and not independent positive evidence for spatial specificity. Training-loop wall times summed to533.17seconds across48 fits; this excludes source feature preparation/model setup and includes the documented concurrent workload. Full per-configuration and standardized inference costs are reported above and in the benchmark CSV.
