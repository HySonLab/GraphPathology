# Artifact distribution review

Reviewed on 2026-11-01. This release does not add a blanket license to the existing software repository or change any manuscript availability statement.

## Public-eligible artifacts

The author explicitly authorized distribution of verified, author-trained artifacts. The three newly trained replacement heads and 96 recovered ResNet18-based exploratory best states contain only author-trained head parameters, not ImageNet/torchvision encoder weights. Forty-eight optimizer/resume files are additional states, not 48 additional trained models. Their source fitting used HEST NCBI783/NCBI785, with TENX95 validation. The official HEST project states that HEST-1k/library/benchmark use CC-BY-NC-SA 4.0. Apply that noncommercial attribution/share-alike license conservatively to the author-owned release head artifacts and permitted derived numerical exports; retain the upstream attribution, license text and dataset citations. This artifact-specific grant does not relicense any separately supplied third-party software or the frozen manuscript.

HEST source: https://github.com/mahmoodlab/HEST ; https://huggingface.co/datasets/MahmoodLab/hest . Original ResNet18 encoder: https://download.pytorch.org/models/resnet18-f37072fd.pth ; torchvision code: https://github.com/pytorch/vision . The encoder is obtained separately from its official distributor and is not an asset here.

Wu evaluation data: https://doi.org/10.5281/zenodo.4739739 (Wu and Swarbrick; record identifies the original study). Bassiouni evaluation data: GEO GSE210616, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210616 ; Bassiouni et al., Cancer Research83:34–48(2023). Raw counts, tissue images and clinical tables are excluded. The public bundles contain first-party model parameters, configs, aggregate per-gene/panel metrics and verification records, not a redistribution of the GEO submission or article. Public availability of GEO is not presented as a blanket raw-data redistribution license.

## Phikon: distribution blocked pending clarification

Pinned encoder: owkin/phikon revision 057cc0295895c2df3dd7681a89680da6015cbefe. Its pinned model card expressly links the Owkin noncommercial license, not MIT/Apache or ordinary CC-BY-NC. The retrieved license hash matches the original local acquisition. The license limits use to nonprofit entities, excludes transfers of material/derivatives/results by or to for-profit organizations, and requires derivative/results licensing on the same terms. HEST has separate CC-BY-NC-SA obligations. Compatible onward distribution of trained Phikon-feature heads/results under both sets of obligations has not been established; no separate permission was found in the authorized local materials.

Therefore the 36 Phikon head binaries and Phikon-derived numerical export archives are NOT uploaded publicly or to the GitHub-hosted draft. They are retained locally, hash-verified and indexed for recovery. A metadata-only draft explains the blocker. This is not a claim that all such distribution is prohibited: clarification or a compatible permission from the applicable rights holders is needed before publication. Do not replace the custom license with the public ResNet artifact license. The draft contains no raw data, pretrained encoder, or uncleared model binary.

Official terms: https://github.com/owkin/HistoSSLscaling/blob/main/LICENSE.txt . Pinned model card and exact fetched license/legal-code bytes, hashes and URLs are recorded in licenses/SOURCES.json. Verification of already recovered Phikon heads and factual comparison to the already published paper are distinguished from permission to distribute those binary/result archives.

## GitHub constraints and integrity

Repository administration/push permissions were confirmed; no releases or tags existed at task start. GitHub documentation states each release asset must be below 2 GiB: https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases . Every upload must have its fetched size and SHA256 checked. Draft upload/download/clean-checkout verification precedes publication of public-eligible assets. No existing release asset is silently replaced. The versioned research tag identifies the tested code commit on the unmerged PR branch; main is unchanged.
