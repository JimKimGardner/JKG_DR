# ESBL/non-ESBL DAA Pipeline

This pipeline reproduces the ESBL/non-ESBL differential abundance analysis workflow used in the project.

## Entry point

Run:

```bash
/usr/local/bin/Rscript Desktop/Dr/scripts/build_esbl_non_esbl_repo_pipeline.R \
  --image "/full/path/to/your/esbl_status_screenshot.png"
```

You can also provide the screenshot path via:

```bash
ESBL_SCREENSHOT_PATH="/full/path/to/your/esbl_status_screenshot.png" \
/usr/local/bin/Rscript Desktop/Dr/scripts/build_esbl_non_esbl_repo_pipeline.R
```

## What the wrapper does

The wrapper runs these three project scripts in sequence:

1. `Desktop/Dr/scripts/extract_esbl_status_from_screenshot.py`
2. `Desktop/Dr/scripts/build_esbl_non_esbl_daa_plots.R`
3. `Desktop/Dr/scripts/compute_esbl_species_mean_abundance.py`

## Inputs

Required project inputs:

- `Desktop/Dr/data/processed/taxa/TAXA_feature_metadata.tsv`
- `Desktop/Dr/data/processed/taxa/derived/TAXA_sgb_tss_for_maaslin2.tsv`
- `Desktop/Dr/data/processed/taxa/derived/TAXA_relative_abundance.tsv`

External input:

- one calibrated screenshot image containing the ESBL status matrix

## Step-by-step data flow

### 1. Screenshot extraction

Script:

- `Desktop/Dr/scripts/extract_esbl_status_from_screenshot.py`

Purpose:

- reads the calibrated screenshot
- classifies each cell as `esbl_positive`, `no_esbl`, or `na`
- maps screenshot timepoints to sequencing timepoints
- matches the extracted status calls to the MetaPhlAn sample IDs

Main outputs:

- `Desktop/Dr/results/esbl_status/esbl_status_from_screenshot.tsv`
- `Desktop/Dr/results/esbl_status/esbl_status_metaphlan_matched.tsv`
- `Desktop/Dr/results/esbl_status/esbl_status_from_screenshot_heatmap.png`

### 2. ESBL/non-ESBL DAA

Script:

- `Desktop/Dr/scripts/build_esbl_non_esbl_daa_plots.R`

Purpose:

- uses the matched ESBL status table
- compares `ESBL+` versus `No ESBL` samples
- computes MaAsLin2 and LinDA results
- builds the consensus classification
- writes the main SGB-level ESBL plot

Main outputs:

- `Desktop/Dr/results/esbl_non_esbl_daa/esbl_non_esbl_sgb_consensus_all.tsv`
- `Desktop/Dr/results/esbl_non_esbl_daa/plots/esbl_non_esbl_sgb_consensus.png`

Additional outputs are also produced for genus and family levels.

### 3. Mean-abundance comparison table

Script:

- `Desktop/Dr/scripts/compute_esbl_species_mean_abundance.py`

Purpose:

- computes per-SGB mean relative abundance in `ESBL+` and `No ESBL`
- adds raw differences and absolute differences
- adds Welch t-test and Mann-Whitney p/q values

Main output:

- `Desktop/Dr/results/esbl_non_esbl_daa/esbl_non_esbl_species_mean_abundance.tsv`

## Notes

- The screenshot extraction is calibrated to the provided ESBL status image layout; if the screenshot geometry changes substantially, the extraction script will need recalibration.
- The ESBL DAA currently uses the project’s pooled ESBL-vs-No-ESBL analysis setup and the current consensus thresholds defined in `build_esbl_non_esbl_daa_plots.R`.
