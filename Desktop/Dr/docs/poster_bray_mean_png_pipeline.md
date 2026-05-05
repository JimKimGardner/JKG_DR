# Poster Bray Mean PNG Pipeline

This pipeline reproduces the two poster figures:

- `results/Poster plots/distance_to_baseline_bray_actual_time_mean_ci.png`
- `results/Poster plots/core_accessory_bray_actual_time_mean_ci.png`

from the project-standard MetaPhlAn SGB table and metadata.

## Entry Point

Run:

```bash
/usr/local/bin/Rscript scripts/build_poster_bray_mean_pngs_from_metaphlan.R
```

This wrapper executes the exact three code steps used in the project:

1. `scripts/01_compute_distance_to_baseline_taxa.R`
2. `scripts/core_accessory_time_factor_models.py`
3. `scripts/build_distance_plots_actual_time.R`

## Inputs

The pipeline expects these project-standard inputs:

- `data/processed/taxa/TAXA_matrix.tsv`
- `data/processed/metadata/METADATA_matrix.tsv`

## Step-by-step data flow

### 1. Distance to personal baseline

Script:

- `scripts/01_compute_distance_to_baseline_taxa.R`

Purpose:

- parses MetaPhlAn SGB sample names
- applies the project prevalence filter
- computes Bray-Curtis and Aitchison distance to each person's pre-travel baseline

Main output used later:

- `results/distance_to_baseline_taxa.csv`

### 2. Core/accessory distance table

Script:

- `scripts/core_accessory_time_factor_models.py`

Purpose:

- re-loads the same MetaPhlAn SGB matrix
- defines `core` and `accessory`
- computes compartment-specific Bray-Curtis and Aitchison distance-to-baseline values

Main output used later:

- `results/core_accessory_distance_long.csv`

### 3. Poster plot build

Script:

- `scripts/build_distance_plots_actual_time.R`

Purpose:

- joins the distance outputs to observed metadata dates
- derives the observed time axis used on the poster
- computes empirical means and 95% confidence intervals per nominal timepoint
- writes the final poster PNGs

Target outputs:

- `results/Poster plots/distance_to_baseline_bray_actual_time_mean_ci.png`
- `results/Poster plots/core_accessory_bray_actual_time_mean_ci.png`

Additional helper used by the plotting script:

- `scripts/poster_plot_style.R`

## Notes

- The plotting step uses the current poster notation where the internal baseline `w0` is displayed as poster `w-1`, and the internal return sample `w1` is displayed as poster `w0`.
- The mean Bray-Curtis plots currently use round, timepoint-colored points and 95% t-based confidence intervals.
- The wrapper script is intended as the cleanest reproducible entry point for GitHub; the scientific logic remains in the underlying analysis scripts listed above.
