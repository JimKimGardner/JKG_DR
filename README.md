# JKG_DR

This repository now includes a reproducible pipeline for the two poster figures:

- `distance_to_baseline_bray_actual_time_mean_ci.png`
- `core_accessory_bray_actual_time_mean_ci.png`

from the MetaPhlAn SGB outputs.

## Poster Bray mean PNG pipeline

Main entry point:

```bash
/usr/local/bin/Rscript Desktop/Dr/scripts/build_poster_bray_mean_pngs_from_metaphlan.R
```

This wrapper runs the exact three analysis steps used for these figures:

1. `Desktop/Dr/scripts/01_compute_distance_to_baseline_taxa.R`
2. `Desktop/Dr/scripts/core_accessory_time_factor_models.py`
3. `Desktop/Dr/scripts/build_distance_plots_actual_time.R`

Additional documentation:

- `Desktop/Dr/docs/poster_bray_mean_png_pipeline.md`

## ESBL/non-ESBL DAA pipeline

Main entry point:

```bash
/usr/local/bin/Rscript Desktop/Dr/scripts/build_esbl_non_esbl_repo_pipeline.R \
  --image "/full/path/to/your/esbl_status_screenshot.png"
```

This wrapper runs:

1. `Desktop/Dr/scripts/extract_esbl_status_from_screenshot.py`
2. `Desktop/Dr/scripts/build_esbl_non_esbl_daa_plots.R`
3. `Desktop/Dr/scripts/compute_esbl_species_mean_abundance.py`

Additional documentation:

- `Desktop/Dr/docs/esbl_non_esbl_daa_pipeline.md`
# JKG_DR
