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
