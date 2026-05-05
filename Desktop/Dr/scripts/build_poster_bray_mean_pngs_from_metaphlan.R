#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tibble)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (basename(project_root) == "analysis") {
  project_root <- dirname(project_root)
}

meta_input_path <- file.path(project_root, "data", "processed", "taxa", "TAXA_matrix.tsv")
metadata_path <- file.path(project_root, "data", "processed", "metadata", "METADATA_matrix.tsv")

distance_script <- file.path(project_root, "scripts", "01_compute_distance_to_baseline_taxa.R")
core_accessory_script <- file.path(project_root, "scripts", "core_accessory_time_factor_models.py")
poster_plot_script <- file.path(project_root, "scripts", "build_distance_plots_actual_time.R")

expected_outputs <- c(
  file.path(project_root, "results", "distance_to_baseline_taxa.csv"),
  file.path(project_root, "results", "core_accessory_distance_long.csv"),
  file.path(project_root, "results", "Poster plots", "distance_to_baseline_bray_actual_time_mean_ci.png"),
  file.path(project_root, "results", "Poster plots", "core_accessory_bray_actual_time_mean_ci.png")
)

required_paths <- c(meta_input_path, metadata_path, distance_script, core_accessory_script, poster_plot_script)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0) {
  stop(
    "Required files are missing:\n",
    paste0("- ", missing_paths, collapse = "\n"),
    call. = FALSE
  )
}

run_command <- function(command, args, env = character()) {
  cat("\n> ", paste(c(command, args), collapse = " "), "\n", sep = "")
  status <- system2(command, args = args, stdout = "", stderr = "", env = env)
  if (!identical(status, 0L)) {
    stop("Command failed with exit code ", status, ": ", command, call. = FALSE)
  }
}

run_command(
  command = "/usr/local/bin/Rscript",
  args = c(
    distance_script,
    paste0("--input=", meta_input_path),
    paste0("--output-dir=", file.path(project_root, "results"))
  )
)

python_env <- c("PYTHONHOME=", "PYTHONPATH=")
run_command(
  command = "/usr/bin/arch",
  args = c("-arm64", "/usr/bin/python3", core_accessory_script),
  env = python_env
)

run_command(
  command = "/usr/local/bin/Rscript",
  args = c(poster_plot_script)
)

missing_outputs <- expected_outputs[!file.exists(expected_outputs)]
if (length(missing_outputs) > 0) {
  stop(
    "Pipeline finished but some expected outputs are missing:\n",
    paste0("- ", missing_outputs, collapse = "\n"),
    call. = FALSE
  )
}

summary_tbl <- tibble(
  step = c(
    "Distance to personal baseline",
    "Core/accessory distance table",
    "Poster Bray mean PNG",
    "Poster core/accessory Bray mean PNG"
  ),
  output = expected_outputs
)

cat("\nPoster Bray mean-png pipeline completed.\n")
cat("Inputs used:\n")
cat(paste0("- MetaPhlAn table: ", meta_input_path, "\n"))
cat(paste0("- Metadata table: ", metadata_path, "\n"))
cat("Generated / refreshed outputs:\n")
cat(paste0("- ", summary_tbl$output, "\n"))
