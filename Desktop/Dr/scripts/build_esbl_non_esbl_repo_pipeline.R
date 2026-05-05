#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tools)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (basename(project_root) == "analysis") {
  project_root <- dirname(project_root)
}

default_image <- "/Users/jimmyg/Library/Mobile Documents/com~apple~CloudDocs/Images/Screenshot 2026-05-04 at 11.17.13.png"

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  image_path <- Sys.getenv("ESBL_SCREENSHOT_PATH", unset = default_image)
  idx <- 1L
  while (idx <= length(args)) {
    arg <- args[[idx]]
    if (startsWith(arg, "--image=")) {
      image_path <- sub("^--image=", "", arg)
    } else if (identical(arg, "--image")) {
      idx <- idx + 1L
      if (idx > length(args)) {
        stop("Missing value after --image", call. = FALSE)
      }
      image_path <- args[[idx]]
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
    idx <- idx + 1L
  }

  list(image_path = normalizePath(path.expand(image_path), winslash = "/", mustWork = FALSE))
}

run_command <- function(command, args, env = character()) {
  cat("\n> ", paste(c(command, args), collapse = " "), "\n", sep = "")
  status <- system2(command, args = args, stdout = "", stderr = "", env = env)
  if (!identical(status, 0L)) {
    stop("Command failed with exit code ", status, ": ", command, call. = FALSE)
  }
}

main <- function() {
  config <- parse_args()
  image_path <- config$image_path

  if (!file.exists(image_path)) {
    stop("Screenshot image not found: ", image_path, call. = FALSE)
  }

  image_ext <- file_ext(image_path)
  temp_image_path <- file.path(tempdir(), paste0("esbl_status_input", ifelse(nzchar(image_ext), paste0(".", image_ext), "")))
  ok_copy <- file.copy(image_path, temp_image_path, overwrite = TRUE)
  if (!isTRUE(ok_copy)) {
    stop("Failed to copy screenshot image to temporary path: ", temp_image_path, call. = FALSE)
  }

  extract_script <- file.path(project_root, "scripts", "extract_esbl_status_from_screenshot.py")
  daa_script <- file.path(project_root, "scripts", "build_esbl_non_esbl_daa_plots.R")
  abundance_script <- file.path(project_root, "scripts", "compute_esbl_species_mean_abundance.py")

  required_paths <- c(
    extract_script,
    daa_script,
    abundance_script,
    file.path(project_root, "data", "processed", "taxa", "TAXA_feature_metadata.tsv"),
    file.path(project_root, "data", "processed", "taxa", "derived", "TAXA_sgb_tss_for_maaslin2.tsv"),
    file.path(project_root, "data", "processed", "taxa", "derived", "TAXA_relative_abundance.tsv")
  )

  missing_paths <- required_paths[!file.exists(required_paths)]
  if (length(missing_paths) > 0) {
    stop(
      "Required files are missing:\n",
      paste0("- ", missing_paths, collapse = "\n"),
      call. = FALSE
    )
  }

  run_command(
    command = "/usr/bin/python3",
    args = c(extract_script, "--image", temp_image_path)
  )

  run_command(
    command = "/usr/local/bin/Rscript",
    args = c(daa_script)
  )

  run_command(
    command = "/usr/bin/python3",
    args = c(abundance_script)
  )

  expected_outputs <- c(
    file.path(project_root, "results", "esbl_status", "esbl_status_metaphlan_matched.tsv"),
    file.path(project_root, "results", "esbl_non_esbl_daa", "esbl_non_esbl_sgb_consensus_all.tsv"),
    file.path(project_root, "results", "esbl_non_esbl_daa", "plots", "esbl_non_esbl_sgb_consensus.png"),
    file.path(project_root, "results", "esbl_non_esbl_daa", "esbl_non_esbl_species_mean_abundance.tsv")
  )

  missing_outputs <- expected_outputs[!file.exists(expected_outputs)]
  if (length(missing_outputs) > 0) {
    stop(
      "Pipeline finished but some expected outputs are missing:\n",
      paste0("- ", missing_outputs, collapse = "\n"),
      call. = FALSE
    )
  }

  cat("\nESBL/non-ESBL DAA pipeline completed.\n")
  cat("Inputs used:\n")
  cat(paste0("- Screenshot image: ", image_path, "\n"))
  cat(paste0("- Temporary normalized image path: ", temp_image_path, "\n"))
  cat("Generated / refreshed outputs:\n")
  cat(paste0("- ", expected_outputs, "\n"))
}

main()
