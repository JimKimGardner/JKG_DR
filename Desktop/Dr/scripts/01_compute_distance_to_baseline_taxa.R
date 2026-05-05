#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ggplot2)
})

DEFAULT_INPUT_PATH <- "data/processed/taxa/TAXA_matrix.tsv"
DEFAULT_OUTPUT_DIR <- "results"
DEFAULT_APPLY_PREVALENCE_FILTER <- TRUE
DEFAULT_MIN_PREVALENCE <- 0.10
DEFAULT_PSEUDOCOUNT <- 1e-6
ALLOWED_WEEKS <- c(0, 1, 2, 4, 6, 8, 10, 12, 16, 52)

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  config <- list(
    input = DEFAULT_INPUT_PATH,
    output_dir = DEFAULT_OUTPUT_DIR,
    apply_prevalence_filter = DEFAULT_APPLY_PREVALENCE_FILTER,
    min_prevalence = DEFAULT_MIN_PREVALENCE,
    pseudocount = DEFAULT_PSEUDOCOUNT,
    baseline_week = 0,
    output_stem = "distance_to_baseline_taxa"
  )

  for (arg in args) {
    if (str_starts(arg, "--input=")) {
      config$input <- str_remove(arg, "^--input=")
    } else if (str_starts(arg, "--output-dir=")) {
      config$output_dir <- str_remove(arg, "^--output-dir=")
    } else if (str_starts(arg, "--min-prevalence=")) {
      config$min_prevalence <- as.numeric(str_remove(arg, "^--min-prevalence="))
    } else if (str_starts(arg, "--pseudocount=")) {
      config$pseudocount <- as.numeric(str_remove(arg, "^--pseudocount="))
    } else if (str_starts(arg, "--baseline-week=")) {
      config$baseline_week <- as.numeric(str_remove(arg, "^--baseline-week="))
    } else if (str_starts(arg, "--output-stem=")) {
      config$output_stem <- str_remove(arg, "^--output-stem=")
    } else if (identical(arg, "--no-prevalence-filter")) {
      config$apply_prevalence_filter <- FALSE
    } else if (identical(arg, "--prevalence-filter")) {
      config$apply_prevalence_filter <- TRUE
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }

  if (!is.finite(config$min_prevalence) || config$min_prevalence < 0 || config$min_prevalence > 1) {
    stop("--min-prevalence must be between 0 and 1.", call. = FALSE)
  }

  if (!is.finite(config$pseudocount) || config$pseudocount <= 0) {
    stop("--pseudocount must be > 0.", call. = FALSE)
  }

  if (!is.finite(config$baseline_week) || config$baseline_week < 0 || !config$baseline_week %in% ALLOWED_WEEKS) {
    stop("--baseline-week must be one of: ", paste(ALLOWED_WEEKS, collapse = ", "), call. = FALSE)
  }

  config
}

ensure_output_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

trim_delimiters <- function(x) {
  x |>
    str_replace_all("^[[:punct:]_[:space:]]+", "") |>
    str_replace_all("[[:punct:]_[:space:]]+$", "") |>
    str_squish()
}

is_sample_like <- function(x) {
  str_detect(x, regex("w0*\\d+", ignore_case = TRUE))
}

load_taxa <- function(path) {
  if (!file.exists(path)) {
    stop("Input file not found: ", path, call. = FALSE)
  }

  raw_df <- read.delim(
    file = path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (ncol(raw_df) < 2 || nrow(raw_df) < 1) {
    stop("Input table must contain at least one ID column and one data column.", call. = FALSE)
  }

  id_values <- as.character(raw_df[[1]])
  data_df <- raw_df[-1]

  numeric_df <- as.data.frame(
    lapply(data_df, function(x) suppressWarnings(as.numeric(x))),
    check.names = FALSE
  )

  if (anyNA(as.matrix(numeric_df))) {
    warning("NA values detected after numeric conversion. Replacing them with 0.")
    numeric_df[is.na(numeric_df)] <- 0
  }

  col_sample_hits <- sum(is_sample_like(names(data_df)))
  row_sample_hits <- sum(is_sample_like(id_values))

  numeric_mat <- as.matrix(numeric_df)

  if (col_sample_hits >= row_sample_hits) {
    rownames(numeric_mat) <- id_values
    taxa_mat <- numeric_mat
    orientation <- "features_in_rows_samples_in_columns"
  } else {
    rownames(numeric_mat) <- id_values
    taxa_mat <- t(numeric_mat)
    orientation <- "features_in_columns_samples_in_rows"
  }

  if (is.null(rownames(taxa_mat)) || is.null(colnames(taxa_mat))) {
    stop("Failed to recover feature and sample identifiers from the input table.", call. = FALSE)
  }

  if (!all(str_detect(rownames(taxa_mat), "(^|\\|)t__SGB"))) {
    sgb_hits <- sum(str_detect(rownames(taxa_mat), "(^|\\|)t__SGB"))
    warning(
      "Not all features look like SGB-level features. Keeping only explicit SGB entries. ",
      "Detected SGB features: ", sgb_hits, "/", nrow(taxa_mat)
    )
    taxa_mat <- taxa_mat[str_detect(rownames(taxa_mat), "(^|\\|)t__SGB"), , drop = FALSE]
  }

  if (nrow(taxa_mat) == 0) {
    stop("No SGB-level features remained after filtering.", call. = FALSE)
  }

  list(
    taxa_mat = taxa_mat,
    orientation = orientation,
    n_features_input = nrow(taxa_mat),
    n_samples_input = ncol(taxa_mat)
  )
}

parse_sample_names <- function(sample_ids) {
  parse_one <- function(sample_id, sample_order) {
    sample_id <- as.character(sample_id)
    timepoint_loc <- str_locate(sample_id, regex("w0*\\d+", ignore_case = TRUE))

    if (any(is.na(timepoint_loc))) {
      return(tibble(
        sample_id = sample_id,
        sample_order = sample_order,
        patient_id = NA_character_,
        timepoint_raw = NA_character_,
        timepoint = NA_character_,
        weeks_post_return = NA_real_,
        parse_status = "unparsed",
        exclusion_reason = "timepoint_not_recognized"
      ))
    }

    timepoint_raw <- str_sub(sample_id, timepoint_loc[1], timepoint_loc[2]) |> str_to_lower()
    week_digits <- str_remove(timepoint_raw, "^w0*")
    weeks_post_return <- ifelse(identical(week_digits, ""), 0, suppressWarnings(as.numeric(week_digits)))

    patient_prefix <- if (timepoint_loc[1] > 1) {
      str_sub(sample_id, 1, timepoint_loc[1] - 1)
    } else {
      ""
    }
    patient_suffix <- if (timepoint_loc[2] < str_length(sample_id)) {
      str_sub(sample_id, timepoint_loc[2] + 1, str_length(sample_id))
    } else {
      ""
    }

    patient_id <- trim_delimiters(patient_prefix)
    if (patient_id == "") {
      patient_id <- trim_delimiters(patient_suffix)
    }

    if (!is.finite(weeks_post_return) || patient_id == "") {
      return(tibble(
        sample_id = sample_id,
        sample_order = sample_order,
        patient_id = ifelse(patient_id == "", NA_character_, patient_id),
        timepoint_raw = timepoint_raw,
        timepoint = NA_character_,
        weeks_post_return = NA_real_,
        parse_status = "unparsed",
        exclusion_reason = "patient_or_timepoint_not_reliably_parsed"
      ))
    }

    canonical_timepoint <- paste0("w", as.integer(weeks_post_return))

    if (!weeks_post_return %in% ALLOWED_WEEKS) {
      return(tibble(
        sample_id = sample_id,
        sample_order = sample_order,
        patient_id = patient_id,
        timepoint_raw = timepoint_raw,
        timepoint = canonical_timepoint,
        weeks_post_return = weeks_post_return,
        parse_status = "unsupported_timepoint",
        exclusion_reason = "timepoint_not_in_analysis_set"
      ))
    }

    tibble(
      sample_id = sample_id,
      sample_order = sample_order,
      patient_id = patient_id,
      timepoint_raw = timepoint_raw,
      timepoint = canonical_timepoint,
      weeks_post_return = weeks_post_return,
      parse_status = "ok",
      exclusion_reason = NA_character_
    )
  }

  map2_dfr(sample_ids, seq_along(sample_ids), parse_one)
}

annotate_eligibility <- function(sample_info, baseline_week = 0) {
  ok_samples <- sample_info |>
    filter(parse_status == "ok") |>
    arrange(sample_order)

  patient_summary <- ok_samples |>
    group_by(patient_id) |>
    summarise(
      n_baseline = sum(weeks_post_return == baseline_week),
      n_post = sum(weeks_post_return > baseline_week),
      baseline_sample_id = first(sample_id[weeks_post_return == baseline_week]),
      baseline_candidates = paste(sample_id[weeks_post_return == baseline_week], collapse = "; "),
      .groups = "drop"
    ) |>
    mutate(
      eligible = n_baseline > 0 & n_post > 0,
      multiple_baseline = n_baseline > 1
    )

  annotated <- sample_info |>
    left_join(patient_summary, by = "patient_id") |>
    mutate(
      eligibility_status = case_when(
        parse_status != "ok" ~ exclusion_reason,
        is.na(eligible) | n_baseline == 0 ~ paste0("excluded_no_w", baseline_week),
        n_post == 0 ~ "excluded_no_post",
        eligible ~ "included",
        TRUE ~ "excluded_other"
      )
    )

  list(sample_info = annotated, patient_summary = patient_summary)
}

apply_prevalence_filter <- function(feature_mat, apply_filter = TRUE, min_prevalence = 0.10) {
  prevalence_tbl <- tibble(
    feature_id = rownames(feature_mat),
    prevalence = rowMeans(feature_mat > 0),
    n_present_samples = rowSums(feature_mat > 0)
  )

  filtered_mat <- feature_mat

  if (apply_filter) {
    keep_features <- prevalence_tbl$prevalence >= min_prevalence
    filtered_mat <- filtered_mat[keep_features, , drop = FALSE]
    prevalence_tbl <- prevalence_tbl[keep_features, , drop = FALSE]
  }

  nonzero_keep <- rowSums(filtered_mat > 0) > 0
  filtered_mat <- filtered_mat[nonzero_keep, , drop = FALSE]
  prevalence_tbl <- prevalence_tbl[nonzero_keep, , drop = FALSE]

  list(
    matrix = filtered_mat,
    prevalence = prevalence_tbl
  )
}

make_relative_abundance <- function(feature_mat) {
  sample_totals <- colSums(feature_mat)

  if (any(sample_totals <= 0)) {
    zero_samples <- names(sample_totals)[sample_totals <= 0]
    stop(
      "Cannot compute relative abundance: samples with zero total abundance detected: ",
      paste(zero_samples, collapse = ", "),
      call. = FALSE
    )
  }

  sweep(feature_mat, 2, sample_totals, "/")
}

compute_clr <- function(relative_mat, pseudocount = 1e-6) {
  adjusted <- relative_mat + pseudocount
  adjusted <- sweep(adjusted, 2, colSums(adjusted), "/")
  log_adjusted <- log(adjusted)
  sweep(log_adjusted, 2, colMeans(log_adjusted), "-")
}

compute_distance_to_baseline <- function(sample_info, bray_mat, aitchison_mat, n_features_after_filter, baseline_week = 0) {
  post_samples <- sample_info |>
    filter(eligibility_status == "included", weeks_post_return > baseline_week) |>
    arrange(patient_id, weeks_post_return, sample_order)

  if (nrow(post_samples) == 0) {
    stop("No post-travel samples remained after parsing and eligibility filtering.", call. = FALSE)
  }

  result_tbl <- post_samples |>
    mutate(
      distance_bray_to_baseline = map2_dbl(
        sample_id,
        baseline_sample_id,
        ~ bray_mat[.x, .y]
      ),
      distance_aitchison_to_baseline = map2_dbl(
        sample_id,
        baseline_sample_id,
        ~ aitchison_mat[.x, .y]
      ),
      n_features_after_filter = n_features_after_filter
    ) |>
    select(
      sample_id,
      patient_id,
      baseline_sample_id,
      timepoint,
      weeks_post_return,
      distance_bray_to_baseline,
      distance_aitchison_to_baseline,
      n_features_after_filter
    )

  list(
    result_tbl = result_tbl
  )
}

make_plots <- function(result_tbl, output_dir, output_stem = "distance_to_baseline_taxa") {
  make_spaghetti_plot <- function(data, y_var, y_label, output_path) {
    plot_obj <- ggplot(data, aes(x = weeks_post_return, y = .data[[y_var]], group = patient_id, color = patient_id)) +
      geom_line(alpha = 0.45, linewidth = 0.6) +
      geom_point(size = 1.8, alpha = 0.8) +
      stat_summary(
        aes(group = 1),
        fun = mean,
        geom = "line",
        color = "black",
        linewidth = 1.1
      ) +
      stat_summary(
        aes(group = 1),
        fun = mean,
        geom = "point",
        color = "black",
        size = 2.3
      ) +
      scale_x_continuous(breaks = sort(unique(data$weeks_post_return))) +
      labs(
        x = "Weeks post return",
        y = y_label,
        color = "Patient"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "right"
      )

    ggsave(output_path, plot = plot_obj, width = 9, height = 5.5)
  }

  bray_plot_path <- file.path(output_dir, paste0("plot_", output_stem, "_bray_spaghetti.pdf"))
  aitchison_plot_path <- file.path(output_dir, paste0("plot_", output_stem, "_aitchison_spaghetti.pdf"))

  make_spaghetti_plot(
    data = result_tbl,
    y_var = "distance_bray_to_baseline",
    y_label = "Bray-Curtis distance to individual baseline",
    output_path = bray_plot_path
  )

  make_spaghetti_plot(
    data = result_tbl,
    y_var = "distance_aitchison_to_baseline",
    y_label = "Aitchison distance to individual baseline",
    output_path = aitchison_plot_path
  )

  c(bray_plot_path, aitchison_plot_path)
}

save_outputs <- function(
  result_tbl,
  sample_info,
  patient_summary,
  prevalence_tbl,
  summary_lines,
  output_dir,
  output_stem = "distance_to_baseline_taxa"
) {
  result_path <- file.path(output_dir, paste0(output_stem, ".csv"))
  summary_path <- file.path(output_dir, paste0(output_stem, "_summary.txt"))
  exclusion_path <- file.path(output_dir, paste0(output_stem, "_exclusions.csv"))
  prevalence_path <- file.path(output_dir, paste0(output_stem, "_prevalence_table.csv"))

  readr::write_csv(result_tbl, result_path)
  readr::write_csv(sample_info |> filter(eligibility_status != "included"), exclusion_path)
  readr::write_csv(prevalence_tbl, prevalence_path)
  writeLines(summary_lines, con = summary_path)

  c(result_path, summary_path, exclusion_path, prevalence_path)
}

main <- function() {
  config <- parse_args()
  ensure_output_dir(config$output_dir)

  taxa_input <- load_taxa(config$input)
  sample_info <- parse_sample_names(colnames(taxa_input$taxa_mat))
  eligibility <- annotate_eligibility(sample_info, baseline_week = config$baseline_week)
  sample_info <- eligibility$sample_info
  patient_summary <- eligibility$patient_summary

  included_samples <- sample_info |>
    filter(eligibility_status == "included") |>
    arrange(sample_order) |>
    pull(sample_id)

  if (length(included_samples) < 2) {
    stop("Fewer than two eligible samples remained after filtering.", call. = FALSE)
  }

  taxa_mat_included <- taxa_input$taxa_mat[, included_samples, drop = FALSE]

  prevalence_res <- apply_prevalence_filter(
    feature_mat = taxa_mat_included,
    apply_filter = config$apply_prevalence_filter,
    min_prevalence = config$min_prevalence
  )

  taxa_mat_filtered <- prevalence_res$matrix

  if (nrow(taxa_mat_filtered) == 0) {
    stop("No features remained after prevalence and zero-only filtering.", call. = FALSE)
  }

  relative_mat <- make_relative_abundance(taxa_mat_filtered)
  clr_mat <- compute_clr(relative_mat, pseudocount = config$pseudocount)

  bray_mat <- as.matrix(vegan::vegdist(t(relative_mat), method = "bray"))
  aitchison_mat <- as.matrix(dist(t(clr_mat), method = "euclidean"))

  baseline_res <- compute_distance_to_baseline(
    sample_info = sample_info,
    bray_mat = bray_mat,
    aitchison_mat = aitchison_mat,
    n_features_after_filter = nrow(taxa_mat_filtered),
    baseline_week = config$baseline_week
  )

  result_tbl <- baseline_res$result_tbl
  plot_paths <- make_plots(result_tbl, config$output_dir, output_stem = config$output_stem)

  unparsed_samples <- sample_info |>
    filter(parse_status == "unparsed") |>
    pull(sample_id)

  unsupported_timepoints <- sample_info |>
    filter(parse_status == "unsupported_timepoint") |>
    pull(sample_id)

  multiple_baseline_patients <- patient_summary |>
    filter(multiple_baseline) |>
    transmute(label = paste0(patient_id, " [", baseline_candidates, "]")) |>
    pull(label)

  excluded_samples <- sample_info |>
    filter(eligibility_status != "included")

  included_patients <- patient_summary |>
    filter(eligible) |>
    pull(patient_id)

  summary_lines <- c(
    "Distance to individual within-person baseline: TAXA / MetaPhlAn4",
    paste0("Run timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("Input file: ", normalizePath(config$input, winslash = "/")),
    paste0("Detected orientation: ", taxa_input$orientation),
    paste0("Allowed timepoints: ", paste0("w", ALLOWED_WEEKS, collapse = ", ")),
    paste0("Baseline week: w", config$baseline_week),
    paste0("Prevalence filter enabled: ", config$apply_prevalence_filter),
    paste0("Minimum prevalence: ", config$min_prevalence),
    paste0("CLR pseudocount: ", config$pseudocount),
    "",
    paste0("Input samples: ", taxa_input$n_samples_input),
    paste0("Input SGB features: ", taxa_input$n_features_input),
    paste0("Included patients: ", length(included_patients)),
    paste0("Included samples (baseline + post): ", length(included_samples)),
    paste0("Post-travel samples in outcome table: ", nrow(result_tbl)),
    paste0("Samples excluded: ", nrow(excluded_samples)),
    paste0("Features after filtering: ", nrow(taxa_mat_filtered)),
    "",
    paste0("Patients with multiple baseline samples: ", length(multiple_baseline_patients)),
    if (length(multiple_baseline_patients) > 0) {
      paste0("Multiple baseline details: ", paste(multiple_baseline_patients, collapse = "; "))
    } else {
      "Multiple baseline details: none"
    },
    if (length(unparsed_samples) > 0) {
      paste0("Unparsed sample names: ", paste(unparsed_samples, collapse = ", "))
    } else {
      "Unparsed sample names: none"
    },
    if (length(unsupported_timepoints) > 0) {
      paste0("Unsupported timepoints excluded: ", paste(unsupported_timepoints, collapse = ", "))
    } else {
      "Unsupported timepoints excluded: none"
    },
    "",
    "Exclusion counts by reason:",
    excluded_samples |>
      count(eligibility_status, name = "n") |>
      mutate(line = paste0("- ", eligibility_status, ": ", n)) |>
      pull(line),
    "",
    "Generated files:"
  )

  output_paths <- save_outputs(
    result_tbl = result_tbl,
    sample_info = sample_info,
    patient_summary = patient_summary,
    prevalence_tbl = prevalence_res$prevalence,
    summary_lines = summary_lines,
    output_dir = config$output_dir,
    output_stem = config$output_stem
  )

  summary_lines <- c(summary_lines, paste0("- ", normalizePath(plot_paths, winslash = "/")))
  summary_lines <- c(summary_lines, paste0("- ", normalizePath(output_paths, winslash = "/")))
  writeLines(summary_lines, con = file.path(config$output_dir, paste0(config$output_stem, "_summary.txt")))

  cat("\nDistance-to-baseline analysis completed.\n")
  cat("Baseline week:", config$baseline_week, "\n")
  cat("Included patients:", length(included_patients), "\n")
  cat("Excluded samples:", nrow(excluded_samples), "\n")
  cat(
    "Unparsed sample names:",
    if (length(unparsed_samples) > 0) paste(unparsed_samples, collapse = ", ") else "none",
    "\n"
  )
  cat("Patients with multiple baseline samples:", length(multiple_baseline_patients), "\n")
  cat("Generated files:\n")
  for (path in c(output_paths, plot_paths)) {
    cat(" -", normalizePath(path, winslash = "/"), "\n")
  }
}

main()
