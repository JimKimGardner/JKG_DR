#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(scales)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "scripts", "poster_plot_style.R"))

metadata_path <- file.path(project_root, "data", "processed", "metadata", "METADATA_matrix.tsv")
distance_path <- file.path(project_root, "results", "distance_to_baseline_taxa.csv")
distance_w4_path <- file.path(project_root, "results", "distance_to_w4_baseline_taxa.csv")
core_accessory_distance_long_path <- file.path(project_root, "results", "core_accessory_distance_long.csv")
aitchison_contributors_path <- file.path(project_root, "results", "aitchison_taxon_contributors_sample_level.csv")
output_dir <- file.path(project_root, "results", "Poster plots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

parse_metadata_date <- function(x) {
  if (is.na(x)) {
    return(as.Date(NA))
  }

  value <- trimws(as.character(x))
  if (identical(value, "") || tolower(value) %in% c("na", "nan", "none", "dropout", "*dropout")) {
    return(as.Date(NA))
  }

  value <- gsub("\\*", "", value)
  value <- gsub(",", ".", value, fixed = TRUE)
  value <- gsub("[^0-9.]", "", value)

  parts <- Filter(nzchar, strsplit(value, "\\.", fixed = FALSE)[[1]])
  if (length(parts) != 3) {
    return(as.Date(NA))
  }

  day <- suppressWarnings(as.integer(parts[[1]]))
  month <- suppressWarnings(as.integer(parts[[2]]))
  year <- suppressWarnings(as.integer(parts[[3]]))

  if (is.na(day) || is.na(month) || is.na(year)) {
    return(as.Date(NA))
  }

  if (nchar(parts[[3]]) == 2) {
    year <- ifelse(year <= 50, 2000 + year, 1900 + year)
  }

  parsed <- suppressWarnings(as.Date(sprintf("%04d-%02d-%02d", year, month, day)))
  parsed
}

format_patient_id <- function(id_value) {
  numeric_id <- suppressWarnings(as.integer(id_value))
  ifelse(
    is.na(numeric_id),
    NA_character_,
    sprintf("p%03d", numeric_id)
  )
}

metadata_time_columns <- c(
  w0 = "before__Date_DD_MM_YY",
  w1 = "after__DateQuestionnaire_DD_MM_YY",
  w2 = "followup__Date_Week_2_DD_MM_YY",
  w4 = "followup__Date_Week_4_DD_MM_YY",
  w6 = "followup__Date_Week_6_DD_MM_YY",
  w8 = "followup__Date_Week_8_DD_MM_YY",
  w10 = "followup__Date_Week_10_DD_MM_YY",
  w12 = "followup__Date_Week_12_DD_MM_YY",
  w16 = "followup__Date_Week_16_DD_MM_YY",
  w20 = "followup__Date_Week_20_DD_MM_YY",
  w52 = "followup__Date_Week_52_DD_MM_YY"
)

nominal_days <- c(
  w1 = 7,
  w2 = 14,
  w4 = 28,
  w6 = 42,
  w8 = 56,
  w10 = 70,
  w12 = 84,
  w16 = 112,
  w20 = 140,
  w52 = 364
)

timepoint_levels <- c("w0", "w1", "w2", "w4", "w6", "w8", "w10", "w12", "w16", "w20", "w52")
timepoint_palette <- c(
  w0 = "#475569",
  w1 = "#9f1239",
  w2 = "#ea580c",
  w4 = "#ca8a04",
  w6 = "#65a30d",
  w8 = "#0f766e",
  w10 = "#0891b2",
  w12 = "#2563eb",
  w16 = "#7c3aed",
  w20 = "#c026d3",
  w52 = "#374151"
)
timepoint_display_labels <- c(
  w0 = "w-1",
  w1 = "w0",
  w2 = "w2",
  w4 = "w4",
  w6 = "w6",
  w8 = "w8",
  w10 = "w10",
  w12 = "w12",
  w16 = "w16",
  w20 = "w20",
  w52 = "w52"
)
phase_levels <- c("pre", "return", "early", "late")
phase_shape_values <- c(
  pre = 16,
  return = 2,
  early = 8,
  late = 15
)
phase_labels <- c(
  pre = "Pre-travel (w-1)",
  return = "Return (w0)",
  early = "Early (w2-w8)",
  late = "Late (w10-w52)"
)

phase_from_timepoint <- function(timepoint) {
  dplyr::case_when(
    timepoint == "w0" ~ "pre",
    timepoint == "w1" ~ "return",
    timepoint %in% c("w2", "w4", "w6", "w8") ~ "early",
    timepoint %in% c("w10", "w12", "w16", "w20", "w52") ~ "late",
    TRUE ~ NA_character_
  )
}

mean_plot_patient_time_limits <- c(
  p009 = "w6",
  p014 = "w8"
)

apply_mean_plot_patient_overrides <- function(data) {
  data |>
    mutate(
      patient_time_limit = unname(mean_plot_patient_time_limits[patient_id]),
      timepoint_index = match(as.character(timepoint), timepoint_levels),
      time_limit_index = match(patient_time_limit, timepoint_levels),
      keep_for_mean_plot = ifelse(
        is.na(patient_time_limit),
        TRUE,
        !is.na(timepoint_index) & !is.na(time_limit_index) & timepoint_index <= time_limit_index
      )
    ) |>
    filter(keep_for_mean_plot) |>
    select(-patient_time_limit, -timepoint_index, -time_limit_index, -keep_for_mean_plot)
}

build_patient_palette <- function(patient_ids) {
  unique_ids <- sort(unique(as.character(patient_ids)))
  n_ids <- length(unique_ids)
  reference_palette <- c(
    "#F04E45", "#4F8FD9", "#63B35D", "#A56BC8", "#FF8A1D",
    "#1A1A1A", "#A87445", "#F38AC1", "#A8A8A8", "#72D1B3",
    "#FF9A63", "#9AB4E8", "#18BCCB", "#95CF5A", "#FF7B69"
  )

  if (n_ids <= length(reference_palette)) {
    palette_values <- reference_palette[seq_len(n_ids)]
  } else {
    extra_colors <- hue_pal(h = c(15, 375), c = 70, l = 60)(n_ids - length(reference_palette))
    palette_values <- c(reference_palette, extra_colors)
  }

  names(palette_values) <- unique_ids
  palette_values
}

phase_breaks_from_data <- function(data) {
  phase_levels[phase_levels %in% unique(as.character(data$phase))]
}

metadata_df <- read_tsv(metadata_path, show_col_types = FALSE) |>
  mutate(patient_id = format_patient_id(ID))

distance_df <- read_csv(distance_path, show_col_types = FALSE)
distance_w4_df <- read_csv(distance_w4_path, show_col_types = FALSE)
core_accessory_distance_long_df <- read_csv(core_accessory_distance_long_path, show_col_types = FALSE)
aitchison_contributors_df <- read_csv(aitchison_contributors_path, show_col_types = FALSE)

included_patients <- sort(unique(distance_df$patient_id))

metadata_df <- metadata_df |>
  filter(patient_id %in% included_patients)

build_time_rows <- function(row_df) {
  patient_id <- row_df$patient_id[[1]]
  observed_dates <- lapply(metadata_time_columns, function(column_name) {
    parse_metadata_date(row_df[[column_name]][[1]])
  })
  baseline_w0_date <- observed_dates[["w0"]]
  baseline_w1_date <- observed_dates[["w1"]]

  inferred_candidates <- lapply(names(nominal_days), function(timepoint) {
    observed_date <- observed_dates[[timepoint]]
    if (is.na(observed_date)) {
      return(as.Date(NA))
    }
    observed_date - nominal_days[[timepoint]]
  })

  inferred_candidates <- as.Date(unlist(inferred_candidates), origin = "1970-01-01")
  inferred_candidates <- inferred_candidates[!is.na(inferred_candidates)]

  inferred_return_date <- if (length(inferred_candidates) == 0) {
    as.Date(NA)
  } else {
    sort(inferred_candidates)[ceiling(length(inferred_candidates) / 2)]
  }

  tibble(
    patient_id = patient_id,
    timepoint = names(metadata_time_columns),
    source_metadata_column = unname(metadata_time_columns),
    observed_date_raw = vapply(metadata_time_columns, function(column_name) row_df[[column_name]][[1]], character(1)),
    observed_date = as.Date(unlist(observed_dates), origin = "1970-01-01"),
    inferred_return_date = inferred_return_date,
    baseline_w0_date = baseline_w0_date
  ) |>
    mutate(
      days_post_return_inferred = as.numeric(observed_date - inferred_return_date),
      days_since_w0 = as.numeric(observed_date - baseline_w0_date),
      days_since_w1 = as.numeric(observed_date - baseline_w1_date)
    )
}

time_lookup_df <- bind_rows(lapply(seq_len(nrow(metadata_df)), function(idx) {
  build_time_rows(metadata_df[idx, , drop = FALSE])
}))

prepare_distance_plot_data <- function(distance_input_df, baseline_timepoint, baseline_weeks_post_return, baseline_sample_label) {
  joined_df <- distance_input_df |>
    left_join(
      time_lookup_df |>
        select(patient_id, timepoint, source_metadata_column, observed_date_raw, observed_date, inferred_return_date, baseline_w0_date, days_post_return_inferred, days_since_w0, days_since_w1),
      by = c("patient_id", "timepoint")
    ) |>
    mutate(
      nominal_days_post_return = weeks_post_return * 7,
      day_deviation_from_nominal = days_post_return_inferred - nominal_days_post_return
    ) |>
    mutate(
      use_observed_time = !is.na(days_post_return_inferred) & abs(day_deviation_from_nominal) <= 60,
      actual_time_days = ifelse(use_observed_time, days_post_return_inferred, NA_real_)
    )

  baseline_df <- time_lookup_df |>
    filter(timepoint == baseline_timepoint, !is.na(observed_date), !is.na(days_post_return_inferred)) |>
    semi_join(joined_df, by = "patient_id") |>
    transmute(
      sample_id = paste0(patient_id, "-", baseline_sample_label, "-baseline"),
      patient_id = patient_id,
      baseline_sample_id = NA_character_,
      timepoint = baseline_timepoint,
      weeks_post_return = baseline_weeks_post_return,
      distance_bray_to_baseline = 0,
      distance_aitchison_to_baseline = 0,
      n_features_after_filter = NA_real_,
      source_metadata_column = source_metadata_column,
      observed_date_raw = observed_date_raw,
      observed_date = observed_date,
      inferred_return_date = inferred_return_date,
      baseline_w0_date = baseline_w0_date,
      days_post_return_inferred = days_post_return_inferred,
      days_since_w0 = days_since_w0,
      days_since_w1 = days_since_w1,
      nominal_days_post_return = baseline_weeks_post_return * 7,
      day_deviation_from_nominal = NA_real_,
      use_observed_time = TRUE,
      actual_time_days = days_post_return_inferred
    )

  plot_df <- bind_rows(joined_df, baseline_df) |>
    mutate(
      timepoint = factor(timepoint, levels = timepoint_levels, ordered = TRUE),
      phase = factor(phase_from_timepoint(as.character(timepoint)), levels = phase_levels)
    ) |>
    filter(use_observed_time) |>
    arrange(patient_id, actual_time_days, timepoint)

  list(
    joined_df = joined_df,
    baseline_df = baseline_df,
    plot_df = plot_df
  )
}

pretravel_data <- prepare_distance_plot_data(
  distance_input_df = distance_df,
  baseline_timepoint = "w0",
  baseline_weeks_post_return = 0,
  baseline_sample_label = "w00"
)

w4_reference_data <- prepare_distance_plot_data(
  distance_input_df = distance_w4_df,
  baseline_timepoint = "w4",
  baseline_weeks_post_return = 4,
  baseline_sample_label = "w04"
)

prepare_core_accessory_plot_data <- function(distance_input_df) {
  joined_df <- distance_input_df |>
    filter(distance_type == "bray", compartment %in% c("core", "accessory")) |>
    left_join(
      time_lookup_df |>
        select(patient_id, timepoint, source_metadata_column, observed_date_raw, observed_date, inferred_return_date, baseline_w0_date, days_post_return_inferred, days_since_w0, days_since_w1),
      by = c("patient_id", "timepoint")
    ) |>
    mutate(
      nominal_days_post_return = weeks_post_return * 7,
      day_deviation_from_nominal = days_post_return_inferred - nominal_days_post_return,
      use_observed_time = !is.na(days_post_return_inferred) & abs(day_deviation_from_nominal) <= 60,
      actual_time_days = ifelse(use_observed_time, days_post_return_inferred, NA_real_),
      compartment = factor(compartment, levels = c("core", "accessory"))
    )

  baseline_df <- time_lookup_df |>
    filter(timepoint == "w0", !is.na(observed_date), !is.na(days_post_return_inferred)) |>
    semi_join(joined_df, by = "patient_id") |>
    tidyr::crossing(compartment = factor(c("core", "accessory"), levels = c("core", "accessory"))) |>
    transmute(
      sample_id = paste0(patient_id, "-w00-", as.character(compartment), "-baseline"),
      patient_id = patient_id,
      baseline_sample_id = NA_character_,
      timepoint = "w0",
      weeks_post_return = 0,
      compartment = compartment,
      distance_type = "bray",
      distance_to_baseline = 0,
      n_features = NA_real_,
      source_metadata_column = source_metadata_column,
      observed_date_raw = observed_date_raw,
      observed_date = observed_date,
      inferred_return_date = inferred_return_date,
      baseline_w0_date = baseline_w0_date,
      days_post_return_inferred = days_post_return_inferred,
      days_since_w0 = days_since_w0,
      days_since_w1 = days_since_w1,
      nominal_days_post_return = 0,
      day_deviation_from_nominal = NA_real_,
      use_observed_time = TRUE,
      actual_time_days = days_post_return_inferred
    )

  plot_df <- bind_rows(joined_df, baseline_df) |>
    mutate(
      timepoint = factor(timepoint, levels = timepoint_levels, ordered = TRUE),
      phase = factor(phase_from_timepoint(as.character(timepoint)), levels = phase_levels)
    ) |>
    filter(use_observed_time) |>
    arrange(compartment, patient_id, actual_time_days, timepoint)

  list(
    joined_df = joined_df,
    baseline_df = baseline_df,
    plot_df = plot_df
  )
}

core_accessory_bray_data <- prepare_core_accessory_plot_data(core_accessory_distance_long_df)
pretravel_bray_mean_data <- apply_mean_plot_patient_overrides(pretravel_data$plot_df)
core_accessory_bray_mean_data <- apply_mean_plot_patient_overrides(core_accessory_bray_data$plot_df)

patient_palette_df <- tibble(
  patient_id = sort(unique(pretravel_data$plot_df$patient_id)),
  color_hex = unname(build_patient_palette(pretravel_data$plot_df$patient_id))
)

plot_data_path <- file.path(output_dir, "distance_to_baseline_actual_time_data.tsv")
write_tsv(pretravel_data$plot_df, plot_data_path, na = "")

core_accessory_plot_data_path <- file.path(output_dir, "core_accessory_bray_actual_time_data.tsv")
write_tsv(core_accessory_bray_data$plot_df, core_accessory_plot_data_path, na = "")

bray_mean_plot_data_path <- file.path(output_dir, "distance_to_baseline_bray_mean_plot_data.tsv")
write_tsv(pretravel_bray_mean_data, bray_mean_plot_data_path, na = "")

core_accessory_mean_plot_data_path <- file.path(output_dir, "core_accessory_bray_mean_plot_data.tsv")
write_tsv(core_accessory_bray_mean_data, core_accessory_mean_plot_data_path, na = "")

plot_data_w4_path <- file.path(output_dir, "distance_to_w4_baseline_actual_time_data.tsv")
write_tsv(w4_reference_data$plot_df, plot_data_w4_path, na = "")

patient_palette_path <- file.path(output_dir, "distance_to_baseline_patient_palette.tsv")
write_tsv(patient_palette_df, patient_palette_path, na = "")

build_info_lines <- c(
  paste("Metadata file:", metadata_path),
  paste("Metadata time anchor:", "Observed sample/questionnaire dates from METADATA_matrix.tsv"),
  paste("Returned time variable:", "days_post_return_inferred"),
  paste("Return date handling:", "Inferred per patient as the median of observed_date - nominal_days across available post-travel timepoints because ReturnDate was empty in the raw after-travel metadata"),
  paste("Distance table:", distance_path),
  paste("Patients in distance table:", length(unique(distance_df$patient_id))),
  paste("Post-travel rows in distance table:", nrow(distance_df)),
  paste("Post-travel rows with usable observed time:", sum(pretravel_data$joined_df$use_observed_time, na.rm = TRUE)),
  paste("Dropped post-travel rows with implausible observed time:", sum(!pretravel_data$joined_df$use_observed_time, na.rm = TRUE)),
  paste("Poster mean-plot patient overrides:", "p009 up to w6; p014 up to w8"),
  paste("Baseline rows added with observed pre-travel date:", nrow(pretravel_data$baseline_df)),
  paste("Baseline rows missing pre-travel date and therefore not plotted:", length(unique(distance_df$patient_id)) - nrow(pretravel_data$baseline_df))
)

build_info_path <- file.path(output_dir, "distance_to_baseline_actual_time_build_info.txt")
writeLines(build_info_lines, con = build_info_path)

build_info_w4_lines <- c(
  paste("Metadata file:", metadata_path),
  paste("Metadata time anchor:", "Observed sample/questionnaire dates from METADATA_matrix.tsv"),
  paste("Returned time variable:", "days_post_return_inferred"),
  paste("Return date handling:", "Inferred per patient as the median of observed_date - nominal_days across available post-travel timepoints because ReturnDate was empty in the raw after-travel metadata"),
  paste("Distance table:", distance_w4_path),
  paste("Patients in distance table:", length(unique(distance_w4_df$patient_id))),
  paste("Post-w4 rows in distance table:", nrow(distance_w4_df)),
  paste("Post-w4 rows with usable observed time:", sum(w4_reference_data$joined_df$use_observed_time, na.rm = TRUE)),
  paste("Dropped post-w4 rows with implausible observed time:", sum(!w4_reference_data$joined_df$use_observed_time, na.rm = TRUE)),
  paste("w4 reference rows added at distance 0:", nrow(w4_reference_data$baseline_df)),
  paste("Patients missing usable w4 metadata date and therefore not plotted as reference:", length(unique(distance_w4_df$patient_id)) - nrow(w4_reference_data$baseline_df))
)

build_info_w4_path <- file.path(output_dir, "distance_to_w4_baseline_actual_time_build_info.txt")
writeLines(build_info_w4_lines, con = build_info_w4_path)

make_publication_plot <- function(data, y_column, y_label, title, pdf_path, png_path) {
  phase_breaks <- phase_breaks_from_data(data)
  patient_palette <- build_patient_palette(data$patient_id)

  plot <- ggplot(data, aes(x = actual_time_days, y = .data[[y_column]], group = patient_id)) +
    geom_vline(xintercept = 0, color = "#6b7280", linetype = "dashed", linewidth = 0.6) +
    geom_line(aes(color = patient_id), linewidth = 0.55, alpha = 0.6, show.legend = FALSE) +
    geom_point(aes(color = patient_id, shape = phase), alpha = 0.95, size = 2.8, stroke = 0.35) +
    scale_color_manual(
      values = patient_palette,
      breaks = names(patient_palette),
      drop = FALSE,
      name = "Person"
    ) +
    scale_shape_manual(
      values = phase_shape_values,
      breaks = phase_breaks,
      labels = phase_labels[phase_breaks],
      drop = FALSE,
      name = "Phase"
    ) +
    labs(
      title = title,
      subtitle = "Colors identify persons; point shapes mark the broad temporal phases at the observed sampling days",
      x = "Observed sampling time (days after inferred return)",
      y = y_label,
      caption = "Pre-travel baseline samples are shown at distance 0 when a metadata-based pre-travel date was available."
    ) +
    poster_theme(base_size = 15, legend_position = "right") +
    theme(
      legend.text = element_text(size = 11)
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(shape = 16, linewidth = 0, size = 3.2, alpha = 1)),
      shape = guide_legend(order = 2, override.aes = list(color = "#1f2933", size = 3.4, alpha = 1))
    )

  save_poster_plot(plot, pdf_path, png_path, width = 13.5, height = 8)
}

make_individual_patient_plot <- function(data, y_column, y_label, title, pdf_path, png_path) {
  patient_palette <- build_patient_palette(data$patient_id)

  plot <- ggplot(
    data,
    aes(
      x = actual_time_days,
      y = .data[[y_column]],
      group = patient_id,
      color = patient_id
    )
  ) +
    geom_vline(xintercept = 0, color = "#6b7280", linetype = "dashed", linewidth = 0.5) +
    geom_line(linewidth = 0.6, alpha = 0.72, show.legend = FALSE) +
    geom_point(aes(shape = phase), size = 2.35, alpha = 0.95, stroke = 0.35) +
    facet_wrap(~ patient_id, ncol = 3, scales = "fixed") +
    scale_color_manual(values = patient_palette, drop = FALSE) +
    scale_shape_manual(
      values = phase_shape_values,
      breaks = phase_breaks_from_data(data),
      labels = phase_labels[phase_breaks_from_data(data)],
      drop = FALSE,
      name = "Phase"
    ) +
    labs(
      title = title,
      subtitle = "Each panel shows one person; point shapes indicate the temporal phases",
      x = "Observed sampling time (days after inferred return)",
      y = y_label,
      caption = "Dashed vertical line marks the inferred day of return; baseline samples are plotted at distance 0."
    ) +
    poster_theme(base_size = 14) +
    theme(
      strip.text = element_text(size = 11.5)
    ) +
    guides(
      color = "none",
      shape = guide_legend(override.aes = list(color = "#1f2933", size = 3.2, alpha = 1))
    )

  save_poster_plot(plot, pdf_path, png_path, width = 13.5, height = 15)
}

make_individual_patient_date_plot <- function(data, y_column, y_label, title, pdf_path, png_path) {
  patient_palette <- build_patient_palette(data$patient_id)

  plot <- ggplot(
    data,
    aes(
      x = observed_date,
      y = .data[[y_column]],
      group = patient_id,
      color = patient_id
    )
  ) +
    geom_line(linewidth = 0.6, alpha = 0.72, show.legend = FALSE) +
    geom_point(aes(shape = phase), size = 2.35, alpha = 0.95, stroke = 0.35) +
    facet_wrap(~ patient_id, ncol = 3, scales = "free_x") +
    scale_color_manual(values = patient_palette, drop = FALSE) +
    scale_shape_manual(
      values = phase_shape_values,
      breaks = phase_breaks_from_data(data),
      labels = phase_labels[phase_breaks_from_data(data)],
      drop = FALSE,
      name = "Phase"
    ) +
    scale_x_date(date_breaks = "2 months", date_labels = "%b\n%Y") +
    labs(
      title = title,
      subtitle = "Each panel shows one person on the recorded metadata dates; point shapes indicate the temporal phases",
      x = "Observed sampling date",
      y = y_label,
      caption = "X-axis uses the actual metadata dates of the corresponding samples."
    ) +
    poster_theme(base_size = 14) +
    theme(
      strip.text = element_text(size = 11.5),
      axis.text.x = element_text(size = 8.5),
    ) +
    guides(
      color = "none",
      shape = guide_legend(override.aes = list(color = "#1f2933", size = 3.2, alpha = 1))
    )

  save_poster_plot(plot, pdf_path, png_path, width = 13.5, height = 15)
}

make_mean_error_plot <- function(data, y_column, y_label, title, pdf_path, png_path) {
  summary_df <- data |>
    mutate(timepoint = factor(as.character(timepoint), levels = timepoint_levels, ordered = TRUE)) |>
    group_by(timepoint) |>
    summarise(
      n_patients = n_distinct(patient_id),
      mean_time_from_w1_days = mean(days_since_w1, na.rm = TRUE),
      mean_distance = mean(.data[[y_column]], na.rm = TRUE),
      sd_distance = sd(.data[[y_column]], na.rm = TRUE),
      se_distance = ifelse(n_patients > 1, sd_distance / sqrt(n_patients), NA_real_),
      t_critical = ifelse(n_patients > 1, qt(0.975, df = n_patients - 1), NA_real_),
      lower_ci = mean_distance - t_critical * se_distance,
      upper_ci = mean_distance + t_critical * se_distance,
      .groups = "drop"
    ) |>
    filter(!is.na(mean_time_from_w1_days), !is.na(mean_distance))

  plot <- ggplot(summary_df, aes(x = mean_time_from_w1_days, y = mean_distance, group = 1)) +
    geom_vline(xintercept = 0, color = "#6b7280", linetype = "dashed", linewidth = 0.6) +
    geom_errorbar(
      aes(ymin = lower_ci, ymax = upper_ci),
      width = 4,
      color = "#8d99a6",
      linewidth = 0.8,
      na.rm = TRUE
    ) +
    geom_line(color = "#4f5d75", linewidth = 0.95) +
    geom_point(aes(fill = timepoint), shape = 21, color = "#1f2937", size = 3.8, stroke = 0.6) +
    scale_fill_manual(
      values = timepoint_palette,
      breaks = timepoint_levels[timepoint_levels %in% unique(as.character(summary_df$timepoint))],
      labels = timepoint_display_labels[timepoint_levels[timepoint_levels %in% unique(as.character(summary_df$timepoint))]],
      drop = FALSE,
      name = "Nominal timepoint"
    ) +
    labs(
      title = title,
      subtitle = "Points show empirical means across persons; return (w0) is fixed at x = 0 and error bars show 95 percent confidence intervals",
      x = "Mean observed sampling time per timepoint (days since return / w0 sample)",
      y = y_label,
      caption = "Means and 95% CIs are computed across patients within each nominal timepoint; legend labels follow the poster notation w-1, w0, w2, ... and x-positions use the mean observed day relative to each patient's return / w0 sample."
    ) +
    poster_theme(base_size = 15, legend_position = "right") +
    theme(
      legend.text = element_text(size = 11)
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(shape = 21, color = "#1f2937", size = 3.8, alpha = 1, stroke = 0.6)
      )
    )

  save_poster_plot(plot, pdf_path, png_path, width = 13.5, height = 8)
}

make_aitchison_contributor_plot <- function(data, patient_id, timepoint, pdf_path, png_path, top_n_each = 10) {
  sample_df <- data |>
    filter(patient_id == !!patient_id, timepoint == !!timepoint) |>
    mutate(
      clr_delta = as.numeric(clr_delta),
      squared_contribution = as.numeric(squared_contribution),
      contribution_share = as.numeric(contribution_share)
    ) |>
    filter(!is.na(clr_delta), !is.na(squared_contribution))

  if (nrow(sample_df) == 0) {
    stop(sprintf("No contributor rows found for %s at %s.", patient_id, timepoint))
  }

  top_positive <- sample_df |>
    filter(clr_delta > 0) |>
    arrange(desc(squared_contribution)) |>
    slice_head(n = top_n_each) |>
    mutate(direction = "Higher than baseline")

  top_negative <- sample_df |>
    filter(clr_delta < 0) |>
    arrange(desc(squared_contribution)) |>
    slice_head(n = top_n_each) |>
    mutate(direction = "Lower than baseline")

  plot_df <- bind_rows(top_positive, top_negative) |>
    mutate(
      direction = factor(direction, levels = c("Higher than baseline", "Lower than baseline")),
      feature_short = ifelse(
        nchar(feature_label) > 58,
        paste0(substr(feature_label, 1, 55), "..."),
        feature_label
      ),
      signed_contribution = ifelse(direction == "Lower than baseline", -squared_contribution, squared_contribution)
    ) |>
    arrange(direction, signed_contribution) |>
    mutate(feature_short = factor(feature_short, levels = feature_short))

  direction_colors <- c(
    "Higher than baseline" = "#c2410c",
    "Lower than baseline" = "#1d4ed8"
  )

  total_share <- sum(plot_df$contribution_share, na.rm = TRUE)

  plot <- ggplot(plot_df, aes(x = signed_contribution, y = feature_short, fill = direction)) +
    geom_col(width = 0.75, alpha = 0.95, show.legend = FALSE) +
    facet_grid(direction ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = direction_colors, drop = FALSE) +
    scale_x_continuous(
      labels = function(x) format(abs(x), trim = TRUE, scientific = FALSE, digits = 3),
      expand = expansion(mult = c(0.05, 0.08))
    ) +
    labs(
      title = sprintf("Top taxon contributors for %s at %s", patient_id, timepoint),
      subtitle = "Largest positive and negative CLR shifts contributing to the Aitchison distance from the personal pre-travel baseline",
      x = "Squared contribution to Aitchison distance",
      y = NULL,
      caption = sprintf(
        "Top %d positive and top %d negative contributors are shown. Together they explain %.1f%% of the reported Aitchison distance signal.",
        top_n_each,
        top_n_each,
        100 * total_share
      )
    ) +
    poster_theme(base_size = 14, major_y = FALSE) +
    theme(
      axis.text.y = element_text(size = 10.5),
      strip.text = element_text(size = 12)
    )

  save_poster_plot(plot, pdf_path, png_path, width = 13.5, height = 11)
}

make_core_accessory_spaghetti_plot <- function(data, pdf_path, png_path) {
  patient_palette <- build_patient_palette(data$patient_id)
  compartment_labels <- c(core = "Core", accessory = "Accessory")

  plot <- ggplot(
    data,
    aes(
      x = actual_time_days,
      y = distance_to_baseline,
      group = patient_id
    )
  ) +
    geom_vline(xintercept = 0, color = "#6b7280", linetype = "dashed", linewidth = 0.5) +
    geom_line(aes(color = patient_id), linewidth = 0.55, alpha = 0.58, show.legend = FALSE) +
    geom_point(aes(color = patient_id, shape = phase), alpha = 0.95, size = 2.5, stroke = 0.35) +
    facet_wrap(~ compartment, ncol = 1, labeller = as_labeller(compartment_labels)) +
    scale_color_manual(
      values = patient_palette,
      breaks = names(patient_palette),
      drop = FALSE,
      name = "Person"
    ) +
    scale_shape_manual(
      values = phase_shape_values,
      breaks = phase_breaks_from_data(data),
      labels = phase_labels[phase_breaks_from_data(data)],
      drop = FALSE,
      name = "Phase"
    ) +
    labs(
      title = "Bray-Curtis abundance trajectories for core and accessory compartments",
      subtitle = "Colors identify persons; point shapes mark the broad temporal phases at the observed sampling days",
      x = "Observed sampling time (days after inferred return)",
      y = "Bray-Curtis distance to pre-travel baseline",
      caption = "Baseline samples are shown at distance 0 when a metadata-based pre-travel date was available."
    ) +
    poster_theme(base_size = 15, legend_position = "right") +
    theme(
      strip.text = element_text(size = 12),
      legend.text = element_text(size = 11)
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(shape = 16, linewidth = 0, size = 3.2, alpha = 1)),
      shape = guide_legend(order = 2, override.aes = list(color = "#1f2933", size = 3.3, alpha = 1))
    )

  save_poster_plot(plot, pdf_path, png_path, width = 13.5, height = 10.5)
}

make_core_accessory_mean_plot <- function(data, pdf_path, png_path) {
  compartment_labels <- c(core = "Core", accessory = "Accessory")
  compartment_colors <- c(core = "#4F8FD9", accessory = "#FF8A1D")

  summary_df <- data |>
    group_by(compartment, timepoint) |>
    summarise(
      n_patients = n_distinct(patient_id),
      mean_time_from_w1_days = mean(days_since_w1, na.rm = TRUE),
      mean_distance = mean(distance_to_baseline, na.rm = TRUE),
      sd_distance = sd(distance_to_baseline, na.rm = TRUE),
      se_distance = ifelse(n_patients > 1, sd_distance / sqrt(n_patients), NA_real_),
      t_critical = ifelse(n_patients > 1, qt(0.975, df = n_patients - 1), NA_real_),
      lower_ci = mean_distance - t_critical * se_distance,
      upper_ci = mean_distance + t_critical * se_distance,
      .groups = "drop"
    ) |>
    filter(!is.na(mean_time_from_w1_days), !is.na(mean_distance)) |>
    mutate(
      compartment = factor(compartment, levels = c("core", "accessory"))
    )

  plot <- ggplot(summary_df, aes(x = mean_time_from_w1_days, y = mean_distance, group = 1)) +
    geom_vline(xintercept = 0, color = "#6b7280", linetype = "dashed", linewidth = 0.6) +
    geom_errorbar(
      aes(ymin = lower_ci, ymax = upper_ci, color = compartment),
      width = 4,
      linewidth = 0.8,
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    geom_line(aes(color = compartment), linewidth = 1.05, show.legend = FALSE) +
    geom_point(aes(fill = timepoint), shape = 21, color = "#1f2937", size = 3.8, stroke = 0.6) +
    facet_wrap(~ compartment, ncol = 1, labeller = as_labeller(compartment_labels)) +
    scale_color_manual(values = compartment_colors, drop = FALSE) +
    scale_fill_manual(
      values = timepoint_palette,
      breaks = timepoint_levels[timepoint_levels %in% unique(as.character(summary_df$timepoint))],
      labels = timepoint_display_labels[timepoint_levels[timepoint_levels %in% unique(as.character(summary_df$timepoint))]],
      drop = FALSE,
      name = "Nominal timepoint"
    ) +
    labs(
      title = "Mean Bray-Curtis distance for core and accessory compartments",
      subtitle = "Points show empirical means across persons; return (w0) is fixed at x = 0 and error bars show 95 percent confidence intervals",
      x = "Mean observed sampling time per timepoint (days since return / w0 sample)",
      y = "Mean Bray-Curtis distance to pre-travel baseline",
      caption = "Means and 95% CIs are computed across patients within each nominal timepoint; legend labels follow the poster notation w-1, w0, w2, ... and x-positions are centered on each patient's return / w0 sample."
    ) +
    poster_theme(base_size = 15, legend_position = "right") +
    theme(
      strip.text = element_text(size = 12),
      legend.text = element_text(size = 11)
    ) +
    guides(
      color = "none",
      fill = guide_legend(
        override.aes = list(shape = 21, color = "#1f2937", size = 3.8, alpha = 1, stroke = 0.6)
      )
    )

  save_poster_plot(plot, pdf_path, png_path, width = 13.5, height = 10.5)
}

bray_pdf <- file.path(output_dir, "distance_to_baseline_bray_actual_time.pdf")
bray_png <- file.path(output_dir, "distance_to_baseline_bray_actual_time.png")
aitchison_pdf <- file.path(output_dir, "distance_to_baseline_aitchison_actual_time.pdf")
aitchison_png <- file.path(output_dir, "distance_to_baseline_aitchison_actual_time.png")
bray_w4_pdf <- file.path(output_dir, "distance_to_w4_baseline_bray_actual_time.pdf")
bray_w4_png <- file.path(output_dir, "distance_to_w4_baseline_bray_actual_time.png")
aitchison_w4_pdf <- file.path(output_dir, "distance_to_w4_baseline_aitchison_actual_time.pdf")
aitchison_w4_png <- file.path(output_dir, "distance_to_w4_baseline_aitchison_actual_time.png")
bray_mean_pdf <- file.path(output_dir, "distance_to_baseline_bray_actual_time_mean_ci.pdf")
bray_mean_png <- file.path(output_dir, "distance_to_baseline_bray_actual_time_mean_ci.png")
core_accessory_bray_spaghetti_pdf <- file.path(output_dir, "core_accessory_bray_actual_time_spaghetti.pdf")
core_accessory_bray_spaghetti_png <- file.path(output_dir, "core_accessory_bray_actual_time_spaghetti.png")
core_accessory_bray_mean_pdf <- file.path(output_dir, "core_accessory_bray_actual_time_mean_ci.pdf")
core_accessory_bray_mean_png <- file.path(output_dir, "core_accessory_bray_actual_time_mean_ci.png")
aitchison_individual_pdf <- file.path(output_dir, "distance_to_baseline_aitchison_actual_time_individuals.pdf")
aitchison_individual_png <- file.path(output_dir, "distance_to_baseline_aitchison_actual_time_individuals.png")
aitchison_individual_dates_pdf <- file.path(output_dir, "distance_to_baseline_aitchison_metadata_dates_individuals.pdf")
aitchison_individual_dates_png <- file.path(output_dir, "distance_to_baseline_aitchison_metadata_dates_individuals.png")
aitchison_mean_pdf <- file.path(output_dir, "distance_to_baseline_aitchison_actual_time_mean_ci.pdf")
aitchison_mean_png <- file.path(output_dir, "distance_to_baseline_aitchison_actual_time_mean_ci.png")
p015_w16_contributors_pdf <- file.path(output_dir, "p015_w16_aitchison_top_taxa_contributors.pdf")
p015_w16_contributors_png <- file.path(output_dir, "p015_w16_aitchison_top_taxa_contributors.png")

make_publication_plot(
  data = pretravel_data$plot_df,
  y_column = "distance_bray_to_baseline",
  y_label = "Bray-Curtis distance to pre-travel baseline",
  title = "Bray-Curtis distance to pre-travel baseline over observed sampling times",
  pdf_path = bray_pdf,
  png_path = bray_png
)

make_publication_plot(
  data = pretravel_data$plot_df,
  y_column = "distance_aitchison_to_baseline",
  y_label = "Aitchison distance to pre-travel baseline",
  title = "Aitchison distance to pre-travel baseline over observed sampling times",
  pdf_path = aitchison_pdf,
  png_path = aitchison_png
)

make_publication_plot(
  data = w4_reference_data$plot_df,
  y_column = "distance_bray_to_baseline",
  y_label = "Bray-Curtis distance to w4 reference",
  title = "Bray-Curtis distance to w4 reference over observed sampling times",
  pdf_path = bray_w4_pdf,
  png_path = bray_w4_png
)

make_publication_plot(
  data = w4_reference_data$plot_df,
  y_column = "distance_aitchison_to_baseline",
  y_label = "Aitchison distance to w4 reference",
  title = "Aitchison distance to w4 reference over observed sampling times",
  pdf_path = aitchison_w4_pdf,
  png_path = aitchison_w4_png
)

make_mean_error_plot(
  data = pretravel_bray_mean_data,
  y_column = "distance_bray_to_baseline",
  y_label = "Mean Bray-Curtis distance to pre-travel baseline",
  title = "Mean Bray-Curtis distance to pre-travel baseline over observed sampling times",
  pdf_path = bray_mean_pdf,
  png_path = bray_mean_png
)

make_core_accessory_spaghetti_plot(
  data = core_accessory_bray_data$plot_df,
  pdf_path = core_accessory_bray_spaghetti_pdf,
  png_path = core_accessory_bray_spaghetti_png
)

make_core_accessory_mean_plot(
  data = core_accessory_bray_mean_data,
  pdf_path = core_accessory_bray_mean_pdf,
  png_path = core_accessory_bray_mean_png
)

make_individual_patient_plot(
  data = pretravel_data$plot_df,
  y_column = "distance_aitchison_to_baseline",
  y_label = "Aitchison distance to pre-travel baseline",
  title = "Aitchison distance to pre-travel baseline for individual patients",
  pdf_path = aitchison_individual_pdf,
  png_path = aitchison_individual_png
)

make_individual_patient_date_plot(
  data = pretravel_data$plot_df,
  y_column = "distance_aitchison_to_baseline",
  y_label = "Aitchison distance to pre-travel baseline",
  title = "Aitchison distance to pre-travel baseline for individual patients on metadata dates",
  pdf_path = aitchison_individual_dates_pdf,
  png_path = aitchison_individual_dates_png
)

make_mean_error_plot(
  data = pretravel_data$plot_df,
  y_column = "distance_aitchison_to_baseline",
  y_label = "Mean Aitchison distance to pre-travel baseline",
  title = "Mean Aitchison distance to pre-travel baseline over observed sampling times",
  pdf_path = aitchison_mean_pdf,
  png_path = aitchison_mean_png
)

make_aitchison_contributor_plot(
  data = aitchison_contributors_df,
  patient_id = "p015",
  timepoint = "w16",
  pdf_path = p015_w16_contributors_pdf,
  png_path = p015_w16_contributors_png,
  top_n_each = 10
)

cat(paste0("Metadata file used: ", metadata_path, "\n"))
cat("Time variable used: days_post_return_inferred\n")
cat("Time variable derivation: observed sample dates from metadata anchored to a patient-specific inferred return date because the raw ReturnDate field was empty\n")
cat(paste0("Distance table used: ", distance_path, "\n"))
cat(paste0("Distance table used for w4 reference: ", distance_w4_path, "\n"))
cat(paste0("Aitchison contributor table used: ", aitchison_contributors_path, "\n"))
cat(paste0("Plot data written: ", plot_data_path, "\n"))
cat(paste0("Core/accessory plot data written: ", core_accessory_plot_data_path, "\n"))
cat(paste0("Bray mean-plot data written: ", bray_mean_plot_data_path, "\n"))
cat(paste0("Core/accessory Bray mean-plot data written: ", core_accessory_mean_plot_data_path, "\n"))
cat(paste0("Plot data written for w4 reference: ", plot_data_w4_path, "\n"))
cat(paste0("Patient palette written: ", patient_palette_path, "\n"))
cat(paste0("Build info written: ", build_info_path, "\n"))
cat(paste0("Build info written for w4 reference: ", build_info_w4_path, "\n"))
cat("Plot files created:\n")
cat(paste0("- ", bray_pdf, "\n"))
cat(paste0("- ", bray_png, "\n"))
cat(paste0("- ", aitchison_pdf, "\n"))
cat(paste0("- ", aitchison_png, "\n"))
cat(paste0("- ", bray_w4_pdf, "\n"))
cat(paste0("- ", bray_w4_png, "\n"))
cat(paste0("- ", aitchison_w4_pdf, "\n"))
cat(paste0("- ", aitchison_w4_png, "\n"))
cat(paste0("- ", bray_mean_pdf, "\n"))
cat(paste0("- ", bray_mean_png, "\n"))
cat(paste0("- ", core_accessory_bray_spaghetti_pdf, "\n"))
cat(paste0("- ", core_accessory_bray_spaghetti_png, "\n"))
cat(paste0("- ", core_accessory_bray_mean_pdf, "\n"))
cat(paste0("- ", core_accessory_bray_mean_png, "\n"))
cat(paste0("- ", aitchison_individual_pdf, "\n"))
cat(paste0("- ", aitchison_individual_png, "\n"))
cat(paste0("- ", aitchison_individual_dates_pdf, "\n"))
cat(paste0("- ", aitchison_individual_dates_png, "\n"))
cat(paste0("- ", aitchison_mean_pdf, "\n"))
cat(paste0("- ", aitchison_mean_png, "\n"))
cat(paste0("- ", p015_w16_contributors_pdf, "\n"))
cat(paste0("- ", p015_w16_contributors_png, "\n"))
cat(
  paste0(
    "Pre-travel samples included: yes; ",
    nrow(pretravel_data$baseline_df),
    " baseline rows were plotted at distance 0 using observed pre-travel dates, and ",
    length(unique(distance_df$patient_id)) - nrow(pretravel_data$baseline_df),
    " patients had no usable pre-travel metadata date.\n"
  )
)
cat(
  paste0(
    "w4 reference samples included: yes; ",
    nrow(w4_reference_data$baseline_df),
    " w4 rows were plotted at distance 0 using observed w4 dates, and ",
    length(unique(distance_w4_df$patient_id)) - nrow(w4_reference_data$baseline_df),
    " patients had no usable w4 metadata date.\n"
  )
)
