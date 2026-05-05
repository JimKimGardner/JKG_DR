#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(forcats)
  library(ggplot2)
  library(Maaslin2)
  library(MicrobiomeStat)
})

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(project_root, "scripts", "poster_plot_style.R"))

esbl_status_path <- file.path(project_root, "results", "esbl_status", "esbl_status_metaphlan_matched.tsv")
sample_meta_path <- file.path(project_root, "data", "processed", "taxa", "TAXA_sample_metadata.tsv")
feature_meta_path <- file.path(project_root, "data", "processed", "taxa", "TAXA_feature_metadata.tsv")
tss_path <- file.path(project_root, "data", "processed", "taxa", "derived", "TAXA_sgb_tss_for_maaslin2.tsv")
relative_path <- file.path(project_root, "data", "processed", "taxa", "derived", "TAXA_relative_abundance.tsv")

results_dir <- file.path(project_root, "results", "esbl_non_esbl_daa")
plot_dir <- file.path(results_dir, "plots")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

daa_consensus_fdr <- 0.05
daa_suggestive_fdr <- 0.05
min_group_n <- 3L

make_short_label <- function(species_label, sgb_label, display_label, feature_id) {
  species_clean <- ifelse(is.na(species_label), "", species_label)
  sgb_clean <- ifelse(is.na(sgb_label), "", sgb_label)
  display_clean <- ifelse(is.na(display_label), "", display_label)

  ifelse(
    nzchar(species_clean) & nzchar(sgb_clean),
    paste0(species_clean, " / ", sgb_clean),
    ifelse(nzchar(display_clean), display_clean, feature_id)
  )
}

extract_genus_label <- function(feature_id) {
  genus_match <- str_extract(feature_id, "g__[^|]+")
  ifelse(is.na(genus_match) | genus_match == "g__", "g__Unclassified", genus_match)
}

extract_family_label <- function(feature_id) {
  family_match <- str_extract(feature_id, "f__[^|]+")
  ifelse(is.na(family_match) | family_match == "f__", "f__Unclassified", family_match)
}

make_aggregate_short_label <- function(rank_label, n_sgbs) {
  ifelse(is.na(n_sgbs), rank_label, paste0(rank_label, " (", n_sgbs, " SGBs)"))
}

direction_palette <- c(
  "Higher abundance in ESBL+" = "#b91c1c",
  "Lower abundance in ESBL+" = "#1d4ed8"
)

direction_palette_clean <- c(
  "Higher abundance in ESBL+" = "#a33b2e",
  "Lower abundance in ESBL+" = "#5a8f80"
)

phylum_palette_clean <- c(
  "Firmicutes" = "#7a8fc2",
  "Actinobacteria" = "#b989bc",
  "Bacteroidetes" = "#b38b57",
  "Proteobacteria" = "#e6d55c",
  "Other" = "#cbd5e1"
)

plot_placeholder <- function(title, subtitle, body_text, pdf_path, png_path) {
  plot <- poster_placeholder_plot(title = title, subtitle = subtitle, body_text = body_text)
  save_poster_plot(plot, pdf_path, png_path, width = 13.2, height = 8.4)
}

extract_phylum_label <- function(feature_id) {
  phylum_match <- str_match(feature_id, "p__([^|]+)")[, 2]
  phylum_raw <- ifelse(is.na(phylum_match) | !nzchar(phylum_match), "Other", phylum_match)
  dplyr::recode(
    phylum_raw,
    "Actinomycetota" = "Actinobacteria",
    "Bacteroidota" = "Bacteroidetes",
    "Pseudomonadota" = "Proteobacteria",
    .default = phylum_raw
  )
}

make_species_plot_label <- function(species_label, feature_short, feature_id) {
  label <- ifelse(
    !is.na(species_label) & nzchar(species_label),
    species_label,
    ifelse(!is.na(feature_short) & nzchar(feature_short), feature_short, feature_id)
  )
  label |>
    str_replace("^s__", "") |>
    str_replace(" /.*$", "") |>
    str_replace_all("_", " ")
}

run_maaslin2_group <- function(tss_sub, meta_sub, output_dir) {
  maaslin_dir <- file.path(output_dir, "maaslin2")
  dir.create(maaslin_dir, recursive = TRUE, showWarnings = FALSE)

  input_data <- as.data.frame(tss_sub)
  rownames(input_data) <- input_data$sample_id
  input_data$sample_id <- NULL

  input_meta <- as.data.frame(meta_sub)
  rownames(input_meta) <- input_meta$sample_id
  input_meta$sample_id <- NULL

  Maaslin2(
    input_data = input_data,
    input_metadata = input_meta,
    output = maaslin_dir,
    fixed_effects = c("esbl_group"),
    random_effects = NULL,
    normalization = "NONE",
    transform = "LOG",
    standardize = FALSE,
    min_prevalence = 0.10,
    min_abundance = 0,
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    max_significance = daa_suggestive_fdr
  )

  read_tsv(file.path(maaslin_dir, "all_results.tsv"), show_col_types = FALSE) |>
    filter(metadata == "esbl_group", value == "esbl_positive") |>
    transmute(
      feature_maaslin = feature,
      effect_maaslin = coef,
      stderr_maaslin = stderr,
      pval_maaslin = pval,
      qval_maaslin = qval,
      n_nonzero_maaslin = `N.not.0`
    )
}

run_linda_group <- function(relative_sub, meta_sub) {
  meta_df <- as.data.frame(meta_sub)
  rownames(meta_df) <- meta_df$sample_id
  meta_df$sample_id <- NULL

  relative_mat <- as.data.frame(relative_sub)
  rownames(relative_mat) <- relative_mat$feature_id
  relative_mat$feature_id <- NULL
  relative_mat <- relative_mat[, meta_sub$sample_id, drop = FALSE]

  fit <- linda(
    feature.dat = as.matrix(relative_mat),
    meta.dat = meta_df,
    formula = "~ esbl_group",
    feature.dat.type = "proportion",
    prev.filter = 0.10,
    mean.abund.filter = 0,
    max.abund.filter = 0,
    zero.handling = "imputation",
    adaptive = FALSE,
    is.winsor = TRUE,
    p.adj.method = "BH",
    alpha = daa_consensus_fdr,
    n.cores = 1,
    verbose = FALSE
  )

  coef_name <- "esbl_groupesbl_positive"
  if (!coef_name %in% names(fit$output)) {
    stop("Expected LinDA coefficient not found: ", coef_name, call. = FALSE)
  }

  fit$output[[coef_name]] |>
    rownames_to_column("feature_id") |>
    as_tibble() |>
    transmute(
      feature_id = feature_id,
      effect_linda = log2FoldChange,
      stderr_linda = lfcSE,
      pval_linda = pvalue,
      qval_linda = padj,
      base_mean_linda = baseMean
    )
}

build_consensus <- function(feature_tbl, maaslin_res, linda_res, level_label, sample_n, participant_n) {
  feature_tbl |>
    left_join(maaslin_res, by = "feature_maaslin") |>
    left_join(linda_res, by = "feature_id") |>
    mutate(
      analysis_label = "ESBL+ vs No ESBL",
      level_label = level_label,
      sample_n = sample_n,
      participant_n = participant_n,
      sign_maaslin = sign(effect_maaslin),
      sign_linda = sign(effect_linda),
      same_direction = !is.na(sign_maaslin) & !is.na(sign_linda) & sign_maaslin == sign_linda & sign_maaslin != 0,
      robust_consensus = same_direction & !is.na(qval_maaslin) & !is.na(qval_linda) &
        qval_maaslin < daa_consensus_fdr & qval_linda < daa_consensus_fdr,
      suggestive_consensus = same_direction & !robust_consensus &
        (
          (!is.na(qval_maaslin) & !is.na(pval_linda) & qval_maaslin < daa_consensus_fdr & pval_linda < 0.05) |
            (!is.na(qval_linda) & !is.na(pval_maaslin) & qval_linda < daa_consensus_fdr & pval_maaslin < 0.05) |
            (!is.na(qval_maaslin) & !is.na(qval_linda) & pmin(qval_maaslin, qval_linda) < daa_suggestive_fdr)
        ),
      consensus_class = case_when(
        robust_consensus ~ "Robust consensus",
        suggestive_consensus ~ "Suggestive consensus",
        same_direction ~ "Same direction only",
        TRUE ~ "No consensus"
      ),
      direction_consensus = case_when(
        same_direction & sign_maaslin > 0 ~ "Higher abundance in ESBL+",
        same_direction & sign_maaslin < 0 ~ "Lower abundance in ESBL+",
        TRUE ~ NA_character_
      ),
      consensus_score = case_when(
        same_direction ~ -log10(pmax(qval_maaslin, qval_linda, 1e-300)),
        TRUE ~ NA_real_
      )
    )
}

make_consensus_plot <- function(data, title, subtitle, pdf_path, png_path) {
  plot_df <- data |>
    filter(consensus_class %in% c("Robust consensus", "Suggestive consensus"), !is.na(consensus_score)) |>
    arrange(desc(consensus_score), desc(abs(effect_maaslin) + abs(effect_linda)), feature_short) |>
    slice_head(n = 20)

  if (nrow(plot_df) == 0) {
    plot_placeholder(
      title = title,
      subtitle = subtitle,
      body_text = "No suggestive or robust consensus hits under the current ESBL vs No-ESBL model.",
      pdf_path = pdf_path,
      png_path = png_path
    )
    return(invisible(NULL))
  }

  plot_df <- plot_df |>
    mutate(
      feature_short = fct_reorder(feature_short, consensus_score),
      consensus_class = factor(consensus_class, levels = c("Suggestive consensus", "Robust consensus"))
    )

  plot <- ggplot(plot_df, aes(x = consensus_score, y = feature_short)) +
    geom_segment(
      aes(x = 0, xend = consensus_score, y = feature_short, yend = feature_short, color = direction_consensus),
      linewidth = 0.8,
      alpha = 0.7
    ) +
    geom_point(
      aes(fill = direction_consensus, shape = consensus_class),
      size = 4.7,
      color = "#111827",
      stroke = 0.35
    ) +
    scale_color_manual(values = direction_palette, drop = FALSE, guide = "none") +
    scale_fill_manual(values = direction_palette, drop = FALSE) +
    scale_shape_manual(values = c("Suggestive consensus" = 22, "Robust consensus" = 23), drop = TRUE) +
    labs(
      title = title,
      subtitle = subtitle,
      x = expression(-log[10]("max FDR")),
      y = NULL,
      fill = "Direction",
      shape = "Consensus class",
      caption = paste0(
        "Exploratory model: all samples pooled across timepoints and participants. ",
        "Robust = FDR < ", formatC(daa_consensus_fdr, format = "f", digits = 2), " in both methods."
      )
    ) +
    poster_theme(base_size = 14, major_y = TRUE, minor_y = FALSE) +
    theme(legend.text = element_text(size = 10.5))

  save_poster_plot(plot, pdf_path, png_path, width = 13.2, height = 9.4)
}

make_sgb_consensus_plot_clean <- function(data, title, subtitle, pdf_path, png_path) {
  plot_df <- data |>
    filter(consensus_class %in% c("Robust consensus", "Suggestive consensus"), !is.na(consensus_score)) |>
    mutate(
      effect_display = (effect_maaslin + effect_linda) / 2,
      phylum_display = extract_phylum_label(feature_id),
      species_plot_label = make_species_plot_label(species_label, feature_short, feature_id)
    ) |>
    arrange(desc(consensus_score), desc(abs(effect_display)), species_plot_label) |>
    slice_head(n = 20)

  if (nrow(plot_df) == 0) {
    plot_placeholder(
      title = title,
      subtitle = subtitle,
      body_text = "No suggestive or robust consensus hits under the current ESBL vs No-ESBL model.",
      pdf_path = pdf_path,
      png_path = png_path
    )
    return(invisible(NULL))
  }

  plot_df <- plot_df |>
    arrange(desc(effect_display), species_plot_label) |>
    mutate(
      species_plot_label = factor(species_plot_label, levels = rev(species_plot_label))
    )

  max_effect <- max(abs(plot_df$effect_display), na.rm = TRUE)
  max_effect <- max(max_effect, 0.04)
  tile_x <- -max_effect * 1.26

  plot <- ggplot(plot_df, aes(x = effect_display, y = species_plot_label)) +
    geom_vline(xintercept = 0, linewidth = 0.45, color = "#9ca3af") +
    geom_col(
      aes(fill = direction_consensus),
      width = 0.72,
      color = NA
    ) +
    geom_point(
      aes(x = tile_x, color = phylum_display),
      shape = 15,
      size = 3.4,
      stroke = 0,
      show.legend = TRUE
    ) +
    annotate(
      "text",
      x = -max_effect * 0.62,
      y = nrow(plot_df) + 1.1,
      label = "Lower in ESBL+",
      color = direction_palette_clean[["Lower abundance in ESBL+"]],
      fontface = "bold",
      size = 3.6
    ) +
    annotate(
      "text",
      x = max_effect * 0.62,
      y = nrow(plot_df) + 1.1,
      label = "Higher in ESBL+",
      color = direction_palette_clean[["Higher abundance in ESBL+"]],
      fontface = "bold",
      size = 3.6
    ) +
    scale_fill_manual(values = direction_palette_clean, drop = FALSE) +
    scale_color_manual(values = phylum_palette_clean, drop = FALSE) +
    scale_x_continuous(
      labels = scales::number_format(accuracy = 0.01),
      expand = expansion(mult = c(0.02, 0.06))
    ) +
    coord_cartesian(
      xlim = c(tile_x - max_effect * 0.14, max_effect * 1.12),
      clip = "off"
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Average signed effect estimate (MaAsLin2 / LinDA)",
      y = NULL,
      fill = "Direction",
      color = "Phylum",
      caption = paste0(
        "Top 20 robust/suggestive consensus SGB hits. Exploratory pooled model; robust = FDR < ",
        formatC(daa_consensus_fdr, format = "f", digits = 2),
        " in both methods."
      )
    ) +
    poster_theme(base_size = 13.5, major_y = FALSE, minor_y = FALSE) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(size = 9.8, color = "#2b2b2b"),
      axis.text.x = element_text(size = 10.5, color = "#4b5563"),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = element_text(size = 10.5),
      legend.text = element_text(size = 9.8),
      plot.margin = margin(t = 14, r = 18, b = 8, l = 18)
    )

  save_poster_plot(plot, pdf_path, png_path, width = 13.2, height = 9.8)
}

make_same_direction_plot <- function(data, title, subtitle, pdf_path, png_path) {
  plot_df <- data |>
    filter(same_direction, !is.na(consensus_score), !is.na(direction_consensus)) |>
    arrange(desc(consensus_score), desc(abs(effect_maaslin) + abs(effect_linda)), feature_short) |>
    slice_head(n = 20)

  if (nrow(plot_df) == 0) {
    plot_placeholder(
      title = title,
      subtitle = subtitle,
      body_text = "No same-direction hits under the current ESBL vs No-ESBL model.",
      pdf_path = pdf_path,
      png_path = png_path
    )
    return(invisible(NULL))
  }

  plot_df <- plot_df |>
    mutate(feature_short = fct_reorder(feature_short, consensus_score))

  plot <- ggplot(plot_df, aes(x = consensus_score, y = feature_short)) +
    geom_segment(
      aes(x = 0, xend = consensus_score, y = feature_short, yend = feature_short, color = direction_consensus),
      linewidth = 0.8,
      alpha = 0.7
    ) +
    geom_point(
      aes(fill = direction_consensus),
      shape = 21,
      size = 4.7,
      color = "#111827",
      stroke = 0.35
    ) +
    scale_color_manual(values = direction_palette, drop = FALSE, guide = "none") +
    scale_fill_manual(values = direction_palette, drop = FALSE) +
    labs(
      title = title,
      subtitle = subtitle,
      x = expression(-log[10]("max FDR")),
      y = NULL,
      fill = "Direction",
      caption = "Exploratory layer: same direction in MaAsLin2 and LinDA, without requiring suggestive or robust consensus."
    ) +
    poster_theme(base_size = 14, major_y = TRUE, minor_y = FALSE) +
    theme(legend.text = element_text(size = 10.5))

  save_poster_plot(plot, pdf_path, png_path, width = 13.2, height = 9.4)
}

esbl_status <- read_tsv(esbl_status_path, show_col_types = FALSE) |>
  filter(esbl_status %in% c("esbl_positive", "no_esbl")) |>
  distinct(sample_id, participant_id, sequencing_timepoint, week, esbl_status) |>
  mutate(
    esbl_group = factor(esbl_status, levels = c("no_esbl", "esbl_positive"))
  )

overall_group_timepoint_summary <- esbl_status |>
  count(sequencing_timepoint, week, esbl_status, name = "n_samples") |>
  tidyr::pivot_wider(names_from = esbl_status, values_from = n_samples, values_fill = 0) |>
  arrange(week)

eligible_timepoints <- overall_group_timepoint_summary |>
  filter(esbl_positive >= min_group_n, no_esbl >= min_group_n) |>
  pull(sequencing_timepoint)

eligible_samples <- esbl_status |>
  filter(sequencing_timepoint %in% eligible_timepoints) |>
  arrange(sample_id)

group_summary <- eligible_samples |>
  group_by(esbl_status) |>
  summarise(
    n_samples = n(),
    n_participants = n_distinct(participant_id),
    .groups = "drop"
  )

group_timepoint_summary <- eligible_samples |>
  count(sequencing_timepoint, week, esbl_status, name = "n_samples") |>
  tidyr::pivot_wider(names_from = esbl_status, values_from = n_samples, values_fill = 0) |>
  arrange(week)

write_tsv(group_summary, file.path(results_dir, "esbl_group_summary.tsv"), na = "")
write_tsv(group_timepoint_summary, file.path(results_dir, "esbl_group_timepoint_summary.tsv"), na = "")
write_tsv(eligible_samples, file.path(results_dir, "esbl_eligible_samples.tsv"), na = "")

if (nrow(eligible_samples) == 0) {
  stop("No eligible ESBL / No-ESBL samples remained after minimum group-size filtering.", call. = FALSE)
}

feature_meta <- read_tsv(feature_meta_path, show_col_types = FALSE)
tss_df <- read_tsv(tss_path, show_col_types = FALSE)
relative_wide <- read_tsv(relative_path, show_col_types = FALSE)

selected_sample_ids <- eligible_samples$sample_id
sample_n <- length(selected_sample_ids)
participant_n <- n_distinct(eligible_samples$participant_id)

meta_sub <- eligible_samples |>
  select(sample_id, participant_id, sequencing_timepoint, week, esbl_group)

sgb_meta <- feature_meta |>
  filter(terminal_rank == "SGB") |>
  mutate(
    feature_id = as.character(feature_id),
    feature_short = make_short_label(species_label, sgb_label, display_label, feature_id),
    feature_maaslin = make.names(feature_id)
  ) |>
  select(feature_id, feature_maaslin, feature_short, species_label, sgb_label, prevalence, n_present_samples)

common_sgb_ids <- intersect(sgb_meta$feature_id, setdiff(names(tss_df), "sample_id"))
sgb_meta <- sgb_meta |>
  filter(feature_id %in% common_sgb_ids)

sgb_tss <- tss_df |>
  select(sample_id, all_of(sgb_meta$feature_id)) |>
  filter(sample_id %in% selected_sample_ids) |>
  mutate(sample_id = factor(sample_id, levels = selected_sample_ids)) |>
  arrange(sample_id) |>
  mutate(sample_id = as.character(sample_id))

sgb_relative <- relative_wide |>
  filter(feature_id %in% sgb_meta$feature_id) |>
  select(feature_id, all_of(selected_sample_ids))

genus_relative <- relative_wide |>
  mutate(
    feature_id = as.character(feature_id),
    genus_label = extract_genus_label(feature_id)
  ) |>
  select(-feature_id) |>
  group_by(genus_label) |>
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

genus_meta <- relative_wide |>
  transmute(
    feature_id = as.character(feature_id),
    genus_label = extract_genus_label(feature_id)
  ) |>
  count(genus_label, name = "n_sgbs") |>
  mutate(
    feature_id = genus_label,
    feature_maaslin = make.names(feature_id),
    feature_short = make_aggregate_short_label(genus_label, n_sgbs)
  ) |>
  select(feature_id, feature_maaslin, feature_short, genus_label, n_sgbs)

genus_tss <- genus_relative |>
  rename(feature_id = genus_label) |>
  pivot_longer(cols = -feature_id, names_to = "sample_id", values_to = "value") |>
  filter(sample_id %in% selected_sample_ids) |>
  pivot_wider(names_from = feature_id, values_from = value) |>
  mutate(sample_id = factor(sample_id, levels = selected_sample_ids)) |>
  arrange(sample_id) |>
  mutate(sample_id = as.character(sample_id))

genus_relative <- genus_relative |>
  rename(feature_id = genus_label) |>
  select(feature_id, all_of(selected_sample_ids))

family_relative <- relative_wide |>
  mutate(
    feature_id = as.character(feature_id),
    family_label = extract_family_label(feature_id)
  ) |>
  select(-feature_id) |>
  group_by(family_label) |>
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

family_meta <- relative_wide |>
  transmute(
    feature_id = as.character(feature_id),
    family_label = extract_family_label(feature_id)
  ) |>
  count(family_label, name = "n_sgbs") |>
  mutate(
    feature_id = family_label,
    feature_maaslin = make.names(feature_id),
    feature_short = make_aggregate_short_label(family_label, n_sgbs)
  ) |>
  select(feature_id, feature_maaslin, feature_short, family_label, n_sgbs)

family_tss <- family_relative |>
  rename(feature_id = family_label) |>
  pivot_longer(cols = -feature_id, names_to = "sample_id", values_to = "value") |>
  filter(sample_id %in% selected_sample_ids) |>
  pivot_wider(names_from = feature_id, values_from = value) |>
  mutate(sample_id = factor(sample_id, levels = selected_sample_ids)) |>
  arrange(sample_id) |>
  mutate(sample_id = as.character(sample_id))

family_relative <- family_relative |>
  rename(feature_id = family_label) |>
  select(feature_id, all_of(selected_sample_ids))

run_level <- function(level_name, level_label, feature_tbl, tss_sub, relative_sub) {
  level_dir <- file.path(results_dir, level_name)
  dir.create(level_dir, recursive = TRUE, showWarnings = FALSE)

  maaslin_res <- run_maaslin2_group(tss_sub, meta_sub, level_dir)
  linda_res <- run_linda_group(relative_sub, meta_sub)
  consensus_res <- build_consensus(
    feature_tbl = feature_tbl,
    maaslin_res = maaslin_res,
    linda_res = linda_res,
    level_label = level_label,
    sample_n = sample_n,
    participant_n = participant_n
  )

  write_tsv(consensus_res, file.path(level_dir, paste0("esbl_non_esbl_", level_name, "_consensus.tsv")), na = "")

  consensus_pdf <- file.path(plot_dir, paste0("esbl_non_esbl_", level_name, "_consensus.pdf"))
  consensus_png <- file.path(plot_dir, paste0("esbl_non_esbl_", level_name, "_consensus.png"))
  exploratory_pdf <- file.path(plot_dir, paste0("esbl_non_esbl_", level_name, "_same_direction.pdf"))
  exploratory_png <- file.path(plot_dir, paste0("esbl_non_esbl_", level_name, "_same_direction.png"))

  subtitle_text <- paste0(
    "Exploratory unpaired ESBL+ vs No-ESBL comparison pooled across eligible timepoints; ",
    sample_n, " samples from ", participant_n, " participants."
  )

  if (identical(level_name, "sgb")) {
    make_sgb_consensus_plot_clean(
      data = consensus_res,
      title = paste0(level_label, ": ESBL+ vs No ESBL consensus DAA"),
      subtitle = subtitle_text,
      pdf_path = consensus_pdf,
      png_path = consensus_png
    )
  } else {
    make_consensus_plot(
      data = consensus_res,
      title = paste0(level_label, ": ESBL+ vs No ESBL consensus DAA"),
      subtitle = subtitle_text,
      pdf_path = consensus_pdf,
      png_path = consensus_png
    )
  }

  make_same_direction_plot(
    data = consensus_res,
    title = paste0(level_label, ": ESBL+ vs No ESBL same-direction signals"),
    subtitle = subtitle_text,
    pdf_path = exploratory_pdf,
    png_path = exploratory_png
  )

  list(
    consensus = consensus_res,
    consensus_pdf = consensus_pdf,
    consensus_png = consensus_png,
    exploratory_pdf = exploratory_pdf,
    exploratory_png = exploratory_png
  )
}

sgb_res <- run_level("sgb", "SGB-level DAA", sgb_meta, sgb_tss, sgb_relative)
genus_res <- run_level("genus", "Genus-level DAA", genus_meta, genus_tss, genus_relative)
family_res <- run_level("family", "Family-level DAA", family_meta, family_tss, family_relative)

write_tsv(sgb_res$consensus, file.path(results_dir, "esbl_non_esbl_sgb_consensus_all.tsv"), na = "")
write_tsv(genus_res$consensus, file.path(results_dir, "esbl_non_esbl_genus_consensus_all.tsv"), na = "")
write_tsv(family_res$consensus, file.path(results_dir, "esbl_non_esbl_family_consensus_all.tsv"), na = "")

build_info_lines <- c(
  "ESBL+ vs No-ESBL differential abundance analysis",
  "Design: unpaired, time-independent, participant-independent pooled comparison",
  paste0("Minimum samples per ESBL group within timepoint for inclusion: ", min_group_n),
  paste0("Eligible timepoints: ", paste(eligible_timepoints, collapse = ", ")),
  paste0("Included samples: ", sample_n),
  paste0("Included participants: ", participant_n),
  paste0("DAA consensus FDR threshold: ", daa_consensus_fdr),
  paste0("DAA suggestive FDR threshold: ", daa_suggestive_fdr)
)
writeLines(build_info_lines, con = file.path(results_dir, "esbl_non_esbl_build_info.txt"))

cat("Created ESBL+ vs No-ESBL DAA outputs in:\n")
cat(paste0("- ", results_dir, "\n"))
