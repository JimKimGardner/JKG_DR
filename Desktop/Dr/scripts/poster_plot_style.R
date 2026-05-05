poster_text_color <- "#1f2933"
poster_muted_text_color <- "#5b6570"
poster_grid_major_color <- "#dfe4ea"
poster_grid_minor_color <- "#eef2f6"
poster_strip_fill <- "#fafbfd"
poster_strip_border <- "#d6dde6"

poster_element_line_or_blank <- function(enabled, color, linewidth) {
  if (isTRUE(enabled)) {
    ggplot2::element_line(color = color, linewidth = linewidth)
  } else {
    ggplot2::element_blank()
  }
}

poster_theme <- function(
  base_size = 14,
  legend_position = "right",
  x_text_angle = 0,
  x_text_hjust = NULL,
  major_x = TRUE,
  major_y = TRUE,
  minor_x = FALSE,
  minor_y = FALSE
) {
  if (is.null(x_text_hjust)) {
    x_text_hjust <- if (x_text_angle == 0) 0.5 else 1
  }

  x_text_vjust <- if (x_text_angle == 0) 0.5 else 1

  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.major.x = poster_element_line_or_blank(major_x, poster_grid_major_color, 0.45),
      panel.grid.major.y = poster_element_line_or_blank(major_y, poster_grid_major_color, 0.45),
      panel.grid.minor.x = poster_element_line_or_blank(minor_x, poster_grid_minor_color, 0.3),
      panel.grid.minor.y = poster_element_line_or_blank(minor_y, poster_grid_minor_color, 0.3),
      axis.title = ggplot2::element_text(face = "plain", color = poster_text_color),
      axis.text.x = ggplot2::element_text(
        angle = x_text_angle,
        hjust = x_text_hjust,
        vjust = x_text_vjust,
        color = poster_muted_text_color
      ),
      axis.text.y = ggplot2::element_text(color = poster_muted_text_color),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.title = ggplot2::element_text(face = "plain", size = ggplot2::rel(1.35), color = poster_text_color, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(1.02), color = poster_text_color, hjust = 0),
      plot.caption = ggplot2::element_text(size = ggplot2::rel(0.82), color = poster_muted_text_color, hjust = 1),
      legend.position = legend_position,
      legend.title = ggplot2::element_text(face = "plain", color = poster_text_color),
      legend.text = ggplot2::element_text(color = poster_text_color),
      strip.text = ggplot2::element_text(face = "plain", color = poster_text_color),
      strip.background = ggplot2::element_rect(fill = poster_strip_fill, color = poster_strip_border, linewidth = 0.4),
      plot.margin = ggplot2::margin(t = 10, r = 12, b = 10, l = 10)
    )
}

poster_placeholder_plot <- function(title, subtitle, body_text) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 1, y = 1, label = body_text, size = 5.6, color = poster_text_color) +
    ggplot2::annotate("text", x = 1, y = 0.84, label = subtitle, size = 4.2, color = poster_muted_text_color) +
    ggplot2::xlim(0.5, 1.5) +
    ggplot2::ylim(0.55, 1.1) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 15) +
    ggplot2::theme(
      plot.title.position = "plot",
      plot.title = ggplot2::element_text(face = "plain", size = ggplot2::rel(1.35), color = poster_text_color, hjust = 0),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(1.02), color = poster_text_color, hjust = 0)
    )
}

save_poster_plot <- function(plot, pdf_path, png_path, width, height, dpi = 320) {
  ggplot2::ggsave(pdf_path, plot = plot, width = width, height = height, units = "in", device = grDevices::pdf)
  ggplot2::ggsave(png_path, plot = plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
}
