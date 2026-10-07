#' Plot within-cluster heterogeneity per cognitive domain
#'
#' For each domain variable, produces a bar chart showing every individual's
#' z-score ordered from lowest to highest, with bars coloured by their cluster
#' assignment. This makes within- and between-cluster heterogeneity immediately
#' visible across domains.
#'
#' @param data `data.frame` containing the domain columns, a `groups` column
#'   (integer cluster assignments), and optionally an id column. Typically
#'   `res$data` from [run_cluster_analysis()].
#' @param vars Character vector of domain column names to plot. If `NULL`
#'   (default), all numeric columns except `groups` and `id_var` are used.
#' @param id_var Name of the subject/id column. Default `"subj"`.
#' @param colors Named or unnamed character vector of colours, one per cluster.
#' @param x_title Label for the x-axis. Default `"Individual"`.
#' @param y_title Label for the y-axis. Default `"Score"`.
#' @param save_figures Logical; write one PDF per domain to `figures_dir`.
#'   Default `FALSE`.
#' @param figures_dir Directory for saved figures. Created if absent.
#'   Default `"figures"`.
#' @param file_tag Prefix used in output filenames. Default `"heterogeneity"`.
#'
#' @return A named list of `ggplot` objects, one per domain.
#'
#' @examples
#' \dontrun{
#' res <- run_cluster_analysis(data = data_T1, vars = domain_vars, k = 4)
#' plots <- plot_domain_heterogeneity(res$data)
#' plots$d_mem
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_hline scale_fill_manual
#'   ggtitle xlab ylab theme element_text element_blank element_line
#' @importFrom ggpubr theme_pubr
#'
#' @export
plot_domain_heterogeneity <- function(data,
                                       vars         = NULL,
                                       id_var       = "subj",
                                       colors       = c("1" = "skyblue",
                                                        "2" = "palegreen",
                                                        "3" = "orange",
                                                        "4" = "tomato",
                                                        "5" = "purple",
                                                        "6" = "gold"),
                                       x_title      = "Individual",
                                       y_title      = "Score",
                                       save_figures = FALSE,
                                       figures_dir  = "figures",
                                       file_tag     = "heterogeneity") {

  if (!"groups" %in% colnames(data)) {
    stop("`data` must contain a `groups` column. ",
         "Use the `$data` output of `run_cluster_analysis()`.")
  }

  if (is.null(vars)) {
    exclude <- c(id_var, "groups")
    vars <- colnames(data)[sapply(data, is.numeric) &
                             !colnames(data) %in% exclude]
    if (length(vars) == 0) stop("No numeric domain columns found in `data`.")
  }

  if (!all(vars %in% colnames(data))) {
    stop("Some `vars` not found in `data`: ",
         paste(setdiff(vars, colnames(data)), collapse = ", "))
  }

  if (save_figures && !dir.exists(figures_dir)) {
    dir.create(figures_dir, recursive = TRUE)
  }

  grp_levels <- as.character(sort(unique(data$groups)))
  col_values <- colors[seq_along(grp_levels)]
  names(col_values) <- grp_levels

  plots <- lapply(vars, function(domain) {

    plot_data <- data.frame(score  = data[[domain]],
                             groups = factor(data$groups))
    plot_data <- plot_data[order(plot_data$score), ]
    plot_data$individual <- seq_len(nrow(plot_data))

    p <- ggplot2::ggplot(plot_data,
                          ggplot2::aes(x = individual, y = score,
                                        fill = groups)) +
      ggplot2::geom_bar(stat = "identity", width = 1) +
      ggplot2::geom_hline(yintercept = 0, colour = "black",
                           linewidth = 0.4) +
      ggplot2::scale_fill_manual(values = col_values, name = "Cluster",
                                  drop = FALSE) +
      ggplot2::ggtitle(domain) +
      ggplot2::xlab(x_title) +
      ggplot2::ylab(y_title) +
      ggpubr::theme_pubr() +
      ggplot2::theme(
        plot.title     = ggplot2::element_text(hjust = 0.5, size = 12,
                                                face = "bold"),
        axis.text.x    = ggplot2::element_blank(),
        axis.ticks.x   = ggplot2::element_blank(),
        axis.title.x   = ggplot2::element_text(size = 10),
        axis.title.y   = ggplot2::element_text(size = 10),
        legend.position = "right"
      )

    if (save_figures) {
      grDevices::cairo_pdf(
        file.path(figures_dir, paste0(file_tag, "_", domain, ".pdf")))
      print(p)
      grDevices::dev.off()
    }
    p
  })

  names(plots) <- vars
  plots
}
