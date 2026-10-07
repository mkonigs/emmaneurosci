#' Run cluster analysis on cognitive domain profiles
#'
#' Performs k-means clustering on a set of cognitive domain variables using
#' either raw variables (`method = "regular"`) or a UMAP embedding
#' (`method = "umap"`), followed by group-wise comparison plots, individual-
#' level heterogeneity plots, and optionally cluster validity indices and
#' bootstrap stability analysis.
#'
#' @param data `data.frame` containing the id column and all clustering
#'   variables. Should already be filtered to the desired time-point / group
#'   (e.g. baseline only).
#' @param vars Character vector of column names to cluster on (e.g.
#'   `c("d_proc", "d_att_ctrl", "d_mem", "d_ver_wm", "d_vis_wm", "d_visuom")`).
#' @param id_var Name of the subject/id column. Default `"subj"`.
#' @param scaling `"unscaled"` (default) leaves scores as-is; `"scaled"`
#'   row-centres each participant's scores by subtracting the row mean.
#' @param method `"umap"` (default) embeds the data with UMAP before
#'   clustering; `"regular"` clusters directly on `vars`.
#' @param k Number of clusters (2–6). Default `4`.
#' @param seed Integer random seed for reproducibility. Default `123`.
#' @param colors Character vector of fill colours, one per cluster.
#'   Defaults to `c("skyblue", "palegreen", "orange", "tomato", "purple", "gold")`.
#' @param save_figures Logical; write per-cluster and combined PDFs to
#'   `figures_dir`. Default `FALSE`.
#' @param figures_dir Directory for saved figures. Created if absent.
#'   Default `"figures"`.
#' @param save_data Logical; write the clustered dataset as a `.csv2` file
#'   to `data_dir`. Default `FALSE`.
#' @param data_dir Directory for saved data. Created if absent.
#'   Default `"databases"`.
#' @param file_tag Prefix used in output filenames. Default `"cluster_analysis"`.
#' @param run_evaluation Logical; run [evaluate_clusters()] after clustering
#'   to compute validity indices and bootstrap stability. Default `TRUE`.
#' @param n_boot Number of bootstrap replicates for stability analysis.
#'   Default `100`. Set to `0` to skip stability (validity only).
#' @param nstart_boot Number of k-means random starts per bootstrap replicate.
#'   Default `25`.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{`data`}{`data.frame` with id column, all `vars`, and a `groups`
#'     column containing the cluster assignment (integer, 1–k).}
#'   \item{`cluster_model`}{The fitted `kmeans` object.}
#'   \item{`nbclust_plot`}{WSS elbow plot (`ggplot`) to help choose `k`.}
#'   \item{`dendrogram`}{Hierarchical clustering dendrogram (`ggplot`).}
#'   \item{`group_plots`}{Named list (`"group_1"`, `"group_2"`, …) of
#'     per-cluster bar plots comparing each cluster against all others.}
#'   \item{`combined_plot`}{A `gtable`/grob of all group plots arranged in a
#'     grid. Draw with `grid::grid.draw(res$combined_plot)`.}
#'   \item{`heterogeneity_plots`}{Named list of per-domain individual-level
#'     bar plots (one per variable in `vars`), coloured by cluster.}
#'   \item{`evaluation`}{`NULL` if `run_evaluation = FALSE`; otherwise the
#'     list returned by [evaluate_clusters()], containing `$validity` (a
#'     data frame of internal validity indices) and `$stability` (ARI,
#'     consensus matrix, and patient-level consensus scores).}
#'   \item{`stats`}{List with `p_values` and `cohen_d`, each a named list
#'     (one entry per cluster) of per-variable statistics.}
#' }
#'
#' @examples
#' \dontrun{
#' library(lubridate)
#'
#' dat <- read.csv2("et_domains_normed_peds_COGCOR_age_2026-05-12.csv",
#'                   stringsAsFactors = FALSE, sep = ";")
#' dat$test <- ymd(dat$test)
#' data_T1 <- dat[dat$group == 1, ]
#'
#' res <- run_cluster_analysis(
#'   data         = data_T1,
#'   vars         = c("d_proc", "d_att_ctrl", "d_mem",
#'                     "d_ver_wm", "d_vis_wm", "d_visuom"),
#'   scaling      = "unscaled",
#'   method       = "umap",
#'   k            = 4,
#'   save_figures = TRUE,
#'   save_data    = TRUE,
#'   n_boot       = 100
#' )
#'
#' res$data                               # clustered dataset
#' res$evaluation$validity                # validity indices table
#' res$evaluation$stability$mean_ari      # mean bootstrap ARI
#' res$evaluation$stability$patient_consensus  # patient-level consensus
#' grid::grid.draw(res$combined_plot)     # all group plots
#' }
#'
#' @importFrom dplyr select all_of
#' @importFrom ggplot2 ggplot aes geom_bar geom_errorbar scale_fill_manual
#'   ggtitle xlab ylab annotate theme element_text element_blank
#'   position_dodge guide_axis scale_x_discrete
#' @importFrom ggpubr theme_pubr
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid grid.draw
#' @importFrom reshape2 melt
#' @importFrom Rmisc summarySE
#' @importFrom psych cohen.d
#' @importFrom factoextra fviz_nbclust fviz_dend
#' @importFrom umap umap umap.defaults
#'
#' @export
run_cluster_analysis <- function(data,
                                  vars,
                                  id_var         = "subj",
                                  scaling        = c("unscaled", "scaled"),
                                  method         = c("umap", "regular"),
                                  k              = 4,
                                  seed           = 123,
                                  colors         = c("skyblue", "palegreen", "orange",
                                                     "tomato", "purple", "gold"),
                                  save_figures   = FALSE,
                                  figures_dir    = "figures",
                                  save_data      = FALSE,
                                  data_dir       = "databases",
                                  file_tag       = "cluster_analysis",
                                  run_evaluation = TRUE,
                                  n_boot         = 100,
                                  nstart_boot    = 25) {

  scaling <- match.arg(scaling)
  method  <- match.arg(method)

  # ---- input validation ----
  if (!all(vars %in% colnames(data))) {
    stop("Some `vars` are missing from `data`: ",
         paste(setdiff(vars, colnames(data)), collapse = ", "))
  }
  if (!id_var %in% colnames(data)) {
    stop("`id_var` ('", id_var, "') not found in `data`.")
  }
  if (k < 2 || k > 6) {
    stop("`k` must be between 2 and 6.")
  }

  subj     <- data[[id_var]]
  data_sel <- dplyr::select(data, dplyr::all_of(vars))

  # ---- optional within-subject scaling ----
  if (scaling == "scaled") {
    data_sel <- data_sel - rowMeans(data_sel)
  }

  # ============================ CLUSTERING ============================

  if (method == "regular") {

    set.seed(seed)
    nbclust_plot <- factoextra::fviz_nbclust(data_sel, kmeans, method = "wss")

    set.seed(seed)
    hc <- data_sel |>
      scale() |>
      dist(method = "manhattan") |>
      hclust(method = "ward.D2")

    dend_plot <- factoextra::fviz_dend(hc, k = k, cex = 0.5,
                                        color_labels_by_k = TRUE, rect = TRUE)

    set.seed(seed)
    clustering_model <- kmeans(as.matrix(data_sel), k,
                                nstart = 100, iter.max = 300)

  } else {

    set.seed(seed)
    umap_result    <- umap::umap(data_sel, config = umap::umap.defaults)
    umap_embedding <- umap_result$layout

    set.seed(seed)
    nbclust_plot <- factoextra::fviz_nbclust(umap_embedding, kmeans,
                                              method = "wss")

    set.seed(seed)
    hc <- umap_embedding |>
      scale() |>
      dist(method = "manhattan") |>
      hclust(method = "ward.D2")

    dend_plot <- factoextra::fviz_dend(hc, k = k, cex = 0.5,
                                        color_labels_by_k = TRUE, rect = TRUE)

    set.seed(seed)
    clustering_model <- kmeans(umap_embedding, k,
                                nstart = 100, iter.max = 300)
  }

  data_sel$groups <- clustering_model$cluster

  # ======================= GROUP-WISE COMPARISON =======================

  if (save_figures && !dir.exists(figures_dir)) {
    dir.create(figures_dir, recursive = TRUE)
  }

  group_plots <- list()
  stats_p     <- list()
  stats_d     <- list()

  for (z in seq_len(k)) {

    grp_data        <- data_sel
    grp_data$groups <- as.numeric(grp_data$groups == z)

    if (sum(grp_data$groups) <= 1) next

    data_long   <- reshape2::melt(grp_data, id = "groups")
    data_long_m <- Rmisc::summarySE(data_long, measurevar = "value",
                                     groupvars = c("variable", "groups"))

    p_vals <- sapply(vars, function(v) {
      round(t.test(grp_data[grp_data$groups == 0, v],
                    grp_data[grp_data$groups == 1, v],
                    alternative = "two.sided")$p.value, 3)
    })

    d_vals <- sapply(vars, function(v) {
      round(psych::cohen.d(grp_data[[v]], group = grp_data$groups)$cohen.d[2], 2)
    })

    stats_p[[paste0("group_", z)]] <- p_vals
    stats_d[[paste0("group_", z)]] <- d_vals

    p_coding <- ifelse(p_vals < .001, "***",
                 ifelse(p_vals < .01,  "**",
                  ifelse(p_vals < .05,  "*", "")))

    col_z <- colors[((z - 1L) %% length(colors)) + 1L]

    plot_z <- ggplot2::ggplot(data_long_m,
                               ggplot2::aes(x    = variable, y = value,
                                             fill = factor(groups))) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = value - se, ymax = value + se, width = 0.5),
        position = ggplot2::position_dodge(width = 0.9), colour = "grey") +
      ggplot2::geom_bar(stat = "identity", position = "dodge",
                         show.legend = FALSE) +
      ggplot2::scale_fill_manual(values = c("grey", col_z)) +
      ggplot2::ggtitle(
        paste0("Group ", z, " (n = ", sum(grp_data$groups == 1), ")")) +
      ggplot2::xlab("Domain") +
      ggplot2::ylab("Performance (z-score)") +
      ggplot2::annotate("text", x = seq_along(vars), y = -0.5,
                         label = paste0("d = ", d_vals), size = 2) +
      ggplot2::annotate("text", x = seq_along(vars), y = 0,
                         label = p_coding, size = 5) +
      ggpubr::theme_pubr() +
      ggplot2::theme(axis.title.x  = ggplot2::element_text(size = 10),
                      axis.title.y  = ggplot2::element_text(size = 10),
                      axis.text.x   = ggplot2::element_text(size = 8),
                      axis.text.y   = ggplot2::element_text(size = 9),
                      plot.title    = ggplot2::element_text(hjust = 0.5),
                      legend.title  = ggplot2::element_blank()) +
      ggplot2::scale_x_discrete(
        guide = ggplot2::guide_axis(n.dodge = 2))

    group_plots[[paste0("group_", z)]] <- plot_z

    if (save_figures) {
      grDevices::cairo_pdf(
        file.path(figures_dir, paste0("plot_group_", z, ".pdf")))
      print(plot_z)
      grDevices::dev.off()
    }
  }

  # ---- combined grid ----
  nrow_lookup   <- c(`2` = 2, `3` = 2, `4` = 2, `5` = 3, `6` = 3)
  combined_plot <- gridExtra::arrangeGrob(
    grobs = group_plots,
    nrow  = nrow_lookup[[as.character(k)]])

  if (save_figures) {
    grDevices::cairo_pdf(
      file.path(figures_dir,
                paste0("plot_", scaling, "_", method,
                        "_clusters_", k, "_all.pdf")))
    grid::grid.draw(combined_plot)
    grDevices::dev.off()
  }

  # ======================= FINAL DATASET =======================

  result_data <- data.frame(
    stats::setNames(list(subj), id_var),
    data_sel,
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )

  if (save_data) {
    if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
    utils::write.csv2(
      result_data,
      file.path(data_dir,
                paste0(file_tag, "_", scaling, "_clusters_", k, "_all.csv")),
      row.names = FALSE
    )
  }

  # ======================= HETEROGENEITY PLOTS =======================

  heterogeneity_plots <- plot_domain_heterogeneity(
    data         = result_data,
    vars         = vars,
    id_var       = id_var,
    colors       = stats::setNames(colors[seq_len(k)], as.character(seq_len(k))),
    save_figures = save_figures,
    figures_dir  = figures_dir,
    file_tag     = paste0(file_tag, "_", scaling, "_", method, "_clusters_", k)
  )

  # ======================= VALIDITY & STABILITY =======================

  evaluation <- NULL
  if (run_evaluation) {
    evaluation <- evaluate_clusters(
      cluster_result = list(data = result_data,
                             cluster_model = clustering_model),
      vars           = vars,
      scaling        = scaling,
      method         = method,
      id_var         = id_var,
      n_boot         = n_boot,
      nstart_boot    = nstart_boot,
      seed           = seed
    )
  }

  # =================================================================

  list(
    data                = result_data,
    cluster_model       = clustering_model,
    nbclust_plot        = nbclust_plot,
    dendrogram          = dend_plot,
    group_plots         = group_plots,
    combined_plot       = combined_plot,
    heterogeneity_plots = heterogeneity_plots,
    evaluation          = evaluation,
    stats               = list(p_values = stats_p, cohen_d = stats_d)
  )
}
