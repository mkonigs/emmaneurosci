#' Evaluate cluster validity and stability
#'
#' Computes internal validity indices for a cluster solution and assesses
#' stability via bootstrap resampling, following the approach described in
#' the COGCO-R analysis plan.
#'
#' **Validity indices** (computed on the reference solution):
#' \itemize{
#'   \item Average silhouette width
#'   \item Between-cluster / total variance ratio (BSS/TSS)
#'   \item Calinski–Harabasz index
#'   \item Davies–Bouldin index
#'   \item Coefficient of variation (CV) of cluster sizes
#' }
#'
#' **Stability** (bootstrap resampling, `n_boot` replicates):
#'
#' In each replicate, N patients are sampled with replacement; duplicates are
#' removed so each unique patient appears once. Scaling, UMAP embedding, and
#' k-means are re-estimated on this subsample. Agreement with the reference
#' partition is quantified using the adjusted Rand index (ARI), calculated only
#' for patients present in the replicate. A pair-wise co-clustering consensus
#' matrix is accumulated across replicates; patient-level consensus is the
#' average proportion of jointly-sampled replicates in which each pair of
#' patients from the same reference cluster was assigned to the same cluster.
#'
#' Designed to accept the list returned by [run_cluster_analysis()] directly.
#'
#' @param cluster_result List returned by [run_cluster_analysis()].
#' @param vars Character vector of domain column names used for clustering
#'   (must match those passed to [run_cluster_analysis()]).
#' @param scaling `"unscaled"` or `"scaled"` — must match the value used in
#'   [run_cluster_analysis()].
#' @param method `"umap"` or `"regular"` — must match the value used in
#'   [run_cluster_analysis()].
#' @param id_var Name of the subject/id column. Default `"subj"`.
#' @param n_boot Number of bootstrap replicates for stability analysis.
#'   Default `100`. Set to `0` to skip stability entirely.
#' @param nstart_boot Number of k-means random starts per replicate.
#'   Default `25` (vs 100 for the reference solution).
#' @param seed Random seed. Default `123`.
#'
#' @return A named list:
#' \describe{
#'   \item{`validity`}{`data.frame` with one row per candidate k (just the
#'     reference k here) and columns: `k`, `silhouette`, `bss_tss`,
#'     `calinski_harabasz`, `davies_bouldin`, `cv_cluster_sizes`.}
#'   \item{`stability`}{`NULL` if `n_boot = 0`, otherwise a list with:
#'     \describe{
#'       \item{`ari`}{Numeric vector of per-replicate ARI values.}
#'       \item{`mean_ari`}{Mean ARI across replicates.}
#'       \item{`sd_ari`}{SD of ARI across replicates.}
#'       \item{`consensus_matrix`}{Patient × patient matrix of consensus
#'         proportions (NA for pairs never co-sampled).}
#'       \item{`patient_consensus`}{Named numeric vector of patient-level
#'         consensus scores (average consensus with reference-cluster peers).}
#'     }}
#' }
#'
#' @examples
#' \dontrun{
#' res <- run_cluster_analysis(
#'   data   = data_T1,
#'   vars   = c("d_proc", "d_att_ctrl", "d_mem",
#'               "d_ver_wm", "d_vis_wm", "d_visuom"),
#'   method = "umap",
#'   k      = 4
#' )
#'
#' ev <- evaluate_clusters(
#'   cluster_result = res,
#'   vars           = c("d_proc", "d_att_ctrl", "d_mem",
#'                       "d_ver_wm", "d_vis_wm", "d_visuom"),
#'   scaling        = "unscaled",
#'   method         = "umap",
#'   n_boot         = 100
#' )
#'
#' ev$validity
#' ev$stability$mean_ari
#' ev$stability$patient_consensus
#' }
#'
#' @importFrom cluster silhouette
#' @importFrom mclust adjustedRandIndex
#' @importFrom umap umap umap.defaults
#'
#' @export
evaluate_clusters <- function(cluster_result,
                               vars,
                               scaling      = c("unscaled", "scaled"),
                               method       = c("umap", "regular"),
                               id_var       = "subj",
                               n_boot       = 100,
                               nstart_boot  = 25,
                               seed         = 123) {

  scaling <- match.arg(scaling)
  method  <- match.arg(method)

  data   <- cluster_result$data
  model  <- cluster_result$cluster_model
  groups <- data$groups
  k      <- length(unique(groups))
  n      <- nrow(data)
  ids    <- data[[id_var]]

  data_sel <- as.matrix(dplyr::select(data, dplyr::all_of(vars)))

  # ================================================================
  # 1. VALIDITY INDICES
  # ================================================================

  # -- silhouette (on Euclidean distance of raw vars) --
  dist_mat <- dist(data_sel, method = "euclidean")
  sil      <- cluster::silhouette(groups, dist_mat)
  avg_sil  <- mean(sil[, "sil_width"])

  # -- BSS / TSS (directly from kmeans object) --
  bss_tss <- model$betweenss / model$totss

  # -- Calinski-Harabasz: (BSS / (k-1)) / (WSS / (n-k)) --
  ch <- (model$betweenss / (k - 1)) / (model$tot.withinss / (n - k))

  # -- Davies-Bouldin (computed manually, no extra dependency) --
  centroids <- model$centers                        # k × p matrix
  # average within-cluster distance to centroid
  s_i <- vapply(seq_len(k), function(ci) {
    pts <- data_sel[groups == ci, , drop = FALSE]
    mean(sqrt(rowSums((pts - centroids[ci, ])^2)))
  }, numeric(1))
  # between-centroid distances
  centroid_dist <- as.matrix(dist(centroids, method = "euclidean"))
  diag(centroid_dist) <- NA
  db_values <- vapply(seq_len(k), function(ci) {
    max((s_i[ci] + s_i[-ci]) / centroid_dist[ci, -ci], na.rm = TRUE)
  }, numeric(1))
  db <- mean(db_values)

  # -- CV of cluster sizes --
  sizes <- as.numeric(table(groups))
  cv    <- stats::sd(sizes) / mean(sizes)

  validity <- data.frame(
    k                   = k,
    silhouette          = round(avg_sil, 4),
    bss_tss             = round(bss_tss, 4),
    calinski_harabasz   = round(ch, 2),
    davies_bouldin      = round(db, 4),
    cv_cluster_sizes    = round(cv, 4)
  )

  # ================================================================
  # 2. STABILITY (bootstrap resampling)
  # ================================================================

  if (n_boot == 0) {
    return(list(validity = validity, stability = NULL))
  }

  message("Running stability analysis (", n_boot, " bootstrap replicates)...")

  # accumulate co-clustering counts
  co_sampled   <- matrix(0L, nrow = n, ncol = n)   # times both were in replicate
  co_clustered <- matrix(0L, nrow = n, ncol = n)   # times both in same cluster
  ari_vec      <- numeric(n_boot)

  set.seed(seed)

  for (b in seq_len(n_boot)) {

    # sample with replacement, then deduplicate
    boot_idx    <- sample(seq_len(n), n, replace = TRUE)
    unique_idx  <- sort(unique(boot_idx))
    boot_data   <- data_sel[unique_idx, , drop = FALSE]

    # within-subject scaling if used
    if (scaling == "scaled") {
      boot_data <- boot_data - rowMeans(boot_data)
    }

    # embed and cluster
    tryCatch({
      if (method == "umap") {
        boot_embed  <- umap::umap(boot_data,
                                   config = umap::umap.defaults)$layout
        boot_model  <- kmeans(boot_embed, centers = k,
                               nstart = nstart_boot, iter.max = 300)
      } else {
        boot_model  <- kmeans(boot_data, centers = k,
                               nstart = nstart_boot, iter.max = 300)
      }

      boot_groups <- boot_model$cluster   # assignments for unique_idx patients

      # ARI vs reference (only for patients in this replicate)
      ref_sub     <- groups[unique_idx]
      ari_vec[b]  <- mclust::adjustedRandIndex(ref_sub, boot_groups)

      # update co-clustering matrices (only for pairs in this replicate)
      m <- length(unique_idx)
      for (i in seq_len(m)) {
        for (j in seq_len(m)) {
          ri <- unique_idx[i]
          rj <- unique_idx[j]
          co_sampled[ri, rj]   <- co_sampled[ri, rj]   + 1L
          if (boot_groups[i] == boot_groups[j]) {
            co_clustered[ri, rj] <- co_clustered[ri, rj] + 1L
          }
        }
      }

    }, error = function(e) {
      warning("Replicate ", b, " failed and was skipped: ", conditionMessage(e))
      ari_vec[b] <<- NA_real_
    })
  }

  # consensus matrix (proportion of co-sampled replicates in same cluster)
  consensus_mat <- matrix(NA_real_, nrow = n, ncol = n)
  nonzero       <- co_sampled > 0
  consensus_mat[nonzero] <- co_clustered[nonzero] / co_sampled[nonzero]
  rownames(consensus_mat) <- ids
  colnames(consensus_mat) <- ids

  # patient-level consensus: average with reference-cluster peers (excl. self)
  patient_consensus <- vapply(seq_len(n), function(i) {
    peers <- which(groups == groups[i] & seq_len(n) != i)
    if (length(peers) == 0) return(NA_real_)
    vals <- consensus_mat[i, peers]
    mean(vals, na.rm = TRUE)
  }, numeric(1))
  names(patient_consensus) <- ids

  ari_clean <- stats::na.omit(ari_vec)

  list(
    validity  = validity,
    stability = list(
      ari              = ari_vec,
      mean_ari         = round(mean(ari_clean), 4),
      sd_ari           = round(stats::sd(ari_clean), 4),
      consensus_matrix = consensus_mat,
      patient_consensus = patient_consensus
    )
  )
}
