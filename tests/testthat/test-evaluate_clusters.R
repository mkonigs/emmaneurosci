test_that("evaluate_clusters returns validity indices", {

  set.seed(1)
  n <- 60
  fake_data <- data.frame(
    subj       = paste0("S", seq_len(n)),
    d_proc     = rnorm(n),
    d_att_ctrl = rnorm(n),
    d_mem      = rnorm(n),
    stringsAsFactors = FALSE
  )
  vars <- c("d_proc", "d_att_ctrl", "d_mem")

  res <- run_cluster_analysis(
    data           = fake_data,
    vars           = vars,
    method         = "regular",
    k              = 3,
    run_evaluation = TRUE,
    n_boot         = 0       # skip bootstrap for speed
  )

  ev <- res$evaluation
  expect_type(ev, "list")
  expect_named(ev, c("validity", "stability"))

  # validity data.frame has the right columns
  expect_s3_class(ev$validity, "data.frame")
  expect_true(all(c("silhouette", "bss_tss", "calinski_harabasz",
                     "davies_bouldin", "cv_cluster_sizes") %in%
                     colnames(ev$validity)))

  # sensible ranges
  expect_true(ev$validity$silhouette >= -1 && ev$validity$silhouette <= 1)
  expect_true(ev$validity$bss_tss    >=  0 && ev$validity$bss_tss    <= 1)
  expect_true(ev$validity$calinski_harabasz > 0)
  expect_true(ev$validity$davies_bouldin    > 0)
  expect_true(ev$validity$cv_cluster_sizes  >= 0)

  # stability skipped
  expect_null(ev$stability)
})

test_that("bootstrap stability returns ARI and consensus", {

  set.seed(2)
  n <- 50
  fake_data <- data.frame(
    subj   = paste0("S", seq_len(n)),
    d_proc = c(rnorm(25, -1), rnorm(25, 1)),
    d_mem  = c(rnorm(25,  1), rnorm(25, -1)),
    stringsAsFactors = FALSE
  )

  res <- run_cluster_analysis(
    data           = fake_data,
    vars           = c("d_proc", "d_mem"),
    method         = "regular",
    k              = 2,
    run_evaluation = TRUE,
    n_boot         = 10,    # small for speed
    nstart_boot    = 5
  )

  stab <- res$evaluation$stability
  expect_length(stab$ari, 10)
  expect_true(is.numeric(stab$mean_ari))
  expect_true(stab$mean_ari >= -1 && stab$mean_ari <= 1)

  # consensus matrix is n × n
  expect_equal(dim(stab$consensus_matrix), c(n, n))

  # patient consensus one value per patient
  expect_length(stab$patient_consensus, n)
})
