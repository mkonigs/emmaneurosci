test_that("run_cluster_analysis returns expected structure", {

  set.seed(1)
  n <- 40
  fake_data <- data.frame(
    subj      = paste0("S", seq_len(n)),
    d_proc    = rnorm(n),
    d_att_ctrl = rnorm(n),
    d_mem     = rnorm(n),
    d_ver_wm  = rnorm(n),
    d_vis_wm  = rnorm(n),
    d_visuom  = rnorm(n)
  )

  vars <- c("d_proc", "d_att_ctrl", "d_mem", "d_ver_wm", "d_vis_wm", "d_visuom")

  res <- run_cluster_analysis(
    data   = fake_data,
    vars   = vars,
    method = "regular",
    k      = 3
  )

  # Return structure
  expect_type(res, "list")
  expect_named(res, c("data", "cluster_model", "nbclust_plot", "dendrogram",
                        "group_plots", "combined_plot", "stats"))

  # Data dimensions
  expect_equal(nrow(res$data), n)
  expect_true("groups" %in% colnames(res$data))
  expect_true(all(res$data$groups %in% 1:3))

  # Cluster model
  expect_s3_class(res$cluster_model, "kmeans")

  # Plots exist
  expect_s3_class(res$nbclust_plot, "ggplot")
  expect_length(res$group_plots, 3)
})

test_that("scaled method row-centres the data", {
  set.seed(2)
  n <- 30
  fake_data <- data.frame(
    subj    = paste0("S", seq_len(n)),
    v1 = rnorm(n, mean = 10),
    v2 = rnorm(n, mean = 10),
    v3 = rnorm(n, mean = 10)
  )

  res <- run_cluster_analysis(
    data    = fake_data,
    vars    = c("v1", "v2", "v3"),
    method  = "regular",
    scaling = "scaled",
    k       = 2
  )

  # Row means of the three domain columns should be ~0 (allowing float error)
  domain_cols <- res$data[, c("v1", "v2", "v3")]
  row_means <- rowMeans(domain_cols)
  expect_true(all(abs(row_means) < 1e-10))
})

test_that("input validation catches bad arguments", {
  fake <- data.frame(subj = 1:10, x = rnorm(10))

  expect_error(
    run_cluster_analysis(fake, vars = c("x", "missing_col")),
    "missing from"
  )
  expect_error(
    run_cluster_analysis(fake, vars = "x", id_var = "nope"),
    "id_var"
  )
  expect_error(
    run_cluster_analysis(fake, vars = "x", k = 7),
    "between 2 and 6"
  )
})
