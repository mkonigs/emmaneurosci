test_that("plot_domain_heterogeneity returns one ggplot per domain", {

  set.seed(1)
  n <- 40
  fake_data <- data.frame(
    subj      = paste0("S", seq_len(n)),
    d_proc    = rnorm(n),
    d_mem     = rnorm(n),
    d_vis_wm  = rnorm(n),
    groups    = sample(1:4, n, replace = TRUE)
  )

  vars  <- c("d_proc", "d_mem", "d_vis_wm")
  plots <- plot_domain_heterogeneity(fake_data, vars = vars)

  expect_type(plots, "list")
  expect_named(plots, vars)
  expect_s3_class(plots$d_proc, "ggplot")
  expect_s3_class(plots$d_mem,  "ggplot")
})

test_that("plot_domain_heterogeneity errors without groups column", {
  fake <- data.frame(subj = 1:10, d_proc = rnorm(10))
  expect_error(plot_domain_heterogeneity(fake), "`groups` column")
})

test_that("auto-detection of vars skips id and groups columns", {
  set.seed(2)
  n <- 20
  fake_data <- data.frame(
    subj   = paste0("S", seq_len(n)),
    d_proc = rnorm(n),
    d_mem  = rnorm(n),
    groups = sample(1:3, n, replace = TRUE)
  )

  plots <- plot_domain_heterogeneity(fake_data)   # no vars supplied
  expect_named(plots, c("d_proc", "d_mem"))       # subj and groups excluded
})
