# emmaneurosci <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)
<!-- badges: end -->

**emmaneurosci** is a private R package containing analysis tools developed for the
Emma Neurocience Group. The package is developed for the Emma Neuroscience Group.

---

## Installation

Because this is a **private** GitHub repository, you need a
[Personal Access Token (PAT)](https://github.com/settings/tokens) with at least
`repo` scope.

```r
# Install remotes if you don't have it
install.packages("remotes")

# Install cogcor from the private repo
remotes::install_github("YOUR_GITHUB_USERNAME/cogcor",
                         auth_token = "YOUR_GITHUB_PAT")
```

> **Tip:** store your PAT in your `.Renviron` as `GITHUB_PAT=...` and use
> `auth_token = Sys.getenv("GITHUB_PAT")` so it never lives in your scripts.

---

## Functions

### `run_cluster_analysis()`

Cluster cognitive domain profiles using either raw variables (`method = "regular"`)
or a UMAP embedding (`method = "umap"`), followed by k-means clustering.

**Returns** a list containing:

| Element | Description |
|---|---|
| `data` | Data frame with cluster assignments added as `groups` |
| `cluster_model` | The fitted `kmeans` object |
| `nbclust_plot` | WSS elbow plot for choosing k |
| `dendrogram` | Hierarchical clustering dendrogram |
| `group_plots` | Named list of per-cluster bar plots |
| `combined_plot` | Grid of all group plots (draw with `grid::grid.draw()`) |
| `stats` | Per-cluster t-test p-values and Cohen's d |

**Example:**

```r
library(emmaneurosci)
library(lubridate)

dat <- read.csv2("et_domains_normed_peds_COGCOR_age_2026-05-12.csv",
                  stringsAsFactors = FALSE, sep = ";")
dat$test <- ymd(dat$test)
data_T1 <- dat[dat$group == 1, ]

res <- run_cluster_analysis(
  data         = data_T1,
  vars         = c("d_proc", "d_att_ctrl", "d_mem",
                    "d_ver_wm", "d_vis_wm", "d_visuom"),
  scaling      = "unscaled",
  method       = "umap",
  k            = 4,
  save_figures = TRUE,
  save_data    = TRUE
)

res$data                           # clustered dataset
res$nbclust_plot                   # elbow plot
res$dendrogram                     # dendrogram
res$group_plots$group_1            # bar plot: cluster 1 vs rest
grid::grid.draw(res$combined_plot) # all clusters combined
```

---

## Development

```r
# Install development dependencies
install.packages(c("devtools", "roxygen2", "usethis"))

# Clone the repo, open cogcor.Rproj, then:
devtools::document()   # regenerate NAMESPACE + man/ from roxygen2 tags
devtools::check()      # run R CMD check
devtools::load_all()   # load the package locally without installing
```

---

## Adding new functions

1. Create `R/your_function.R` with proper `#' @export` roxygen2 tags.
2. Run `devtools::document()` to update `NAMESPACE` and `man/`.
3. Commit and push — collaborators can then re-install from GitHub.
