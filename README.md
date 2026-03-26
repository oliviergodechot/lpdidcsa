# lpdidcsa

**Local Projections DiD and Callaway & Sant'Anna Event Studies in R.**

## License and Credits

- The code in this package is licensed under the [MIT License](LICENSE).
- The `mpdta` dataset is sourced from the [`did` package](https://cloud.r-project.org/web/packages/did/index.html) and is licensed under [GPL-2](LICENSE.did).

## Description

`lpdidcsa` provides two functions that together cover the full workflow:

1.  **`lpdidcsa_data()`** — reshapes a panel into the format required for estimation
2.  **`lpdidcsa()`** — estimates event study coefficients and produces a plot

The package implements the LP-DiD estimator of Dube, Girardi, Jordà & Taylor (2025) and a local-projection version of the Callaway & Sant'Anna (2021) estimator, with support for absorbing and non-absorbing treatments. It enables to run LP-DiD and reweighted LP-DiD with or without control variables, adjusted LP-DiD with control variables, inverse probability weighting LP-DiD, LP-CSA with or without control variables and inverse probability weighting LP-CSA.

------------------------------------------------------------------------

## Installation

``` r
# Install from GitHub (requires remotes)
remotes::install_github("oliviergodechot/lpdidcsa")
```

------------------------------------------------------------------------

## Workflow

```         
Raw panel data
      │
      ▼
lpdidcsa_data()   →   Wide or long horizon dataset
      │
      ▼
lpdidcsa()        →   Estimates + event study plot
```

------------------------------------------------------------------------

## Step 1 — `lpdidcsa_data()`

Prepares the panel by computing treatment timing variables and creating leads/lags of the outcome at each horizon.

### Arguments

| Argument | Type | Default | Description |
|------------------|------------------|------------------|------------------|
| `data` | `data.frame` / `data.table` | — | Input panel |
| `unit` | `character` | — | Unit identifier column |
| `time` | `character` | `"year"` | Time period column (integer-valued) |
| `dependent` | `character` | — | Outcome variable column |
| `treat` | `character` | `"treat"` | Binary treatment indicator (0/1) |
| `absorbing` | `logical` | `TRUE` | Whether treatment is absorbing |
| `nb_pre` | `integer` / `NULL` | `NULL` | Pre-treatment periods (all if `NULL`) |
| `nb_post` | `integer` / `NULL` | `NULL` | Post-treatment periods (all if `NULL`) |
| `type_horizon` | `character` | `"wide"` | `"wide"` or `"long"` output format |
| `h_variables` | `character` / `NULL` | `NULL` | Controls carried at every horizon |
| `p_variables` | named `list` | `list()` | Horizon-specific controls |

### `p_variables` naming convention

List names follow `"t0"`, `"t1"`, `"tm1"`, `"tm2"`, … where `"tm"` stands for *t minus*:

``` r
p_variables <- list(
  tm2 = c("x1", "x2"),   # controls measured at t-2
  tm1 = c("x1", "x2"),   # controls measured at t-1
  t0  = c("x1")          # control measured at t
)
```

### Output columns added

| Column | Description |
|------------------------------------|------------------------------------|
| `first.treat` | First treatment period per unit (`NA` for never-treated). Absorbing only. |
| `last.treat` | Most recent treatment period as of *t* |
| `next.treat` | Next upcoming treatment period |
| `dtreat` | Treatment change: *D_it − D_i,t−1* |
| `<dep>_tm2`, `<dep>_t0`, … | Outcome at each horizon (wide format) |
| `horizon` | Horizon integer *h = t′ − t* (long format) |


### Example

``` r
library(data.table)
data(mpdta_r, package = "lpdidcsa")

p_variables <- list(
  tm2 = c("lpop", "lpop2"),
  tm1 = c("lpop", "lpop2")
)

# Wide format (one column per horizon)
df_wide <- lpdidcsa_data(
  data        = mpdta_r,
  unit        = "countyreal",
  time        = "year",
  dependent   = "lemp",
  treat       = "treat",
  absorbing   = TRUE,
  nb_pre      = 3,
  nb_post     = 3,
  p_variables = p_variables
)

# Long format (one row per unit × reference-period × horizon)
df_long <- lpdidcsa_data(
  data         = mpdta_r,
  unit         = "countyreal",
  time         = "year",
  dependent    = "lemp",
  treat        = "treat",
  absorbing    = TRUE,
  nb_pre       = 3,
  nb_post      = 3,
  type_horizon = "long",
  p_variables  = p_variables
)
```

------------------------------------------------------------------------

## Step 2 — `lpdidcsa()`

Estimates event study coefficients horizon by horizon (or in a single stacked regression) and returns a `ggplot2` plot.

### Arguments

| Argument | Type | Default | Description |
|------------------|------------------|------------------|------------------|
| `data` | `data.frame` / `data.table` | — | Output of `lpdidcsa_data()` |
| `unit` | `character` | — | Unit identifier column |
| `time` | `character` | `"year"` | Time period column |
| `dependent` | `character` | — | Outcome variable column |
| `dtreat` | `character` | `"dtreat"` | Treatment change indicator column |
| `nb_pre` | `integer` / `NULL` | `NULL` | Pre-treatment horizons (all if `NULL`) |
| `nb_post` | `integer` / `NULL` | `NULL` | Post-treatment horizons (all if `NULL`) |
| `controls` | `character` / `NULL` | `NULL` | Control variables (fixest-style strings accepted) |
| `FE` | `character` / `NULL` | `NULL` | Additional fixed effects |
| `cluster` | `character` | `unit` | Clustering column |
| `weight` | `character` / `NULL` | `NULL` | Sampling weight column |
| `meth` | `character` | `"lpdid_rw"` | Estimator (see table below) |
| `absorbing` | `logical` | `TRUE` | Absorbing treatment assumption |
| `nonabs_reentry` | `integer` / `NULL` | `NULL` | Min. periods since last treatment for re-entry (non-absorbing only) |
| `type_horizon` | `character` | `"wide"` | Format of the input data: `"wide"` or `"long"` |
| `horizon` | `character` | `"horizon"` | Horizon column name (long format only) |
| `one_reg` | `logical` | `FALSE` | Run all horizons in one stacked regression (long format only) |

### Estimators

| `meth` | Description | ATT type |
|------------------------|------------------------|------------------------|
| `"lpdid"` | LP-DiD with control variables | Variance-weighted |
| `"lpdid_rw"` | Reweighted LP-DiD with control variables | Equally-weighted |
| `"lpdid_adj"` | LP-DiD with adjusted regression (`avg_comparisons`) | Equally-weighted |
| `"lpdid_ipw"` | LP-DiD with inverse probability weighting | Equally-weighted |
| `"lpcsa_ipw"` | LP-CSA with inverse probability weighting | Equally-weighted |
| `"lpcsa"` | LP-CSA with control variables | Equally-weighted |

### Return value

A named list with four elements:

| Element | Description |
|------------------------------------|------------------------------------|
| `est` | Main results: `data.table` with columns `h`, `variable`, `estimate`, `se`, `T`, `pvalue`, `nb_obs`, `formula` |
| `est_det` | Cohort-level estimates (CSA and adjusted methods only) |
| `ps` | Propensity score estimates (IPW methods only) |
| `plot` | `ggplot2` event study plot with 95% confidence intervals |

### Example

``` r
library(data.table)
library(fixest)
library(ggplot2)
data(mpdta_r, package = "lpdidcsa")

# Step 1: prepare data
df <- lpdidcsa_data(
  data      = mpdta_r,
  unit      = "countyreal",
  time      = "year",
  dependent = "lemp",
  treat     = "treat"
)

# Step 2: estimate
res <- lpdidcsa(
  data      = df,
  unit      = "countyreal",
  time      = "year",
  dependent = "lemp",
  meth      = "lpdid_rw",
  nb_pre    = 3,
  nb_post   = 3
)

res$plot   # event study plot
res$est    # coefficient table
```

------------------------------------------------------------------------

## Notes

-   `nb_pre` and `nb_post` are capped at the maximum lags/leads available in the data in both functions.
-   Horizon **−1** is the reference period and is always normalised to zero in the plot.
-   `lpdid_adj` is the slowest method as it calls `avg_comparisons()` for each horizon.

------------------------------------------------------------------------

## References

Dube, A., Girardi, D., Jordà, O., & Taylor, A. M. (2025). A local projections approach to difference-in-differences. *Journal of Applied Econometrics*, 40(7), 741-758.

Callaway, B., & Sant'Anna, P. H. C. (2021). Difference-in-differences with multiple time periods. *Journal of Econometrics*, 225(2), 200–230.
