# lpdidcsa

**Local Projections DiD and Callaway & Sant'Anna Event Studies in R.**

## License and Credits

-   The code in this package is licensed under the [MIT License](LICENSE).
-   The `mpdta` dataset is sourced from the [`did` package](https://cloud.r-project.org/web/packages/did/index.html) and is licensed under [GPL-2](LICENSE.did).
-   The lpdidcsa is inspired by the following packages:\
    Brantly Callaway's and Pedro H. C. Sant'Anna's [`did`](https://cloud.r-project.org/web/packages/did/index.html), Laurent Bergé's [`fixest`](https://cloud.r-project.org/web/packages/did/index.html), Daniele Girardi's [`lpdid`](https://github.com/danielegirardi/lpdid/) and Alex Cardazzi's [`lpdid`](https://github.com/alexCardazzi/lpdid)
-   We thank Daniele Girardi and Pedro Sant'Anna for conversations which helped us to clarify the different methods.

## Description

`lpdidcsa` provides two functions that together cover the full workflow:

1.  **`lpdidcsa_data()`** — reshapes a panel into the format required for estimation
2.  **`lpdidcsa()`** — estimates event study coefficients and produces a plot

The package implements the LP-DiD estimator of Dube, Girardi, Jordà & Taylor (2025) and a local-projection version of the Callaway & Sant'Anna (2021) estimator, with support for absorbing and non-absorbing treatments. It enables to run LP-DiD and reweighted LP-DiD with or without control variables, adjusted LP-DiD with control variables, inverse probability weighting LP-DiD, LP-CSA with or without control variables and inverse probability weighting LP-CSA.

In this [methodological note](docs/LPDiD_methodological_note.pdf), we clarify the links between CSA's, LP-CSA and LP-DID estimates and introduce the inverse probability weighting strategy for LP-DiD.

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
| `n_pre` | `integer` / `NULL` | `NULL` | Pre-treatment periods (all if `NULL`) |
| `n_post` | `integer` / `NULL` | `NULL` | Post-treatment periods (all if `NULL`) |
| `p_variables` | named `list` | `list()` | Horizon-specific variables |
| `h_variables` | `character` / `NULL` | `NULL` | Variables carried at every horizon |
| `type_horizon` | `character` | `"wide"` | `"wide"` or `"long"` output format |

### `p_variables` naming convention

List names follow `"t0"`, `"t1"`, `"tm1"`, `"tm2"`, … where `"tm"` stands for *t minus*:

``` r
p_variables <- list(
  tm2 = c("x1", "x2"),   # controls measured at t-2
  tm1 = c("x1", "x2"),   # controls measured at t-1
  t1  = c("x1")          # control measured at t+1
)
```

### Output columns added

| Column | Description |
|------------------------------------|------------------------------------|
| `first.treat` | First treatment period per unit (`NA` for never-treated). Absorbing only. |
| `last.treat` | Most recent treatment period as of *t*. Non-absorbing only. |
| `next.treat` | Next upcoming treatment period. Non-absorbing only. |
| `dtreat` | Treatment change: *D_it − D_i,t−1* |
| `<dep>_tm2`, `<dep>_t0`, … | Outcome at each horizon (wide format) |
| `horizon` | Horizon integer *h = t′ − t* (long format) |

### Example

``` r
library(data.table)
data(mpdta_r, package = "lpdidcsa")

p_variables <- list(
  tm2 = c("lpop"),
  tm1 = c("lpop")
)

# Wide format (one column per horizon)
df_wide <- lpdidcsa_data(
  data        = mpdta_r,
  unit        = "countyreal",
  time        = "year",
  dependent   = "lemp",
  treat       = "treat",
  absorbing   = TRUE,
  p_variables = p_variables
)

# Long format (one row per unit × horizon)
df_long <- lpdidcsa_data(
  data         = mpdta_r,
  unit         = "countyreal",
  time         = "year",
  dependent    = "lemp",
  treat        = "treat",
  absorbing    = TRUE,
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
| `n_pre` | `integer` / `NULL` | `NULL` | Pre-treatment horizons (all if `NULL`) |
| `n_post` | `integer` / `NULL` | `NULL` | Post-treatment horizons (all if `NULL`) |
| `controls` | `character` / `NULL` | `NULL` | Control variables (fixest-style strings accepted) |
| `FE` | `character` / `NULL` | `NULL` | Additional fixed effects |
| `clusters` | `character` | `unit` | Clustering column |
| `weight` | `character` / `NULL` | `NULL` | Sampling weight column |
| `controls_h` | `character` / `NULL` | `NULL` | Horizon-specific control variables |
| `FE_h` | `character` / `NULL` | `NULL` | Horizon-specific additional fixed effects |
| `clusters_h` | `character` | `unit` | Horizon-specific clustering column |
| `weight_h` | `character` / `NULL` | `NULL` | Horizon-specific sampling weight column |
| `meth` | `character` | `"lpdid_ipw"` | Estimator (see table below) |
| `anticipation` | `integer` | `0` | Number of time periods ahead of the treatment where participants can anticipate the treatment and adapt their behavior |
| `absorbing` | `logical` | `TRUE` | Absorbing treatment assumption |
| `reentry` | `integer` / `NULL` | `NULL` | Min. periods since last treatment for re-entry in the analysis (only when absorbing = FALSE) |
| `type_horizon` | `character` | `"wide"` | Format of the input data: `"wide"` or `"long"` |
| `horizon` | `character` | `"horizon"` | Horizon column name (long format only) |
| `one_reg` | `logical` | `FALSE` | Run all horizons in one stacked regression (long format only) |

### Estimators

#### Without control variables

| `meth` | Description | ATT type |
|------------------------|------------------------|------------------------|
| `"lpdid"` | LP-DiD | Variance-weighted |
| `"lpdid_rw"` | Reweighted LP-DiD | Equally-weighted |
| `"lpdid_adj"` | LP-DiD with adjusted regression (`avg_comparisons`) | Equally-weighted |
| `"lpcsa"` | LP-CSA | Equally-weighted |
| `"lpdid_ipw"` | runs `"lpdid_rw"` | Equally-weighted |
| `"lpcsa_ipw"` | runs `"lpcsa"` | Equally-weighted |

#### With control variables

| `meth` | Description | ATT type | Hypothesis |
|------------------|------------------|------------------|------------------|
| `"lpdid"` | LP-DiD with control variables | Variance-weighted | Treatment effects do not vary with values of covariates |
| `"lpdid_rw"` | Reweighted LP-DiD with control variables | Equally-weighted | Treatment effects do not vary with values of covariates |
| `"lpdid_adj"` | LP-DiD with adjusted regression (`avg_comparisons`) | Equally-weighted |  |
| `"lpdid_ipw"` | LP-DiD with inverse probability weighting | Equally-weighted |  |
| `"lpcsa"` | LP-CSA with control variables | Equally-weighted | Treatment effect do not vary with values of covariates |
| `"lpcsa_ipw"` | LP-CSA with inverse probability weighting | Equally-weighted |  |

### Return value

A named list with four elements:

| Element | Description |
|------------------------------------|------------------------------------|
| `est` | Main results: `data.table` with columns `h`, `variable`, `estimate`, `se`, `T`, `pvalue`, `n_obs`, `formula` |
| `est_det` | Cohort-level estimates (CSA and adjusted methods only) |
| `ps` | Propensity score estimates (IPW methods only) |
| `plot` | `ggplot2` event study plot with 95% confidence intervals |

### Example

``` r
library("lpdidcsa")

# Example 1 with mpdta dataset ----
data(mpdta_r, package = "lpdidcsa")

# Step 1: prepare data
df <- lpdidcsa_data(
  data      = mpdta_r,
  unit      = "countyreal",
  time      = "year",
  dependent = "lemp",
  treat     = "treat"
)

# Step 2: estimate (method lpdid_ipw with one covariate)
res <- lpdidcsa(
  data      = df,
  unit      = "countyreal",
  time      = "year",
  dependent = "lemp",
  controls = "lpop",
  meth      = "lpdid_ipw"
)

res$plot   # event study plot
res$est    # coefficient table
res$ps    # coefficient table for the first stage propensity score

# Example 2 with toydataset mimicking a linked employer-employee panel dataset ----
toydata <- sim_staggered_panel(n=10000,n_firm=500)
nrow(toydata)

df_w  <- lpdidcsa_data(toydata, 
                     unit = "id", 
                     time = "t",
                     dependent = "log_earnings",
                     treat = "treat", # treatment dummy variable
                     n_pre = 12, # number of pre-treatment periods
                     n_post= 12, # number of post-treatment periods
                     h_variables="firm_id" # other variables for which all horizons are computed
                     )
colnames(df_w)

res <- lpdidcsa(df_w, 
                unit = "id", 
                time = "t",
                dependent = "log_earnings", 
                meth = "lpdid_ipw", 
                n_pre=12, # number of pre-treatment periods
                n_post=12, # number of post-treatment periods
                controls=c("female"), # horizon invariant control variables
                clusters="firm_id_tm1", # horizon invariant clustering variable
                clusters_h="firm_id", # horizon specific clustering variable
                FE="firm_id_tm1" # horizon invariant fixed effects
                )
res$plot   # event study plot
res$est    # coefficient table
res$ps    # coefficient table for the first stage propensity score
```

------------------------------------------------------------------------

## Notes

-   `n_pre` and `n_post` are capped at the maximum lags/leads available in the data in both functions.
-   Horizon **−1** is the reference period and is always normalised to zero in the plot.
-   `lpdid_adj` is the slowest method as it calls `avg_comparisons()` for each horizon.

------------------------------------------------------------------------

## References

Dube, A., Girardi, D., Jordà, O., & Taylor, A. M. (2025). A local projections approach to difference-in-differences. *Journal of Applied Econometrics*, 40(7), 741-758.

Callaway, B., & Sant'Anna, P. H. C. (2021). Difference-in-differences with multiple time periods. *Journal of Econometrics*, 225(2), 200–230.
