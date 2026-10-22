# lpdidcsa

**Local Projections DiD and Callaway & Sant'Anna Event Studies in R.**

## License and Credits

-   The code in this package is licensed under the [MIT License](LICENSE).
-   The `mpdta` dataset is sourced from the [`did` package](https://cloud.r-project.org/web/packages/did/index.html) and is licensed under [GPL-2](LICENSE.did).
-   The lpdidcsa is inspired by the following packages:\
    Brantly Callaway's and Pedro H. C. Sant'Anna's [`did`](https://cloud.r-project.org/web/packages/did/index.html) and [`DRDID`](https://cloud.r-project.org/web/packages/DRDID/index.html), Laurent Bergé's [`fixest`](https://cloud.r-project.org/web/packages/fixest/index.html), Daniele Girardi's [`lpdid`](https://github.com/danielegirardi/lpdid/) and Alex Cardazzi's [`lpdid`](https://github.com/alexCardazzi/lpdid)
-   We thank Daniele Girardi and Pedro Sant'Anna for conversations which helped us to clarify the different methods.

## Description

`lpdidcsa` provides two functions that together cover the full workflow:

1.  **`lpdidcsa_data()`** — reshapes a panel into the format required for estimation
2.  **`lpdidcsa()`** — estimates event study coefficients, produces an event study plot and, optionally, predicted outcome paths of the treated with and without treatment

The package implements the LP-DiD estimator of Dube, Girardi, Jordà & Taylor (2025) and a local-projection version of the Callaway & Sant'Anna (2021) estimator, with support for absorbing and non-absorbing treatments. Seven estimators are available: LP-DiD and reweighted LP-DiD with or without control variables, LP-DiD with regression adjustment, inverse probability weighting LP-DiD, doubly robust LP-DiD (Sant'Anna & Zhao, 2020), LP-CSA with or without control variables and inverse probability weighting LP-CSA.

For the regression-adjustment and doubly robust estimators (and optionally for the reweighted and IPW ones), standard errors can be corrected for the estimation of the first-stage nuisance models (outcome regression, propensity score, reweighting shares) through cluster-level influence functions.

In this [methodological note](docs/LPDiD_methodological_note.pdf), we clarify the links between CSA's, LP-CSA and LP-DID estimates and introduce the inverse probability weighting strategy for LP-DiD.

------------------------------------------------------------------------

## Installation

``` r
# Install from GitHub (requires remotes)
remotes::install_github("oliviergodechot/lpdidcsa")
```

The package depends on `data.table`, `fixest`, `ggplot2` and `sandwich`.

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
                      (+ predicted paths and figures if prediction = TRUE)
```

------------------------------------------------------------------------

## Step 1 — `lpdidcsa_data()`

Prepares the panel by computing treatment timing variables and creating leads/lags of the outcome (and of selected variables) at each horizon.

### Arguments

| Argument | Type | Default | Description |
|------------------|------------------|------------------|------------------|
| `data` | `data.frame` / `data.table` | — | Input panel (left untouched: the function works on a copy) |
| `unit` | `character` | — | Unit identifier variable |
| `time` | `character` | `"year"` | Time period variable (integer-valued) |
| `dependent` | `character` | — | Outcome variable |
| `treat` | `character` | `"treat"` | Binary treatment indicator (0/1) |
| `absorbing` | `logical` | `TRUE` | Whether treatment is absorbing |
| `n_pre` | `integer` / `NULL` | `NULL` | Pre-treatment periods (all if `NULL`) |
| `n_post` | `integer` / `NULL` | `NULL` | Post-treatment periods (all if `NULL`) |
| `p_variables` | named `list` | `list()` | Variables to measure at specific horizons |
| `h_variables` | `character` / `NULL` | `NULL` | Variables carried at every horizon |
| `type_horizon` | `character` | `"wide"` | `"wide"` or `"long"` output format |

### `p_variables` naming convention

List names follow `"t0"`, `"t1"`, `"tm1"`, `"tm2"`, … where `"tm"` stands for *t minus* (`"tm0"` is accepted as a synonym of `"t0"`):

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
| `dtreat` | Treatment switch: 1 if *D_it − D_i,t−1 = 1*, 0 otherwise |
| `<dep>_tm1` | Outcome at the reference period *t−1* (both formats) |
| `<dep>_tm2`, `<dep>_t0`, … | Outcome at each horizon (wide format) |
| `<var>_tm2`, `<var>_t1`, … | `p_variables` at their horizons, `h_variables` at every horizon (wide format) |
| `<dep>_th`, `<var>_th` | Outcome and `h_variables` at the horizon of the row (long format) |
| `horizon` | Horizon integer *h = t′ − t* (long format) |

Horizon −1 is the reference period and is never computed as a horizon. If some units have gaps in their time series, the missing intermediate periods are added (with a message) so that leads and lags refer to calendar time; duplicated unit × period observations stop the function.

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

# Long format (one row per unit × reference period × horizon)
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

Estimates event study coefficients horizon by horizon (or in a single stacked regression for some methods) and returns a `ggplot2` plot.

### Arguments

| Argument | Type | Default | Description |
|------------------|------------------|------------------|------------------|
| `data` | `data.frame` / `data.table` | — | Output of `lpdidcsa_data()` |
| `unit` | `character` | — | Unit identifier variable |
| `time` | `character` | `"year"` | Time period variable |
| `dependent` | `character` | — | Outcome variable (initial name, without suffix) |
| `dtreat` | `character` | `"dtreat"` | Treatment change indicator variable |
| `n_pre` | `integer` / `NULL` | `NULL` | Pre-treatment horizons (all if `NULL`) |
| `n_post` | `integer` / `NULL` | `NULL` | Post-treatment horizons (all if `NULL`) |
| `controls` | `character` / `NULL` | `NULL` | Horizon-invariant control variables (fixest-style strings accepted, e.g. `"i(x1) + log(x2_tm1)"`) |
| `controls_h` | `character` / `NULL` | `NULL` | Horizon-variant control variables (initial names; suffixes added automatically) |
| `FE` | `character` / `NULL` | `NULL` | Horizon-invariant fixed effects (fixest-style strings accepted, e.g. `"x1^x2"`) |
| `FE_h` | `character` / `NULL` | `NULL` | Horizon-variant fixed effects (initial names) |
| `clusters` | `character` | `unit` | Horizon-invariant clustering variables |
| `clusters_h` | `character` / `NULL` | `NULL` | Horizon-variant clustering variables (initial names) |
| `weight` | `character` / `NULL` | `NULL` | Horizon-invariant sampling weight variable |
| `weight_h` | `character` / `NULL` | `NULL` | Horizon-variant sampling weight variable (initial name) |
| `meth` | `character` | `"lpdid_ipw"` | Estimator (see tables below) |
| `anticipation` | `integer` | `0` | Number of periods before treatment in which units can anticipate it; the reference period becomes *t − 1 − anticipation* |
| `absorbing` | `logical` | `TRUE` | Absorbing treatment assumption |
| `reentry` | `integer` / `NULL` | `NULL` | Min. periods since last treatment for re-entry in the analysis (only when `absorbing = FALSE`) |
| `type_horizon` | `character` | `"wide"` | Format of the input data: `"wide"` or `"long"` (detected and corrected if inconsistent with the data) |
| `horizon` | `character` | `"horizon"` | Horizon variable name (long format only) |
| `one_reg` | `logical` | `FALSE` | Run all horizons in one stacked regression (long format only; not available for `lpdid_adj` and `lpdid_dr`) |
| `se_correction` | `logical` | `TRUE` | For `lpdid_adj` and `lpdid_dr`: add a `se_correct` column propagating the first-stage estimation uncertainty |
| `se_correction_rw_ipw` | `logical` | `FALSE` | For `lpdid_rw` and `lpdid_ipw` (`one_reg = FALSE`): add a `se_correct` column propagating the uncertainty of the reweighting shares and, for IPW, of the propensity score |
| `ps_trim` | `numeric` | `0.995` | Control observations with propensity score ≥ `ps_trim` are dropped (`lpdid_ipw`, `lpcsa_ipw`, `lpdid_dr`); `1` disables trimming |
| `time_x_controls` | `logical` | `FALSE` | Let the effect of every control (and absorbed FE) vary by period, in the regressions where controls enter (`one_reg = FALSE` only) |
| `prediction` | `logical` | `FALSE` | Compute predicted outcome paths of the treated with and without treatment, and draw three figures (`one_reg = FALSE` only) |
| `prediction_se` | `logical` | `FALSE` | Add standard errors and 95% confidence intervals to the predictions (implies `prediction = TRUE`) |

### Estimators

#### Without control variables

| `meth` | Description | ATT type |
|------------------------|------------------------|------------------------|
| `"lpdid"` | LP-DiD | Variance-weighted across cohorts |
| `"lpdid_rw"` | Reweighted LP-DiD | Equally-weighted |
| `"lpdid_adj"` | LP-DiD with regression adjustment on period | Equally-weighted |
| `"lpdid_dr"` | Doubly robust LP-DiD (period only) | Equally-weighted |
| `"lpcsa"` | LP-CSA | Equally-weighted |
| `"lpdid_ipw"` | runs `"lpdid_rw"` | Equally-weighted |
| `"lpcsa_ipw"` | runs `"lpcsa"` | Equally-weighted |

Without control variables, all equally-weighted estimators yield the same ATT; only their standard errors may differ. `"lpdid"` yields a convex, variance-weighted average of cohort-specific ATTs; it equals the equally-weighted ATT when treatment effects do not vary across cohorts.

*Cohort* refers here to the period in which treatment starts (with non-absorbing treatments, the period of the treatment switch). At a given horizon *h*, it coincides with the calendar period of the treated observation.

#### With control variables

| `meth` | Description | ATT type | Hypothesis |
|------------------|------------------|------------------|------------------|
| `"lpdid"` | LP-DiD with control variables | Variance-weighted across cohorts | Treatment effects do not vary with covariates. If they also do not vary across cohorts, the variance-weighted ATT equals the equally-weighted ATT |
| `"lpdid_rw"` | Reweighted LP-DiD with control variables | Equally-weighted | Treatment effects do not vary with covariates (they may vary across cohorts) |
| `"lpdid_adj"` | LP-DiD with regression adjustment: outcome regression fitted on clean controls, predicted for the treated | Equally-weighted | Effects may vary with covariates and cohorts; correct outcome model for controls |
| `"lpdid_ipw"` | LP-DiD with inverse probability weighting | Equally-weighted | Effects may vary with covariates and cohorts; correct propensity score model |
| `"lpdid_dr"` | Doubly robust LP-DiD (outcome regression + IPW, AIPW-style score) | Equally-weighted | Effects may vary with covariates and cohorts; correct outcome model *or* correct propensity score model |
| `"lpcsa"` | LP-CSA with control variables | Equally-weighted | Treatment effects do not vary with covariates (they may vary across cohorts) |
| `"lpcsa_ipw"` | LP-CSA with inverse probability weighting | Equally-weighted | Effects may vary with covariates and cohorts; correct propensity score model |

All estimators rely on (conditional) parallel trends between treated and clean controls. The cohort dimension is handled by reweighting (`lpdid_rw`, `lpdid_ipw`) or by cohort-specific coefficients (`lpcsa`, `lpcsa_ipw`); `lpdid` keeps the regression's variance weighting, which is always non-negative with clean controls. Within a cohort, linear controls weight treated units by the residual variance of the treatment given the covariates (Słoczyński, 2022), hence the homogeneity condition for `lpdid`, `lpdid_rw` and `lpcsa`. The regression-adjustment, IPW and doubly robust estimators avoid it, at the cost of a model for the outcome or the treatment.

### Standard errors

| `meth` | `se` | `se_correct` |
|------------------|------------------------|------------------------|
| `"lpdid"`, `"lpcsa"`, `"lpcsa_ipw"` | Clustered (fixest) | — |
| `"lpdid_rw"`, `"lpdid_ipw"` | Clustered, weights treated as known | With `se_correction_rw_ipw = TRUE`: propagates the reweighting shares and, for IPW, the propensity score |
| `"lpdid_adj"` | Clustered, first stage treated as known | With `se_correction = TRUE` (default): propagates the outcome regression; reproduces `lpdid.ado`'s RA standard error |
| `"lpdid_dr"` | Clustered, first stages treated as known | With `se_correction = TRUE` (default): propagates the outcome regression and the propensity score |

The event study plot uses `se_correct` where it has been computed and `se` otherwise. The corrections use single-way clustering (first element of `clusters`, or `unit`). With absorbed fixed effects (`FE`) in a first-stage model, the correction is partial: the fixed effects' own sampling uncertainty is not propagated. Enter low-dimensional fixed effects as `factor()` terms in `controls` to have them fully propagated. Trimming is treated as fixed.

### Return value

A named list with eight elements:

| Element | Description |
|------------------------------------|------------------------------------|
| `est` | Main results: `data.table` with columns `h`, `variable`, `estimate`, `se`, (`se_correct`), `T`, `pvalue`, `n_obs`, `n_treat`, `n_ctrl`, (`sw_treat`, `sw_ctrl`, `weight` with weights), `formula`, `clusters` (non-CSA methods) |
| `est_det` | Cohort-level estimates (CSA methods only) |
| `ps` | First-stage propensity score estimates (IPW and DR methods only) |
| `or` | First-stage outcome regression on the control group (`lpdid_adj` and `lpdid_dr` only) |
| `plot` | `ggplot2` event study plot with 95% confidence intervals |
| `pred` | Predicted paths (if `prediction = TRUE`, see below) |
| `pred_anchor` | Baseline levels used as anchors, with their standard errors (if `prediction = TRUE`) |
| `pred_plot` | List of three `ggplot2` figures: `change`, `common`, `group` (if `prediction = TRUE`) |

### Predictions

With `prediction = TRUE`, `lpdidcsa()` also returns, for each horizon, the predicted outcome change since the reference period for the treated:

-   **with treatment**: μ₁(h), the weighted mean of the observed change among treated observations of the estimation sample;
-   **without treatment** (counterfactual, labelled `Control`): μ₀(h) = μ₁(h) − ATT(h), with the ATT reported in `est`.

For `lpdid_adj`, μ₀ is the mean of the control-model predictions over the treated; for `lpdid_dr`, it adds the IPW-weighted control residual term; for `lpdid`, it equals the average prediction over the treated with `dtreat` set to 0 (it inherits the variance-weighted estimand). For CSA methods, cohort means are aggregated with the same cohort weights as `fixest::aggregate()`.

Two anchors (reference-period levels) are computed on the *h* = 0 sample: the treated's weighted mean of the reference-period outcome, and the controls' mean adjusted to the treated's period mix (not to covariates). Three figures are drawn:

| Figure | Content | Gap between curves |
|------------------|------------------------|------------------------|
| `pred_plot$change` | Predicted changes, both curves at 0 at the reference period | ATT |
| `pred_plot$common` | Levels: treated anchor + predicted changes for both curves | ATT |
| `pred_plot$group` | Levels: each group's own anchor + predicted changes | Period-adjusted baseline gap + ATT |

With `prediction_se = TRUE`, standard errors come from cluster-level influence functions: IF(μ₁) for the treated, IF(μ₁) − IF(ATT) for the counterfactual (keeping their covariance), using the corrected IF(ATT) where it is available. In the common-anchor figure the anchor is treated as fixed (it shifts both curves jointly), so its confidence intervals are those of the changes; the group-anchor figure includes the anchors' sampling uncertainty. Overlapping intervals should not be read as a test of the gap: its standard error is that of the ATT in `est`.

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
  controls  = "lpop",
  meth      = "lpdid_ipw"
)

res$plot   # event study plot
res$est    # coefficient table
res$ps     # coefficient table for the first stage propensity score

# Doubly robust estimation, with corrected standard errors and predictions
res <- lpdidcsa(
  data          = df,
  unit          = "countyreal",
  time          = "year",
  dependent     = "lemp",
  controls      = "lpop",
  meth          = "lpdid_dr",
  prediction_se = TRUE
)

res$est              # se (naive) and se_correct (influence-function corrected)
res$or               # first-stage outcome regression on the control group
res$pred             # predicted paths of the treated with and without treatment
res$pred_plot$change # figure without anchor
res$pred_plot$common # levels, common anchor
res$pred_plot$group  # levels, group anchors

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
res$ps     # coefficient table for the first stage propensity score
```

------------------------------------------------------------------------

## Notes

-   `n_pre` and `n_post` are capped at the maximum lags/leads available in the data in both functions.
-   `lpdid_adj` and `lpdid_dr` are only available with `one_reg = FALSE`; `prediction`, `time_x_controls` and `se_correction_rw_ipw` are ignored (with a warning) when `one_reg = TRUE`.
-   With `time_x_controls = TRUE`, the generated formulas use base-R `factor()` terms instead of fixest's `i()`; the treatment coefficient is then named `dtreat` instead of `dtreat::1`, and controls written `i(var)` must not carry further arguments (use `factor(var)`).
-   `se_correction = FALSE` skips the influence-function computations of `lpdid_adj` and `lpdid_dr`, which is faster on large datasets (they are still computed if `prediction_se = TRUE`).
-   Prediction standard errors use the *G/(G−1)* small-sample factor and can differ slightly from the ATT standard errors in `est`, which follow each method's own convention.

------------------------------------------------------------------------

## References

Callaway, B., & Sant'Anna, P. H. C. (2021). Difference-in-differences with multiple time periods. *Journal of Econometrics*, 225(2), 200–230.

Dube, A., Girardi, D., Jordà, O., & Taylor, A. M. (2025). A local projections approach to difference-in-differences. *Journal of Applied Econometrics*, 40(7), 741-758.

Sant'Anna, P. H. C., & Zhao, J. (2020). Doubly robust difference-in-differences estimators. *Journal of Econometrics*, 219(1), 101–122.

Słoczyński, T. (2022). Interpreting OLS estimands when treatment effects are heterogeneous: Smaller groups get larger weights. *Review of Economics and Statistics*, 104(3), 501–509.
