#' LP-DiD and LP-CSA Event Study Estimation
#'
#' Estimates event study coefficients using Local Projections
#' Difference-in-Differences (LP-DiD). Seven estimators are supported:
#'
#' \describe{
#'   \item{\code{"lpdid"}}{LP-DiD with control variables. Yields a
#'     variance-weighted ATT under the hypothesis that treatment effects do not
#'     vary with the value of covariates.}
#'   \item{\code{"lpdid_rw"}}{Reweighted LP-DiD with control variables.
#'     Yields an equally-weighted ATT under the hypothesis that treatment effects do not
#'     vary with the value of covariates.}
#'   \item{\code{"lpdid_adj"}}{LP-DiD with regression adjustment (RA), via a
#'     control-only outcome regression and a cluster-level combined
#'     influence-function standard error. Yields an equally-weighted ATT.}
#'   \item{\code{"lpdid_dr"}}{LP-DiD doubly robust (RA + IPW combined), via
#'     an AIPW-style score and a cluster-level combined influence-function
#'     standard error. Yields an equally-weighted ATT. \code{one_reg = TRUE}
#'     is not implemented for this method.}
#'   \item{\code{"lpdid_ipw"}}{LP-DiD with inverse probability weighting.
#'     Yields an equally-weighted ATT.}
#'   \item{\code{"lpcsa"}}{LP-CSA with control variables. Yields an
#'     equally-weighted ATT under the hypothesis that treatment effects do not
#'     vary with the value of covariates.}
#'   \item{\code{"lpcsa_ipw"}}{LP-CSA with inverse probability weighting
#'     (local-projection version of Callaway & Sant'Anna).}
#' }
#'
#' @param data A \code{data.frame} or \code{data.table}, typically the output
#'   of \code{lpdidcsa_data()}.
#' @param unit Character. Column name for the unit identifier.
#' @param time Character. Column name for the time period
#'   (default: \code{"year"}).
#' @param dependent Character. Column name for the dependent variable.
#' @param dtreat Character. Column name for the treatment change indicator
#'   (0/1; default: \code{"dtreat"}).
#' @param n_pre Integer >= 0 or \code{NULL}. Number of pre-treatment horizons.
#'   If \code{NULL}, all available pre-periods are used.
#' @param n_post Integer >= 0 or \code{NULL}. Number of post-treatment
#'   horizons. If \code{NULL}, all available post-periods are used.
#' @param controls Character vector of the column names for horizon-invariant control variable, or
#'   \code{NULL}. Accepts \code{fixest}-style formula strings (e.g.
#'   \code{controls=c("i(var1) + log(var2_tm1)","var3_tm2")}).
#' @param controls_h Character vector of the column names for horizon-variant control
#'  variables, or \code{NULL}. Does not accept \code{fixest}-style
#'  formula strings. Enter the initial variable names (e.g., \code{controls_h=c("var1","var2","var3")} ).
#'  The function adds the suffix (e.g. '_tm5',...,'_t5') if \code{type_horizon = "wide"} and
#'  the suffix '_th' if \code{type_horizon = "long"}.
#' @param FE Character vector of column names for high-dimensional horizon-invariant
#'   fixed-effect variables , or \code{NULL}. Accepts \code{fixest}-style
#'   FE formula strings (e.g., \code{FE=c("var1 + var2_tm1","var3_tm2^var4_tm2","var5") }).
#' @param FE_h Character vector of column names for for high-dimensional horizon-variant
#'   fixed-effect variables, or \code{NULL}. Enter the initial variable names (e.g.
#'    \code{FE_h=c("var1","var2")} ). Does not accept \code{fixest}-style FE formula strings.
#'   The function adds the suffix (e.g. '_tm5',...,'_t5') if \code{type_horizon = "wide"}
#'   and the suffix '_th' if \code{type_horizon = "long"}. The function
#'   adds the suffix (e.g. '_tm5',...,'_t5') if \code{type_horizon = "wide"} and the
#'   suffix '_th' if \code{type_horizon = "long"}.
#' @param clusters Character. Character vector of column names for horizon-invariant variables used for
#'   clustering standard errors (e.g. \code{clusters=c("var1","var2_tm1")} ). Defaults to \code{unit}.
#' @param clusters_h Character vector of column names for horizon-variant variables
#'   used for clustering standard errors (e.g. \code{clusters_h=c("var1","var2")} ). Defaults to
#'   \code{NULL}. Enter the initial variable names. The function
#'   adds the suffix (e.g. '_tm5',...,'_t5') if \code{type_horizon = "wide"} and the
#'   suffix '_th' if \code{type_horizon = "long"}.
#' @param weight Character. Column name for horizon-invariant survey/sampling
#'   weight variable, or \code{NULL}.
#' @param weight_h Character. Column name for horizon-variant survey/sampling
#'  weight variable, or \code{NULL}. Enter the initial variable name. The function
#'   adds the suffix (e.g. '_tm5',...,'_t5') if \code{type_horizon = "wide"} and the
#'   suffix '_th' if \code{type_horizon = "long"}.
#' @param meth Character. Estimator to use. One of \code{"lpdid"},
#'   \code{"lpdid_rw"}, \code{"lpdid_adj"}, \code{"lpdid_dr"}, \code{"lpdid_ipw"},
#'   \code{"lpcsa_ipw"}, \code{"lpcsa"} (default: \code{"lpdid_ipw"}).
#' @param anticipation Integer. Number of time periods ahead of the treatment
#' where participants can anticipate the treatment and adapt their behavior
#' (default: \code{0L}).
#' @param absorbing Logical. If \code{TRUE} (default), treatment is assumed
#'   absorbing and not-yet-treated units serve as controls. If \code{FALSE},
#'   a non-absorbing design is used.
#' @param reentry Integer or \code{NULL}. If \code{absorbing = FALSE }, minimum
#'   number of periods  since last treatment for a unit to re-enter the DiD comparison
#'   set as a control unit or for a new treatment.
#'   \code{NULL} forbids reentry of already treated units (default: \code{NULL}).
#' @param type_horizon Character. Whether the input \code{data} is in
#'   \code{"wide"} (default) or \code{"long"} format, as produced by
#'   \code{lpdidcsa_data()}.
#' @param horizon Character. Column name for the horizon variable when
#'   \code{type_horizon = "long"} (default: \code{"horizon"}).
#' @param one_reg Logical. If \code{TRUE} and \code{type_horizon = "long"},
#'   runs all horizons in a single stacked regression instead of one per
#'   horizon (default: \code{FALSE}). Not supported for \code{meth = "lpdid_dr"},
#'   and no longer supported for \code{meth = "lpdid_adj"} (see Details).
#' @param se_correction Logical. If \code{TRUE} (default) and
#'   \code{meth \%in\% c("lpdid_adj","lpdid_dr")}, computes a corrected
#'   standard error (\code{se_correct} column) via a cluster-level combined
#'   influence function, in addition to the naive standard error
#'   (\code{se} column). The correction propagates the estimation
#'   uncertainty of the first-stage nuisance model(s) (outcome regression
#'   for RA; outcome regression and propensity score for DR) into the final
#'   ATT's standard error, and has been validated to reproduce
#'   \code{lpdid.ado}'s \code{listreg}-based RA standard error exactly on a
#'   test panel. If \code{FALSE}, only the naive standard error is computed,
#'   which is faster on large datasets since it skips the additional
#'   \code{sandwich::estfun()}/\code{bread()} calls (no extra model refits
#'   are involved either way). Ignored for all other methods.
#' @param se_correction_rw_ipw Logical. If \code{TRUE} and
#'   \code{meth \%in\% c("lpdid_rw","lpdid_ipw")} (with \code{one_reg = FALSE}),
#'   adds a \code{se_correct} column whose standard error propagates the
#'   estimation uncertainty of the reweighting stage (the per-period treated
#'   shares behind \code{rewgt_h}) and, for \code{"lpdid_ipw"}, of the
#'   propensity-score logit, via a cluster-level combined influence function
#'   (see Details). Default \code{FALSE}: the correction is usually small and
#'   requires extra \code{estfun()}/\code{bread()}/\code{model.matrix()} calls.
#'   Ignored for all other methods and for \code{one_reg = TRUE}.
#' @param ps_trim Numeric in (0, 1]. Control observations whose estimated
#'   propensity score is \code{>= ps_trim} are dropped before weighting
#'   (treated observations are always kept). Applies to \code{"lpdid_ipw"},
#'   \code{"lpcsa_ipw"} and \code{"lpdid_dr"} (default: \code{0.995});
#'   \code{ps_trim = 1} disables trimming.
#' @param time_x_controls Logical. If \code{TRUE}, the effect of every
#'   control is allowed to vary by period and every absorbed fixed effect is
#'   interacted with time (\code{FE^time}) in the regressions where controls
#'   enter: the final regression for \code{"lpdid"}, \code{"lpdid_rw"} and
#'   \code{"lpcsa"}; the control-group outcome regression for
#'   \code{"lpdid_adj"} and \code{"lpdid_dr"}; the propensity-score logit for
#'   \code{"lpdid_ipw"}, \code{"lpcsa_ipw"} and \code{"lpdid_dr"}. Controls
#'   are interacted with \code{factor(time)} using base-R formula operators,
#'   so compound terms (\code{log(x)}, \code{a*b}, \code{poly(x, 2)},
#'   \code{factor(v)}) are interacted as a whole: the RHS becomes
#'   \code{factor(time)*(controls)} when time is not otherwise in the model
#'   (PS logit, outcome regression without \code{FE}), and
#'   \code{(controls) + factor(time):(controls)} when time is absorbed (final
#'   regressions, or any regression with \code{FE}). Factor controls written
#'   \code{i(var)} are rewritten as \code{factor(var)}; \code{i()} terms with
#'   further arguments are not accepted (use \code{factor()}). The
#'   interaction is only added when the estimation sample of that regression
#'   contains at least two periods; otherwise the controls enter
#'   uninteracted. Note that a unit-level FE interacted with time is an
#'   observation-level FE in a per-horizon regression. With this option the
#'   generated formulas contain no fixest \code{i()} term (mixing \code{i()}
#'   with base-R \code{factor()} or \code{:} terms fails when fixest is not
#'   attached), so the treatment coefficient is named \code{dtreat} instead
#'   of \code{dtreat::1}, and the \code{"lpcsa"} cohort coefficients
#'   \code{dtreat:factor(year)2004} instead of \code{dtreat::1:year::2004}.
#'   Only for \code{one_reg = FALSE} (default: \code{FALSE}).
#'
#' @return A named list with five elements:
#'   \describe{
#'     \item{\code{est}}{Main estimation results as a \code{data.table} with
#'       columns \code{h}, \code{variable}, \code{estimate}, \code{se},
#'       \code{T}, \code{pvalue}, \code{n_obs}, \code{formula}, \code{clusters},
#'       and, for \code{meth \%in\% c("lpdid_adj","lpdid_dr")} (and for
#'       \code{"lpdid_rw"}/\code{"lpdid_ipw"} when
#'       \code{se_correction_rw_ipw = TRUE}), \code{se_correct}. The
#'       \code{formula} column is the final-stage regression; for
#'       \code{"lpdid_adj"} and \code{"lpdid_dr"} the first-stage formulas
#'       are in \code{or} and \code{ps}.}
#'     \item{\code{est_det}}{Cohort-level (disaggregated) estimates as a
#'       \code{data.table}. Non-\code{NULL} for CSA methods only.}
#'     \item{\code{ps}}{Propensity score model estimates as a
#'       \code{data.table}. Non-\code{NULL} for IPW and DR methods only.}
#'     \item{\code{or}}{First-step outcome regression on the control group
#'       (one per horizon) as a \code{data.table}, in the same format as
#'       \code{ps}. Non-\code{NULL} for \code{"lpdid_adj"} and
#'       \code{"lpdid_dr"} only.}
#'     \item{\code{plot}}{A \code{ggplot2} event study plot with 95\%
#'       confidence intervals, built from \code{se_correct} where it has
#'       been computed and from the naive \code{se} otherwise.}
#'   }
#'
#' @details
#' As of this version, \code{meth = "lpdid_adj"} uses a control-only outcome
#' regression (fit on untreated observations) plus a cluster-level combined
#' influence function for its standard error, replacing the previous
#' \code{marginaleffects::avg_comparisons()}-based implementation. This
#' correctly accounts for the cross-cluster covariance between the
#' control-fitting stage and the treated-outcome stage, which a naive
#' additive \code{V_naive + V_extra} formula would omit whenever the same
#' cluster (e.g. the same unit) contributes both a clean-control observation
#' and its own treated observation within a horizon's estimation sample --
#' the typical case in a staggered-adoption panel. \code{one_reg = TRUE} is
#' no longer supported for \code{"lpdid_adj"}, since the corresponding
#' influence-function derivation for the stacked one-regression case has not
#' been validated; use \code{one_reg = FALSE}.
#'
#' \code{meth = "lpdid_dr"} combines this same outcome regression with an
#' inverse-probability-weighted control-arm and treated-arm score (an
#' AIPW-style estimator), together with a combined influence function that
#' propagates uncertainty from both the outcome regression and the
#' propensity-score model. It is only implemented for \code{one_reg = FALSE}.
#'
#' In both cases, if \code{FE} is supplied (fixed effects absorbed via
#' \code{fixest}'s \code{|} syntax) in the outcome-regression formula, the
#' standard-error correction is a partial correction only: the absorbed
#' fixed effects' own sampling uncertainty is not propagated, since
#' \code{fixest}'s \code{coef()} for an absorbed-FE model excludes them.
#' Clustering for the correction uses only the first element of
#' \code{clusters} (or \code{unit} if \code{clusters} is \code{NULL}); no
#' multi-way clustering is supported by the correction itself, matching the
#' single-cluster-variable limitation of \code{lpdid.ado} that this
#' correction was validated against.
#'
#' \code{meth \%in\% c("lpdid_rw","lpdid_ipw")}: the final WLS uses weights
#' \eqn{w_i = b_i / (1 - q_{t(i)})} (up to a global constant, which does not
#' affect the estimate), where \eqn{q_t} is the \eqn{b}-weighted share of
#' treated observations in period \eqn{t} and \eqn{b_i} is the survey weight
#' (\code{"lpdid_rw"}) or the survey weight times the ATT-IPW weight
#' (1 for treated, \eqn{\hat p/(1-\hat p)} for controls; \code{"lpdid_ipw"}).
#' The naive \code{se} treats \eqn{w_i} as known. With
#' \code{se_correction_rw_ipw = TRUE}, \code{se_correct} adds to each
#' observation's WLS score \eqn{s_i} the first-order effect of estimating
#' the nuisance parameters: \eqn{G_t\, b_j (D_j - q_t)/\sum_{t} b} with
#' \eqn{G_t = \sum_{i \in t} s_i / (1-q_t)} for the period shares, and, for
#' IPW, \eqn{J\, \mathrm{IF}_\gamma(j)} for the logit coefficients, where
#' \eqn{J = \partial \sum_i s_i / \partial \gamma} (through both the control
#' odds \eqn{e^{x'\gamma}} and the shares \eqn{q_t}). The combined terms are
#' summed by cluster and passed through the WLS inverse Hessian. To keep
#' \code{fixest}'s small-sample adjustment, \code{se_correct} is reported as
#' \code{se * sqrt(V_corrected / V_uncorrected)}, both variances being
#' computed with the same raw cluster-sum formula (the uncorrected one
#' reproduces \code{fixest}'s clustered SE up to that adjustment). Trimming
#' is treated as fixed. The same partial-correction caveats apply for
#' absorbed \code{FE} in the PS model; in the PS, time now enters as explicit
#' dummies (\code{i(time)}) rather than an absorbed fixed effect, so that its
#' uncertainty is propagated (identical point estimates).
#'
#' @references
#'   Dube, A., Girardi, D., Jorda, O., & Taylor, A. M. (2023).
#'   A local projections approach to difference-in-differences.
#'   \emph{NBER Working Paper 31184}.
#'
#'   Callaway, B., & Sant'Anna, P. H. C. (2021).
#'   Difference-in-differences with multiple time periods.
#'   \emph{Journal of Econometrics}, 225(2), 200-230.
#'
#'   Sant'Anna, P. H. C., & Zhao, J. (2020).
#'   Doubly robust difference-in-differences estimators.
#'   \emph{Journal of Econometrics}, 219(1), 101-122.
#'
#' @examples
#' \donttest{
#' data(mpdta_r, package = "lpdidcsa")
#'
#' # Data preparation (with wide format horizons)
#' suppressWarnings({df  <- lpdidcsa_data(mpdta_r, unit = "countyreal", time = "year",
#'                      dependent = "lemp", treat = "treat")})
#'
#'
#' # 1. Estimation with lpdid_rw and no covariates
#' res <- lpdidcsa(df, unit = "countyreal", time = "year",
#'                 dependent = "lemp", meth = "lpdid_rw")
#'
#' # Horizon estimates
#' print(res$plot)
#' print(res$est)
#'
#'
#' # 2. Estimation with lpcsa and no covariates
#' res <- lpdidcsa(df, unit = "countyreal", time = "year",
#'                 dependent = "lemp", meth = "lpcsa")
#'
#' # Horizon estimates
#' print(res$plot)
#' print(res$est)
#'
#' # Cohort * horizon estimates
#' print(res$est_det)
#'
#'
#' # 3. Estimation with lpdid_ipw and covariates
#' res <- lpdidcsa(df, unit = "countyreal", time = "year",
#'                 dependent = "lemp", meth = "lpdid_ipw",controls="lpop")
#'
#' # Horizon estimates
#' print(res$plot)
#' print(res$est)
#'
#' # First stage propensity score estimates
#' print(res$ps)
#'
#' # 4. Estimation with lpcsa_ipw and covariates
#' res <- lpdidcsa(df, unit = "countyreal", time = "year",
#'                 dependent = "lemp", meth = "lpcsa_ipw",controls="lpop")
#'
#' # Horizon estimates
#' print(res$plot)
#' print(res$est)
#'
#' # Cohort * horizon estimates
#' print(res$est_det)
#'
#' # First stage propensity score estimates
#' print(res$ps)
#'
#' # 5. Estimation with lpdid_adj (regression adjustment) and covariates,
#' #    with the corrected standard error
#' res <- lpdidcsa(df, unit = "countyreal", time = "year",
#'                 dependent = "lemp", meth = "lpdid_adj", controls = "lpop")
#' print(res$est)   # includes both se (naive) and se_correct
#' print(res$or)    # first-step outcome regression on the control group
#'
#' # 6. Same, but skipping the SE correction for speed on large data
#' res <- lpdidcsa(df, unit = "countyreal", time = "year",
#'                 dependent = "lemp", meth = "lpdid_adj", controls = "lpop",
#'                 se_correction = FALSE)
#' print(res$est)   # se_correct is NA
#'
#' # 7. Estimation with lpdid_dr (doubly robust) and covariates
#' res <- lpdidcsa(df, unit = "countyreal", time = "year",
#'                 dependent = "lemp", meth = "lpdid_dr", controls = "lpop")
#' print(res$est)
#' print(res$ps)
#' print(res$or)
#'
#' # 8. Same, with covariate effects allowed to vary by period
#' res <- lpdidcsa(df, unit = "countyreal", time = "year",
#'                 dependent = "lemp", meth = "lpdid_dr", controls = "lpop",
#'                 time_x_controls = TRUE)
#' print(res$ps)    # dtreat ~ factor(year)*(lpop)
#' print(res$or)
#'
#' }
#'
#' @import data.table fixest ggplot2
#' @export
lpdidcsa <- function(data,
                     dtreat         = "dtreat",
                     time           = "year",
                     unit           = NULL,
                     dependent      = NULL,
                     n_pre         = NULL,
                     n_post        = NULL,
                     controls       = NULL,
                     controls_h     = NULL,
                     FE             = NULL,
                     FE_h           = NULL,
                     clusters       = NULL,
                     clusters_h     = NULL,
                     weight         = NULL,
                     weight_h       = NULL,
                     meth           = "lpdid_ipw",
                     anticipation   = 0,
                     absorbing      = TRUE,
                     reentry        = NULL,
                     type_horizon   = "wide",
                     horizon        = "horizon",
                     one_reg        = FALSE,
                     se_correction  = TRUE,
                     se_correction_rw_ipw = FALSE,
                     ps_trim        = 0.995,
                     time_x_controls = FALSE)
{
  # ── Input validation ───────────────────────────────────────────────────────
  valid_meths <- c("lpdid", "lpdid_rw", "lpdid_adj", "lpdid_dr",
                   "lpdid_ipw", "lpcsa_ipw", "lpcsa")
  stopifnot(
    is.data.frame(data),
    is.character(unit)      && length(unit)      == 1,
    is.character(time)      && length(time)      == 1,
    is.character(dependent) && length(dependent) == 1,
    is.character(dtreat)    && length(dtreat)    == 1,
    is.logical(absorbing)   && length(absorbing) == 1,
    is.logical(se_correction) && length(se_correction) == 1,
    is.logical(se_correction_rw_ipw) && length(se_correction_rw_ipw) == 1,
    is.numeric(ps_trim) && length(ps_trim) == 1 && ps_trim > 0 && ps_trim <= 1,
    is.logical(time_x_controls) && length(time_x_controls) == 1,
    meth %in% valid_meths
  )

  if (isTRUE(time_x_controls) && isTRUE(one_reg)) {
    warning("time_x_controls is only implemented for one_reg = FALSE; ignored.")
    time_x_controls <- FALSE
  }

  if (isTRUE(se_correction_rw_ipw) && isTRUE(one_reg) &&
      meth %in% c("lpdid_rw", "lpdid_ipw"))
    warning("se_correction_rw_ipw is only implemented for one_reg = FALSE; ",
            "ignored.")

  if (meth == "lpdid_dr" && isTRUE(one_reg))
    stop("meth = 'lpdid_dr' with one_reg = TRUE is not implemented. ",
         "The stacked-horizon influence-function derivation for DR has not ",
         "been validated; use one_reg = FALSE.")

  if (meth == "lpdid_adj" && isTRUE(one_reg))
    stop("meth = 'lpdid_adj' with one_reg = TRUE is no longer supported. ",
         "The regression-adjustment engine was replaced with a validated ",
         "influence-function correction that has only been derived for the ",
         "per-horizon (one_reg = FALSE) case. Use one_reg = FALSE.")

  if (meth=="lpdid_ipw" &
      is.null(FE) & is.null(FE_h) & is.null(controls) & is.null(controls_h)) {meth <- "lpdid_rw"}
  if (meth=="lpcsa_ipw" &
      is.null(FE) & is.null(FE_h) & is.null(controls) & is.null(controls_h)) {meth <- "lpcsa"}


  # ── Helpers ────────────────────────────────────────────────────────────────

  # Extract bare column names from fixest-style formula strings
  # e.g. "i(x) + log(y)" -> c("x", "y")
  parse_col_names <- function(x) {
    x <- gsub("[a-zA-Z][a-zA-Z0-9_]*\\(([^)]+)\\)", "\\1", x)
    x <- unlist(strsplit(x, "[*+:^|[:space:]]+"))  # [:space:]: "\\s" is a literal "s" inside [...]
    x <- x[grepl("^[a-zA-Z\\.][a-zA-Z0-9\\._]*$", x)]
    unique(x)
  }

  required_cols <- c(unit, time, paste0(dependent, "_tm1"), dtreat,
                     parse_col_names(controls),
                     parse_col_names(FE),
                     clusters, weight)

  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0)
    stop("Columns missing from 'data': ", paste(missing_cols, collapse = ", "))

  if (paste0(dependent,"_th") %in% colnames(data) & type_horizon=="wide") {
    warning("Horizons seem to be in the long format. type_horizon changed to 'long'")
    type_horizon="long"
  }
  if ((paste0(dependent,"_th") %in% colnames(data))==F & type_horizon=="long") {
    warning("Horizons seem to be in the wide format. type_horizon changed to 'wide'")
    type_horizon="wide"
  }
  # ── Setup ──────────────────────────────────────────────────────────────────
  df <- setDT(data)

  col_unit    <- unit
  col_time    <- time
  col_dep     <- dependent
  col_dtreat  <- dtreat
  col_clusters <- if (!is.null(clusters)) clusters else unit
  col_weight  <- weight
  col_dep_ref <- paste0(col_dep, "_tm",1+anticipation)
  col_horizon <- horizon

  # Backup of variables for managing h_variables
  b_controls <- controls
  b_col_clusters <- col_clusters
  b_FE <- FE

  # Determine the available horizon range from the data
  if (type_horizon == "wide") {
    col_dep_all <- colnames(df)[substr(colnames(df), 1, nchar(dependent) + 1) ==
                                  paste0(dependent, "_")]
    time_dep    <- as.numeric(gsub("m", "-",
                                   substr(col_dep_all,
                                          nchar(dependent) + 3,
                                          nchar(col_dep_all))))
    max_n_pre  <- -min(time_dep)
    max_n_post <- max(time_dep)

  } else if (type_horizon == "long") {
    max_n_pre  <- -df[, min(get(col_horizon), na.rm = TRUE)]
    max_n_post <-  df[, max(get(col_horizon), na.rm = TRUE)]
  }

  # Cap n_pre / n_post at what the data supports
  if (is.null(n_post) | !is.numeric(n_post)) {
    n_post <- max_n_post
  } else if (n_post == floor(n_post)) {
    n_post <- pmin(n_post, max_n_post)
  } else {
    n_post <- max_n_post
  }

  if (is.null(n_pre) | !is.numeric(n_pre)) {
    n_pre <- max_n_pre
  } else if (n_pre == floor(n_pre)) {
    n_pre <- pmin(max_n_pre, abs(n_pre))
  } else {
    n_pre <- max_n_pre
  }

  horizons <- setdiff(-abs(n_pre):n_post, -1L-anticipation)

  # ── Inner helpers ──────────────────────────────────────────────────────────

  # Outcome column name at horizon h
  # Wide: e.g. h = -2 -> "lemp_tm2", h = 3 -> "lemp_t3"
  # Long: always "lemp_th" (the stacked horizon column)
  dep_at <- function(h) {
    if (type_horizon == "wide") {
      paste0(col_dep, "_t", gsub("-", "m", h))
    } else {
      paste0(col_dep, "_th")
    }
  }

  var_at <- function(h,col_var=NULL) {
    if (type_horizon == "wide") {
      paste0(col_var, "_t", gsub("-", "m", h))
    } else {
      paste0(col_var, "_th")
    }
  }

  # LHS of regression: outcome change relative to t-1
  lhs <- function(h) {
    paste0("(", dep_at(h), " - ", col_dep_ref, ")")
  }

  # Filter rows to the clean DiD comparison set for horizon h.
  # Treated units are always kept; control units are kept only if they are
  # not yet treated (absorbing) or satisfy the non-absorbing entry conditions.
  clean_sample <- function(df, h) {
    col_dep_th  <- dep_at(h)
    base_filter <- !is.na(df[[col_dep_ref]]) & !is.na(df[[col_dep_th]]) & !is.na(df[[col_dtreat]])

    if (absorbing) {
      keep <- (df[[col_dtreat]] == 1 |
                 is.na(df$first.treat) | #never treated
                 (!is.na(df$first.treat)
                  & df[[col_time]] +  (h + anticipation) * (h + anticipation > 0 ) <  df$first.treat )
      )
    } else if (is.null(reentry)) {
      keep <- (is.na(df$last.treat) &
                 (is.na(df$next.treat) |
                    (!is.na(df$next.treat) &
                       df[[col_time]] + (h + anticipation) * (h + anticipation > 0) < df$next.treat )))
    } else {
      keep <- ((is.na(df$last.treat) |
                  df[[col_time]] - df$last.treat > reentry) &
                 (is.na(df$next.treat) |
                    (!is.na(df$next.treat) &
                       df[[col_time]] + (h + anticipation) * (h + anticipation > 0) < df$next.treat )))
    }

    if (type_horizon == "long") {
      df[keep & base_filter & get(col_horizon) == h]
    } else {
      df[keep & base_filter]
    }
  }

  # Same filtering logic for long format when one_reg = TRUE (all horizons at once)
  clean_sample_long <- function(df) {

    col_dep_th  <- dep_at(h)
    base_filter <- !is.na(df[[col_dep_ref]]) & !is.na(df[[col_dep_th]]) & !is.na(df[[col_dtreat]])

    if (absorbing) {
      keep <- (df[[col_dtreat]] == 1 |
                 is.na(df$first.treat) |
                 (!is.na(df$first.treat)
                  & df[[col_time]] + (df[[col_horizon]] + anticipation) * (df[[col_horizon]] + anticipation > 0) < df$first.treat))
    } else if (is.null(reentry)) {
      keep <- (is.na(df$last.treat) &
                 (is.na(df$next.treat) |
                    (!is.na(df$next.treat) &
                       df[[col_time]] + (df[[col_horizon]] + anticipation) * (df[[col_horizon]] + anticipation > 0) < df$next.treat)))
    } else {
      keep <- ((is.na(df$last.treat) |
                  df[[col_time]] - df$last.treat > reentry) &
                 (is.na(df$next.treat) |
                    (!is.na(df$next.treat) &
                       df[[col_time]] + (df[[col_horizon]] + anticipation) * (df[[col_horizon]] + anticipation > 0 ) < df$next.treat)))
    }
    df[keep & base_filter]
  }

  # Drop missing observations from the data
  drop_missing <- function(df,model){
    if (length(model$obs_selection$obsRemoved)>0) {
      df[model$obs_selection$obsRemoved,]
    } else {df}
  }

  # Build a fixest formula string from its components
  build_formula <- function(lhs, rhs_main, controls, FE, time_fe = col_time) {
    rhs <- if (!is.null(controls))
      paste(rhs_main, "+", paste(controls, collapse = " + "))
    else
      rhs_main
    fe_part <- if (!is.null(FE))
      paste(FE, collapse = " + ")
    else
      time_fe
    if (!is.null(fe_part))
      paste0(lhs, " ~ ", rhs, " | ", fe_part)
    else
      paste0(lhs, " ~ ", rhs)
  }

  # Split a character vector of formula pieces into top-level '+' terms
  # (a '+' inside parentheses, e.g. "log(a + b)", is not a split point).
  split_terms <- function(x) {
    chars <- strsplit(paste(x, collapse = " + "), "")[[1]]
    depth <- 0L; cur <- ""; out <- character(0)
    for (ch in chars) {
      if (ch == "(") depth <- depth + 1L
      if (ch == ")") depth <- depth - 1L
      if (ch == "+" && depth == 0L) { out <- c(out, trimws(cur)); cur <- "" }
      else cur <- paste0(cur, ch)
    }
    out <- c(out, trimws(cur))
    out[nzchar(out)]
  }

  # time_x_controls: interact controls and absorbed FE with time, for the
  # regression whose estimation sample is `data`. Controls are interacted
  # with factor(time) using base-R formula operators, so any term R can
  # expand (log(x), poly(x, 2), a*b, factor(v), ...) is interacted as a whole:
  #   time absorbed elsewhere   -> (controls) + factor(time):(controls)
  #       (final regressions with "| time" or "FE^time"; since the controls'
  #        main effects are present, R contrast-codes factor(time) in the
  #        interaction: T-1 slope shifts per control, nothing collinear with
  #        the absorbed time FE)
  #   time not otherwise in the model -> factor(time)*(controls)
  #       (PS logit and control-group outcome regression without FE; the
  #        factor(time) main effect then replaces the separate time term,
  #        signalled by time_in_controls = TRUE)
  #   FE term f -> f^time  (fixest FE syntax; time itself is left as is)
  # fixest's i() is not used here: inside base-R interaction operators it is
  # fragile, and its two-argument form fails when fixest is imported but not
  # attached. Simple factor terms i(v) are rewritten as factor(v); i() terms
  # with further arguments (ref = ..., two variables) are not translatable
  # and raise an error asking for factor() syntax.
  # Only when `data` has >= 2 periods; otherwise everything is returned
  # unchanged.
  #
  #   time_absorbed : TRUE when the calling regression already absorbs time
  #                   (final regressions). For the PS / OR builders pass
  #                   FALSE: time is then absorbed only if FE are present
  #                   (FE^time).
  # Returns list(controls, FE, time_in_fe, time_in_controls, active),
  # controls being a single formula string (or NULL); active = TRUE when the
  # interaction was applied.
  tx_apply <- function(controls, FE, data, time_absorbed = TRUE) {
    res <- list(controls = controls, FE = FE, time_in_fe = FALSE,
                time_in_controls = FALSE, active = FALSE)
    if (!isTRUE(time_x_controls)) return(res)
    # i(v) -> factor(v) whenever time_x_controls is requested, also for
    # single-period samples, so that coefficient names and the i()-free
    # formulas are the same at every horizon.
    if (!is.null(controls)) {
      ctrl <- paste(controls, collapse = " + ")
      ctrl <- gsub("\\bi\\(\\s*([a-zA-Z.][a-zA-Z0-9._]*)\\s*\\)", "factor(\\1)", ctrl, perl = TRUE)
      if (grepl("\\bi\\(", ctrl, perl = TRUE))
        stop("time_x_controls: controls contain an i() term with arguments ('",
             ctrl, "'), which cannot be interacted with ", col_time,
             ". Write factor terms as factor(var) instead.")
      res$controls <- ctrl
    }
    if (length(unique(data[[col_time]])) < 2) return(res)
    res$active <- TRUE
    if (!is.null(FE)) {
      fe_terms <- split_terms(FE)
      res$FE <- ifelse(fe_terms == col_time, fe_terms,
                       paste0(fe_terms, "^", col_time))
      res$time_in_fe <- TRUE
    }
    if (!is.null(controls)) {
      ctrl <- res$controls
      ft <- paste0("factor(", col_time, ")")
      if (time_absorbed || res$time_in_fe) {
        res$controls <- paste0("(", ctrl, ") + ", ft, ":(", ctrl, ")")
      } else {
        res$controls <- paste0(ft, "*(", ctrl, ")")
        res$time_in_controls <- TRUE
      }
    }
    res
  }

  # Propensity-score logit formula (lpdid_ipw, lpcsa_ipw, lpdid_dr).
  # Time enters as explicit dummies i(time) -- not absorbed -- so that its
  # uncertainty can be propagated by the SE corrections; dropped when the
  # sample has a single period, when time_x_controls interacts the absorbed
  # FE with time (FE^time then absorbs it), or when it already enters as the
  # factor(time) main effect of factor(time)*(controls).
  build_ps_formula <- function(controls, FE, data) {
    multi_time <- length(unique(data[[col_time]])) > 1
    tx <- tx_apply(controls, FE, data, time_absorbed = FALSE)
    time_term <- if (multi_time && !tx$time_in_fe && !tx$time_in_controls)
      paste0("i(", col_time, ")") else NULL   # never combined with tx terms
    rhs <- c(tx$controls, time_term)
    f <- paste0(col_dtreat, " ~ ",
                if (length(rhs) > 0) paste(rhs, collapse = " + ") else "1")
    if (!is.null(tx$FE)) f <- paste0(f, " | ", paste(tx$FE, collapse = " + "))
    f
  }

  # Build the control-only outcome-regression formula used by lpdid_adj / lpdid_dr.
  # No treatment indicator on the RHS -- this is fit on untreated observations only.
  #
  # IMPORTANT: when FE is NULL, time enters as an EXPLICIT covariate
  # (factor(col_time)), not as an absorbed fixed effect. This matches the
  # validated standalone RA script this correction was checked against
  # (which used "factor(year) + lpop", not "lpop | year"). Absorbing time by
  # default would silently exclude its sampling uncertainty from
  # ra_combined_se/dr_combined_se's correction, since coef() on a model with
  # absorbed FE does not include the FE's own coefficients.
  #
  # Because time is now explicit rather than absorbed, a horizon whose clean
  # sample spans only one time period must drop the time term entirely
  # (a single-level factor is either degenerate or causes an error) -- same
  # guard as the user's own "uniqueN(clean_df$year) > 1" check in the
  # validated standalone script. Absorption-based FE (the FE-supplied
  # branch) doesn't need this: demeaning by a single group is a no-op.
  #
  # `ctrl_df` is the sample the regression is fit on (control observations);
  # the single-period guard and the time_x_controls check use it. With
  # time_x_controls, the RHS becomes factor(time)*(controls) when FE is NULL,
  # and (controls) + factor(time):(controls) | FE^time otherwise.
  build_control_formula <- function(h, controls, FE, ctrl_df) {
    multi_time <- length(unique(ctrl_df[[col_time]])) > 1
    tx <- tx_apply(controls, FE, ctrl_df, time_absorbed = FALSE)
    controls <- tx$controls
    FE <- tx$FE
    rhs_terms <- if (is.null(FE) && multi_time && !tx$time_in_controls)
      paste0("factor(", col_time, ")") else character(0)
    if (!is.null(controls)) rhs_terms <- c(rhs_terms, controls)
    rhs <- if (length(rhs_terms) > 0) paste(rhs_terms, collapse = " + ") else "1"
    if (!is.null(FE))
      paste0(lhs(h), " ~ ", rhs, " | ", paste(FE, collapse = " + "))
    else
      paste0(lhs(h), " ~ ", rhs)
  }

  # Model matrix of a fitted fixest model on new data, with exactly the
  # columns `cols` (the coefficient names) in that order. fixest's
  # model.matrix(newdata) omits an i(time, x) column when that level is
  # absent from `data` (e.g. a year with no treated rows under
  # time_x_controls); such a column is identically zero on `data`, so it is
  # zero-filled here rather than causing a subscript error.
  mm_aligned <- function(model, data, cols) {
    mm <- model.matrix(model, data = data)
    if (nrow(mm) != nrow(data))
      stop("model.matrix() dropped rows of the supplied data (", nrow(mm),
           " vs ", nrow(data), "); check for missing values in the controls.")
    out <- matrix(0, nrow(mm), length(cols), dimnames = list(NULL, cols))
    common <- intersect(cols, colnames(mm))
    out[, common] <- mm[, common, drop = FALSE]
    out
  }

  # Resolve which rows of `data_used` (the exact data.frame/data.table that
  # was passed to a fitted model) actually fed the fit, using the row names
  # sandwich::estfun() attaches. This is more robust than obs_selection$obsRemoved
  # alone: obs_selection$obsRemoved only tracks NA-driven drops, NOT
  # observations feols/feglm silently drops for other reasons (most commonly
  # a singleton fixed-effect group -- a level with a single observation,
  # perfectly explained by its own FE dummy). Relying on obs_selection alone
  # can leave the influence function shorter than the clustering vector it's
  # paired with, which surfaces as a tapply() "arguments must have same
  # length" error with no indication of the actual cause.
  resolve_used_rows <- function(estfun_mat, data_used) {
    kept <- rownames(estfun_mat)
    if (is.null(kept)) {
      if (nrow(estfun_mat) != nrow(data_used))
        stop("Cannot align the model's estimation sample with the supplied ",
             "data: estfun() returned no row names and row counts differ (",
             nrow(estfun_mat), " vs ", nrow(data_used), "). This can happen ",
             "when feols/feglm silently drops observations (e.g. singleton ",
             "fixed-effect groups).")
      return(seq_len(nrow(data_used)))
    }
    as.integer(kept)
  }

  # Combined cluster-level influence-function SE for regression adjustment (RA).
  # Validated to reproduce lpdid.ado's listreg-based RA standard error exactly
  # on a test panel (mpdta). Cluster-combines the control-fitting stage's
  # propagated uncertainty with the treated-arm's own residual noise, which
  # correctly captures the cross-cluster covariance a naive additive
  # V_naive + Xbar1'V_beta*Xbar1 formula would omit.
  #
  # ctrl_df_h must be the EXACT data.table object that was passed to
  # fit_feols() to produce control_reg (same row order), so that estfun()'s
  # row names can be used to recover exactly which rows were kept.
  #
  # NOTE: assumes control_reg has NO absorbed fixed effects contributing
  # uncertainty beyond what coef(control_reg) captures -- if FE is supplied,
  # this is a partial correction only (absorbed FE's own sampling
  # uncertainty is not propagated, since fixest's coef() excludes it).
  # Single-way clustering only (cl_col), matching lpdid.ado's own limitation.
  ra_combined_se <- function(control_reg, ctrl_df_h, trt_df, att_hat, e1, wgt_trt, cl_col) {
    bread_c <- sandwich::bread(control_reg) / control_reg$nobs   # undo sandwich's n-scaling
    score_beta_ctrl <- sandwich::estfun(control_reg)
    # sandwich::estfun() does not reliably attach column names for fixest
    # objects -- set them explicitly from coef(), as already done for ps_mod
    # elsewhere in this file. Without this, colnames(score_beta_ctrl) is
    # NULL, and matrix[, NULL, drop=FALSE] silently returns a matrix with
    # ZERO columns rather than "all columns" -- collapsing Xbar1 to
    # numeric(0) and breaking the correction downstream.
    stopifnot(ncol(score_beta_ctrl) == length(coef(control_reg)))
    colnames(score_beta_ctrl) <- names(coef(control_reg))
    IF_beta_c <- score_beta_ctrl %*% bread_c

    used_local <- resolve_used_rows(score_beta_ctrl, ctrl_df_h)
    ctrl_df_used <- ctrl_df_h[used_local, ]
    stopifnot(nrow(ctrl_df_used) == nrow(score_beta_ctrl))

    Xbeta_names <- colnames(score_beta_ctrl)
    X1 <- mm_aligned(control_reg, trt_df, Xbeta_names)
    N1w <- sum(wgt_trt)
    Xbar1 <- colSums(X1 * wgt_trt) / N1w

    psi_ctrl <- -as.numeric(IF_beta_c %*% Xbar1)
    psi_trt  <- as.numeric(wgt_trt * (e1 - att_hat) / N1w)

    S <- tapply(c(psi_ctrl, psi_trt),
                c(ctrl_df_used[[cl_col]], trt_df[[cl_col]]), sum)
    sqrt(sum(S^2))
  }

  # Combined cluster-level standard error for the doubly robust (DR)
  # estimator. Propagates uncertainty from BOTH the outcome regression
  # (control_reg) and the propensity score model (ps_mod) into the final
  # AIPW-style ATT's standard error, following the standard asymptotic-
  # linear-representation (ALR) convention (Sant'Anna & Zhao 2020): the
  # combined per-observation term att_if is treated as the estimator's ALR,
  # and its cluster-robust standard error is obtained the standard way --
  # via feols(att_if ~ 1, cluster = ...)'s own SE, i.e. the SE of att_if's
  # sample MEAN. This matches the user's validated standalone DR script
  # exactly.
  #
  # NOTE: this is a DIFFERENT aggregation convention from ra_combined_se,
  # which sums (not means) the per-observation terms within cluster with no
  # further N-normalization -- that convention is specific to replicating
  # lpdid.ado's listreg implementation (validated bit-for-bit against
  # Stata) and does not generalize to DR, which has no Stata analog to
  # replicate. bread() is likewise used RAW here (not divided by nobs),
  # matching the reference script -- unlike ra_combined_se's bread, which
  # IS divided by nobs to match listreg's specific internal construction.
  #
  # NOTE: same absorbed-FE and single-way-clustering caveats as
  # ra_combined_se above. ps_mod is assumed fit on clean_df directly (no
  # further remapping needed there under normal circumstances), but a guard
  # is still applied in case ps_mod also silently drops rows.
  dr_combined_se <- function(control_reg, ctrl_full_idx, ps_mod, clean_df,
                             w_treat, w_cont, dy, mu0, eta0, eta1, cl_col) {
    n <- nrow(clean_df)

    bread_beta <- sandwich::bread(control_reg)
    score_beta_ctrl <- sandwich::estfun(control_reg)
    # Same fix as in ra_combined_se: estfun() doesn't reliably name columns
    # for fixest objects, so set them explicitly from coef().
    stopifnot(ncol(score_beta_ctrl) == length(coef(control_reg)))
    colnames(score_beta_ctrl) <- names(coef(control_reg))
    common_or <- colnames(score_beta_ctrl)

    used_local <- resolve_used_rows(score_beta_ctrl, data.frame(row = seq_along(ctrl_full_idx)))
    ctrl_used_idx <- ctrl_full_idx[used_local]
    stopifnot(length(ctrl_used_idx) == nrow(score_beta_ctrl))

    score_beta <- matrix(0, n, ncol(score_beta_ctrl))
    score_beta[ctrl_used_idx, ] <- score_beta_ctrl
    asy_lin_rep_beta <- score_beta %*% bread_beta

    bread_ps <- sandwich::bread(ps_mod)
    score_ps_raw <- sandwich::estfun(ps_mod)
    colnames(score_ps_raw) <- names(coef(ps_mod))
    common_ps <- colnames(score_ps_raw)

    ps_used_local <- resolve_used_rows(score_ps_raw, clean_df)
    if (length(ps_used_local) == n && identical(ps_used_local, seq_len(n))) {
      asy_lin_rep_ps <- score_ps_raw %*% bread_ps
    } else {
      asy_lin_rep_ps <- matrix(0, n, ncol(score_ps_raw))
      asy_lin_rep_ps[ps_used_local, ] <- score_ps_raw %*% bread_ps
    }

    X_full_or <- mm_aligned(control_reg, clean_df, common_or)
    X_full_ps <- mm_aligned(ps_mod, clean_df, common_ps)

    Xbar1 <- colSums(w_treat * X_full_or) / sum(w_treat)
    Xbar0 <- colSums(w_cont  * X_full_or) / sum(w_cont)
    M2    <- colSums(w_cont * ((dy - mu0) - eta0) * X_full_ps) / n

    or_correction <- as.numeric(asy_lin_rep_beta %*% (Xbar1 - Xbar0))
    ps_correction <- as.numeric(asy_lin_rep_ps %*% M2) / mean(w_cont)

    inf_treat <- (w_treat * (dy - mu0) - w_treat * eta1) / mean(w_treat)
    inf_cont  <- (w_cont  * (dy - mu0) - w_cont  * eta0) / mean(w_cont)

    att_if <- inf_treat - inf_cont - or_correction - ps_correction

    se_df <- data.frame(att_if = att_if, .cl = clean_df[[cl_col]])
    se_mod <- fixest::feols(att_if ~ 1, cluster = ~.cl, data = se_df)
    unname(se_mod$coeftable[1, "Std. Error"])
  }



  # Combined cluster-level influence-function SE for the reweighted (rw) and
  # IPW LP-DiD estimators.
  #
  # Both estimators run a WLS of the outcome change on dtreat (+ controls,
  # FE) with weights w_i = b_i * S / (1 - q_t(i)), where q_t is the
  # b-weighted treated share in period t (first stage of add_rewgt) and
  #   b_i = survey weight (or 1)                       for lpdid_rw
  #   b_i = survey weight * (1 if treated, exp(x'g))   for lpdid_ipw
  # (exp(x'g) = ps/(1-ps) under the logit). The global constant S does not
  # affect the estimate. With the WLS score sum S_beta(eta) = sum_i s_i:
  #   beta_hat - beta ~ H^-1 [ sum_i s_i + dS_beta/deta (eta_hat - eta) ].
  # Since dS_beta/dw_i = s_i / w_i (the derivatives through FE demeaning
  # and the FE coefficients vanish by the FE first-order conditions):
  #   * period shares: dS/dq_t = G_t = sum_{i in t} s_i / (1 - q_t), and
  #     q_hat_t - q_t ~ sum_{j in t} b_j (D_j - q_t) / B_t;
  #   * logit (ipw only): J = dS/dgamma
  #       = sum_i s_i (1-D_i) x_i'  +  sum_t G_t dq_t/dgamma,
  #     dq_t/dgamma = -(q_t / B_t) sum_{i in t, D=0} b_i x_i,
  #     with gamma_hat - gamma ~ sum_j H_ps^-1 sc_ps_j.
  # Per-observation terms are summed by cluster and premultiplied by
  # H^-1 = bread()/nobs (fixest convention, checked numerically). Returns the
  # raw corrected and uncorrected variances of the dtreat coefficient; the
  # caller rescales fixest's SE by sqrt(Vc / Vu) so the small-sample
  # adjustment is kept. Trimming is treated as fixed.
  #
  #   mod    : final feols fit on est_df (the same data add_rewgt was run on)
  #   bw     : base weight b over the rows of est_df
  #   ps_mod : logit (ipw only), fit on ps_df (rows after drop_missing,
  #            before trimming)
  rw_ipw_combined_se <- function(mod, est_df, bw, cl_col,
                                 ps_mod = NULL, ps_df = NULL) {
    cn <- names(coef(mod))
    k  <- grep(col_dtreat, cn, fixed = TRUE)[1]
    if (is.na(k)) stop("treatment coefficient not found in the final model")

    sc <- sandwich::estfun(mod)
    stopifnot(ncol(sc) == length(cn))
    Hinv <- sandwich::bread(mod) / mod$nobs

    n_e   <- nrow(est_df)
    used  <- fixest::obs(mod)            # rows of est_df used (0-weight rows dropped)
    stopifnot(length(used) == nrow(sc))
    score <- matrix(0, n_e, ncol(sc))
    score[used, ] <- sc

    D  <- est_df[[col_dtreat]]
    tt <- as.character(est_df[[col_time]])
    BW <- tapply(bw, tt, sum)
    Q  <- tapply(bw * D, tt, sum) / BW
    q_i  <- as.numeric(Q[tt])
    BW_i <- as.numeric(BW[tt])
    inv1q <- ifelse(q_i < 1 - 1e-10, 1 / (1 - q_i), 0)

    G   <- rowsum(score * inv1q, tt)     # periods x K
    G_i <- G[tt, , drop = FALSE]
    phi_e <- score + G_i * (bw * (D - q_i) / BW_i)
    phi <- phi_e
    cl  <- est_df[[cl_col]]

    if (!is.null(ps_mod)) {
      gn <- names(coef(ps_mod))
      Xe <- mm_aligned(ps_mod, est_df, gn)
      ctrl <- 1 - D
      J  <- crossprod(score, ctrl * Xe)  # direct effect through exp(x'g)
      dC <- rowsum(bw * ctrl * Xe, tt)   # periods x Kg
      dq <- -(as.numeric(Q[rownames(dC)]) / as.numeric(BW[rownames(dC)])) * dC
      J  <- J + crossprod(G[rownames(dC), , drop = FALSE], dq)

      sc_ps <- sandwich::estfun(ps_mod)
      # ps_df is already restricted to the PS estimation sample (drop_missing)
      stopifnot(nrow(sc_ps) == nrow(ps_df))
      IF_g  <- sc_ps %*% (sandwich::bread(ps_mod) / ps_mod$nobs)
      phi   <- rbind(phi, IF_g %*% t(J))
      cl    <- c(cl, ps_df[[cl_col]])
    }

    Vc <- sum((rowsum(phi, cl) %*% Hinv[, k])^2)
    Vu <- sum((rowsum(score, est_df[[cl_col]]) %*% Hinv[, k])^2)
    c(Vc = Vc, Vu = Vu)
  }

  # Coerce a feols/feglm coeftable to a standardised data.table (one horizon)
  export_coeftable <- function(model, h, formula_str) {
    out <- data.table(
      h        = h,
      variable = row.names(model$coeftable),
      model$coeftable,
      n_obs   = model$nobs,
      formula  = formula_str,
      clusters = attr(model$coeftable,"vcov_type")
    )
    out[,clusters:=gsub("Clustered ","",clusters)]
    setnames(out, c("h", "variable", "estimate", "se", "T", "pvalue",
                    "n_obs", "formula","clusters"))
    out
  }

  # Same as above for the one_reg = TRUE case: extracts h from the variable name
  export_coeftable_long <- function(model, formula_str) {
    out <- data.table(
      variable = row.names(model$coeftable),
      model$coeftable,
      n_obs   = model$nobs,
      formula  = formula_str,
      clusters = attr(model$coeftable,"vcov_type")
    )
    out[,clusters:=gsub("Clustered ","",clusters)]
    h_part <- ifelse(
      regexpr(paste0(col_horizon,"::"), out$variable) > 0,
      substr(out$variable,
             regexpr(paste0(col_horizon,"::"), out$variable) + nchar(col_horizon) + 2,
             nchar(out$variable))
      ,NA)
    h <- as.numeric(ifelse(regexpr(":", h_part) > 0,
                           substr(h_part,1,regexpr(":", h_part)-1),
                           h_part))

    out <- data.frame(h, out)
    setDT(out)
    setnames(out, c("h", "variable", "estimate", "se", "T", "pvalue",
                    "n_obs", "formula","clusters"))
    out
  }

  # Coerce aggregate() output to a standardised data.table (one horizon)
  export_agg <- function(model, model_det, h, formula_str) {
    out <- data.table(
      h        = h,
      variable = row.names(model),
      model,
      n_obs   = model_det$nobs,
      formula  = formula_str
    )
    setnames(out, c("h", "variable", "estimate", "se", "T", "pvalue",
                    "n_obs", "formula"))
    out
  }

  # Same as above for the one_reg = TRUE case
  export_agg_long <- function(model, model_det, formula_str) {
    out <- data.table(
      variable = row.names(model),
      model,
      n_obs   = model_det$nobs,
      formula  = formula_str
    )
    h_part <- ifelse(
      regexpr(paste0(col_horizon,"::"), out$variable) > 0,
      substr(out$variable,
             regexpr(paste0(col_horizon,"::"), out$variable) + nchar(col_horizon) + 2,
             nchar(out$variable))
      ,NA)
    h <- as.numeric(ifelse(regexpr(":", h_part) > 0,
                           substr(h_part,1,regexpr(":", h_part)-1),
                           h_part))

    out <- data.frame(h, out)
    setDT(out)
    setnames(out, c("h", "variable", "estimate", "se", "T", "pvalue",
                    "n_obs", "formula"))
    out
  }

  # Compute reweighting residuals for equally-weighted ATT (wide / one horizon)
  # Regresses dtreat on time FE; residuals are used as importance weights.
  add_rewgt <- function(df, h, col_w = NULL) {
    formula_rw <- as.formula(paste0(col_dtreat, " ~ 1 | ", col_time))
    wgtreg <- if (is.null(col_w))
      feols(formula_rw, data = df)
    else
      feols(formula_rw, weights = df[[col_w]], data = df)

    df[, num_wgt_h := wgtreg$residuals]
    df[, num_wgt_h := fifelse(is.na(get(col_dtreat)) | get(col_dtreat) == 0,
                              NA_real_, num_wgt_h)]
    df[, wgt_h     := num_wgt_h / sum(num_wgt_h, na.rm = TRUE)]
    df[, gwgt_h    := max(wgt_h, na.rm = TRUE), by = col_time]
    df[, wgt_h     := fifelse(is.na(wgt_h), gwgt_h, wgt_h)]
    df[, rewgt_h   := if (is.null(col_w)) 1 / wgt_h else get(col_w) / wgt_h]
    df[, c("num_wgt_h", "wgt_h", "gwgt_h") := NULL]
    df
  }

  # Same reweighting for long format / one_reg = TRUE
  # Uses horizon x time FE interaction in the first-stage regression.
  add_rewgt_long <- function(df, col_w = NULL) {
    formula_rw <- as.formula(paste0(col_dtreat, " ~ 1 | ",col_horizon,"^", col_time))
    wgtreg <- if (is.null(col_w))
      feols(formula_rw, data = df)
    else
      feols(formula_rw, weights = df[[col_w]], data = df)

    df[, num_wgt_h := wgtreg$residuals]
    df[, num_wgt_h := fifelse(is.na(get(col_dtreat)) | get(col_dtreat) == 0,
                              NA_real_, num_wgt_h)]
    df[, wgt_h     := num_wgt_h / sum(num_wgt_h, na.rm = TRUE)]
    df[, gwgt_h    := max(wgt_h, na.rm = TRUE), by = list(get(col_time), get(col_horizon))]
    df[, wgt_h     := fifelse(is.na(wgt_h), gwgt_h, wgt_h)]
    df[, rewgt_h   := if (is.null(col_w)) 1 / wgt_h else get(col_w) / wgt_h]
    df[, c("num_wgt_h", "wgt_h", "gwgt_h") := NULL]
    df
  }

  # Fit feols with optional weights and fixed clustering
  fit_feols <- function(formula, data, cl, col_w = NULL, ...) {
    # Select clustering columns based on data type
    if (inherits(data, "data.table")) {
      cluster_cols <- data[, ..cl]  # Use ..cl for data.table
    } else {
      cluster_cols <- data[, cl, drop = FALSE]  # Use standard syntax for data.frame
    }

    args <- list(
      fml     = as.formula(formula),
      data    = data,
      cluster = cluster_cols,
      ...
    )
    if (!is.null(col_w)) args$weights <- data[[col_w]]

    model <- do.call(feols, args)

    model
  }
  # Fit feglm (logit) with optional weights — used for propensity score models
  fit_feglm <- function(formula, data, col_w = NULL) {
    args <- list(
      fml      = as.formula(formula),
      data     = data,
      fixef.rm = "none",
      family   = binomial
    )
    if (!is.null(col_w)) args$weights <- data[[col_w]]
    do.call(feglm, args)
  }


  # Attach treatment year and number of treated units to CSA cohort-level output
  add_meta <- function(dt, clean_df) {

    n_treat <- clean_df[get(col_dtreat) == 1, .N]
    n_ctrl <- clean_df[get(col_dtreat) == 0, .N]
    n_t <- data.frame(n_treat, n_ctrl)

    if (!is.null(col_weight)){
      sw_treat <- clean_df[get(col_dtreat) == 1, sum(get(col_weight))]
      sw_ctrl <- clean_df[get(col_dtreat) == 0, sum(get(col_weight))]
      n_t <- data.frame(n_t, sw_treat, sw_ctrl)
    }
    dt <- setDT(data.frame(dt, n_t))   # setDT: ':=' below needs a data.table
    move_cols(dt,c("n_treat","n_ctrl"),before="formula")
    if (!is.null(col_weight)){
      move_cols(dt,c("sw_treat","sw_ctrl"),before="formula")
      dt[,weight:=col_weight]
    }
    dt
  }

  add_meta_long <- function(dt, clean_df) {

    n_t <- clean_df[get(col_dtreat) == 1, .N, by = c(col_horizon)]
    n_c <- clean_df[get(col_dtreat) == 0, .N, by = c(col_horizon)]
    setnames(n_t, c(col_horizon, "n_treat"))
    setnames(n_c, c(col_horizon, "n_ctrl"))
    n_t <- merge(n_t, n_c, by = c(col_horizon), all.x=T)

    if (!is.null(col_weight)){
      sw_t <- clean_df[get(col_dtreat) == 1, sum(get(col_weight)), by = c(col_horizon)]
      sw_c <- clean_df[get(col_dtreat) == 0, sum(get(col_weight)), by = c(col_horizon)]
      setnames(sw_t, c(col_horizon, "sw_treat"))
      setnames(sw_c, c(col_horizon, "sw_ctrl"))
      sw_t <- merge(sw_t, sw_c, by = c(col_horizon), all.x=T)
      n_t <- merge(n_t, sw_t, by = c(col_horizon), all.x=T)
    }

    dt <- merge(dt, n_t, by.x = c("h"), by.y = c(col_horizon), all.x=T)
    move_cols(dt,c("n_treat","n_ctrl"),before="formula")
    if (!is.null(col_weight)){
      move_cols(dt,c("sw_treat","sw_ctrl"),before="formula")
      dt[,weight:=col_weight]
    }

    dt

  }

  # Attach treatment year and number of treated units to CSA cohort-level output
  add_csa_meta <- function(dt, clean_df) {
    # Cohort year from "dtreat::1:time::2004" (i() syntax) or
    # "dtreat:factor(time)2004" (time_x_controls syntax); NA otherwise.
    pat <- paste0("(?:", col_time, "::|factor\\(", col_time, "\\))(-?[0-9]+)")
    m <- regmatches(dt$variable, regexec(pat, dt$variable, perl = TRUE))
    dt[, (col_time) := vapply(m, function(x) if (length(x) == 2) as.integer(x[2]) else NA_integer_,
                              integer(1))]
    # Plain "dtreat" (single-period sample under time_x_controls): the
    # cohort is the only period with treated rows
    trt_years <- unique(clean_df[get(col_dtreat) == 1][[col_time]])
    if (length(trt_years) == 1)
      dt[variable == col_dtreat & is.na(get(col_time)), (col_time) := as.integer(trt_years)]

    n_t <- clean_df[get(col_dtreat) == 1, .N, by = col_time]
    n_c <- clean_df[get(col_dtreat) == 0, .N, by = col_time]
    setnames(n_t, c(col_time, "n_treat"))
    setnames(n_c, c(col_time, "n_ctrl"))
    n_t <- merge(n_t, n_c, by = col_time, all.x=T)

    if (!is.null(col_weight)){
      sw_t <- clean_df[get(col_dtreat) == 1, sum(get(col_weight)), by = col_time]
      sw_c <- clean_df[get(col_dtreat) == 0, sum(get(col_weight)), by = col_time]
      setnames(sw_t, c(col_time, "sw_treat"))
      setnames(sw_c, c(col_time, "sw_ctrl"))
      sw_t <- merge(sw_t, sw_c, by = col_time, all.x=T)
      n_t <- merge(n_t, sw_t, by = col_time, all.x=T)
    }

    dt <- merge(dt, n_t, by = col_time, all.x=T)
    move_cols(dt,c("n_treat","n_ctrl"),before="formula")

    if (!is.null(col_weight)){
      move_cols(dt,c("sw_treat","sw_ctrl"),before="formula")
      dt[,weight:=col_weight]
    }

    dt <- dt[order(get(col_time),h),]
    dt
  }

  add_csa_meta_long <- function(dt, clean_df) {
    dt[,t_part:=ifelse(
      regexpr(paste0(":",col_time,"::"), variable) > 0,
      substr(variable,
             regexpr(paste0(":",col_time,"::"),variable) + nchar(col_time) + 3,
             nchar(variable))
      ,NA)]
    dt[,(col_time) := as.integer(ifelse(regexpr(":", t_part) > 0,
                                        substr(h_part,1,regexpr(":", t_part)-1),
                                        t_part))]
    dt[,t_part:=NULL]

    n_t <- clean_df[get(col_dtreat) == 1, .N, by = c(col_time,col_horizon)]
    n_c <- clean_df[get(col_dtreat) == 0, .N, by = c(col_time,col_horizon)]
    setnames(n_t, c(col_time,col_horizon, "n_treat"))
    setnames(n_c, c(col_time,col_horizon, "n_ctrl"))
    n_t <- merge(n_t, n_c, by = c(col_time,col_horizon), all.x=T)

    if (!is.null(col_weight)){
      sw_t <- clean_df[get(col_dtreat) == 1, sum(get(col_weight)), by = c(col_time,col_horizon)]
      sw_c <- clean_df[get(col_dtreat) == 0, sum(get(col_weight)), by = c(col_time,col_horizon)]
      setnames(sw_t, c(col_time,col_horizon, "sw_treat"))
      setnames(sw_c, c(col_time,col_horizon, "sw_ctrl"))
      sw_t <- merge(sw_t, sw_c, by = c(col_time,col_horizon), all.x=T)
      n_t <- merge(n_t, sw_t, by = c(col_time,col_horizon), all.x=T)
    }

    dt <- merge(dt, n_t, by.x = c(col_time,"h"), by.y = c(col_time,col_horizon), all.x=T)
    move_cols(dt,c("n_treat","n_ctrl"),before="formula")
    if (!is.null(col_weight)){
      move_cols(dt,c("sw_treat","sw_ctrl"),before="formula")
      dt[,weight:=col_weight]
    }
    dt <- dt[order(get(col_time),h),]
    dt


  }

  move_cols <- function(dt, cols, before = NULL, after = NULL) {
    others <- setdiff(names(dt), cols)
    if (!is.null(before)) {
      pos <- match(before, others) - 1L
    } else if (!is.null(after)) {
      pos <- match(after, others)
    } else {
      pos <- 0L
    }
    if (is.na(pos)) stop("colonne de référence introuvable ou incluse dans cols")
    setcolorder(dt, c(head(others, pos), cols, tail(others, length(others) - pos)))
    invisible(dt)
  }
  # ── Method dispatch ────────────────────────────────────────────────────────

  did_est     <- NULL
  did_est_ps  <- NULL
  did_est_det <- NULL
  did_est_or  <- NULL

  # ── lpdid: LP-DiD with control variables ──────────────────────────────────
  if (meth == "lpdid" & one_reg == FALSE) {
    results <- rbindlist(lapply(horizons, function(h) {
      print(paste ("Horizon",h))
      clean_df <- clean_sample(df, h)
      if (nrow(clean_df)>0 & clean_df[,sum(get(col_dtreat))]>0) {
        if (!is.null(clusters_h)){
          b_clusters_h <- var_at(h,col_var=clusters_h)
          col_clusters <- c(b_col_clusters,b_clusters_h)}
        if (!is.null(controls_h)){
          b_controls_h <- var_at(h,col_var=controls_h)
          controls <- c(b_controls,b_controls_h)    }
        if (!is.null(FE_h)){
          b_FE_h <- var_at(h,col_var=FE_h)
          FE <- c(b_FE,b_FE_h)}
        if (!is.null(weight_h) & is.null(weight)){
          b_weight_h <- var_at(h,col_var=weight_h)
          col_weight <- b_weight_h}

        # time_x_controls: controls -> i(time, x), FE -> FE^time (which then
        # absorbs time, so the explicit i(time) is dropped)
        tx <- tx_apply(controls, FE, clean_df)

        # With time_x_controls, plain dtreat rather than fixest's i(dtreat):
        # mixing i() with base-R factor() or ':' terms in the controls fails
        # when fixest is imported but not attached. Time is then absorbed by
        # "| time" or FE^time (or the sample has a single period).
        rhs_main <- if (isTRUE(time_x_controls))
          col_dtreat
        else if (!is.null(tx$FE) && !tx$time_in_fe)
          paste0("i(", col_dtreat, ") + i(", col_time, ")")
        else
          paste0("i(", col_dtreat, ")")

        fe_part <- if (!is.null(tx$FE)) tx$FE else col_time

        formula <- build_formula(lhs(h), rhs_main, tx$controls, fe_part)

        print("- LPDiD Regression")

        mod     <- fit_feols(formula, clean_df, col_clusters, col_weight)

        clean_df <- drop_missing(clean_df,mod)

        est <- export_coeftable(mod, h, formula)

        add_meta(est,clean_df)

      }
    }))
    did_est <- results
  }

  if (meth == "lpdid" & one_reg == TRUE) {
    clean_df    <- clean_sample_long(df)

    # add h variables
    if (!is.null(clusters_h)){
      b_clusters_h <- var_at(h,col_var=clusters_h)
      col_clusters <- c(b_col_clusters,b_clusters_h)}
    if (!is.null(controls_h)){
      b_controls_h <- var_at(h,col_var=controls_h)
      controls <- c(b_controls,b_controls_h)    }
    if (!is.null(FE_h)){
      b_FE_h <- var_at(h,col_var=FE_h)
      FE <- c(b_FE,b_FE_h)}
    if (!is.null(weight_h) & is.null(weight)){
      b_weight_h <- var_at(h,col_var=weight_h)
      col_weight <- b_weight_h}

    rhs_main    <- if (!is.null(FE))
      paste0(col_dtreat, ":i(",col_horizon,") + i(", col_time, "i.",col_horizon,")")
    else
      paste0(col_dtreat, ":i(",col_horizon,")")
    if (!is.null(controls)) controls <- paste0("i(",col_horizon,")*(", controls, ")")
    fe_part     <- if (!is.null(FE)) paste(FE, "^",col_horizon) else paste0(col_horizon,"^", col_time)
    formula     <- build_formula(lhs(h), rhs_main, controls, fe_part)

    print("- LPDiD Regression")

    mod         <- fit_feols(formula, clean_df, col_clusters, col_weight)

    clean_df <- drop_missing(clean_df,mod)

    est <- export_coeftable_long(mod, formula)

    did_est <- add_meta_long(est,clean_df)
  }

  # ── lpdid_rw: Reweighted LP-DiD with control variables ────────────────────
  if (meth == "lpdid_rw" & one_reg == FALSE) {
    results <- rbindlist(lapply(horizons, function(h) {
      print(paste ("Horizon",h))

      clean_df  <- clean_sample(df, h)
      if (nrow(clean_df)>0 & clean_df[,sum(get(col_dtreat))]>0) {
        # add h variables
        if (!is.null(clusters_h)){
          b_clusters_h <- var_at(h,col_var=clusters_h)
          col_clusters <- c(b_col_clusters,b_clusters_h)}
        if (!is.null(controls_h)){
          b_controls_h <- var_at(h,col_var=controls_h)
          controls <- c(b_controls,b_controls_h)    }
        if (!is.null(FE_h)){
          b_FE_h <- var_at(h,col_var=FE_h)
          FE <- c(b_FE,b_FE_h)}
        if (!is.null(weight_h) & is.null(weight)){
          b_weight_h <- var_at(h,col_var=weight_h)
          col_weight <- b_weight_h}

        keep_cols <- unique(c(col_dtreat, col_dep_ref, dep_at(h),
                              col_unit, col_time, col_clusters, col_weight,
                              parse_col_names(controls),
                              parse_col_names(FE)))

        clean_df  <- na.omit(clean_df[, .SD,
                                      .SDcols = intersect(keep_cols, names(clean_df))])

        clean_df  <- add_rewgt(clean_df, h, col_weight)

        # time_x_controls: controls -> i(time, x), FE -> FE^time (which then
        # absorbs time, so the explicit i(time) is dropped)
        tx <- tx_apply(controls, FE, clean_df)

        # With time_x_controls, plain dtreat rather than fixest's i(dtreat):
        # mixing i() with base-R factor() or ':' terms in the controls fails
        # when fixest is imported but not attached. Time is then absorbed by
        # "| time" or FE^time (or the sample has a single period).
        rhs_main <- if (isTRUE(time_x_controls))
          col_dtreat
        else if (!is.null(tx$FE) && !tx$time_in_fe)
          paste0("i(", col_dtreat, ") + i(", col_time, ")")
        else
          paste0("i(", col_dtreat, ")")

        fe_part <- if (!is.null(tx$FE)) tx$FE else col_time

        formula <- build_formula(lhs(h), rhs_main, tx$controls, fe_part)

        print("- LPDiD Regression")

        mod <- feols(as.formula(formula),
                     cluster = clean_df[,..col_clusters],
                     weights = clean_df$rewgt_h,
                     data    = clean_df)

        rw_df    <- clean_df            # sample of add_rewgt and of mod
        clean_df <- drop_missing(clean_df,mod)

        est <- export_coeftable(mod, h, formula)

        if (isTRUE(se_correction_rw_ipw)) {

          print("- Standard error correction")

          bw <- if (!is.null(col_weight)) rw_df[[col_weight]] else rep(1, nrow(rw_df))
          V  <- tryCatch(
            rw_ipw_combined_se(mod, rw_df, bw, col_clusters[1]),
            error = function(e) {
              warning("RW SE correction failed at h=", h, ": ", conditionMessage(e))
              c(Vc = NA_real_, Vu = NA_real_)
            })
          est[, se_correct := fifelse(grepl(col_dtreat, variable, fixed = TRUE),
                                      se * sqrt(V[["Vc"]] / V[["Vu"]]), NA_real_)]
        }

        add_meta(est,clean_df)
      }
    }), fill = TRUE)
    did_est <- results
  }

  if (meth == "lpdid_rw" & one_reg == TRUE) {
    clean_df    <- clean_sample_long(df)

    # add h variables
    if (!is.null(clusters_h)){
      b_clusters_h <- var_at(h,col_var=clusters_h)
      col_clusters <- c(b_col_clusters,b_clusters_h)}
    if (!is.null(controls_h)){
      b_controls_h <- var_at(h,col_var=controls_h)
      controls <- c(b_controls,b_controls_h)    }
    if (!is.null(FE_h)){
      b_FE_h <- var_at(h,col_var=FE_h)
      FE <- c(b_FE,b_FE_h)}
    if (!is.null(weight_h) & is.null(weight)){
      b_weight_h <- var_at(h,col_var=weight_h)
      col_weight <- b_weight_h}

    keep_cols   <- unique(c(col_dtreat, col_dep_ref, dep_at(h),
                            col_unit, col_time, col_clusters, col_weight,
                            parse_col_names(controls),
                            parse_col_names(FE)))
    clean_df    <- add_rewgt_long(clean_df, col_weight)

    rhs_main    <- if (!is.null(FE))
      paste0(col_dtreat, ":i(",col_horizon,") + i(", col_time, ",i.",col_horizon,")")
    else
      paste0(col_dtreat, ":i(",col_horizon,")")

    fe_part     <- if (!is.null(FE)) paste0(FE, "^",col_horizon) else paste0(col_horizon,"^", col_time)
    if (!is.null(controls)) controls <- paste0("i(",col_horizon,")*(", controls, ")")

    formula     <- build_formula(lhs(h), rhs_main, controls, fe_part)

    print("- LPDiD Regression")

    mod <- feols(as.formula(formula),
                 cluster = clean_df[,..col_clusters],
                 weights = clean_df$rewgt_h,
                 data    = clean_df)

    clean_df <- drop_missing(clean_df,mod)

    est <- export_coeftable_long(mod, formula)

    did_est <- add_meta_long(est,clean_df)


  }

  # ── lpdid_adj: LP-DiD with regression adjustment (RA) ─────────────────────
  # Control-only outcome regression + cluster-level combined influence-
  # function standard error, validated against lpdid.ado's listreg-based RA.
  if (meth == "lpdid_adj" & one_reg == FALSE) {
    agg_list <- vector("list", length(horizons))
    or_list  <- vector("list", length(horizons))

    for (i in seq_along(horizons)) {
      h        <- horizons[i]
      print(paste ("Horizon",h))
      clean_df <- clean_sample(df, h)
      if (nrow(clean_df) > 0 & clean_df[, sum(get(col_dtreat))] > 0) {

        # add h variables
        if (!is.null(clusters_h)) {
          b_clusters_h <- var_at(h, col_var = clusters_h)
          col_clusters <- c(b_col_clusters, b_clusters_h)
        }
        if (!is.null(controls_h)) {
          b_controls_h <- var_at(h, col_var = controls_h)
          controls <- c(b_controls, b_controls_h)
        }
        if (!is.null(FE_h)) {
          b_FE_h <- var_at(h, col_var = FE_h)
          FE <- c(b_FE, b_FE_h)
        }
        if (!is.null(weight_h) & is.null(weight)) {
          b_weight_h <- var_at(h, col_var = weight_h)
          col_weight <- b_weight_h
        }

        ctrl_df_h <- clean_df[get(col_dtreat) == 0]
        formula_ctrl <- build_control_formula(h, controls, FE, ctrl_df_h)

        print("- Outcome regression on the control group")

        control_reg <- fit_feols(formula_ctrl, ctrl_df_h, col_clusters, col_weight)
        or_list[[i]] <- export_coeftable(control_reg, h, formula_ctrl)

        trt_df <- clean_df[get(col_dtreat) == 1]
        trt_df[, predicted_no_treat := predict(control_reg, newdata = trt_df)]

        wgt_trt <- if (!is.null(col_weight)) trt_df[[col_weight]] else rep(1, nrow(trt_df))
        e1 <- trt_df[[dep_at(h)]] - trt_df[[col_dep_ref]] - trt_df$predicted_no_treat
        att_hat <- sum(wgt_trt * e1) / sum(wgt_trt)

        print("- LPDiD regression")

        formula_adj <- paste0("(", dep_at(h), "-", col_dep_ref, "-predicted_no_treat) ~ 1")
        did_adj <- fit_feols(formula_adj, trt_df, col_clusters, col_weight)
        se_naive <- did_adj$coeftable[1, "Std. Error"]

        se_correct <- if (isTRUE(se_correction)) {
          print("- Standard error correction")
          tryCatch(
            ra_combined_se(control_reg, ctrl_df_h, trt_df,
                           att_hat, e1, wgt_trt, col_clusters[1]),
            error = function(e) {
              warning("RA SE correction failed at h=", h, ": ", conditionMessage(e))
              NA_real_
            }
          )
        } else NA_real_

        est <- data.table(
          h = h, variable = col_dtreat, estimate = att_hat,
          se = se_naive, se_correct = se_correct,
          T = att_hat / se_naive,
          pvalue = 2 * pnorm(-abs(att_hat / se_naive)),
          n_obs = did_adj$nobs, formula = formula_adj,   # final stage; first stage in $or
          clusters = col_clusters[1]
        )
        agg_list[[i]] <- add_meta(est, clean_df)
      }
    }
    did_est    <- rbindlist(agg_list, fill = TRUE)
    did_est_or <- rbindlist(or_list, fill = TRUE)
  }

  # ── lpdid_dr: LP-DiD doubly robust (RA + IPW combined) ────────────────────
  # AIPW-style score (outcome regression + propensity-score weighting) with a
  # cluster-level combined influence-function standard error propagating
  # uncertainty from both nuisance models. one_reg = TRUE not implemented.
  if (meth == "lpdid_dr" & one_reg == FALSE) {
    ps_list  <- vector("list", length(horizons))
    or_list  <- vector("list", length(horizons))
    agg_list <- vector("list", length(horizons))

    for (i in seq_along(horizons)) {
      h        <- horizons[i]
      print(paste ("Horizon",h))
      clean_df <- clean_sample(df, h)
      if (nrow(clean_df) > 0 & clean_df[, sum(get(col_dtreat))] > 0) {

        # add h variables
        if (!is.null(clusters_h)) {
          b_clusters_h <- var_at(h, col_var = clusters_h)
          col_clusters <- c(b_col_clusters, b_clusters_h)
        }
        if (!is.null(controls_h)) {
          b_controls_h <- var_at(h, col_var = controls_h)
          controls <- c(b_controls, b_controls_h)
        }
        if (!is.null(FE_h)) {
          b_FE_h <- var_at(h, col_var = FE_h)
          FE <- c(b_FE, b_FE_h)
        }
        if (!is.null(weight_h) & is.null(weight)) {
          b_weight_h <- var_at(h, col_var = weight_h)
          col_weight <- b_weight_h
        }

        # Propensity score: controls + explicit time dummies (see
        # build_ps_formula; time_x_controls interactions applied there).
        formula_ps <- build_ps_formula(controls, FE, clean_df)

        print("- Propensity score regression")

        ps_mod       <- fit_feglm(formula_ps, clean_df, col_weight)
        ps_list[[i]] <- export_coeftable(ps_mod, h, formula_ps)
        clean_df     <- drop_missing(clean_df, ps_mod)

        clean_df[, ps := pmin(predict(ps_mod, newdata = clean_df, type = "response"), 1 - 1e-6)]
        clean_df[, trim := (get(col_dtreat) == 1 | ps < ps_trim) * 1]

        ctrl_full_idx <- which(clean_df[[col_dtreat]] == 0)
        ctrl_df_h <- clean_df[ctrl_full_idx, ]
        formula_ctrl <- build_control_formula(h, controls, FE, ctrl_df_h)

        print("- Outcome regression on the control group")

        control_reg <- fit_feols(formula_ctrl, ctrl_df_h, col_clusters, col_weight)
        or_list[[i]] <- export_coeftable(control_reg, h, formula_ctrl)
        clean_df[, mu0 := predict(control_reg, newdata = clean_df)]

        iw <- if (!is.null(col_weight)) clean_df[[col_weight]] / mean(clean_df[[col_weight]]) else rep(1, nrow(clean_df))
        clean_df[, dy := get(dep_at(h)) - get(col_dep_ref)]
        clean_df[, w_treat := trim * get(col_dtreat) * iw]
        clean_df[, w_cont  := trim * (1 - get(col_dtreat)) * (ps / (1 - ps)) * iw]

        eta1 <- clean_df[, sum(w_treat * (dy - mu0)) / sum(w_treat)]
        eta0 <- clean_df[, sum(w_cont  * (dy - mu0)) / sum(w_cont)]
        att_hat <- eta1 - eta0

        clean_df[, treat_hat := (w_treat * (dy - mu0)) / mean(w_treat)]
        clean_df[, cont_hat  := (w_cont  * (dy - mu0)) / mean(w_cont)]
        did_dr <- fit_feols("(treat_hat - cont_hat) ~ 1", clean_df, col_clusters, NULL)
        se_naive <- did_dr$coeftable[1, "Std. Error"]

        se_correct <- if (isTRUE(se_correction)) {
          print("- Standard error correction")

          tryCatch(
            dr_combined_se(control_reg, ctrl_full_idx, ps_mod, clean_df,
                           clean_df$w_treat, clean_df$w_cont, clean_df$dy, clean_df$mu0,
                           eta0, eta1, col_clusters[1]),
            error = function(e) {
              warning("DR SE correction failed at h=", h, ": ", conditionMessage(e))
              NA_real_
            }
          )
        } else NA_real_

        est <- data.table(
          h = h, variable = col_dtreat, estimate = att_hat,
          se = se_naive, se_correct = se_correct,
          T = att_hat / se_naive,
          pvalue = 2 * pnorm(-abs(att_hat / se_naive)),
          n_obs = did_dr$nobs,
          formula = "(treat_hat - cont_hat) ~ 1",   # final stage; first stages in $or and $ps
          clusters = col_clusters[1]
        )
        agg_list[[i]] <- add_meta(est, clean_df)
      }
    }
    did_est    <- rbindlist(agg_list, fill = TRUE)
    did_est_ps <- rbindlist(ps_list, fill = TRUE)
    did_est_or <- rbindlist(or_list, fill = TRUE)
  }

  # ── lpdid_ipw: LP-DiD with inverse probability weighting ──────────────────
  if (meth == "lpdid_ipw" & one_reg == FALSE) {
    ps_list  <- vector("list", length(horizons))
    est_list <- vector("list", length(horizons))

    for (i in seq_along(horizons)) {
      h        <- horizons[i]
      print(paste ("Horizon",h))
      clean_df <- clean_sample(df, h)
      if (nrow(clean_df)>0 & clean_df[,sum(get(col_dtreat))]>0) {

        # add h variables
        if (!is.null(clusters_h)){
          b_clusters_h <- var_at(h,col_var=clusters_h)
          col_clusters <- c(b_col_clusters,b_clusters_h)}
        if (!is.null(controls_h)){
          b_controls_h <- var_at(h,col_var=controls_h)
          controls <- c(b_controls,b_controls_h)    }
        if (!is.null(FE_h)){
          b_FE_h <- var_at(h,col_var=FE_h)
          FE <- c(b_FE,b_FE_h)}
        if (!is.null(weight_h) & is.null(weight)){
          b_weight_h <- var_at(h,col_var=weight_h)
          col_weight <- b_weight_h}

        # Propensity score: logit of dtreat on controls and explicit time
        # dummies (see build_ps_formula; identical MLE to an absorbed time FE,
        # but its uncertainty can be propagated by the SE correction).
        formula_ps <- build_ps_formula(controls, FE, clean_df)

        print("- Propensity score regression")
        ps_mod       <- fit_feglm(formula_ps, clean_df, col_weight)
        ps_list[[i]] <- export_coeftable(ps_mod, h, formula_ps)

        clean_df <- drop_missing(clean_df,ps_mod)
        # ATT IPW weights: 1 for treated, ps/(1-ps) for controls
        clean_df[, ps := ps_mod$fitted.values]
        clean_df[, attwgt := fifelse(get(col_dtreat) == 1, 1, ps / (1 - ps))]
        if (!is.null(col_weight))
          clean_df[, attwgt := get(col_weight) * attwgt]

        # Trimming: drop controls with ps >= ps_trim (treated always kept)
        ps_df    <- clean_df           # PS estimation sample (for the SE correction)
        clean_df <- clean_df[get(col_dtreat) == 1 | ps < ps_trim]

        clean_df <- add_rewgt(clean_df, h, col_w = "attwgt")

        print("- LPDiD regression")
        formula_reg <- paste0(lhs(h), " ~ i(", col_dtreat, ") | ", col_time)
        mod <- feols(as.formula(formula_reg),
                     cluster = clean_df[,..col_clusters],
                     weights = clean_df$rewgt_h,
                     data    = clean_df)

        est <- export_coeftable(mod, h, paste0(formula_reg))

        if (isTRUE(se_correction_rw_ipw)) {
          print("- Standard error correction")
          V <- tryCatch(
            rw_ipw_combined_se(mod, clean_df, clean_df$attwgt, col_clusters[1],
                               ps_mod = ps_mod, ps_df = ps_df),
            error = function(e) {
              warning("IPW SE correction failed at h=", h, ": ", conditionMessage(e))
              c(Vc = NA_real_, Vu = NA_real_)
            })
          est[, se_correct := fifelse(grepl(col_dtreat, variable, fixed = TRUE),
                                      se * sqrt(V[["Vc"]] / V[["Vu"]]), NA_real_)]
        }

        est_list[[i]] <- add_meta(est,clean_df)


      }
    }

    did_est    <- rbindlist(est_list, fill = TRUE)
    did_est_ps <- rbindlist(ps_list, fill = TRUE)
  }

  if (meth == "lpdid_ipw" & one_reg == TRUE) {
    clean_df <- clean_sample_long(df)

    # add h variables
    if (!is.null(clusters_h)){
      b_clusters_h <- var_at(h,col_var=clusters_h)
      col_clusters <- c(b_col_clusters,b_clusters_h)}
    if (!is.null(controls_h)){
      b_controls_h <- var_at(h,col_var=controls_h)
      controls <- c(b_controls,b_controls_h)    }
    if (!is.null(FE_h)){
      b_FE_h <- var_at(h,col_var=FE_h)
      FE <- c(b_FE,b_FE_h)}
    if (!is.null(weight_h) & is.null(weight)){
      b_weight_h <- var_at(h,col_var=weight_h)
      col_weight <- b_weight_h}

    formula_ps <- if (!is.null(controls) & is.null(FE)) {
      paste0(col_dtreat, " ~ (", paste(controls, collapse = " + "),
             "):i(",col_horizon,") | ", col_time, "^",col_horizon)
    } else if (is.null(controls) & !is.null(FE)) {
      paste0(col_dtreat, " ~ i(", col_time, ",i.",col_horizon,") | ",
             paste(paste0(FE, "^",col_horizon), collapse = " + "))
    } else if (!is.null(controls) & !is.null(FE)) {
      paste0(col_dtreat, " ~ (", paste(controls, collapse = " + "),
             "):i(",col_horizon,") + i(", col_time, ",i.",col_horizon,") | ",
             paste(paste0(FE, "^",col_horizon), collapse = " + "))
    } else if (is.null(controls) & is.null(FE)) {
      paste0(col_dtreat, " ~ i(", col_time, ",i.",col_horizon,")")
    }

    print("- Propensity score regression")
    ps_mod  <- fit_feglm(formula_ps, clean_df, col_weight)
    ps_list <- export_coeftable_long(ps_mod, formula_ps)

    clean_df <- drop_missing(clean_df,ps_mod)
    ps_vals  <- ps_mod$fitted.values
    clean_df[, attwgt := fifelse(get(col_dtreat) == 1, 1,
                                 ps_vals / (1 - ps_vals))]
    if (!is.null(col_weight))
      clean_df[, attwgt := get(col_weight) * attwgt]
    # Trimming: drop controls with ps >= ps_trim (treated always kept)
    clean_df <- clean_df[get(col_dtreat) == 1 | ps_vals < ps_trim]

    clean_df <- add_rewgt_long(clean_df, col_w = "attwgt")

    formula_reg <- paste0(lhs(h), " ~ ", col_dtreat,
                          ":i(",col_horizon,") | ", col_time, "^",col_horizon)

    print("- LPDiD regression")
    mod <- feols(as.formula(formula_reg),
                 cluster = clean_df[,..col_clusters],
                 weights = clean_df$rewgt_h,
                 data    = clean_df)

    est_list <- export_coeftable_long(mod, formula_reg)

    est_list <- add_meta_long(est_list,clean_df)


    did_est    <- est_list
    did_est_ps <- ps_list
  }

  # ── lpcsa_ipw: LP-CSA with inverse probability weighting ──────────────────
  if (meth == "lpcsa_ipw" & one_reg == FALSE) {
    ps_list  <- vector("list", length(horizons))
    det_list <- vector("list", length(horizons))
    agg_list <- vector("list", length(horizons))

    for (i in seq_along(horizons)) {
      h        <- horizons[i]
      print(paste ("Horizon",h))
      clean_df <- clean_sample(df, h)
      if (nrow(clean_df)>0 & clean_df[,sum(get(col_dtreat))]>0) {

        # add h variables
        if (!is.null(clusters_h)){
          b_clusters_h <- var_at(h,col_var=clusters_h)
          col_clusters <- c(b_col_clusters,b_clusters_h)}
        if (!is.null(controls_h)){
          b_controls_h <- var_at(h,col_var=controls_h)
          controls <- c(b_controls,b_controls_h)    }
        if (!is.null(FE_h)){
          b_FE_h <- var_at(h,col_var=FE_h)
          FE <- c(b_FE,b_FE_h)}
        if (!is.null(weight_h) & is.null(weight)){
          b_weight_h <- var_at(h,col_var=weight_h)
          col_weight <- b_weight_h}

        # Propensity score (see build_ps_formula; time as explicit dummies,
        # identical MLE to the former absorbed time FE)
        formula_ps <- build_ps_formula(controls, FE, clean_df)

        print("- Propensity score regression")
        ps_mod       <- fit_feglm(formula_ps, clean_df, col_weight)
        ps_list[[i]] <- export_coeftable(ps_mod, h, formula_ps)

        clean_df <- drop_missing(clean_df,ps_mod)
        ps_vals <- ps_mod$fitted.values
        clean_df[, attwgt := fifelse(get(col_dtreat) == 1, 1,
                                     ps_vals / (1 - ps_vals))]
        if (!is.null(col_weight))
          clean_df[, attwgt := get(col_weight) * attwgt]
        # Trimming: drop controls with ps >= ps_trim (treated always kept)
        clean_df <- clean_df[get(col_dtreat) == 1 | ps_vals < ps_trim]

        # CSA cohort regression: treatment interacted with treatment cohort (time)
        formula_reg <- paste0(lhs(h), " ~ i(", col_dtreat, ",i.", col_time,
                              ",ref=c(0)) | ", col_time)

        print("- LPDiD regression")
        mod <- feols(as.formula(formula_reg),
                     data    = clean_df,
                     weights = clean_df$attwgt,
                     cluster = clean_df[,..col_clusters])

        det_h         <- export_coeftable(mod, h, formula_reg)
        det_h         <- add_csa_meta(det_h, clean_df)
        det_list[[i]] <- det_h

        # Aggregate cohort estimates to a single ATT per horizon
        mod_agg       <- aggregate(mod, paste0("(", col_dtreat, ")"))
        est  <- export_agg(model     = mod_agg,
                           model_det = mod,
                           h         = h,
                           formula_str = paste0("aggregate(mod,(",
                                                col_dtreat, "))"))


        agg_list[[i]] <- add_meta(est,clean_df)
      }
    }

    did_est     <- rbindlist(agg_list)
    did_est_det <- rbindlist(det_list, fill = TRUE)
    did_est_ps  <- rbindlist(ps_list)
  }

  if (meth == "lpcsa_ipw" & one_reg == TRUE) {
    clean_df <- clean_sample_long(df)

    # add h variables
    if (!is.null(clusters_h)){
      b_clusters_h <- var_at(h,col_var=clusters_h)
      col_clusters <- c(b_col_clusters,b_clusters_h)}
    if (!is.null(controls_h)){
      b_controls_h <- var_at(h,col_var=controls_h)
      controls <- c(b_controls,b_controls_h)    }
    if (!is.null(FE_h)){
      b_FE_h <- var_at(h,col_var=FE_h)
      FE <- c(b_FE,b_FE_h)}
    if (!is.null(weight_h) & is.null(weight)){
      b_weight_h <- var_at(h,col_var=weight_h)
      col_weight <- b_weight_h}

    formula_ps <- if (!is.null(controls) & is.null(FE)) {
      paste0(col_dtreat, " ~ (", paste(controls, collapse = " + "),
             "):i(",col_horizon,") | ", col_time, "^",col_horizon)
    } else if (is.null(controls) & !is.null(FE)) {
      paste0(col_dtreat, " ~ i(", col_time, ",i.",col_horizon,") | ",
             paste(paste0(FE, "^",col_horizon), collapse = " + "))
    } else if (!is.null(controls) & !is.null(FE)) {
      paste0(col_dtreat, " ~ (", paste(controls, collapse = " + "),
             "):i(",col_horizon,") + i(", col_time, ",i.",col_horizon,") | ",
             paste(paste0(FE, "^",col_horizon), collapse = " + "))
    } else if (is.null(controls) & is.null(FE)) {
      paste0(col_dtreat, " ~ i(", col_time, ",i.",col_horizon,")")
    }

    print("- Propensity score regression")
    ps_mod  <- fit_feglm(formula_ps, clean_df, col_weight)
    ps_list <- export_coeftable_long(ps_mod, formula_ps)

    clean_df <- drop_missing(clean_df,ps_mod)
    ps_vals <- ps_mod$fitted.values
    clean_df[, attwgt := fifelse(get(col_dtreat) == 1, 1,
                                 ps_vals / (1 - ps_vals))]
    if (!is.null(col_weight))
      clean_df[, attwgt := get(col_weight) * attwgt]
    # Trimming: drop controls with ps >= ps_trim (treated always kept)
    clean_df <- clean_df[get(col_dtreat) == 1 | ps_vals < ps_trim]

    formula_reg <- paste0(lhs(h), " ~  ", col_dtreat,
                          ":i(",col_horizon,",i.", col_time, ") | ",
                          col_time, "^",col_horizon)

    print("- LPDiD regression")
    mod <- feols(as.formula(formula_reg),
                 data    = clean_df,
                 weights = clean_df$attwgt,
                 cluster = clean_df[,..col_clusters])

    det_h   <- export_coeftable_long(mod, formula_reg)
    det_h   <- add_csa_meta_long(det_h, clean_df)
    det_list <- det_h

    mod_agg  <- aggregate(mod,
                          paste0("(", col_dtreat,
                                 ":",col_horizon,")::(-?[[:digit:]]+)"))

    agg_h    <- export_agg_long(
      model       = mod_agg,
      model_det   = mod,
      formula_str = paste0("aggregate(mod,(", col_dtreat,
                           ":",col_horizon,")::(-?[[:digit:]]+)")
    )
    agg_list <- add_meta_long(agg_h,clean_df)

    did_est     <- agg_list
    did_est_det <- det_list
    did_est_ps  <- ps_list
  }

  # ── lpcsa: LP-CSA with outcome control variables ───────────────────────────
  if (meth == "lpcsa" & one_reg == FALSE) {
    det_list <- vector("list", length(horizons))
    agg_list <- vector("list", length(horizons))

    for (i in seq_along(horizons)) {
      h        <- horizons[i]
      print(paste ("Horizon",h))
      clean_df <- clean_sample(df, h)

      if (nrow(clean_df)>0 & clean_df[,sum(get(col_dtreat))]>0) {

        # add h variables
        if (!is.null(clusters_h)){
          b_clusters_h <- var_at(h,col_var=clusters_h)
          col_clusters <- c(b_col_clusters,b_clusters_h)}
        if (!is.null(controls_h)){
          b_controls_h <- var_at(h,col_var=controls_h)
          controls <- c(b_controls,b_controls_h)    }
        if (!is.null(FE_h)){
          b_FE_h <- var_at(h,col_var=FE_h)
          FE <- c(b_FE,b_FE_h)}
        if (!is.null(weight_h) & is.null(weight)){
          b_weight_h <- var_at(h,col_var=weight_h)
          col_weight <- b_weight_h}

        # CSA cohort regression: treatment interacted with treatment cohort (time)
        tx       <- tx_apply(controls, FE, clean_df)   # time_x_controls
        # With time_x_controls, the cohort term is written dtreat:factor(time)
        # -- same columns as i(dtreat, i.time, ref = 0) -- because fixest's
        # i() cannot be mixed with the base-R factor(time):(controls) terms
        # when fixest is imported but not attached. Years without treated
        # rows then give all-zero columns, which fixest drops (notes
        # silenced in that case only).
        # (single-period sample: plain dtreat, the cohort year being filled
        # in by add_csa_meta)
        rhs_main <- if (!isTRUE(time_x_controls))
          paste0("i(", col_dtreat, ",i.", col_time, ",ref=c(0))")
        else if (tx$active)
          paste0(col_dtreat, ":factor(", col_time, ")")
        else
          col_dtreat
        fe_part  <- if (!is.null(tx$FE)) tx$FE else col_time
        formula  <- build_formula(lhs(h), rhs_main, tx$controls, fe_part)

        print("- LPDiD regression")
        mod <- fit_feols(formula, clean_df, col_clusters, col_weight,
                         notes = !isTRUE(time_x_controls))

        det_h         <- export_coeftable(mod, h, formula)

        clean_df <- drop_missing(clean_df,mod)
        det_h         <- add_csa_meta(det_h, clean_df)
        det_list[[i]] <- det_h

        # Aggregate cohort estimates to a single ATT per horizon
        mod_agg       <- aggregate(mod, paste0("(", col_dtreat, ")"))
        est <- export_agg(model       = mod_agg,
                          model_det   = mod,
                          h           = h,
                          formula_str = paste0("aggregate(mod,(",
                                               col_dtreat, "))"))


        agg_list[[i]] <- add_meta(est,clean_df)
      }
    }

    did_est     <- rbindlist(agg_list)
    did_est_det <- rbindlist(det_list, fill = TRUE)
  }

  if (meth == "lpcsa" & one_reg == TRUE) {
    clean_df <- clean_sample_long(df)

    # add h variables
    if (!is.null(clusters_h)){
      b_clusters_h <- var_at(h,col_var=clusters_h)
      col_clusters <- c(b_col_clusters,b_clusters_h)}
    if (!is.null(controls_h)){
      b_controls_h <- var_at(h,col_var=controls_h)
      controls <- c(b_controls,b_controls_h)    }
    if (!is.null(FE_h)){
      b_FE_h <- var_at(h,col_var=FE_h)
      FE <- c(b_FE,b_FE_h)}
    if (!is.null(weight_h) & is.null(weight)){
      b_weight_h <- var_at(h,col_var=weight_h)
      col_weight <- b_weight_h}

    rhs_main         <- paste0(col_dtreat, ":i(",col_horizon,",i.", col_time, ")")
    fe_part          <- if (!is.null(FE)) paste0(FE, "^",col_horizon) else paste0(col_time, "^",col_horizon)
    if (!is.null(controls)) p_variables_part <- paste0("i(",col_horizon,"):(", paste(controls, collapse = " + "), ")")
    else p_variables_part <- NULL
    formula          <- build_formula(lhs(h), rhs_main, p_variables_part, fe_part)

    print("- LPDiD regression")
    mod <- fit_feols(formula, clean_df, col_clusters, col_weight)

    det_h    <- export_coeftable_long(mod, formula)

    clean_df <- drop_missing(clean_df,mod)
    det_h    <- add_csa_meta_long(det_h, clean_df)

    det_list <- det_h


    mod_agg  <- aggregate(mod,
                          paste0("(", col_dtreat,
                                 ":",col_horizon,")::(-?[[:digit:]]+)"))
    agg_h    <- export_agg_long(
      model       = mod_agg,
      model_det   = mod,
      formula_str = paste0("aggregate(mod,(", col_dtreat,
                           ":",col_horizon,")::(-?[[:digit:]]+)")
    )

    agg_h    <- add_meta_long(agg_h, clean_df)
    agg_list <- agg_h

    did_est     <- agg_list
    did_est_det <- det_list
  }

  if (!is.null(did_est) && "se_correct" %in% names(did_est))
    move_cols(did_est, "se_correct", after = "se")

  # ── Event study plot ───────────────────────────────────────────────────────

  # Keep only the treatment coefficient rows, add a zero bar at h = -1
  # (the reference period, normalised to zero by construction)
  ddplot <- did_est[grepl(col_dtreat, variable, fixed = TRUE)]

  ref_row          <- ddplot[h == 0][1]
  ref_row$h        <- -1L - anticipation
  ref_row$estimate <- 0
  ref_row$se       <- 0
  if ("se_correct" %in% names(ref_row)) ref_row$se_correct <- 0
  ddplot <- rbind(ddplot, ref_row, fill = TRUE)

  # Confidence intervals from the corrected SE where it was estimated,
  # from the naive (default) SE otherwise
  if (!"se_correct" %in% names(ddplot)) ddplot[, se_correct := NA_real_]
  ddplot[, se_plot := fifelse(!is.na(se_correct), se_correct, se)]
  plot_caption <- if (ddplot[h != -1L - anticipation, any(!is.na(se_correct))])
    "95% CI from corrected standard errors" else "95% CI from naive standard errors"

  plotdid <- ggplot(ddplot, aes(x = h, y = estimate)) +
    geom_hline(yintercept = 0,    colour = "grey70") +
    geom_vline(xintercept = -0.5, colour = "red", linetype = "dashed") +
    geom_errorbar(aes(ymin = estimate - 1.96 * se_plot,
                      ymax = estimate + 1.96 * se_plot,
                      width = 0)) +
    geom_line()  +
    geom_point() +
    scale_x_continuous(breaks = -abs(n_pre):n_post) +
    labs(x = "Horizon", y = "Estimate", caption = plot_caption) +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank())

  # ── Return ─────────────────────────────────────────────────────────────────
  list(
    est     = did_est,
    est_det = did_est_det,
    ps      = did_est_ps,
    or      = did_est_or,
    plot    = plotdid
  )
}
