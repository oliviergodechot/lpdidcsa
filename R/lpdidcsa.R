# Estimation: lpdidcsa function ####

#' LP-DiD and LP-CSA Event Study Estimation
#'
#' Estimates event study coefficients using Local Projections
#' Difference-in-Differences (LP-DiD). Six estimators are supported:
#'
#' \describe{
#'   \item{\code{"lpdid"}}{LP-DiD with control variables. Yields a
#'     variance-weighted ATT under the hypothesis that treatment effects do not
#'     vary with the value of covariates.}
#'   \item{\code{"lpdid_rw"}}{Reweighted LP-DiD with control variables.
#'     Yields an equally-weighted ATT under the hypothesis that treatment effects do not
#'     vary with the value of covariates.}
#'   \item{\code{"lpdid_adj"}}{LP-DiD with adjusted regression via
#'     \code{avg_comparisons}. Yields an equally-weighted ATT (slower).}
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
#'   \code{"lpdid_rw"}, \code{"lpdid_adj"}, \code{"lpdid_ipw"},
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
#'   horizon (default: \code{FALSE}).
#'
#' @return A named list with four elements:
#'   \describe{
#'     \item{\code{est}}{Main estimation results as a \code{data.table} with
#'       columns \code{h}, \code{variable}, \code{estimate}, \code{se},
#'       \code{T}, \code{pvalue}, \code{n_obs}, \code{formula}, \code{clusters}.}
#'     \item{\code{est_det}}{Cohort-level (disaggregated) estimates as a
#'       \code{data.table}. Non-\code{NULL} for CSA and adjusted methods
#'       only.}
#'     \item{\code{ps}}{Propensity score model estimates as a
#'       \code{data.table}. Non-\code{NULL} for IPW methods only.}
#'     \item{\code{plot}}{A \code{ggplot2} event study plot with 95\%
#'       confidence intervals.}
#'   }
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
#' }
#'
#' @import data.table fixest ggplot2
#' @importFrom marginaleffects avg_comparisons
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
                     one_reg        = FALSE)
{
  # ── Input validation ───────────────────────────────────────────────────────
  valid_meths <- c("lpdid", "lpdid_rw", "lpdid_adj",
                   "lpdid_ipw", "lpcsa_ipw", "lpcsa")
  stopifnot(
    is.data.frame(data),
    is.character(unit)      && length(unit)      == 1,
    is.character(time)      && length(time)      == 1,
    is.character(dependent) && length(dependent) == 1,
    is.character(dtreat)    && length(dtreat)    == 1,
    is.logical(absorbing)   && length(absorbing) == 1,
    meth %in% valid_meths
  )
  
  
  if (meth=="lpdid_ipw" & 
      is.null(FE) & is.null(FE_h) & is.null(controls) & is.null(controls_h)) {meth <- "lpdid_rw"}
  if (meth=="lpcsa_ipw" & 
      is.null(FE) & is.null(FE_h) & is.null(controls) & is.null(controls_h)) {meth <- "lpcsa"}

  
    # ── Helpers ────────────────────────────────────────────────────────────────
  
  # Extract bare column names from fixest-style formula strings
  # e.g. "i(x) + log(y)" -> c("x", "y")
  parse_col_names <- function(x) {
    x <- gsub("[a-zA-Z][a-zA-Z0-9_]*\\(([^)]+)\\)", "\\1", x)
    x <- unlist(strsplit(x, "[\\*\\+\\:\\^\\|\\s]+"))
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
    else if (meth != "lpdid_adj")
      time_fe
    else
      NULL
    if (!is.null(fe_part))
      paste0(lhs, " ~ ", rhs, " | ", fe_part)
    else
      paste0(lhs, " ~ ", rhs)
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
                              NA, num_wgt_h)]
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
                              NA, num_wgt_h)]
    df[, wgt_h     := num_wgt_h / sum(num_wgt_h, na.rm = TRUE)]
    df[, gwgt_h    := max(wgt_h, na.rm = TRUE), by = list(get(col_time), get(col_horizon))]
    df[, wgt_h     := fifelse(is.na(wgt_h), gwgt_h, wgt_h)]
    df[, rewgt_h   := if (is.null(col_w)) 1 / wgt_h else get(col_w) / wgt_h]
    df[, c("num_wgt_h", "wgt_h", "gwgt_h") := NULL]
    df
  }
  
  # Fit feols with optional weights and fixed clustering
   fit_feols <- function(formula, data, cl, col_w = NULL) {
    # Select clustering columns based on data type
    if (inherits(data, "data.table")) {
      cluster_cols <- data[, ..cl]  # Use ..cl for data.table
    } else {
      cluster_cols <- data[, cl, drop = FALSE]  # Use standard syntax for data.frame
    }
    
    args <- list(
      fml     = as.formula(formula),
      data    = data,
      cluster = cluster_cols
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
    dt <- data.frame(dt, n_t)
    setcolorder(dt,c("n_treat","n_ctrl"),before="formula")
    if (!is.null(col_weight)){
      setcolorder(dt,c("sw_treat","sw_ctrl"),before="formula")
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
    setcolorder(dt,c("n_treat","n_ctrl"),before="formula")
    if (!is.null(col_weight)){
      setcolorder(dt,c("sw_treat","sw_ctrl"),before="formula")
      dt[,weight:=col_weight]
      }
    
    dt
    
  }
  
  # Attach treatment year and number of treated units to CSA cohort-level output
  add_csa_meta <- function(dt, clean_df) {
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
    setcolorder(dt,c("n_treat","n_ctrl"),before="formula")
    if (!is.null(col_weight)){
      setcolorder(dt,c("sw_treat","sw_ctrl"),before="formula")
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
    setcolorder(dt,c("n_treat","n_ctrl"),before="formula")
    if (!is.null(col_weight)){
      setcolorder(dt,c("sw_treat","sw_ctrl"),before="formula")
      dt[,weight:=col_weight]
      }
    dt <- dt[order(get(col_time),h),]
    dt
    
    
  }
  
  # ── Method dispatch ────────────────────────────────────────────────────────
  
  did_est     <- NULL
  did_est_ps  <- NULL
  did_est_det <- NULL
  
  # ── lpdid: LP-DiD with control variables ──────────────────────────────────
  if (meth == "lpdid" & one_reg == FALSE) {
    results <- rbindlist(lapply(horizons, function(h) {
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
        
        rhs_main <- if (!is.null(FE))
          paste0("i(", col_dtreat, ") + i(", col_time, ")")
        else
          paste0("i(", col_dtreat, ")")
        
          fe_part <- if (!is.null(FE)) FE else col_time
          
          formula <- build_formula(lhs(h), rhs_main, controls, fe_part)
          
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
    
    mod         <- fit_feols(formula, clean_df, col_clusters, col_weight)
    
    clean_df <- drop_missing(clean_df,mod)
    
    est <- export_coeftable_long(mod, formula)
    
    did_est <- add_meta_long(est,clean_df)
  }
  
  # ── lpdid_rw: Reweighted LP-DiD with control variables ────────────────────
  if (meth == "lpdid_rw" & one_reg == FALSE) {
    results <- rbindlist(lapply(horizons, function(h) {
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
        
        rhs_main <- if (!is.null(FE))
          paste0("i(", col_dtreat, ") + i(", col_time, ")")
        else
          paste0("i(", col_dtreat, ")")
        
        fe_part <- if (!is.null(FE)) FE else col_time
        
        formula <- build_formula(lhs(h), rhs_main, controls, fe_part)
        
        mod <- feols(as.formula(formula),
                     cluster = clean_df[,..col_clusters],
                     weights = clean_df$rewgt_h,
                     data    = clean_df)
        
        clean_df <- drop_missing(clean_df,mod)
        
        est <- export_coeftable(mod, h, formula)
        
        add_meta(est,clean_df)
      }
    }))
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
    
    mod <- feols(as.formula(formula),
                 cluster = clean_df[,..col_clusters],
                 weights = clean_df$rewgt_h,
                 data    = clean_df)
    
    clean_df <- drop_missing(clean_df,mod)
    
    est <- export_coeftable_long(mod, formula)
    
    did_est <- add_meta_long(est,clean_df)
    
    
  }
  
  # ── lpdid_adj: LP-DiD with adjusted regression via avg_comparisons ─────────
  if (meth == "lpdid_adj" & one_reg == FALSE) {
    det_list <- vector("list", length(horizons))
    agg_list <- vector("list", length(horizons))
    
    for (i in seq_along(horizons)) {
      h        <- horizons[i]
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
        
        
        # Fully interacted formula: treatment x time, so avg_comparisons
        # can marginalise over the time distribution of treated units
        if (is.null(controls) & is.null(FE)) {
          formula <- paste0("(", lhs(h), ")~ i(", col_dtreat, ",i.", col_time,
                            ",ref=c(0)) + i(", col_time, ")")
        } else if (is.null(controls) & !is.null(FE)) {
          formula <- paste0("(", lhs(h), ")~ i(", col_dtreat, ",i.", col_time,
                            ",ref=c(0)) + i(", col_time, ")  | ",
                            paste(paste0(FE, "^", col_dtreat), collapse = " + "))
        } else if (!is.null(controls) & is.null(FE)) {
          formula <- paste0("(", lhs(h), ")~ i(", col_dtreat, ",i.", col_time,
                            ",ref=c(0)) + i(", col_time, ") + ", col_dtreat,
                            " * (", paste(controls, collapse = " + "), ")")
        } else {
          formula <- paste0("(", lhs(h), ")~ i(", col_dtreat, ",i.", col_time,
                            ",ref=c(0))  + i(", col_time, ") + ", col_dtreat,
                            " * (", paste(controls, collapse = " + "), ")  | ",
                            paste(paste0(FE, "^", col_dtreat), collapse = " + "))
        }
        
        clean_df <- as.data.frame(clean_df)
        mod      <- fit_feols(formula, clean_df, col_clusters, col_weight)
        
        out_h <- avg_comparisons(
          mod,
          variables = col_dtreat,
          type      = "response",
          newdata   = clean_df[clean_df[[col_dtreat]] == 1, ],
          vcov      = "HC1"
        )
        
        clean_df <- drop_missing(clean_df,mod)
        
        est <- as.data.table(data.frame(h, out_h[, 1:6], mod$nobs))
        # est <- add_meta(est,clean_df)
        
        agg_list[[i]] <- est
        det_list[[i]] <- export_coeftable(mod, h, formula)
        
         
        
        
      }
      
      did_est     <- rbindlist(agg_list, fill = TRUE)
      did_est_det <- rbindlist(det_list)
      
      # Harmonise column names from avg_comparisons output
      setnames(did_est,
               intersect(c("term", "std.error"), names(did_est)),
               c("variable", "se")[seq_along(
                 intersect(c("term", "std.error"), names(did_est)))])
    }
  }
  
  if (meth == "lpdid_adj" & one_reg == TRUE) {
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
    
    rhs_main <- paste0(col_dtreat, "*factor(",col_horizon,")*factor(", col_time, ")")
    
    if (is.null(controls) & is.null(FE)) {
      formula <- paste0("(", lhs(h), ")~ ", rhs_main)
    } else if (is.null(controls) & !is.null(FE)) {
      formula <- paste0("(", lhs(h), ")~ ", rhs_main, "   | ",
                        paste(paste0(FE, "^", col_dtreat, "^",col_horizon),
                              collapse = " + "))
    } else if (!is.null(controls) & is.null(FE)) {
      formula <- paste0("(", lhs(h), ")~", rhs_main, "  + ", col_dtreat,
                        " * factor(",col_horizon,") * (",
                        paste(controls, collapse = " + "), ")")
    } else {
      formula <- paste0("(", lhs(h), ")~ ", rhs_main, "  + ", col_dtreat,
                        " * factor(",col_horizon,") * (",
                        paste(controls, collapse = " + "), ") | ",
                        paste(paste0(FE, "^", col_dtreat, "^",col_horizon),
                              collapse = " + "))
    }
    
    clean_df <- as.data.frame(clean_df)
    mod      <- fit_feols(formula, clean_df, col_clusters, col_weight)

    
    out_h <- avg_comparisons(
      mod,
      variables = col_dtreat,
      type      = "response",
      newdata   = clean_df[clean_df[[col_dtreat]] == 1, ],
      by        = col_horizon,
      vcov      = "HC1"
    )

    did_est_det <- export_coeftable_long(mod, formula)
    
    did_est     <- as.data.table(data.frame(out_h[, 1:6], mod$nobs))
    setnames(did_est, col_horizon, "h")
    
    # Harmonise column names from avg_comparisons output
    setnames(did_est,
             intersect(c("term", "std.error"), names(did_est)),
             c("variable", "se")[seq_along(
               intersect(c("term", "std.error"), names(did_est)))])
    
    
    clean_df <- drop_missing(clean_df,mod)
    
    # did_est <- add_meta_long(did_est,clean_df)
    
  }
  
  # ── lpdid_ipw: LP-DiD with inverse probability weighting ──────────────────
  if (meth == "lpdid_ipw" & one_reg == FALSE) {
    ps_list  <- vector("list", length(horizons))
    est_list <- vector("list", length(horizons))
    
    for (i in seq_along(horizons)) {
      h        <- horizons[i]
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
        
        # Propensity score: logit of dtreat on controls, with time FE
        formula_ps <- if (!is.null(controls) & is.null(FE)) 
          {
            paste0(col_dtreat, " ~ ", paste(controls, collapse = " + "),
                   " | ", col_time) 
          } else if (is.null(controls) & !is.null(FE)) {
            paste0(col_dtreat, " ~ i(", col_time, ") | ",
                   paste(FE, collapse = " + "))
          } else if (!is.null(controls) & !is.null(FE)) {
            paste0(col_dtreat, " ~ ", paste(controls, collapse = " + "),
                   " + i(", col_time, ") | ",paste(FE, collapse = " + "))  
          } else if (is.null(controls) & is.null(FE)) {
            paste0(col_dtreat, " ~ i(", col_time, ")")
          }
        
        
        ps_mod       <- fit_feglm(formula_ps, clean_df, col_weight)
        ps_list[[i]] <- export_coeftable(ps_mod, h, formula_ps)
        
        clean_df <- drop_missing(clean_df,ps_mod)
        # ATT IPW weights: 1 for treated, ps/(1-ps) for controls
        ps_vals  <- ps_mod$fitted.values
        clean_df[, attwgt := fifelse(get(col_dtreat) == 1, 1,
                                     ps_vals / (1 - ps_vals))]
        if (!is.null(col_weight))
          clean_df[, attwgt := get(col_weight) * attwgt]
        
        clean_df <- add_rewgt(clean_df, h, col_w = "attwgt")
        
        formula_reg <- paste0(lhs(h), " ~ i(", col_dtreat, ") | ", col_time)
        mod <- feols(as.formula(formula_reg),
                     cluster = clean_df[,..col_clusters],
                     weights = clean_df$rewgt_h,
                     data    = clean_df)
        
        est <- export_coeftable(mod, h, paste0(formula_reg))
        
        est_list[[i]] <- add_meta(est,clean_df)
        
        
      }
    }
    
    did_est    <- rbindlist(est_list)
    did_est_ps <- rbindlist(ps_list)
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
    
    ps_mod  <- fit_feglm(formula_ps, clean_df, col_weight)
    ps_list <- export_coeftable_long(ps_mod, formula_ps)
    
    clean_df <- drop_missing(clean_df,ps_mod)
    ps_vals  <- ps_mod$fitted.values
    clean_df[, attwgt := fifelse(get(col_dtreat) == 1, 1,
                                 ps_vals / (1 - ps_vals))]
    if (!is.null(col_weight))
      clean_df[, attwgt := get(col_weight) * attwgt]
    
    clean_df <- add_rewgt_long(clean_df, col_w = "attwgt")
    
    formula_reg <- paste0(lhs(h), " ~ ", col_dtreat,
                          ":i(",col_horizon,") | ", col_time, "^",col_horizon)
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
        
        formula_ps <- if (!is.null(controls) && !is.null(FE))
          paste0(col_dtreat, " ~ i(", col_time, ") + ",
                 paste(controls, collapse = " + "), " | ",
                 paste(FE, collapse = " + "))
        else if (!is.null(controls))
          paste0(col_dtreat, " ~ ", paste(controls, collapse = " + "),
                 " | ", col_time)
        else if (!is.null(FE))
          paste0(col_dtreat, " ~ i(", col_time, ") | ",
                 paste(FE, collapse = " + "))
        else
          paste0(col_dtreat, " ~ i(", col_time, ")")
        
        ps_mod       <- fit_feglm(formula_ps, clean_df, col_weight)
        ps_list[[i]] <- export_coeftable(ps_mod, h, formula_ps)
        
        clean_df <- drop_missing(clean_df,ps_mod)
        ps_vals <- ps_mod$fitted.values
        clean_df[, attwgt := fifelse(get(col_dtreat) == 1, 1,
                                     ps_vals / (1 - ps_vals))]
        if (!is.null(col_weight))
          clean_df[, attwgt := get(col_weight) * attwgt]
        
        # CSA cohort regression: treatment interacted with treatment cohort (time)
        formula_reg <- paste0(lhs(h), " ~ i(", col_dtreat, ",i.", col_time,
                              ",ref=c(0)) | ", col_time)
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
    
    
    ps_mod  <- fit_feglm(formula_ps, clean_df, col_weight)
    ps_list <- export_coeftable_long(ps_mod, formula_ps)
    
    clean_df <- drop_missing(clean_df,ps_mod)
    ps_vals <- ps_mod$fitted.values
    clean_df[, attwgt := fifelse(get(col_dtreat) == 1, 1,
                                 ps_vals / (1 - ps_vals))]
    if (!is.null(col_weight))
      clean_df[, attwgt := get(col_weight) * attwgt]
    
    formula_reg <- paste0(lhs(h), " ~  ", col_dtreat,
                          ":i(",col_horizon,",i.", col_time, ") | ",
                          col_time, "^",col_horizon)

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
        rhs_main <- paste0("i(", col_dtreat, ",i.", col_time, ",ref=c(0))")
        fe_part  <- if (!is.null(FE)) FE else col_time
        formula  <- build_formula(lhs(h), rhs_main, controls, fe_part)
        
        mod <- fit_feols(formula, clean_df, col_clusters, col_weight)
        
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
  
  # ── Event study plot ───────────────────────────────────────────────────────
  
  # Keep only the treatment coefficient rows, add a zero bar at h = -1
  # (the reference period, normalised to zero by construction)
  ddplot <- did_est[grepl(col_dtreat, variable, fixed = TRUE)]
  
  ref_row          <- ddplot[h == 0][1]
  ref_row$h        <- -1L - anticipation
  ref_row$estimate <- 0
  ref_row$se       <- 0
  ddplot <- rbind(ddplot, ref_row)
  
  plotdid <- ggplot(ddplot, aes(x = h, y = estimate)) +
    geom_hline(yintercept = 0,    colour = "grey70") +
    geom_vline(xintercept = -0.5, colour = "red", linetype = "dashed") +
    geom_errorbar(aes(ymin = estimate - 1.96 * se,
                      ymax = estimate + 1.96 * se,
                      width = 0)) +
    geom_line()  +
    geom_point() +
    scale_x_continuous(breaks = -abs(n_pre):n_post) +
    labs(x = "Horizon", y = "Estimate") +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank())
  
  # ── Return ─────────────────────────────────────────────────────────────────
  list(
    est     = did_est,
    est_det = did_est_det,
    ps      = did_est_ps,
    plot    = plotdid
  )
}