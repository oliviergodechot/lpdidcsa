# Estimation: lpdidcsa function ####

#' LP-DiD and LP-CSA Event Study Estimation
#'
#' Estimates event study coefficients using Local Projections
#' Difference-in-Differences (LP-DiD). Six estimators are supported:
#'
#' \describe{
#'   \item{\code{"lpdid"}}{LP-DiD with control variables. Yields a
#'     variance-weighted ATT under the control-variable homogeneity
#'     hypothesis.}
#'   \item{\code{"lpdid_rw"}}{Reweighted LP-DiD with control variables.
#'     Yields an equally-weighted ATT under the control-variable homogeneity
#'     hypothesis.}
#'   \item{\code{"lpdid_adj"}}{LP-DiD with adjusted regression via
#'     \code{avg_comparisons}. Yields an equally-weighted ATT (slower).}
#'   \item{\code{"lpdid_ipw"}}{LP-DiD with inverse probability weighting.
#'     Yields an equally-weighted ATT.}
#'   \item{\code{"lpcsa_ipw"}}{LP-CSA with inverse probability weighting
#'     (local-projection version of Callaway & Sant'Anna).}
#'   \item{\code{"lpcsa"}}{LP-CSA with control variables. Yields an
#'     equally-weighted ATT under the control-variable homogeneity
#'     hypothesis.}
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
#' @param nb_pre Integer >= 0 or \code{NULL}. Number of pre-treatment horizons.
#'   If \code{NULL}, all available pre-periods are used.
#' @param nb_post Integer >= 0 or \code{NULL}. Number of post-treatment
#'   horizons. If \code{NULL}, all available post-periods are used.
#' @param controls Character vector of control variable column names, or
#'   \code{NULL}. Accepts \code{fixest}-style formula strings (e.g.
#'   \code{"i(x)"}, \code{"s(x)"}).
#' @param controls_h Character vector of control horizon-specific 
#' variable column names, needed or \code{NULL}. Does not accept \code{fixest}-style 
#' formula strings (e.g. \code{"i(x)"}, \code{"s(x)"}). Enter the initial variable 
#' name. The function adds the suffix (e.g. '_tm5',...,'_t5') if \code{type_horizon = "wide"} 
#' and the suffix '_th' if \code{type_horizon = "long"}.
#' @param FE Character vector of additional fixed-effect column names, or
#'   \code{NULL}.
#' @param FE_h Character vector of additional horizon-specific fixed-effect column names, or
#'   \code{NULL}.
#' @param clusters Character. Column names used for clustering standard errors.
#'   Defaults to \code{unit}.
#' @param clusters_h Character vector of control horizon-specific variable column names, needed or
#'   \code{NULL}. Enter the initial variable names. The function
#'   adds the suffix (e.g. '_tm5',...,'_t5') if \code{type_horizon = "wide"} and the
#'   suffix '_th' if \code{type_horizon = "long"}.
#' @param weight Character. Column name for survey/sampling weights, or
#'   \code{NULL}.
#' @param weight_h Character. Horizon-specific column name for survey/sampling 
#' weights, or \code{NULL}.
#' @param meth Character. Estimator to use. One of \code{"lpdid"},
#'   \code{"lpdid_rw"}, \code{"lpdid_adj"}, \code{"lpdid_ipw"},
#'   \code{"lpcsa_ipw"}, \code{"lpcsa"} (default: \code{"lpdid_rw"}).
#' @param absorbing Logical. If \code{TRUE} (default), treatment is assumed
#'   absorbing and not-yet-treated units serve as controls. If \code{FALSE},
#'   a non-absorbing design is used.
#' @param nonabs_reentry Integer or \code{NULL}. Only used when
#'   \code{absorbing = FALSE}. Minimum number of periods since last treatment
#'   for a unit to re-enter the DiD comparison set. \code{NULL} imposes no
#'   restriction on prior treatment.
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
#'       \code{T}, \code{pvalue}, \code{nb_obs}, \code{formula}.}
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
#' \dontrun{
#' library(data.table); library(fixest); library(ggplot2); library(did)
#' data(mpdta_r, package = "lpdidcsa")
#'
#' # Prepare data, then estimate with reweighted LP-DiD
#' df  <- lpdidcsa_data(mpdta_r, unit = "countyreal", time = "year",
#'                      dependent = "lemp", treat = "treat")
#' res <- lpdidcsa(df, unit = "countyreal", time = "year",
#'                 dependent = "lemp", meth = "lpdid_rw")
#' res$plot
#' res$est
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
                     nb_pre         = NULL,
                     nb_post        = NULL,
                     controls       = NULL,
                     controls_h     = NULL,
                     FE             = NULL,
                     FE_h           = NULL,
                     clusters       = NULL,
                     clusters_h     = NULL,
                     weight         = NULL,
                     weight_h       = NULL,                     
                     meth           = "lpdid_rw",
                     absorbing      = TRUE,
                     nonabs_reentry = NULL,
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
  # e.g. "i(x) + s(y)" -> c("x", "y")
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
  
  # ── Setup ──────────────────────────────────────────────────────────────────
  df <- setDT(data)
  
  col_unit    <- unit
  col_time    <- time
  col_dep     <- dependent
  col_dtreat  <- dtreat
  col_clusters <- if (!is.null(clusters)) clusters else unit
  col_weight  <- weight
  col_dep_tm1 <- paste0(col_dep, "_tm1")
  
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
    max_nb_pre  <- -min(time_dep)
    max_nb_post <- max(time_dep)
    
  } else if (type_horizon == "long") {
    max_nb_pre  <- -df[, min(horizon, na.rm = TRUE)]
    max_nb_post <-  df[, max(horizon, na.rm = TRUE)]
  }
  
  # Cap nb_pre / nb_post at what the data supports
  if (is.null(nb_post) | !is.numeric(nb_post)) {
    nb_post <- max_nb_post
  } else if (nb_post == floor(nb_post)) {
    nb_post <- pmin(nb_post, max_nb_post)
  } else {
    nb_post <- max_nb_post
  }
  
  if (is.null(nb_pre) | !is.numeric(nb_pre)) {
    nb_pre <- max_nb_pre
  } else if (nb_pre == floor(nb_pre)) {
    nb_pre <- pmin(max_nb_pre, abs(nb_pre))
  } else {
    nb_pre <- max_nb_pre
  }
  
  horizons <- setdiff(-abs(nb_pre):nb_post, -1L)
  
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
    paste0("(", dep_at(h), " - ", col_dep_tm1, ")")
  }
  
  # Filter rows to the clean DiD comparison set for horizon h.
  # Treated units are always kept; control units are kept only if they are
  # not yet treated (absorbing) or satisfy the non-absorbing entry conditions.
  clean_sample <- function(df, h) {
    col_dep_th  <- dep_at(h)
    base_filter <- !is.na(df[[col_dep_tm1]]) & !is.na(df[[col_dep_th]]) & !is.na(df[[col_dtreat]])
    
    if (absorbing) {
      keep <- (df[[col_dtreat]] == 1 |
                 is.na(df$first.treat) |
                 (!is.na(df$first.treat) & h <  0 & df[[col_time]] <  df$first.treat) |
                 (!is.na(df$first.treat) & h >= 0 & df[[col_time]] + h < df$first.treat))
    } else if (is.null(nonabs_reentry)) {
      keep <- (is.na(df$last.treat) &
                 (is.na(df$next.treat) |
                    (!is.na(df$next.treat) &
                       df[[col_time]] + h * (h >= 0) < df$next.treat)))
    } else {
      keep <- ((is.na(df$last.treat) |
                  df[[col_time]] - df$last.treat > nonabs_reentry) &
                 (is.na(df$next.treat) |
                    (!is.na(df$next.treat) &
                       df[[col_time]] + h * (h >= 0) < df$next.treat)))
    }
    
    if (type_horizon == "long") {
      df[keep & base_filter & horizon == h]
    } else {
      df[keep & base_filter]
    }
  }
  
  # Same filtering logic for long format when one_reg = TRUE (all horizons at once)
  clean_sample_long <- function(df) {
    
    col_dep_th  <- dep_at(h)
    base_filter <- !is.na(df[[col_dep_tm1]]) & !is.na(df[[col_dep_th]]) & !is.na(df[[col_dtreat]])
    
    if (absorbing) {
      keep <- (df[[col_dtreat]] == 1 |
                 is.na(df$first.treat) |
                 (!is.na(df$first.treat) & df[[horizon]] <  0 & df[[col_time]] <  df$first.treat) |
                 (!is.na(df$first.treat) & df[[horizon]] >= 0 & df[[col_time]] + df[[horizon]] < df$first.treat))
    } else if (is.null(nonabs_reentry)) {
      keep <- (is.na(df$last.treat) &
                 (is.na(df$next.treat) |
                    (!is.na(df$next.treat) &
                       df[[col_time]] + df[[horizon]] * (df[[horizon]] >= 0) < df$next.treat)))
    } else {
      keep <- ((is.na(df$last.treat) |
                  df[[col_time]] - df$last.treat > nonabs_reentry) &
                 (is.na(df$next.treat) |
                    (!is.na(df$next.treat) &
                       df[[col_time]] + df[[horizon]] * (df[[horizon]] >= 0) < df$next.treat)))
    }
    df[keep & base_filter]
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
      nb_obs   = model$nobs,
      formula  = formula_str
    )
    setnames(out, c("h", "variable", "estimate", "se", "T", "pvalue",
                    "nb_obs", "formula"))
    out
  }
  
  # Same as above for the one_reg = TRUE case: extracts h from the variable name
  export_coeftable_long <- function(model, formula_str) {
    out <- data.table(
      variable = row.names(model$coeftable),
      model$coeftable,
      nb_obs   = model$nobs,
      formula  = formula_str
    )
    h <- as.numeric(ifelse(
      regexpr("horizon::", out$variable) > 0,
      substr(out$variable,
             regexpr("horizon::", out$variable) + 9,
             nchar(out$variable)),
      -1
    ))
    out <- data.frame(h, out)
    setDT(out)
    setnames(out, c("h", "variable", "estimate", "se", "T", "pvalue",
                    "nb_obs", "formula"))
    out
  }
  
  # Coerce aggregate() output to a standardised data.table (one horizon)
  export_agg <- function(model, model_det, h, formula_str) {
    out <- data.table(
      h        = h,
      variable = row.names(model),
      model,
      nb_obs   = model_det$nobs,
      formula  = formula_str
    )
    setnames(out, c("h", "variable", "estimate", "se", "T", "pvalue",
                    "nb_obs", "formula"))
    out
  }
  
  # Same as above for the one_reg = TRUE case
  export_agg_long <- function(model, model_det, formula_str) {
    out <- data.table(
      variable = row.names(model),
      model,
      nb_obs   = model_det$nobs,
      formula  = formula_str
    )
    h <- as.numeric(ifelse(
      regexpr("horizon::", out$variable) > 0,
      substr(out$variable,
             regexpr("horizon::", out$variable) + 9,
             nchar(out$variable)),
      -1
    ))
    out <- data.frame(h, out)
    setDT(out)
    setnames(out, c("h", "variable", "estimate", "se", "T", "pvalue",
                    "nb_obs", "formula"))
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
    formula_rw <- as.formula(paste0(col_dtreat, " ~ 1 | horizon^", col_time))
    wgtreg <- if (is.null(col_w))
      feols(formula_rw, data = df)
    else
      feols(formula_rw, weights = df[[col_w]], data = df)
    
    df[, num_wgt_h := wgtreg$residuals]
    df[, num_wgt_h := fifelse(is.na(get(col_dtreat)) | get(col_dtreat) == 0,
                              NA, num_wgt_h)]
    df[, wgt_h     := num_wgt_h / sum(num_wgt_h, na.rm = TRUE)]
    df[, gwgt_h    := max(wgt_h, na.rm = TRUE), by = list(get(col_time), horizon)]
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
  add_csa_meta <- function(dt, clean_df) {
    dt[, (col_time) := as.integer(substr(variable,
                                         nchar(variable) - 3,
                                         nchar(variable)))]
    n_t <- clean_df[get(col_dtreat) == 1, .N, by = col_time]
    setnames(n_t, c(col_time, "nb_treat"))
    merge(dt, n_t, by = col_time)
  }
  
  # ── Method dispatch ────────────────────────────────────────────────────────
  
  did_est     <- NULL
  did_est_ps  <- NULL
  did_est_det <- NULL
  
  # ── lpdid: LP-DiD with control variables ──────────────────────────────────
  if (meth == "lpdid" & one_reg == FALSE) {
    results <- rbindlist(lapply(horizons, function(h) {
      clean_df <- clean_sample(df, h)
      if (nrow(clean_df)>0) {
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
          export_coeftable(mod, h, formula)
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
      paste0(col_dtreat, ":i(horizon) + i(", col_time, "i.horizon)")
    else
      paste0(col_dtreat, ":i(horizon)")
    controls <- paste0("i(horizon)*(", controls, ")")
    fe_part     <- if (!is.null(FE)) paste(FE, "^horizon") else paste0("horizon^", col_time)
    formula     <- build_formula(lhs(h), rhs_main, controls, fe_part)
    mod         <- fit_feols(formula, clean_df, col_clusters, col_weight)
    did_est     <- export_coeftable_long(mod, formula)
  }
  
  # ── lpdid_rw: Reweighted LP-DiD with control variables ────────────────────
  if (meth == "lpdid_rw" & one_reg == FALSE) {
    results <- rbindlist(lapply(horizons, function(h) {
      clean_df  <- clean_sample(df, h)
      if (nrow(clean_df)>0) {
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
        
        keep_cols <- unique(c(col_dtreat, col_dep_tm1, dep_at(h),
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
        export_coeftable(mod, h, formula)
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
    
    keep_cols   <- unique(c(col_dtreat, col_dep_tm1, dep_at(h),
                            col_unit, col_time, col_clusters, col_weight,
                            parse_col_names(controls),
                            parse_col_names(FE)))
    clean_df    <- add_rewgt_long(clean_df, col_weight)
    rhs_main    <- if (!is.null(FE))
      paste0(col_dtreat, ":i(horizon) + i(", col_time, ",i.horizon)")
    else
      paste0(col_dtreat, ":i(horizon)")
    fe_part     <- if (!is.null(FE)) paste0(FE, "^horizon") else paste0("horizon^", col_time)
    controls <- paste0("i(horizon)*(", controls, ")")
    formula     <- build_formula(lhs(h), rhs_main, controls, fe_part)
    
    mod <- feols(as.formula(formula),
                 cluster = clean_df[,..col_clusters],
                 weights = clean_df$rewgt_h,
                 data    = clean_df)
    did_est <- export_coeftable_long(mod, formula)
  }
  
  # ── lpdid_adj: LP-DiD with adjusted regression via avg_comparisons ─────────
  if (meth == "lpdid_adj" & one_reg == FALSE) {
    det_list <- vector("list", length(horizons))
    agg_list <- vector("list", length(horizons))
    
    for (i in seq_along(horizons)) {
      h        <- horizons[i]
      clean_df <- clean_sample(df, h)
      if (nrow(clean_df)>0) {
        
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
        
        agg_list[[i]] <- as.data.table(data.frame(h, out_h[, 1:6], mod$nobs))
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
    
    rhs_main <- paste0(col_dtreat, "*factor(horizon)*factor(", col_time, ")")
    
    if (is.null(controls) & is.null(FE)) {
      formula <- paste0("(", lhs(h), ")~ ", rhs_main)
    } else if (is.null(controls) & !is.null(FE)) {
      formula <- paste0("(", lhs(h), ")~ ", rhs_main, "   | ",
                        paste(paste0(FE, "^", col_dtreat, "^horizon"),
                              collapse = " + "))
    } else if (!is.null(controls) & is.null(FE)) {
      formula <- paste0("(", lhs(h), ")~", rhs_main, "  + ", col_dtreat,
                        " * factor(horizon) * (",
                        paste(controls, collapse = " + "), ")")
    } else {
      formula <- paste0("(", lhs(h), ")~ ", rhs_main, "  + ", col_dtreat,
                        " * factor(horizon) * (",
                        paste(controls, collapse = " + "), ") | ",
                        paste(paste0(FE, "^", col_dtreat, "^horizon"),
                              collapse = " + "))
    }
    
    clean_df <- as.data.frame(clean_df)
    mod      <- fit_feols(formula, clean_df, col_clusters, col_weight)

    
    out_h <- avg_comparisons(
      mod,
      variables = col_dtreat,
      type      = "response",
      newdata   = clean_df[clean_df[[col_dtreat]] == 1, ],
      by        = "horizon",
      vcov      = "HC1"
    )

    
    did_est     <- as.data.table(data.frame(out_h[, 1:6], mod$nobs))
    did_est_det <- export_coeftable_long(mod, formula)
    setnames(did_est, "horizon", "h")
    
    # Harmonise column names from avg_comparisons output
    setnames(did_est,
             intersect(c("term", "std.error"), names(did_est)),
             c("variable", "se")[seq_along(
               intersect(c("term", "std.error"), names(did_est)))])
  }
  
  # ── lpdid_ipw: LP-DiD with inverse probability weighting ──────────────────
  if (meth == "lpdid_ipw" & one_reg == FALSE) {
    ps_list  <- vector("list", length(horizons))
    est_list <- vector("list", length(horizons))
    
    for (i in seq_along(horizons)) {
      h        <- horizons[i]
      clean_df <- clean_sample(df, h)
      if (nrow(clean_df)>0) {
        
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
        formula_ps <- if (!is.null(controls))
          paste0(col_dtreat, " ~ ", paste(controls, collapse = " + "),
                 " | ", col_time)
        else if (!is.null(FE))
          paste0(col_dtreat, " ~ i(", col_time, ") | ",
                 paste(FE, collapse = " + "))
        else
          paste0(col_dtreat, " ~ i(", col_time, ")")
        
        ps_mod       <- fit_feglm(formula_ps, clean_df, col_weight)
        ps_list[[i]] <- export_coeftable(ps_mod, h, formula_ps)
        
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
        est_list[[i]] <- export_coeftable(mod, h, formula_reg)
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
    
    formula_ps <- if (!is.null(controls)) {
      controls <- paste0("i(horizon)*(", controls, ")")
      paste0(col_dtreat, " ~ ", paste(controls, collapse = " + "),
             ":i(horizon) | ", col_time, "^horizon")
    } else if (!is.null(FE)) {
      paste0(col_dtreat, " ~ i(", col_time, ",i.horizon) | ",
             paste(paste0(FE, "^horizon"), collapse = " + "))
    } else {
      paste0(col_dtreat, " ~ i(", col_time, ",i.horizon)")
    }
    
    ps_mod  <- fit_feglm(formula_ps, clean_df, col_weight)
    ps_list <- export_coeftable_long(ps_mod, formula_ps)
    
    ps_vals  <- ps_mod$fitted.values
    clean_df[, attwgt := fifelse(get(col_dtreat) == 1, 1,
                                 ps_vals / (1 - ps_vals))]
    if (!is.null(col_weight))
      clean_df[, attwgt := get(col_weight) * attwgt]
    
    clean_df <- add_rewgt_long(clean_df, col_w = "attwgt")
    
    formula_reg <- paste0(lhs(h), " ~ ", col_dtreat,
                          ":i(horizon) | ", col_time, "^horizon")
    mod <- feols(as.formula(formula_reg),
                 cluster = clean_df[,..col_clusters],
                 weights = clean_df$rewgt_h,
                 data    = clean_df)
    est_list <- export_coeftable_long(mod, formula_reg)
    
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
      if (nrow(clean_df)>0) {
      
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
        agg_list[[i]] <- export_agg(model     = mod_agg,
                                    model_det = mod,
                                    h         = h,
                                    formula_str = paste0("aggregate(mod,(",
                                                         col_dtreat, "))"))
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
    
    formula_ps <- if (!is.null(controls) && !is.null(FE)) {
      controls <- paste0("i(horizon)*(", controls, ")")
      paste0(col_dtreat,
             " ~ factor(horizon)*(i(", col_time, ") + (",
             paste(controls, collapse = " + "), ")) | ",
             paste(paste0(FE, "^horizon"), collapse = " + "))
    } else if (!is.null(controls)) {
      paste0(col_dtreat, " ~ factor(horizon)*(",
             paste(controls, collapse = " + "), ") | ",
             col_time, "^horizon")
    } else if (!is.null(FE)) {
      paste0(col_dtreat, " ~ factor(horizon)*i(", col_time, ") | ",
             paste(paste0(FE, "^horizon"), collapse = " + "))
    } else {
      paste0(col_dtreat, " ~ factor(horizon)*i(", col_time, ")")
    }
    
    ps_mod  <- fit_feglm(formula_ps, clean_df, col_weight)
    ps_list <- export_coeftable_long(ps_mod, formula_ps)
    
    ps_vals <- ps_mod$fitted.values
    clean_df[, attwgt := fifelse(get(col_dtreat) == 1, 1,
                                 ps_vals / (1 - ps_vals))]
    if (!is.null(col_weight))
      clean_df[, attwgt := get(col_weight) * attwgt]
    
    formula_reg <- paste0(lhs(h), " ~  ", col_dtreat,
                          ":i(horizon,i.", col_time, ") | ",
                          col_time, "^horizon")

    mod <- feols(as.formula(formula_reg),
                 data    = clean_df,
                 weights = clean_df$attwgt,
                 cluster = clean_df[,..col_clusters])

    det_h   <- export_coeftable_long(mod, formula_reg)
    det_h   <- add_csa_meta(det_h, clean_df)
    det_list <- det_h
    
    mod_agg  <- aggregate(mod,
                          paste0("(", col_dtreat,
                                 ":horizon)::(-?[[:digit:]]+)"))

    agg_h    <- export_agg_long(
      model       = mod_agg,
      model_det   = mod,
      formula_str = paste0("aggregate(mod,(", col_dtreat,
                           ":horizon)::(-?[[:digit:]]+)")
    )
    agg_list <- agg_h
    
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
      
      if (nrow(clean_df)>0) {
        
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
        det_h         <- add_csa_meta(det_h, clean_df)
        det_list[[i]] <- det_h
        
        # Aggregate cohort estimates to a single ATT per horizon
        mod_agg       <- aggregate(mod, paste0("(", col_dtreat, ")"))
        agg_list[[i]] <- export_agg(model       = mod_agg,
                                    model_det   = mod,
                                    h           = h,
                                    formula_str = paste0("aggregate(mod,(",
                                                         col_dtreat, "))"))
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
    
    rhs_main         <- paste0(col_dtreat, ":i(horizon,i.", col_time, ")")
    fe_part          <- if (!is.null(FE)) paste0(FE, "^horizon") else paste0(col_time, "^horizon")
    p_variables_part <- paste0("i(horizon):(", paste(controls, collapse = " + "), ")")
    formula          <- build_formula(lhs(h), rhs_main, p_variables_part, fe_part)
    
    mod <- fit_feols(formula, clean_df, col_clusters, col_weight)
    
    det_h    <- export_coeftable_long(mod, formula)
    det_h    <- add_csa_meta(det_h, clean_df)
    det_list <- det_h
    
    mod_agg  <- aggregate(mod,
                          paste0("(", col_dtreat,
                                 ":horizon)::(-?[[:digit:]]+)"))
    agg_h    <- export_agg_long(
      model       = mod_agg,
      model_det   = mod,
      formula_str = paste0("aggregate(mod,(", col_dtreat,
                           ":horizon)::(-?[[:digit:]]+)")
    )
    agg_list <- agg_h
    
    did_est     <- agg_list
    did_est_det <- det_list
  }
  
  # ── Event study plot ───────────────────────────────────────────────────────
  
  # Keep only the treatment coefficient rows, add a zero bar at h = -1
  # (the reference period, normalised to zero by construction)
  ddplot <- did_est[grepl(col_dtreat, variable, fixed = TRUE)]
  
  ref_row          <- ddplot[h == 0][1]
  ref_row$h        <- -1L
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
    scale_x_continuous(breaks = -abs(nb_pre):nb_post) +
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
