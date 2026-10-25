# Data preparation : lpdidcsa_data function ####

#' Prepare Data for LP-DiD and LP-CSA Event Studies
#'
#' Constructs a wide or long horizon dataset with leads and lags of the
#' dependent variable and treatment indicator, suitable for Local Projections
#' Difference-in-Differences (LP-DiD) estimation following Dube et al. (2023)
#' and Callaway & Sant'Anna (2021) style heterogeneity-robust event studies.
#'
#' @param data A \code{data.frame} or \code{data.table}.
#' @param unit Character. Name of the unit identifier column.
#' @param time Character. Name of the time period column (default:
#'   \code{"year"}). Must be on a regular grid with unit step (e.g. years).
#' @param dependent Character. Name of the dependent variable column.
#' @param treat Character or \code{NULL} (default). Name of the binary
#'   treatment indicator column (0, 1 or \code{NA}). At least one of
#'   \code{treat} and \code{last_treat} must be supplied; \code{last_treat}
#'   is recommended, as it is robust to missing periods. Treatment onsets are
#'   defined as switches from 0 to 1 between two consecutive periods, so that
#'   \code{treat} can be either a status (1 as long as the unit is treated) or
#'   a one-period event indicator. When \code{last_treat} is also supplied,
#'   \code{treat} is only used for consistency checks and, if
#'   \code{absorbing = FALSE}, to exclude periods of ongoing treatment (see
#'   Details): with a treatment status and \code{absorbing = FALSE}, supply
#'   both.
#' @param last_treat Character or \code{NULL} (default). Name of a column
#'   giving the period of the most recent treatment onset known at \emph{t}
#'   (\code{NA} if the unit has not been treated yet), e.g. the birth year of
#'   the last child. Must be expressed on the same grid as \code{time}.
#'   Because onset dates are read from every observed row, an onset is
#'   recorded even if the unit is not observed in the onset period. When
#'   \code{absorbing = TRUE}, a time-invariant cohort variable (first
#'   treatment period, as in \code{csdid}'s \code{gvar()}) can be used
#'   instead. Never-treated units must have \code{NA}; a value of 0 is read as
#'   \code{NA} (with a message) when all periods are positive.
#' @param absorbing Logical. If \code{TRUE} (default), treatment is assumed
#'   absorbing (once treated, always treated): only the first onset of each
#'   unit is used. If \code{FALSE}, units can be treated several times.
#' @param n_pre Integer or \code{NULL}. Number of pre-treatment periods to
#'   include. Defaults to all available pre-periods if \code{NULL} or
#'   non-integer. Capped at the maximum available in the data.
#' @param n_post Integer or \code{NULL}. Number of post-treatment periods to
#'   include. Defaults to all available post-periods if \code{NULL} or
#'   non-integer. Capped at the maximum available in the data.
#' @param p_variables Named list of character vectors specifying
#'   horizon-specific control variables. List names follow the convention
#'   \code{"t0"}, \code{"t1"}, \code{"tm1"}, \code{"tm2"}, etc., where
#'   \code{"tm"} stands for \emph{t minus}. Use an empty list (default:
#'   \code{list()}) when no horizon-specific controls are needed. With
#'   \code{type_horizon = "long"} and \code{anticipation = a > 0} in
#'   \code{lpdidcsa()}, the reference outcome must be created here, e.g.
#'   \code{p_variables = list(tm2 = "lemp")} for \code{a = 1}.
#' @param h_variables Character vector of column names to be carried through at
#'   every horizon (e.g. time-varying controls). Use \code{NULL} (default) if
#'   not needed.
#' @param type_horizon Character. Output structure: \code{"wide"} (default)
#'   creates one column per horizon; \code{"long"} stacks all horizons into a
#'   single \code{horizon} column.
#' @param keep_filled Logical. If \code{TRUE} (default), periods missing
#'   between the first appearance of a unit and the last period of the panel
#'   are added, and the added
#'   rows are kept as reference rows when at least one lag of the dependent
#'   variable (\code{<dependent>_tm1}, \code{<dependent>_tm2}, ...) is
#'   observed. This allows a unit observed at \emph{t-1} and \emph{t+h} but
#'   not at \emph{t} to contribute to horizon \emph{h} (as in the
#'   pair-balanced samples of Stata's \code{csdid}). If \code{FALSE}, added
#'   rows are only used to compute leads and lags and are dropped from the
#'   output (behaviour of earlier versions).
#'
#' @return A \code{data.table} with the original columns plus:
#'   \describe{
#'     \item{\code{first.treat}}{First treatment onset per unit (\code{NA}
#'       for never-treated units). For an onset hidden in a gap of
#'       observation (treatment coded from \code{treat} only), the earliest
#'       possible onset.}
#'     \item{\code{last.treat}}{Most recent treatment onset strictly before
#'       \emph{t} (\code{NA} if none). For an onset hidden in a gap, the
#'       latest possible onset.}
#'     \item{\code{next.treat}}{Next treatment onset strictly after \emph{t}
#'       (\code{NA} if none). For an onset hidden in a gap, the earliest
#'       possible onset.}
#'     \item{\code{dtreat}}{1 if a treatment onset occurs at \emph{t}, 0
#'       otherwise, and \code{NA} if this cannot be determined (onset hidden
#'       in a gap, first period of a unit already treated, or ongoing
#'       treatment when \code{absorbing = FALSE}).}
#'     \item{Horizon columns}{For \code{type_horizon = "wide"}: outcome and
#'       control columns at each horizon, e.g. \code{lemp_tm2}, \code{lemp_t0}.
#'       Horizon \emph{-1} is always omitted (reference period).}
#'     \item{\code{horizon}}{For \code{type_horizon = "long"}: integer giving
#'       the horizon \eqn{h = t' - t} for each merged row.}
#'   }
#'
#' @details
#' \strong{Treatment onsets.} With \code{treat} only, an onset at \emph{t}
#' requires \code{treat} = 0 at \emph{t-1} and 1 at \emph{t}. When the unit is
#' not observed (or \code{treat} is \code{NA}) between a period with
#' \code{treat} = 0 and a later period with \code{treat} = 1, the onset date is
#' only known to lie in that interval: \code{dtreat} is set to \code{NA} on
#' the whole interval, \code{first.treat} and \code{next.treat} use its lower
#' bound and \code{last.treat} its upper bound, so that such units never enter
#' the control group by mistake. With \code{absorbing = FALSE}, the same
#' applies to a gap between two periods with \code{treat} = 1 (a new onset may
#' have occurred), and a one-period event occurring during a gap cannot be
#' detected at all; a warning reports the number of units concerned.
#' Supplying \code{last_treat} removes these ambiguities.
#'
#' With \code{absorbing = FALSE} and \code{treat} available, periods where
#' \code{treat} = 1 without an onset (ongoing treatment status) get
#' \code{dtreat = NA}, so that they are used neither as treated nor as control
#' observations. Onsets are counted from their start; the time since the end of
#' a treatment is not used.
#'
#' \strong{Missing periods.} Periods missing between the first appearance of
#' a unit and the last period of the panel are added before computing leads
#' and lags (a period after the last appearance can still serve as the
#' reference row of a pre-treatment horizon). Columns
#' that are constant over the observed rows of a unit (e.g. time-invariant
#' controls, fixed effects, clusters or weights) are copied to the added rows
#' of that unit; other columns are left \code{NA}.
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
#' # Initial mpdta_r characteristics
#' print(colnames(mpdta_r))
#' print(nrow(mpdta_r))
#' print(ncol(mpdta_r))
#' print(format(object.size(mpdta_r),units="Mb"))
#'
#' p_variables <- list(
#'   tm1 = c("lpop"),
#'   tm2 = c("lpop")
#' )
#'
#' # Wide format: one column per horizon
#' suppressWarnings({df_wide <- lpdidcsa_data(
#'   data        = mpdta_r,
#'   unit        = "countyreal",
#'   time        = "year",
#'   dependent   = "lemp",
#'   treat       = "treat",
#'   absorbing   = TRUE,
#'   p_variables = p_variables
#' )})
#'
#' print(colnames(df_wide))
#' print(nrow(df_wide))
#' print(ncol(df_wide))
#' print(format(object.size(df_wide),units="Mb"))
#'
#' # Same, with treatment onsets given by the cohort variable (recommended,
#' # robust to missing periods; 0 = never treated is read as NA)
#' suppressWarnings({df_wide2 <- lpdidcsa_data(
#'   data        = mpdta_r,
#'   unit        = "countyreal",
#'   time        = "year",
#'   dependent   = "lemp",
#'   last_treat  = "first.treat",
#'   p_variables = p_variables
#' )})
#'
#' # Long format: one row per unit x reference-period x horizon
#' suppressWarnings({df_long <- lpdidcsa_data(
#'   data         = mpdta_r,
#'   unit         = "countyreal",
#'   time         = "year",
#'   dependent    = "lemp",
#'   treat        = "treat",
#'   absorbing    = TRUE,
#'   p_variables  = p_variables,
#'   type_horizon = "long"
#' )})
#'
#'
#' print(colnames(df_long))
#' print(nrow(df_long))
#' print(ncol(df_long))
#' print(format(object.size(df_long),units="Mb"))
#'
#' }
#' @import data.table
#' @export
lpdidcsa_data <- function(data,
                          time        = "year",
                          treat       = NULL,
                          unit        = NULL,
                          dependent   = NULL,
                          last_treat  = NULL,
                          n_pre       = NULL,
                          n_post      = NULL,
                          absorbing   = TRUE,
                          p_variables = list(),
                          h_variables = NULL,
                          type_horizon= "wide",
                          keep_filled = TRUE)
{
  # ── Input validation ───────────────────────────────────────────────────────
  stopifnot(
    is.data.frame(data),
    is.character(unit)      && length(unit)      == 1,
    is.character(time)      && length(time)      == 1,
    is.character(dependent) && length(dependent) == 1,
    is.null(treat)      || (is.character(treat)      && length(treat)      == 1),
    is.null(last_treat) || (is.character(last_treat) && length(last_treat) == 1),
    is.logical(absorbing)   && length(absorbing)   == 1,
    is.logical(keep_filled) && length(keep_filled) == 1,
    is.list(p_variables),
    is.null(h_variables) || is.character(h_variables),
    type_horizon %in% c("wide", "long")
  )

  if (length(p_variables) > 0 &&
      (is.null(names(p_variables)) ||
       !all(grepl("^(t|tm)[0-9]+$", names(p_variables))))) {
    stop("'p_variables' must be a named list with names such as 't0', 't1', 'tm1', 'tm2'.")
  }

  if (is.null(last_treat) && is.null(treat)) {
    stop("Supply 'last_treat' (recommended: period of the most recent ",
         "treatment onset, or the cohort variable, i.e. first treatment ",
         "period, if absorbing = TRUE) or 'treat' (0/1 treatment indicator).")
  }

  required_cols <- c(unit, time, dependent, treat, last_treat,
                     unlist(p_variables, use.names = FALSE), h_variables)

  missing_cols  <- setdiff(required_cols, names(data))

  if (length(missing_cols) > 0) {
    stop("The following columns are missing from 'data': ",
         paste(missing_cols, collapse = ", "))
  }

  if (!is.null(treat) && !all(data[[treat]] %in% c(0, 1, NA))) {
    stop("'", treat, "' must only contain 0, 1 or NA.")
  }

  # ── Helpers ────────────────────────────────────────────────────────────────

  # Build the column suffix for a given horizon h ("_tm2", "_t0", "_t3", ...)
  horizon_suffix <- function(h) {
    if (h < 0) paste0("_tm", abs(h)) else paste0("_t", h)
  }

  # Return p_variables for horizon h, or NULL if not specified
  # (h = 0 -> "t0"; "tm0" is also accepted for backward compatibility)
  p_variables_at <- function(h) {
    if (h == 0) return(unique(c(p_variables[["t0"]], p_variables[["tm0"]])))
    key <- if (h < 0) paste0("tm", abs(h)) else paste0("t", h)
    p_variables[[key]]
  }

  # Add columns <cols><suffix> shifted by h within unit (lag if h < 0,
  # lead if h > 0, plain copy if h = 0). Does nothing if 'cols' is empty.
  # 'cols' is de-duplicated: duplicated LHS names in := are handled
  # differently (error / silent misassignment) across data.table versions.
  add_shifted <- function(dt, cols, h, suffix) {
    cols <- unique(cols)
    if (length(cols) == 0) return(invisible(dt))
    new_cols <- paste0(cols, suffix)
    if (h == 0) {
      dt[, (new_cols) := .SD, .SDcols = cols]
    } else {
      dt[, (new_cols) := shift(.SD, n = abs(h),
                               type = if (h < 0) "lag" else "lead"),
         .SDcols = cols, by = c(col_unit)]
    }
    invisible(dt)
  }

  # Cast a vector to the storage type of the time column (needed for joins)
  as_time <- function(x) if (is.integer(df[[col_time]])) as.integer(x) else as.numeric(x)

  # Lags of the dependent variable present in dt (<dep>_tm1, <dep>_tm2, ...)
  dep_lag_cols <- function(dt) {
    pre <- paste0(col_dep, "_tm")
    nm  <- names(dt)
    nm[startsWith(nm, pre) & grepl("^[0-9]+$", substring(nm, nchar(pre) + 1))]
  }

  # ── Setup ──────────────────────────────────────────────────────────────────
  df <- as.data.table(data)   # copy: the user's object is left untouched

  col_unit     <- unit
  col_time     <- time
  col_dep      <- dependent
  col_trt      <- treat
  col_ltr      <- last_treat
  col_h_variables <- h_variables
  col_time_th  <- paste0(col_time, "_th")

  generated <- intersect(c("first.treat", "last.treat", "next.treat", "dtreat"),
                         names(df))
  # Inputs named like a generated column (e.g. last_treat = "first.treat")
  # are read from an internal copy, dropped before returning
  internal <- character(0)
  if (!is.null(col_ltr) && col_ltr %in% generated) {
    df[, .last_treat_input := get(col_ltr)]
    col_ltr  <- ".last_treat_input"
    internal <- c(internal, col_ltr)
  }
  if (!is.null(col_trt) && col_trt %in% generated) {
    df[, .treat_input := get(col_trt)]
    col_trt  <- ".treat_input"
    internal <- c(internal, col_trt)
  }
  if (length(generated) > 0) {
    warning("Existing column(s) overwritten: ", paste(generated, collapse = ", "))
    df[, (generated) := NULL]
  }

  setorderv(df, c(col_unit, col_time))
  df[, rownumber := seq_len(.N)]

  # Duplicated unit x time
  df[, dtime := get(col_time) - shift(get(col_time)) - 1, by = c(col_unit)]
  if (df[, any(dtime < 0, na.rm = TRUE)]) {
    print(df[, .N, by = dtime])
    stop("Some units have multiple observations for the same time period! (e.g., dtime = -1,  where dtime = time - lag(time) - 1)")
  }

  # ── 1. Treatment onsets (event table) ──────────────────────────────────────
  # One row per unit x onset, with the interval [lo, hi] in which the onset
  # lies (lo = hi when the onset date is known) and certain = (lo == hi).

  if (!is.null(col_ltr)) {

    ltr <- df[[col_ltr]]
    if (any(ltr %in% 0) && all(df[[col_time]] > 0)) {
      message("'", col_ltr, "' = 0 is read as never treated (NA).")
      df[get(col_ltr) %in% 0, (col_ltr) := NA]
    }

    if (!absorbing) {
      if (df[, any(get(col_ltr) > get(col_time), na.rm = TRUE)])
        stop("'", col_ltr, "' must not be later than '", col_time,
             "' (most recent onset known at t) when absorbing = FALSE.")
      dec <- df[!is.na(get(col_ltr)),
                any(diff(get(col_ltr)) < 0), by = c(col_unit)][, any(V1)]
      if (dec) stop("'", col_ltr, "' must be non-decreasing over time within units.")
    }

    ev <- unique(df[!is.na(get(col_ltr)),
                    .(u = get(col_unit), lo = as.numeric(get(col_ltr)))])
    ev[, `:=`(hi = lo, certain = TRUE)]

    # Consistency with treat, if available
    if (!is.null(col_trt)) {
      n_inc <- df[get(col_trt) %in% 1 & is.na(get(col_ltr)), .N]
      if (absorbing)
        n_inc <- n_inc + df[get(col_trt) %in% 0 & !is.na(get(col_ltr)) &
                              get(col_ltr) <= get(col_time), .N]
      if (n_inc > 0)
        warning(n_inc, " observation(s) where '", col_trt, "' and '", col_ltr,
                "' are inconsistent. Onsets are taken from '", col_ltr, "'.")
    }

  } else {

    # Onsets from the treatment indicator: 0 -> 1 switches between
    # consecutive non-missing observations
    obs <- df[!is.na(get(col_trt)),
              .(u = get(col_unit), t = as.numeric(get(col_time)),
                d = get(col_trt))]
    obs[, `:=`(pd = shift(d), pt = shift(t)), by = u]

    ev <- obs[d == 1 & (is.na(pd) | pd == 0 |
                          (!absorbing & pd == 1 & t - pt > 1))]
    ev[, `:=`(lo      = fifelse(is.na(pd), t, pt + 1),
              hi      = t,
              certain = !is.na(pd) & pd == 0 & t - pt == 1)]

    # Warnings on onsets made ambiguous by missing periods
    # (units already treated at their first observation are not counted:
    # their onset is unknown with or without missing periods)
    amb_units <- ev[!certain & !is.na(pd), uniqueN(u)]
    if (amb_units > 0) {
      warning(amb_units, " unit(s) have a treatment onset hidden in a gap of ",
              "observation (missing periods or missing '", col_trt, "'). ",
              "Their onset date is unknown: dtreat is set to NA over the gap. ",
              "Consider supplying 'last_treat'.")
    }
    if (!absorbing) {
      gap_units <- obs[!is.na(pt) & t - pt > 1, uniqueN(u)]
      if (gap_units > 0) {
        warning(gap_units, " unit(s) have missing periods (or missing '",
                col_trt, "'). With absorbing = FALSE, a one-period event ",
                "occurring in such a gap cannot be detected, and the unit may ",
                "wrongly serve as a control. Consider supplying 'last_treat'.")
      }
    }
    ev <- ev[, .(u, lo, hi, certain)]
    rm(obs)
  }

  setorder(ev, u, lo)

  if (absorbing) {
    n_multi <- ev[, .N, by = u][N > 1, .N]
    if (n_multi > 0)
      warning(n_multi, " unit(s) have several treatment onsets while ",
              "absorbing = TRUE: only the first one is used.")
    ev <- ev[, .SD[1], by = u]
  }

  if (nrow(ev) == 0) stop("No treatment onset found in the data.")

  # ── 2. Completing for missing years ───────────────────────────────────────
  # Periods are added between the first appearance of each unit and the last
  # period of the panel: a row after the last appearance can still be the
  # reference row of a pre-treatment horizon (y(t+h) - y(t-1), h <= -2, only
  # uses earlier periods). Periods before the first appearance are useless
  # (no observed t-1) and are not added.
  t_max <- max(df[[col_time]])
  df[, incomplete := (sum(dtime, na.rm = TRUE) > 0 |
                        max(get(col_time)) < t_max) * 1, by = c(col_unit)]

  if (df[, any(incomplete == 1)]) {

    message("Some units have incomplete time periods.")

    print(df[, .N, by = incomplete])

    message("Missing time periods were added (between the first appearance ",
            "of each unit and the last period of the panel).")

    # All unit x time combinations from first appearance to the panel's end
    df_in   <- df[incomplete == 1]
    df_comp <- df_in[, .(tt = seq(min(get(col_time)), t_max)),
                     by = c(col_unit)]
    setnames(df_comp, c(col_unit, col_time))
    df_in2  <- merge(df_comp, df_in, by = c(col_unit, col_time), all.x = TRUE)

    # Copy unit-constant columns to the added rows
    if (keep_filled) {
      cand <- setdiff(names(df_in2),
                      c(col_unit, col_time, col_dep, col_trt, col_ltr,
                        "rownumber", "dtime", "incomplete"))
      if (length(cand) > 0) {
        fill_const <- function(x) {
          if (!anyNA(x)) return(x)
          u <- unique(x[!is.na(x)])
          if (length(u) == 1L) x[is.na(x)] <- u
          x
        }
        df_in2[, (cand) := lapply(.SD, fill_const),
               by = c(col_unit), .SDcols = cand]
      }
    }

    df <- rbind(df[incomplete == 0], df_in2)
    setorderv(df, c(col_unit, col_time))

    rm(df_in, df_comp, df_in2)
  }
  df[, incomplete := NULL]

  df[, dtime := NULL]

  # ── 3. Treatment timing variables ──────────────────────────────────────────
  # Rolling joins on the event table (time is on a unit-step grid):
  #   last.treat = latest hi  with hi < t   (key hi + 1 <= t, roll forward)
  #   next.treat = earliest lo with lo > t  (key lo - 1 >= t, roll backward)
  join_ev <- function(keys, vals, roll) {
    tab <- data.table(keys$u, as_time(keys$k), vals)
    setnames(tab, c(col_unit, col_time, "val"))
    tab[df, on = c(col_unit, col_time), roll = roll, x.val]
  }

  first <- ev[, .(first.treat = lo[1]), by = u]
  setnames(first, "u", col_unit)
  df[first, first.treat := i.first.treat, on = c(col_unit)]

  df[, last.treat := join_ev(list(u = ev$u, k = ev$hi + 1), ev$hi,  Inf)]
  df[, next.treat := join_ev(list(u = ev$u, k = ev$lo - 1), ev$lo, -Inf)]

  # Onset at t (certain onsets only)
  evc <- ev[certain == TRUE]
  onset <- !is.na(join_ev(list(u = evc$u, k = evc$lo), evc$lo, FALSE))

  # Inside the interval of an uncertain onset
  evu <- ev[certain == FALSE]
  amb_hi <- join_ev(list(u = evu$u, k = evu$lo), evu$hi, Inf)
  in_amb <- !is.na(amb_hi) & df[[col_time]] <= amb_hi

  df[, dtreat := as.numeric(onset)]
  df[in_amb, dtreat := NA_real_]
  if (!absorbing && !is.null(col_trt)) {
    # ongoing treatment status (treat = 1 without onset)
    df[get(col_trt) %in% 1 & !onset, dtreat := NA_real_]
  }

  # ── 4. Available horizons ──────────────────────────────────────────────────
  obs_rows <- !is.na(df$rownumber)
  if (absorbing) {
    t2t <- (df[[col_time]] - df$first.treat)[obs_rows]
    max_n_post <- max(t2t, na.rm = TRUE)
    max_n_pre  <- -min(t2t, na.rm = TRUE)
  } else {
    max_n_post <- max((df[[col_time]] - df$last.treat)[obs_rows], na.rm = TRUE)
    max_n_pre  <- max((df$next.treat - df[[col_time]])[obs_rows], na.rm = TRUE)
  }

  # Cap n_post / n_pre at what the data actually supports
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

  # Rows kept as reference rows: observed rows, plus (keep_filled) added rows
  # with at least one observed lag of the dependent variable
  keep_rows <- function(dt) {
    lags <- dep_lag_cols(dt)
    has_lag <- if (length(lags) > 0)
      rowSums(!is.na(dt[, ..lags])) > 0 else rep(FALSE, nrow(dt))
    if (keep_filled) !is.na(dt$rownumber) | has_lag else !is.na(dt$rownumber)
  }

  # ── 5. Wide horizon format ─────────────────────────────────────────────────
  if (type_horizon == "wide")
  {
    # 5.1 Merge the t-1 slice: reference outcome and reference-period
    #     p_variables
    cols_tm1 <- unique(c(col_dep, p_variables_at(-1), col_h_variables))
    add_shifted(df, cols_tm1, -1, "_tm1")

    # 5.2 Leads & lags of the outcome and horizon-specific controls
    # Horizon -1 is the reference period and is excluded by convention
    horizons <- setdiff(-abs(n_pre):n_post, -1L)
    message("Computing horizons: ", paste(horizons, collapse = ", "))

    for (h in horizons) {
      p_var  <- p_variables_at(h)
      cols_h <- unique(c(col_dep, p_var, col_h_variables))
      add_shifted(df, cols_h, h, horizon_suffix(h))
    }

    df <- df[keep_rows(df)]
    n_added <- df[is.na(rownumber), .N]
    if (n_added > 0)
      message(n_added, " added row(s) kept as reference rows (observed lag of '",
              col_dep, "').")
    df[, rownumber := NULL]
    if (length(internal) > 0) df[, (internal) := NULL]

    return(df[])

    # ── 6. Long horizon format ─────────────────────────────────────────────────
  } else if (type_horizon == "long") {

    # 6.1 Merge the t-1 slice for the reference outcome and reference-period
    #     p_variables
    cols_tm1 <- unique(c(col_dep, p_variables_at(-1), h_variables))
    add_shifted(df, cols_tm1, -1, "_tm1")

    # 6.2 Wide merge: attach horizon-specific p_variables (all horizons != -1)
    w_th <- as.numeric(gsub("t", "", gsub("tm", "-", names(p_variables))))
    w_th <- unique(w_th[w_th != -1])

    for (h in w_th) {
      add_shifted(df, p_variables_at(h), h, horizon_suffix(h))
    }

    # 6.3 Long (cartesian) merge: cross each reference row with all horizon rows
    horizons <- setdiff(-abs(n_pre):n_post, -1L)

    # Horizon-side table (outcome and h_variables at each observed period)
    cols_th <- intersect(unique(c(col_unit, col_time, col_dep, h_variables)),
                         names(df))
    df_th   <- df[!is.na(get(col_dep)), ..cols_th]
    setnames(df_th,
             setdiff(cols_th, col_unit),
             paste0(setdiff(cols_th, col_unit), "_th"))

    # Reference side: rows with a defined dtreat and an observed lag
    lags <- dep_lag_cols(df)
    has_lag <- if (length(lags) > 0) rowSums(!is.na(df[, ..lags])) > 0 else FALSE
    df_ref <- df[keep_rows(df) & has_lag & !is.na(dtreat)]

    # Merge on unit only -> all (reference period, horizon period) combinations
    df_l <- merge(df_ref, df_th,
                  by              = col_unit,
                  all.x           = TRUE,
                  allow.cartesian = TRUE)

    # 6.4 Compute horizon and filter to the desired window
    # horizon = (horizon period) - (reference period); -1 kept as reference
    df_l[, horizon := get(col_time_th) - get(col_time)]
    df_l[, (col_time_th) := NULL]

    df_l <- df_l[horizon %in% c(-1, horizons)]

    df_l <- df_l[order(horizon, get(col_unit), get(col_time))]
    df_l[, rownumber := NULL]
    if (length(internal) > 0) df_l[, (internal) := NULL]

    return(df_l[])
  }
}
