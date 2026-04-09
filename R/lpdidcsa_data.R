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
#' @param time Character. Name of the time period column (default: \code{"year"}).
#' @param dependent Character. Name of the dependent variable column.
#' @param treat Character. Name of the binary treatment indicator column
#'   (must contain only 0 and 1; default: \code{"treat"}).
#' @param absorbing Logical. If \code{TRUE} (default), treatment is assumed
#'   absorbing (once treated, always treated). If \code{FALSE}, treatment can
#'   switch on and off across periods.
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
#'   \code{list()}) when no horizon-specific controls are needed.
#' @param h_variables Character vector of column names to be carried through at
#'   every horizon (e.g. time-varying controls). Use \code{NULL} (default) if
#'   not needed.
#' @param type_horizon Character. Output structure: \code{"wide"} (default)
#'   creates one column per horizon; \code{"long"} stacks all horizons into a
#'   single \code{horizon} column.
#'
#' @return A \code{data.table} with the original columns plus:
#'   \describe{
#'     \item{\code{first.treat}}{First treatment period per unit (\code{NA}
#'       for never-treated). Only present when \code{absorbing = TRUE}.}
#'     \item{\code{last.treat}}{Most recent treatment period as of \emph{t}
#'       (\code{NA} if never treated yet). Only present when \code{absorbing = FALSE}.}
#'     \item{\code{next.treat}}{Next upcoming treatment period (\code{NA} if
#'       no future treatment). Only present when \code{absorbing = FALSE}.}
#'     \item{\code{dtreat}}{Change in treatment status:
#'       \eqn{D_{it} - D_{i,t-1}}.}
#'     \item{Horizon columns}{For \code{type_horizon = "wide"}: outcome and
#'       control columns at each horizon, e.g. \code{lemp_tm2}, \code{lemp_t0}.
#'       Horizon \emph{-1} is always omitted (reference period).}
#'     \item{\code{horizon}}{For \code{type_horizon = "long"}: integer giving
#'       the horizon \eqn{h = t' - t} for each merged row.}
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
                          treat       = "treat",
                          unit        = NULL,
                          dependent   = NULL,
                          n_pre      = NULL,
                          n_post     = NULL,
                          absorbing   = TRUE,
                          p_variables = list(),                          
                          h_variables = NULL,
                          type_horizon= "wide")
{
  # ── Input validation ───────────────────────────────────────────────────────
  stopifnot(
    is.data.frame(data),
    is.character(unit)      && length(unit)      == 1,
    is.character(time)      && length(time)      == 1,
    is.character(dependent) && length(dependent) == 1,
    is.character(treat)     && length(treat)     == 1,
    is.list(p_variables)
  )
  
  required_cols <- c(unit, time, dependent, treat,
                     unlist(p_variables, use.names = FALSE))
  
  missing_cols  <- setdiff(required_cols, names(data))
  
  if (length(missing_cols) > 0) {
    stop("The following columns are missing from 'data': ",
         paste(missing_cols, collapse = ", "))
  }
  
  # ── Helpers ────────────────────────────────────────────────────────────────
  
  # Build the column suffix for a given horizon h ("_tm2", "_t0", "_t3", ...)
  horizon_suffix <- function(h) {
    if (h < 0) paste0("_tm", abs(h)) else paste0("_t", h)
  }
  
  # Return p_variables for horizon h, or NULL if not specified
  p_variables_at <- function(h) {
    key <- if (h <= 0) {paste0("tm", abs(h))} else {paste0("t", h)}
    p_variables[[key]]
  }
  
  # ── Setup ──────────────────────────────────────────────────────────────────
  df <- copy(setDT(data))
  
  col_unit     <- unit
  col_time     <- time
  col_dep      <- dependent
  col_trt      <- treat
  col_h_variables <- h_variables
  col_time_th  <- paste0(col_time, "_th")
  col_dep_th   <- paste0(col_dep,  "_th")
  col_trt_tm1  <- paste0(col_trt,  "_tm1")
  col_dep_tm1  <- paste0(col_dep,  "_tm1")
  
  setorderv(data, c(col_unit, col_time))
  
  # ── 1. Treatment timing variables ──────────────────────────────────────────
  if (absorbing == TRUE)
  {
    # First treatment period per unit (NA for never-treated)
    df[, first.treat := min(get(col_time) * fifelse(get(col_trt) == 0, NA, 1),
                            na.rm = TRUE),
       by = list(get(col_unit))]
    df[, first.treat := fifelse(is.infinite(first.treat), NA, first.treat)]
    
    df[, time2treat  := fifelse(is.na(first.treat), NA, get(col_time) - first.treat)]
    
    # # Last treatment period: most recent t at which treatment switched on
    # df[, last.treat := fifelse(shift(get(col_time), type = "lag") == first.treat,
    #                            shift(get(col_time), type = "lag"), NA_real_),
    #    by = list(get(col_unit))]
    # df[, last.treat := nafill(last.treat, type = "locf"), by = list(get(col_unit))]
    # df[, time2lasttreat  := fifelse(is.na(last.treat), NA, get(col_time) - last.treat)]
    # 
    # # Next treatment period: next t at which treatment will switch on
    # df[, next.treat := fifelse(shift(get(col_time), type = "lead") == first.treat,
    #                            shift(get(col_time), type = "lead"), NA_real_),
    #    by = list(get(col_unit))]
    # df[, next.treat := nafill(next.treat, type = "nocb"), by = list(get(col_unit))]
    # df[, time2nexttreat  := fifelse(is.na(next.treat), NA, get(col_time) - next.treat)]
    
    max_n_post <- df[, max(time2treat, na.rm = TRUE)]
    max_n_pre  <- -df[, min(time2treat, na.rm = TRUE)]
    
  } else {
    
    # Non-absorbing case: last/next treatment defined from the treatment indicator
    df[, last.treat := fifelse(shift(get(col_trt), type = "lag") == 1,
                               shift(get(col_time), type = "lag"), NA_real_),
       by = list(get(col_unit))]
    df[, last.treat := nafill(last.treat, type = "locf"), by = list(get(col_unit))]
    df[, time2lasttreat  := fifelse(is.na(last.treat), NA,
                                    get(col_time) - last.treat)]
    
    df[, next.treat := fifelse(shift(get(col_trt), type = "lead") == 1,
                               shift(get(col_time), type = "lead"), NA_real_),
       by = list(get(col_unit))]
    df[, next.treat := nafill(next.treat, type = "nocb"), by = list(get(col_unit))]
    df[, time2nexttreat  := fifelse(is.na(next.treat), NA,
                                    get(col_time) - next.treat)]
    
    max_n_post <- df[, max(time2lasttreat, na.rm = TRUE)]
    max_n_pre  <- df[, max(abs(time2nexttreat), na.rm = TRUE)]
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
  
  df[, c("time2treat", "time2nexttreat", "time2lasttreat") := NULL]
  
  # ── 2. Completing for missing years ───────────────────────────────────────
  
  setorderv(df, cols=c(col_unit,col_time))
  
  df[,rownumber:=c(1L:nrow(df))]
  
  df[, dtime:=eval(as.name(col_time))-shift(eval(as.name(col_time)))-1,by=col_unit]
  sdtime <- df[,sum(abs(dtime),na.rm=T)]
  
  if (sdtime>0){
    mindtime <- df[,min(dtime,na.rm=T)]
    
    if (mindtime < 0) {
      print(df[,.N,by=dtime])
      stop("Some units have multiple observations for the same time period! (e.g., dtime = -1,  where dtime = time - lag(time) - 1)")
    }
    
    message("Some units have incomplete time periods.")
    
    df[, incomplete:=((sum(dtime,na.rm=T))>0)*1,by=col_unit]
    
    print(df[, .N,by=incomplete])
    
    message("Missing intermediary time periods were added.")
    
    # Create a data.table with all county-year combinations
    df_in <- df[incomplete==1,]
    df_comp <- df[incomplete==1, .(time = seq(min(eval(as.name(col_time))), max(eval(as.name(col_time))))), by = col_unit]
    colnames(df_comp) <- c(col_unit,col_time)
    df_in2 <- merge(df_comp,df_in,by=c(col_unit,col_time),all.x=T)
    
    df <- rbind(df[incomplete==0],df_in2)
    setorderv(df, cols=c(col_unit,col_time))
    
    rm(df_in,df_com,df_in2)
    df[,incomplete:=NULL]
  }
  
  df[,dtime:=NULL]
  
  # ── 3. Wide horizon format ─────────────────────────────────────────────────
  if (type_horizon == "wide")
  {
    # 3.1 Merge the t-1 slice: used to compute dtreat and to attach
    #     reference-period p_variables
    cols_tm1   <- c(col_trt, col_dep,
                    p_variables_at(-1), col_h_variables)
    h <- -1
    
    df[, paste0(cols_tm1, "_tm1") := shift(.SD, n = abs(h), 
                                           type = ifelse(h < 0, "lag", "lead")),
       .SDcols = cols_tm1,
       by = col_unit]
    
    # 3.2 Treatment change indicator: 1 if unit switches from 0 to 1
    col_trt_tm1 <- paste0(col_trt, "_tm1")
    df[, dtreat        := (get(col_trt) - get(col_trt_tm1) == 1) * 1]
    df[, (col_trt_tm1) := NULL]
    
    # 3.3 Leads & lags of the outcome and horizon-specific controls
    # Horizon -1 is the reference period and is excluded by convention
    horizons <- setdiff(-abs(n_pre):n_post, -1L)
    message("Computing horizons: ", paste(horizons, collapse = ", "))
    
    for (h in horizons) {
      
      p_var   <- p_variables_at(h)
      cols_h  <- unique(c(col_dep, p_var, col_h_variables))
      suffix <- horizon_suffix(h)
      
      df[, paste0(cols_h, suffix) := shift(.SD, n = abs(h), 
                                           type = ifelse(h < 0, "lag", "lead")),
         .SDcols = cols_h,
         by = col_unit]
      
    }
    
    df <- df[!is.na(rownumber),]
    df <- df[order(rownumber),]
    df[,rownumber:=NULL]
    
    return(df)
    
    # ── 4. Long horizon format ─────────────────────────────────────────────────
  } else if (type_horizon == "long") {
    
    # 4.1 Merge the t-1 slice for dtreat and reference-period p_variables
    cols_tm1   <- unique(c(col_dep, col_trt,
                           p_variables_at(-1), h_variables))
    h <- -1
    
    df[, paste0(cols_tm1, "_tm1") := shift(.SD, n = abs(h), 
                                           type = ifelse(h < 0, "lag", "lead")),
       .SDcols = cols_tm1,
       by = col_unit]
    
    # 4.2 Wide merge: attach horizon-specific p_variables (all horizons != -1)
    w_th <- as.numeric(gsub("t", "", gsub("tm", "-", names(p_variables))))
    w_th <- w_th[w_th != -1]
    
    if (length(w_th) > 0) {
      for (i in seq_along(w_th)) {
        cols_th   <- unique(c(p_variables_at(w_th[i])))
        h <- w_th[i]
        
        suffix <- paste0("_t",gsub("-", "m", w_th[i]))
        
        df[, paste0(cols_th, suffix) := shift(.SD, n = abs(h), 
                                              type = ifelse(h < 0, "lag", "lead")),
           .SDcols = cols_th,
           by = col_unit]
        
      }
    }
    
    # 4.3 Long (cartesian) merge: cross each reference row with all horizon rows
    # Reference rows are anchored at t-1 relative to treatment (or untreated).
    # The right-hand side contains all periods, used as horizon observations.
    horizons  <- setdiff(-abs(n_pre):n_post, -1L)
    all_p_var <- unique(unlist(p_variables, use.names = FALSE))
    
    cols_ref <- intersect(
      c(col_unit, col_time, "first.treat", "last.treat", "next.treat",
        "dtreat", col_dep_tm1, col_trt_tm1, col_trt, all_p_var),
      names(df)
    )
    
    # Prepare the horizon-side table (outcome and h_variables at each period)
    cols_th  <- c(col_unit, col_time, col_dep, h_variables)
    cols_th  <- intersect(cols_th, names(df))
    df_th    <- df[, ..cols_th]
    setnames(df_th,
             setdiff(cols_th, col_unit),
             paste0(setdiff(cols_th, col_unit), "_th"))
    
    # Merge on unit only -> all (reference period, horizon period) combinations
    df_l <- merge(df, df_th,
                  by             = col_unit,
                  all.x          = TRUE,
                  allow.cartesian = TRUE)
    
    # 4.4 Compute horizon and filter to the desired window
    # horizon = (horizon period) - (reference period); -1 excluded by convention
    df_l[, horizon := get(col_time_th) - get(col_time)]
    df_l[, dtreat  := get(col_trt) - get(col_trt_tm1)]
    df_l[, (col_trt_tm1)  := NULL]
    df_l[, (col_time_th)  := NULL]
    
    df_l <- df_l[horizon %in% c(-1, horizons) & !is.na(dtreat)]
    
    df_l <- df_l[!is.na(rownumber),]
    df_l <- df_l[order(horizon,rownumber),]
    df_l[,rownumber:=NULL]
    
    
    return(df_l)
  }
}
