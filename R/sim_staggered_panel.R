# =============================================================================
# sim_staggered_panel.R
# For use in the lpdidcsa package
# =============================================================================

# -----------------------------------------------------------------------------
# Internal helpers (unexported)
# -----------------------------------------------------------------------------

.gen_firms <- function(n_firms, firm_size_meanlog, firm_size_sdlog) {
  data.table::data.table(
    firm_id = 1L:n_firms,
    size    = stats::rlnorm(n_firms,
                            meanlog = firm_size_meanlog,
                            sdlog   = firm_size_sdlog),
    psi_raw = stats::rnorm(n_firms)
  )
}

.gen_individuals <- function(n, p_treat, t_min_treat, t_max_treat,
                             early_cutoff, sigma_alpha,
                             p_female, female_treat_rr) {
  female  <- stats::rbinom(n, 1L, p_female)
  alpha_i <- stats::rnorm(n, mean = -0.10 * female, sd = sigma_alpha)
  
  # Back-solve gender-specific treatment probs to hit the marginal p_treat
  share_f      <- mean(female)
  p_m          <- p_treat / (female_treat_rr * share_f + (1 - share_f))
  p_f          <- min(female_treat_rr * p_m, 1)
  ever_treated <- stats::rbinom(n, 1L, ifelse(female == 1L, p_f, p_m))
  
  g_i <- rep(NA_integer_, n)
  n_t <- sum(ever_treated)
  if (n_t > 0L)
    g_i[ever_treated == 1L] <- sample(t_min_treat:t_max_treat,
                                      size = n_t, replace = TRUE)
  
  data.table::data.table(
    id           = 1L:n,
    female       = female,
    alpha_i      = alpha_i,
    ever_treated = ever_treated,
    g_i          = g_i,
    early_cohort = data.table::fifelse(
      !is.na(g_i) & g_i <= early_cutoff, 1L, 0L)
  )
}

.gen_firm_spells <- function(n, Tmax, n_firms, prob_j, p_move) {
  firm_mat       <- matrix(0L, nrow = n, ncol = Tmax)
  firm_mat[, 1L] <- sample(1L:n_firms, size = n, replace = TRUE, prob = prob_j)
  
  for (tt in 2L:Tmax) {
    prev   <- firm_mat[, tt - 1L]
    movers <- which(stats::rbinom(n, 1L, p_move) == 1L)
    new_f  <- prev
    if (length(movers) > 0L)
      new_f[movers] <- sample(1L:n_firms, size = length(movers),
                              replace = TRUE, prob = prob_j)
    firm_mat[, tt] <- new_f
  }
  
  data.table::data.table(
    id      = rep(1L:n,    each  = Tmax),
    t       = rep(1L:Tmax, times = n),
    firm_id = as.integer(t(firm_mat))
  )
}

.get_att_vec <- function(event_time, early_cohort, att_profile,
                         early_multiplier) {
  # k=0 convention: event_time == 0 is the first treated period.
  max_h    <- length(att_profile)
  et_idx   <- pmin(pmax(event_time + 1L, 0L), max_h)
  base_att <- c(0, att_profile)[et_idx + 1L]
  base_att[is.na(event_time) | event_time < 0L] <- 0
  base_att * data.table::fifelse(early_cohort == 1L, early_multiplier, 1.0)
}

# -----------------------------------------------------------------------------
#' Simulate a Staggered Adoption Panel Dataset
#'
#' Generates a balanced individual-level panel suitable for evaluating
#' staggered difference-in-differences and event-study estimators. The DGP
#' features individual and firm fixed effects, gender heterogeneity in earnings
#' growth and treatment selection, cohort-heterogeneous treatment effects, and
#' realistic worker mobility across firms of heterogeneous size.
#'
#' @param n Integer. Number of individuals. Default: \code{50000}.
#' @param Tmax Integer. Number of periods (indexed \code{1:Tmax}).
#'   Default: \code{20}.
#' @param p_treat Numeric in (0, 1). Share of individuals ever treated.
#'   Default: \code{0.10}.
#' @param t_min_treat Integer. Earliest possible treatment cohort.
#'   Default: \code{1}.
#' @param t_max_treat Integer. Latest possible treatment cohort.
#'   Must satisfy \code{t_max_treat <= Tmax}. Default: \code{20}.
#' @param att_profile Numeric vector. Average treatment effect (in log-points)
#'   at event-times \eqn{k = 0, 1, 2, \ldots} (\strong{k=0 convention}: k=0
#'   is the first period of treatment, i.e. \eqn{t = g_i}). Element \code{j}
#'   of the vector applies to event-time \eqn{k = j - 1}. The last element
#'   governs all subsequent periods. Negative values indicate earnings
#'   reductions. Default: \code{c(-0.02, -0.03, -0.037, -0.039, -0.04)}.
#' @param early_cutoff Integer. Cohorts with \eqn{g_i \le} \code{early_cutoff}
#'   are labelled "early" and receive a scaled-up ATT. Default: \code{10}.
#' @param early_multiplier Numeric. Multiplicative scaling of \code{att_profile}
#'   applied to early cohorts. Default: \code{2.0}.
#' @param sigma_alpha Numeric. Standard deviation of the individual fixed
#'   effect \eqn{\alpha_i}. Default: \code{0.30}.
#' @param p_female Numeric in (0, 1). Population share of women.
#'   Default: \code{0.50}.
#' @param female_treat_rr Numeric. Relative risk of treatment for women vs.
#'   men; the overall treated share is held fixed at \code{p_treat}.
#'   Default: \code{2.0}.
#' @param male_growth Numeric. Annual log-earnings growth for men
#'   (log-points per period). Default: \code{0.020}.
#' @param female_growth Numeric. Annual log-earnings growth for women.
#'   Default: \code{0.012}.
#' @param n_firms Integer. Number of firms. Default: \code{500}.
#' @param firm_size_meanlog Numeric. \code{meanlog} of the log-normal
#'   distribution used to draw firm sizes (allocation weights).
#'   Default: \code{4.0}.
#' @param firm_size_sdlog Numeric. \code{sdlog} of the log-normal firm-size
#'   distribution. Larger values produce a fatter size tail. Default: \code{1.2}.
#' @param firm_fe_var_share Numeric in (0, 1). Target share of total
#'   \code{log_earnings} variance attributable to firm fixed effects.
#'   The firm FE standard deviation is calibrated analytically to hit this
#'   target. Default: \code{0.10}.
#' @param p_move Numeric in (0, 1). Per-period probability that a worker
#'   changes firm. The destination is drawn proportional to firm size (no
#'   sorting on firm FEs). Default: \code{0.05}.
#' @param base_log_earnings Numeric. Log-earnings intercept at period 1.
#'   Default: \code{10.0}.
#' @param sigma_eps Numeric. Standard deviation of the idiosyncratic iid
#'   shock \eqn{\varepsilon_{it}}. Default: \code{0.15}.
#' @param seed Integer or \code{NULL}. Passed to \code{\link{set.seed}}.
#'   Use \code{NULL} to skip seeding (e.g. inside a Monte Carlo loop that
#'   manages its own seed). Default: \code{42}.
#' @param verbose Logical. If \code{TRUE}, prints sanity-check tables
#'   (treatment shares, cohort distribution, true ATT by event-time, firm
#'   FE variance, mobility rate) and renders four diagnostic ggplot2 figures.
#'   Default: \code{TRUE}.
#'
#' @return A \code{data.table} with one row per (individual \eqn{\times}
#'   period) and the following columns, returned \strong{invisibly}
#'   (so that diagnostic output is visible at the console without printing
#'   the full table):
#' \describe{
#'   \item{\code{id}}{Individual identifier (integer).}
#'   \item{\code{t}}{Period (integer, \code{1:Tmax}).}
#'   \item{\code{firm_id}}{Firm identifier in period \code{t} (integer).}
#'   \item{\code{female}}{Gender indicator: 1 = woman, 0 = man (integer).}
#'   \item{\code{ever_treated}}{1 if the individual is ever treated (integer).}
#'   \item{\code{g_i}}{Treatment cohort (period of first treatment);
#'     \code{NA} for never-treated (integer).}
#'   \item{\code{early_cohort}}{1 if \eqn{g_i \le} \code{early_cutoff}
#'     (integer).}
#'   \item{\code{treat}}{Binary treatment indicator: 1 if \eqn{t \ge g_i}
#'     (integer).}
#'   \item{\code{event_time}}{Periods since treatment (\eqn{t - g_i}); 0 is
#'     the first treated period; \code{NA} for never-treated (integer).}
#'   \item{\code{tau_it}}{True individual-level treatment effect (numeric).}
#'   \item{\code{alpha_i}}{Individual fixed effect (numeric).}
#'   \item{\code{psi}}{Firm fixed effect in period \code{t} (numeric).}
#'   \item{\code{trend}}{Gender-specific linear earnings trend (numeric).}
#'   \item{\code{eps_it}}{Idiosyncratic shock (numeric).}
#'   \item{\code{log_earnings}}{Outcome variable (numeric).}
#' }
#'
#' @details
#' \strong{Outcome equation:}
#' \deqn{
#'   \log w_{it} = \mu + \alpha_i + \psi_{J(i,t)}
#'               + \delta_{\mathrm{gender}(i)} (t - 1)
#'               + \tau_{it} + \varepsilon_{it}
#' }
#' with \eqn{\alpha_i \sim N(-0.1 \cdot \mathbf{1}[\mathrm{female}_i],\,
#' \sigma_\alpha^2)}, gender-specific slopes \eqn{\delta_g}, firm fixed
#' effect \eqn{\psi_j} calibrated so that
#' \eqn{\mathrm{Var}(\psi_{J(i,t)}) / \mathrm{Var}(\log w_{it}) \approx}
#' \code{firm_fe_var_share}, and
#' \eqn{\varepsilon_{it} \overset{iid}{\sim} N(0, \sigma_\varepsilon^2)}.
#'
#' \strong{k=0 convention:} event-time \eqn{k = 0} corresponds to
#' \eqn{t = g_i}. The ATT at \eqn{k \ge 0} equals
#' \code{att_profile[k + 1]} (last element used for \eqn{k} beyond the
#' profile length), multiplied by \code{early_multiplier} for early cohorts.
#' Pre-treatment effects (\eqn{k < 0}) are exactly zero.
#'
#' \strong{Firm mobility:} each period a worker moves to a new firm with
#' probability \code{p_move}; the destination is drawn proportional to firm
#' size, so there is no systematic sorting on firm quality.
#'
#' @examples
#' # Quick run with a small panel
#' panel <- sim_staggered_panel(n = 2000, seed = 1)
#' head(panel)
#'
#' # Change the ATT profile and cohort heterogeneity
#' panel2 <- sim_staggered_panel(
#'   att_profile      = c(-0.05, -0.03, -0.01, 0),
#'   early_multiplier = 1.5,
#'   seed             = 99
#' )
#'
#' @importFrom data.table data.table setorder setcolorder fifelse fcase uniqueN
#' @importFrom stats rnorm rbinom rlnorm var sd
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline geom_vline
#'   geom_histogram geom_bar labs theme_minimal
#' @export
# -----------------------------------------------------------------------------
sim_staggered_panel <- function(
    # --- Panel dimensions ---
  n                 = 50000L,
  Tmax              = 20L,
  # --- Treatment ---
  p_treat           = 0.10,
  t_min_treat       = 1L,
  t_max_treat       = 20L,
  # --- ATT profile (k=0 convention) ---
  att_profile       = c(-0.02, -0.03, -0.037, -0.039, -0.04),
  # --- Cohort heterogeneity ---
  early_cutoff      = 10L,
  early_multiplier  = 2.0,
  # --- Individual heterogeneity ---
  sigma_alpha       = 0.30,
  # --- Gender ---
  p_female          = 0.50,
  female_treat_rr   = 2.0,
  male_growth       = 0.020,
  female_growth     = 0.012,
  # --- Firms ---
  n_firms           = 500L,
  firm_size_meanlog = 4.0,
  firm_size_sdlog   = 1.2,
  firm_fe_var_share = 0.10,
  p_move            = 0.05,
  # --- Earnings ---
  base_log_earnings = 10.0,
  sigma_eps         = 0.15,
  # --- Misc ---
  seed              = 42L,
  verbose           = TRUE
) {
  
  # ---- Input validation -------------------------------------------------------
  stopifnot(
    "`n` must be a positive integer"           = n >= 1L,
    "`Tmax` must be >= 2"                      = Tmax >= 2L,
    "`p_treat` must be in (0, 1)"              = p_treat > 0 & p_treat < 1,
    "`p_female` must be in (0, 1)"             = p_female > 0 & p_female < 1,
    "`t_min_treat` >= 1"                       = t_min_treat >= 1L,
    "`t_max_treat` <= Tmax"                    = t_max_treat <= Tmax,
    "`t_min_treat` <= `t_max_treat`"           = t_min_treat <= t_max_treat,
    "`firm_fe_var_share` must be in (0, 1)"    = firm_fe_var_share > 0 & firm_fe_var_share < 1,
    "`p_move` must be in [0, 1]"               = p_move >= 0 & p_move <= 1,
    "`att_profile` must have length >= 1"      = length(att_profile) >= 1L
  )
  
  if (!is.null(seed)) set.seed(seed)
  
  # ---- Build firms & individuals ----------------------------------------------
  firms <- .gen_firms(n_firms, firm_size_meanlog, firm_size_sdlog)
  indiv <- .gen_individuals(n, p_treat, t_min_treat, t_max_treat,
                            early_cutoff, sigma_alpha,
                            p_female, female_treat_rr)
  
  # ---- Expand to balanced panel -----------------------------------------------
  panel <- indiv[rep(seq_len(nrow(indiv)), each = Tmax)]
  panel[, t := rep(1L:Tmax, times = nrow(indiv))]
  data.table::setorder(panel, id, t)
  
  # ---- Firm spells ------------------------------------------------------------
  prob_j <- firms$size / sum(firms$size)
  spells <- .gen_firm_spells(n, Tmax, n_firms, prob_j, p_move)
  panel  <- merge(panel, spells, by = c("id", "t"), sort = FALSE)
  data.table::setorder(panel, id, t)
  
  # ---- Calibrate firm FEs to hit target variance share -----------------------
  # Var(Y) ≈ Var(alpha) + Var(trend) + Var(eps) + Var(psi)  [covs ≈ 0]
  # => sigma_psi^2 = target / (1 - target) * Var(others) / Var(psi_raw)
  avg_growth <- (male_growth + female_growth) / 2
  var_others <- sigma_alpha^2 +
    avg_growth^2 * mean((0L:(Tmax - 1L))^2) +
    sigma_eps^2
  target_var <- firm_fe_var_share / (1 - firm_fe_var_share) * var_others
  sigma_psi  <- sqrt(target_var / stats::var(firms$psi_raw))
  firms[, psi     := psi_raw * sigma_psi]
  firms[, psi_raw := NULL]
  
  panel <- merge(panel, firms[, .(firm_id, size, psi)],
                 by = "firm_id", sort = FALSE)
  data.table::setorder(panel, id, t)
  
  # ---- Event time, treatment indicator, ATT ----------------------------------
  panel[, event_time := data.table::fifelse(
    ever_treated == 1L, as.integer(t - g_i), NA_integer_)]
  
  panel[, treat := data.table::fifelse(
    !is.na(g_i), as.integer(t >= g_i), 0L)]
  
  panel[, tau_it := .get_att_vec(event_time, early_cohort,
                                 att_profile, early_multiplier)]
  
  # ---- Gender-specific trend & outcome ---------------------------------------
  panel[, trend := data.table::fifelse(
    female == 1L,
    female_growth * (t - 1L),
    male_growth   * (t - 1L))]
  
  panel[, eps_it      := stats::rnorm(.N, 0, sigma_eps)]
  panel[, log_earnings := base_log_earnings + alpha_i + psi +
          trend + tau_it + eps_it]
  
  data.table::setcolorder(
    panel,
    c("id", "t", "firm_id", "female", "ever_treated", "g_i",
      "early_cohort", "treat", "event_time", "tau_it",
      "alpha_i", "psi", "trend", "eps_it", "log_earnings")
  )
  
  # ---- Sanity checks & figures -----------------------------------------------
  if (verbose) {
    
    indiv_sum <- unique(panel[, .(id, female, ever_treated, g_i, early_cohort)])
    
    # -- Text summaries --
    message("=== sim_staggered_panel: DGP summary ===\n")
    
    message("--- Treatment ---")
    cat(sprintf("  N individuals            : %d\n",    nrow(indiv_sum)))
    cat(sprintf("  Share treated            : %.1f%%\n", 100 * indiv_sum[, mean(ever_treated)]))
    cat(sprintf("  Share female             : %.1f%%\n", 100 * indiv_sum[, mean(female)]))
    cat(sprintf("  Share female | treated   : %.1f%%\n",
                100 * indiv_sum[ever_treated == 1L, mean(female)]))
    cat(sprintf("  Share female | untreated : %.1f%%\n",
                100 * indiv_sum[ever_treated == 0L, mean(female)]))
    
    message("\n--- Cohort distribution (treated only) ---")
    print(indiv_sum[ever_treated == 1L,
                    .N, keyby = .(cohort = g_i, early_cohort)])
    
    message("\n--- True ATT by event time (k = -4 to 6) ---")
    print(
      panel[ever_treated == 1L & !is.na(event_time) &
              event_time >= -4L & event_time <= 6L,
            .(mean_tau       = mean(tau_it),
              mean_tau_early = mean(tau_it[early_cohort == 1L]),
              mean_tau_late  = mean(tau_it[early_cohort == 0L])),
            keyby = event_time]
    )
    
    message("\n--- Firms ---")
    cat(sprintf("  N firms                  : %d\n",   nrow(firms)))
    cat(sprintf("  SD of firm FEs (psi)     : %.4f\n", firms[, stats::sd(psi)]))
    cat(sprintf("  Var(psi) / Var(Y)        : %.3f  (target: %.3f)\n",
                panel[, stats::var(psi) / stats::var(log_earnings)],
                firm_fe_var_share))
    
    message("\n--- Worker mobility ---")
    panel[, .moved := as.integer(firm_id != data.table::shift(firm_id, 1L)),
          by = id]
    cat(sprintf("  Annual move rate         : %.1f%%  (target: %.0f%%)\n",
                100 * panel[t > 1L, mean(.moved)],
                100 * p_move))
    cat(sprintf("  Avg distinct firms/worker: %.2f\n",
                panel[, data.table::uniqueN(firm_id), by = id][, mean(V1)]))
    panel[, .moved := NULL]
    
    # -- Figures --
    
    # Figure 1: Average log earnings by group x period
    panel[, .group := data.table::fcase(
      ever_treated == 0L, "Never treated",
      early_cohort == 1L, "Early treated",
      default            = "Late treated"
    )]
    avg_earn <- panel[, .(mean_log_earnings = mean(log_earnings)),
                      keyby = .(t, group = .group)]
    
    fig1 <- ggplot2::ggplot(avg_earn,
                            ggplot2::aes(x = t, y = mean_log_earnings,
                                         color = group)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 1.5) +
      ggplot2::labs(title = "Average log earnings by group and period",
                    x = "Period (t)", y = "Mean log earnings",
                    color = "Group") +
      ggplot2::theme_minimal(base_size = 13)
    
    panel[, .group := NULL]
    
    # Figure 2: True ATT profile by cohort + aggregate
    et_range <- panel[ever_treated == 1L & !is.na(event_time) &
                        event_time >= -5L & event_time <= 20L]
    
    ep_cohort <- et_range[, .(mean_tau = mean(tau_it)),
                          keyby = .(event_time, early_cohort)]
    ep_cohort[, cohort := data.table::fifelse(
      early_cohort == 1L, "Early cohort", "Late cohort")]
    
    ep_agg <- et_range[, .(mean_tau = mean(tau_it), cohort = "Aggregate"),
                       keyby = event_time]
    
    event_plot <- data.table::rbindlist(
      list(ep_cohort[, .(event_time, mean_tau, cohort)],
           ep_agg[,   .(event_time, mean_tau, cohort)]),
      use.names = TRUE
    )
    
    fig2 <- ggplot2::ggplot(event_plot,
                            ggplot2::aes(x = event_time, y = mean_tau,
                                         color = cohort)) +
      ggplot2::geom_hline(yintercept = 0,    linetype = "dashed",
                          color = "grey50") +
      ggplot2::geom_vline(xintercept = -0.5, linetype = "dotted",
                          color = "grey50") +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(title = "True ATT profile by cohort",
                    x     = "Event time (k=0 is first treated period)",
                    y     = "Average treatment effect (log points)",
                    color = "Cohort") +
      ggplot2::theme_minimal(base_size = 13)
    
    # Figure 3: Firm FE distribution
    fig3 <- ggplot2::ggplot(firms, ggplot2::aes(x = psi)) +
      ggplot2::geom_histogram(bins = 40, fill = "steelblue",
                              color = "white", alpha = 0.8) +
      ggplot2::labs(title = "Distribution of firm fixed effects",
                    x = "Firm FE (\u03c8)", y = "Count") +
      ggplot2::theme_minimal(base_size = 13)
    
    # Figure 4: Distinct firms per worker
    n_firms_wkr <- panel[, .(n_firms_seen = data.table::uniqueN(firm_id)),
                         by = id]
    fig4 <- ggplot2::ggplot(n_firms_wkr,
                            ggplot2::aes(x = factor(n_firms_seen))) +
      ggplot2::geom_bar(fill = "steelblue", alpha = 0.8) +
      ggplot2::labs(title = "Number of distinct firms per worker",
                    x = "Distinct firms over the panel",
                    y = "Number of workers") +
      ggplot2::theme_minimal(base_size = 13)
    
    print(fig1)
    print(fig2)
    print(fig3)
    print(fig4)
  }
  
  invisible(panel)
}
