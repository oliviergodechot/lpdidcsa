#' mpdta_r dataset
#'
#' A toy dataset based on the `mpdta` dataset from the `did` package, modified for use in the `lpdidcsa` package.
#'   The original `mpdta` dataset is described in the `did` package documentation:
#'   https://cran.r-project.org/web/packages/did/index.html
#'   As explained on this page, "is a simplified example of the effect of states increasing 
#'   their minimum wages on county-level teen employment rates which comes from 
#'   Callaway and Sant’Anna (2021). The dataset contains 500 observations of county-level 
#'   teen employment rates from 2003-2007. Some states are first treated in 2004, 
#'   some in 2006, and some in 2007 (see the paper for more details)".
#'
#' @format A data frame with 500 rows and the following columns:
#' \describe{
#'   \item{lemp}{Log of county-level teen employment (outcome variable)}
#'   \item{first.treat}{Year of first minimum wage increase (2004, 2006, or 2007)}
#'   \item{year}{Year (time variable)}
#'   \item{countyreal}{County ID (individual identifier)}
#'   \item{treat_back}{Original treatment variable from `did::mpdta`}
#'   \item{treat}{Absorbing treatment dummy (1 if treated at any point)}
#'   \item{lpop}{Log population size}
#'   \item{treat_nabs}{Non-absorbing treatment dummy (1 only in the treatment year)}
#'   \item{lpop2}{Squared log population (for nonlinear models)}
#'   \item{noise}{Random noise for robustness checks}
#'   \item{drop}{Random drop indicator (1 if dropped, 0 otherwise)}
#' }
#' @source Original dataset: \url{https://cran.r-project.org/web/packages/did/index.html}
"mpdta_r"