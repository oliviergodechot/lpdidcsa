#' --- title: "Preparation of the mpdta_r Toy Dataset"
#' --- description: This script modifies the `mpdta` dataset from the `did` package
#'   to create a simplified toy dataset (`mpdta_r`) for the `lpdidcsa` package.
#'   The original `mpdta` dataset is described in the `did` package documentation:
#'   https://cran.r-project.org/web/packages/did/index.html
#'   As explained on this page, "is a simplified example of the effect of states increasing 
#'   their minimum wages on county-level teen employment rates which comes from 
#'   Callaway and Sant’Anna (2021). The dataset contains 500 observations of county-level 
#'   teen employment rates from 2003-2007. Some states are first treated in 2004, 
#'   some in 2006, and some in 2007 (see the paper for more details)".
#'   
#'   Key variables:
#'   - lemp: Log of county-level teen employment (outcome variable)
#'   - first.treat: Year of first minimum wage increase (2004, 2006, or 2007)
#'   - year: Year (time variable)
#'   - countyreal: County ID (individual identifier)
#'   - treat: Treatment status (modified in this script)
#'   - lpop : Log population size
#'   - treat_nabs: Non-absorbing treatment dummy (created in this script)
#'   - lpop2: Squared log population (created in this script)
#'   - drop: Random drop indicator (created in this script)"

# Load required libraries ------------------------------------------------------
library(data.table)
library(did)

# Set working directory to script location -----------------------------------
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Import and copy the original dataset -----------------------------------------
data("mpdta", package = "did")
mpdta_r <- setDT(copy(mpdta))  # Use data.table for efficiency

# Modify treatment variables ---------------------------------------------------
# 1. Backup the original treatment variable
mpdta_r[, treat_back := treat]

# 2. Redefine treatment as an absorbing dummy (1 if treated at any point)
mpdta_r[, treat := fifelse(first.treat == 0, 0, (year - first.treat >= 0) * 1)]

# 3. Create a non-absorbing treatment dummy (1 only in the treatment year)
mpdta_r[, treat_nabs := (first.treat == year) * 1]

# Add derived variables ---------------------------------------------------------
# 1. Squared log population (for nonlinear models)
mpdta_r[, lpop2 := lpop^2]

# 2. Random drop indicator (for robustness checks)
set.seed(444)  # For reproducibility
mpdta_r[, noise := runif(nrow(mpdta_r))]
mpdta_r[, drop := 1 - (noise > 0.1 | first.treat > 0)]

# Save the modified dataset ----------------------------------------------------
save(mpdta_r, file = "../data/mpdta_r.rda")  # Note: Saving mpdta_r, not mpdta
