# This is Abdullah's version of the model
#
library(pacman)
pacman::p_load(
  haven,
  Rcpp,
  RcppArmadillo,
  profvis,
  readr,
  tidyverse,
  here
)


# About
# This script is my re-writing of the model.


# Call functions
invisible(lapply(fs::dir_ls("cpp_functions"), Rcpp::sourceCpp))
invisible(lapply(fs::dir_ls("functions"), source))

# Call datasets
## Coefficients dataset
df_UKPDS_coef <- read.csv("data/ukpds_coef.csv")

## Population dataset
df_ukpds_pop <- read.csv("data/population.csv")

## PSA bootstraps
df_ukpds_boots <- readRDS("data/ukpds_bootstraps.rds")

## Run model ----
# Initialize event variables and update history
events <- c("amp", "blindness", "chf", "esrd", "ihd", "mi", "stroke", "ulcer")

# create event and history column names once and save each group in a vector
v_event_cols <- paste0(events, "_event")
v_history_cols <- paste0(events, "_hist")

## Check if the output is identical --------------------------------------------------------------------------------
set.seed(123)
l_PSA_results <- run_psa(
    df_UKPDS_coef = df_UKPDS_coef,
    df_ukpds_pop = df_ukpds_pop,
    df_bootstraps = df_ukpds_boots[1:2, ],
    cycles = 20,
    parallel = TRUE
)
set.seed(123)
l_PSA_results <- run_psa(
  df_UKPDS_coef = df_UKPDS_coef,
  df_ukpds_pop = df_ukpds_pop,
  df_bootstraps = df_ukpds_boots[1:2, ],
  cycles = 20,
  parallel = FALSE
)
