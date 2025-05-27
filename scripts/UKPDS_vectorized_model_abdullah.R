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

## Run model ----
# Initialize event variables and update history
events <- c("amp", "blindness", "chf", "esrd", "ihd", "mi", "stroke", "ulcer")

# create event and history column names once and save each group in a vector
v_event_cols <- paste0(events, "_event")
v_history_cols <- paste0(events, "_hist")


## Check if the output is identical --------------------------------------------------------------------------------
# We will repeat the seed to test reproducibility
set.seed(123)
model_R <- run_microsim_model(
  df_UKPDS_coef = df_UKPDS_coef,
  df_ukpds_pop = df_ukpds_pop,
  cycles = 20
)


set.seed(123)
model_C <- run_microsim_modelC(
  df_UKPDS_coef = df_UKPDS_coef,
  df_ukpds_pop = df_ukpds_pop,
  cycles = 20
)


death_R <- model_R$a_other_ind_traits[, "death", ]
events_R <- model_R$a_ind_traits[, v_event_cols, ]
history_R <- model_R$a_ind_traits[, v_history_cols, ]
death_C <- model_C$a_other_ind_traits[, "death", ]
events_C <- model_C$a_ind_traits[, v_event_cols, ]
history_C <- model_C$a_ind_traits[, v_history_cols, ]


if (as.logical(mean(death_R == death_C)) &
  as.logical(mean(events_R == events_C)) &
  as.logical(mean(history_R == history_C))) {
  print("The output is identical!")
} else {
  print("The output is not identical! Check the model")
}



# Compare performance ---------------------------------------------------------------------------------------------

microbenchmark::microbenchmark(
  R_version = run_microsim_model(
    df_UKPDS_coef = df_UKPDS_coef,
    df_ukpds_pop = df_ukpds_pop,
    cycles = 20
  ),
  C_version = run_microsim_modelC(
    df_UKPDS_coef = df_UKPDS_coef,
    df_ukpds_pop = df_ukpds_pop,
    cycles = 20
  ),
  times = 10
)
