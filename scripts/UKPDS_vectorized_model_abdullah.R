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



## Check results ----
# df_results_abdullah <- abdullah_model[["a_ind_traits"]] |>
#   as.data.frame.table(responseName = "value", stringsAsFactors = F) |>
#   as_tibble() |>
#   pivot_wider(
#     values_from = value,
#     names_from = Var2
#   ) |>
#   mutate(
#     cycle = parse_number(as.character(Var3)),
#     id = parse_number(as.character(Var1)),
#     .keep = "unused",
#     .before = everything()
#   )
#
#
# # Sample some patients
# my_sample <- sample(1:250e3, 1000)
#
# # Check trends ----
# df_results_abdullah |>
#   filter(id %in% my_sample) |>
#   ggplot(
#     aes(x = cycle, y = wbc, color = as.factor(id))
#   ) +
#   geom_line() +
#   theme(
#     legend.position = "none"
#   )

# ## Compare speed performance between the three models ----
# microbenchmark::microbenchmark(
#   abdullah = run_microsim_model(
#     df_UKPDS_coef = df_UKPDS_coef,
#     df_ukpds_pop = df_ukpds_pop,
#     cycles = 20,
#     model = "abdullah"
#   ),
#   # Aaron's model needs to be wrapped
#   aaron = {
#     a_ind_traits <- update_all_biomarkers(a_ind_traits, a_coef_ukpds_ind_traits,
#       time_step = 1, next_row = 2
#     )
#
#     # create a patient population
#     # for (time_step in 2:num_cycles)
#     for (time_step in 2:20) {
#       a_ind_traits[, "age", time_step] <- a_ind_traits[, "age", max(1, time_step - 1)] + 1
#       a_ind_traits[, "diab_dur", time_step] <- a_ind_traits[, "diab_dur", max(1, time_step - 1)] + 1
#       a_ind_traits[, "diab_dur_log", time_step] <- (log(a_ind_traits[, "diab_dur", time_step]))
#
#       # a_ind_traits<- update_all_biomarkers(a_ind_traits, a_coef_ukpds_ind_traits,
#       #                                      time_step = time_step, next_row = time_step+1)
#
#       # ready to simulate
#       # event prediction at t
#       a_ind_traits <- update_health_events(a_ind_traits, a_coef_ukpds_ind_traits,
#         time_step = time_step
#       )
#       # mortality prediction at t
#
#       a_other_ind_traits <- mortality(a_ind_traits, a_other_ind_traits, a_coef_ukpds_ind_traits, time_step = time_step)
#       # predict the risk factors for the next cycle (t+1)
#
#       a_ind_traits <- update_all_biomarkers(a_ind_traits, a_coef_ukpds_ind_traits, time_step = time_step, next_row = time_step + 1)
#     }
#   },
#   wael = {
#     a_ind_traits[, , 2] <- update_all_biomarkers_w(
#       m_ind_traits_init = m_ind_traits_init,
#       m_ind_traits_nStep = a_ind_traits[, , 2],
#       m_ind_traits_step = a_ind_traits[, , 1],
#       a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
#       a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits
#     )
#
#     # create a patient population
#     for (time_step in 2:20) {
#       a_ind_traits[, "age", time_step] <- a_ind_traits[
#         , "age", max(1, time_step - 1)
#       ] + 1
#       a_ind_traits[, "diab_dur", time_step] <- a_ind_traits[
#         , "diab_dur", max(1, time_step - 1)
#       ] + 1
#       a_ind_traits[, "diab_dur_log", time_step] <- log(
#         a_ind_traits[, "diab_dur", time_step]
#       )
#
#       # ready to simulate
#
#       # event prediction at t
#       a_ind_traits[, , time_step] <- update_health_events_w(
#         m_ind_traits_init = m_ind_traits_init,
#         m_ind_traits_step = a_ind_traits[, , time_step],
#         m_ind_traits_pStep = a_ind_traits[, , time_step - 1],
#         a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
#         a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits
#       )
#
#       # mortality prediction at t  time_step <- 1
#       a_other_ind_traits[, "death", time_step] <- mortality_w(
#         m_ind_traits_step = a_ind_traits[, , time_step],
#         m_ind_traits_pStep = a_ind_traits[, , max(1, time_step - 1)],
#         m_other_ind_traits_step = a_other_ind_traits[, , time_step],
#         m_other_ind_traits_pStep = a_other_ind_traits[, , max(1, time_step - 1)],
#         a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
#         a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits
#       )
#
#       # predict the risk factors for the next cycle (t+1)
#       a_ind_traits[, , time_step + 1] <- update_all_biomarkers_w(
#         m_ind_traits_init = m_ind_traits_init,
#         m_ind_traits_nStep = a_ind_traits[, , time_step + 1],
#         m_ind_traits_step = a_ind_traits[, , time_step],
#         a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
#         a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits
#       )
#
#       # # provide some feedback to the console
#       # if(time_step/2 == round(time_step/2,0)) {
#       #   cat('\r', paste(time_step/num_cycles * 100, "% done", sep = " "))
#       # }
#     }
#   },
#   times = 10
# )
