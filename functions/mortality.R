#' @title Calculate Mortality Events for a Given Time Step
#'
#' @description This function calculates mortality events for a given time step
#' based on new health events and medical history using Gompertz and logistic models.
#'
#' @param a_ind_traits An array containing patient-level data, including health event history.
#' @param a_other_ind_traits An array containing lambda, rhos and death.
#' @param a_coef_ukpds_ind_traits A coefficient matrix used in Gompertz and logistic event calculations.
#' @param time_step An integer specifying the time step at which mortality should be calculated.
#'
#' @return The updated `m_ind_traits` matrix with the mortality status recorded for the specified time step.
#'
#' @examples
#' \dontrun{
#' m_ind_traits <- mortality(m_ind_traits, a_coef_ukpds_ind_traits, time_step = 5)
#' }
#'
#' @export
mortality <- function(a_ind_traits, a_other_ind_traits, a_coef_ukpds_ind_traits, time_step) {
  # Calculate new health event occurrence and prior history

  # Define events of interest
  events <- c("amp", "blindness", "chf", "esrd", "ihd", "mi", "stroke", "ulcer")

  # Create vectors containing events and event-history names:
  v_event_cols <- paste0(events, "_event")
  v_hist_cols <- paste0(events, "_hist")

  # Get the maximum across those columns, for the given time_step
  # Calculate any new health event
  new_event <- max(a_ind_traits[, v_event_cols, time_step])
  # Calculate any prior history of health events
  any_history <- max(a_ind_traits[, v_hist_cols, time_step])


  # Determine event-history combinations
  nhne <- new_event == 0 & any_history == 0 # No history, no event
  yhne <- new_event == 0 & any_history == 1 # Yes history, no event
  nhye <- new_event == 1 & any_history == 0 # No history, new event
  yhye <- new_event == 1 & any_history == 1 # Yes history, new event

  # Mortality calculations using Gompertz and logistic models
  death_nhne <- gompertz_event(a_ind_traits, a_coef_ukpds_ind_traits,
    health_outcome = "death_nhne",
    health_event = "death_nhne", time_step = time_step
  )

  death_yhne <- gompertz_event(a_ind_traits, a_coef_ukpds_ind_traits,
    health_outcome = "death_yhne",
    health_event = "death_yhne", time_step = time_step
  )

  death_nhye <- logistic_event(a_ind_traits, a_coef_ukpds_ind_traits,
    health_outcome = "death_1st_event",
    health_event = "death_nhye", time_step = time_step
  )

  death_yhye <- logistic_event(a_ind_traits, a_coef_ukpds_ind_traits,
    health_outcome = "death_yhye",
    health_event = "death_yhye", time_step = time_step
  )

  # Calculate new mortality status
  new_death <- nhne * death_nhne + yhne * death_yhne + nhye * death_nhye + yhye * death_yhye

  # Update the mortality status in the matrix for the given time step
  a_other_ind_traits[, "death", time_step] <- new_death + a_other_ind_traits[, "death", max(time_step - 1, 1)]

  return(a_other_ind_traits)
}
