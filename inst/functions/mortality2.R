#' Title
#'
#' @param m_ind_traits
#' @param m_other_ind_traits_previous
#' @param m_coef_ukpds_ind_traits
#' @param m_coef_ukpds_other_ind_traits
#' @param m_other_ind_traits
#'
#' @returns
#' @export
#'
#' @examples
mortality2 <- function(
    m_ind_traits,
    m_other_ind_traits,
    m_other_ind_traits_previous,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits) {
  # Calculate new health event occurrence and prior history

  # Define events of interest
  events <- c("amp", "blindness", "chf", "esrd", "ihd", "mi", "stroke", "ulcer")

  # Create vectors containing events and event-history names:
  v_event_cols <- paste0(events, "_event")
  v_hist_cols <- paste0(events, "_hist")

  # Get the maximum across those columns, for the given time_step
  # Calculate any new health event
  new_event <- max(m_ind_traits[, v_event_cols])
  # Calculate any prior history of health events
  any_history <- max(m_ind_traits[, v_hist_cols])

  # Determine event-history combinations
  nhne <- new_event == 0 & any_history == 0 # No history, no event
  yhne <- new_event == 0 & any_history == 1 # Yes history, no event
  nhye <- new_event == 1 & any_history == 0 # No history, new event
  yhye <- new_event == 1 & any_history == 1 # Yes history, new event

  # Mortality calculations using Gompertz and logistic models
  death_nhne <- gompertz_event2(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome = "death_nhne"
  )

  death_yhne <- gompertz_event2(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome = "death_yhne"
  )

  death_nhye <- logistic_event2(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome = "death_1st_event"
  )

  death_yhye <- logistic_event2(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome = "death_yhye"
  )

  # Calculate new mortality status
  new_death <-
    nhne *
    death_nhne +
    yhne * death_yhne +
    nhye * death_nhye +
    yhye * death_yhye

  # Update the mortality status in the matrix for the given time step
  m_other_ind_traits[, "death"] <- new_death + m_other_ind_traits_previous[, "death"]

  return(m_other_ind_traits)
}
