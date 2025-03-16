# Step 6: Initialize event and history variables ####
#' @title Update Health Events Over Time Steps
#' @description This function updates health events in a patient data matrix (`m_ind_traits`) by applying Weibull
#' and logistic event functions in a randomized order across multiple time steps.
#'
#' @param a_ind_traits An array containing patient-level data, including health event history.
#' @param a_coef_ukpds_ind_traits A coefficient matrix used in Weibull and logistic event calculations.
#' @param time_step An integer indicating the current time step to update events.
#'
#' @return Updated `m_ind_traits` matrix with event and history values updated for the given time step.
#'
#' @examples
#' \dontrun{
#' a_ind_traits <- update_health_events(m_ind_traits, a_coef_ukpds_ind_traits, time_step = 1)
#' }
#'
#' @export
update_health_events <- function(a_ind_traits, a_coef_ukpds_ind_traits, time_step) {
  # Ensure m_ind_traits remains a matrix
  if (!is.array(a_ind_traits)) {
    stop("m_ind_traits must be an array.")
  }

  # Initialize event variables and update history
  events <- c("amp", "blindness", "chf", "esrd", "ihd", "mi", "stroke", "ulcer")

  # create event and history column names once and save each group in a vector
  v_event_cols <- paste0(events, "_event")
  v_history_cols <- paste0(events, "_hist")

  # Update history columns in one vectorized call
  a_ind_traits[, v_event_cols, time_step] <- 0
  a_ind_traits[, v_history_cols, time_step] <- pmax(
    a_ind_traits[, v_history_cols, max(1, time_step - 1)],
    a_ind_traits[, v_event_cols, max(1, time_step - 1)]
  )



  a_ind_traits[, "amp_event2", time_step] <- 0
  # Randomize event order
  randomized_events <- sample(events)

  for (events in randomized_events) {
    if (events == "amp") {
      amp1_no_ulcer <- weibull_event(a_ind_traits, a_coef_ukpds_ind_traits,
        health_outcome = "amp1_no_ulcer", health_event = "amp_event", time_step = time_step
      )
      amp1_yes_ulcer <- weibull_event(a_ind_traits, a_coef_ukpds_ind_traits,
        health_outcome = "amp1_yes_ulcer", health_event = "amp_event", time_step = time_step
      )

      a_ind_traits[, "amp_event", time_step] <-
        (amp1_no_ulcer * (a_ind_traits[, "ulcer_hist", time_step] == 0)) +
        (amp1_yes_ulcer * (a_ind_traits[, "ulcer_hist", time_step] == 1))


      # Ensure that this is a new event
      a_ind_traits[, "amp_event", time_step] <-
        a_ind_traits[, "amp_event", time_step] * (a_ind_traits[, "amp_hist", time_step] == 0)

      # Calculate amp2 event
      amp2 <- weibull_event(a_ind_traits, a_coef_ukpds_ind_traits,
        health_outcome = "amp2", health_event = "amp_event2", time_step = time_step
      )

      a_ind_traits[, "amp_event2", time_step] <- 0

      a_ind_traits[, "amp_event2", time_step] <-
        amp2 * (a_ind_traits[, "amp_hist", time_step] == 1)
    } else if (events == "mi") {
      mi1_male <- weibull_event(a_ind_traits, a_coef_ukpds_ind_traits,
        health_outcome = "mi1_male", health_event = "mi_event", time_step = time_step
      )
      mi1_female <- weibull_event(a_ind_traits, a_coef_ukpds_ind_traits,
        health_outcome = "mi1_female", health_event = "mi_event", time_step = time_step
      )

      a_ind_traits[, "mi_event", time_step] <-
        (mi1_male * (a_ind_traits[, "female", time_step] == 0)) +
        (mi1_female * (a_ind_traits[, "female", time_step] == 1))

      a_ind_traits[, "mi_event", time_step] <-
        a_ind_traits[, "mi_event", time_step] * (a_ind_traits[, "mi_hist", time_step] == 0)

      mi2 <- weibull_event(a_ind_traits, a_coef_ukpds_ind_traits,
        health_outcome = "mi2", health_event = "mi_event", time_step = time_step
      )

      a_ind_traits[, "mi_event", time_step] <-
        (a_ind_traits[, "mi_hist", time_step] == 0) * a_ind_traits[, "mi_event", time_step] +
        (a_ind_traits[, "mi_hist", time_step] == 1) * mi2
    } else if (events == "stroke") {
      stroke1 <- weibull_event(a_ind_traits, a_coef_ukpds_ind_traits,
        health_outcome = "stroke_1", health_event = "stroke_event", time_step = time_step
      )
      stroke2 <- weibull_event(a_ind_traits, a_coef_ukpds_ind_traits,
        health_outcome = "stroke_2", health_event = "stroke_event", time_step = time_step
      )

      a_ind_traits[, "stroke_event", time_step] <-
        (stroke1 * (a_ind_traits[, "stroke_hist", time_step] == 0)) +
        (stroke2 * (a_ind_traits[, "stroke_hist", time_step] == 1))

      a_ind_traits[, "stroke_event", time_step] <-
        a_ind_traits[, "stroke_event", time_step] * (a_ind_traits[, "stroke_hist", time_step] == 0)
    } else if (events == "ulcer") {
      a_ind_traits[, "ulcer_event", time_step] <-
        logistic_event(a_ind_traits, a_coef_ukpds_ind_traits,
          health_outcome = "ulcer", health_event = "ulcer_event", time_step = time_step
        )

      a_ind_traits[, "ulcer_event", time_step] <-
        a_ind_traits[, "ulcer_event", time_step] * (a_ind_traits[, "ulcer_hist", time_step] == 0)
    } else {
      a_ind_traits[, paste0(events, "_event"), time_step] <-
        weibull_event(a_ind_traits, a_coef_ukpds_ind_traits,
          health_outcome = events, health_event = paste0(events, "_event"), time_step = time_step
        )

      a_ind_traits[, paste0(events, "_event"), time_step] <-
        a_ind_traits[, paste0(events, "_event"), time_step] *
          (a_ind_traits[, paste0(events, "_hist"), time_step] == 0)
    }
  }

  return(a_ind_traits)
}
