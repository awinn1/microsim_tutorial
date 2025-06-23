#' Title
#'
#' @param m_ind_traits
#' @param m_coef_ukpds_other_ind_traits
#' @param m_coef_ukpds_ind_traits
#'
#' @returns
#' @export
#'
#' @examples
update_health_events2 <- function(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits) {
  # Ensure m_ind_traits remains a matrix
  if (!is.matrix(m_ind_traits)) {
    stop("m_ind_traits must be a matrix")
  }

  # Initialize event variables and update history
  events <- c("amp", "blindness", "chf", "esrd", "ihd", "mi", "stroke", "ulcer")

  # create event and history column names once and save each group in a vector
  v_event_cols <- paste0(events, "_event")
  v_history_cols <- paste0(events, "_hist")

  # Update history columns in one vectorized call
  m_ind_traits[, v_history_cols] <- pmax(
    m_ind_traits[, v_history_cols],
    m_ind_traits[, v_event_cols]
  )

  # Randomize event order
  randomized_events <- sample(events)

  for (event in randomized_events) {
    if (event == "amp") {
      amp1_no_ulcer <- weibull_event2(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        health_outcome = "amp1_no_ulcer"
      )
      amp1_yes_ulcer <- weibull_event2(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        health_outcome = "amp1_yes_ulcer"
      )

      m_ind_traits[, "amp_event"] <-
        (amp1_no_ulcer * (m_ind_traits[, "ulcer_hist"] == 0)) +
        (amp1_yes_ulcer * (m_ind_traits[, "ulcer_hist"] == 1))

      # Ensure that this is a new event
      m_ind_traits[, "amp_event"] <-
        m_ind_traits[, "amp_event"] *
          (m_ind_traits[, "amp_hist"] == 0)

      # Calculate amp2 event
      amp2 <- weibull_event2(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        health_outcome = "amp2"
      )

      m_ind_traits[, "amp_event2"] <-
        amp2 * (m_ind_traits[, "amp_hist"] == 1)
    } else if (event == "mi") {
      mi1_male <- weibull_event2(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        health_outcome = "mi1_male"
      )
      mi1_female <- weibull_event2(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        health_outcome = "mi1_female"
      )

      m_ind_traits[, "mi_event"] <-
        (mi1_male * (m_ind_traits[, "female"] == 0)) +
        (mi1_female * (m_ind_traits[, "female"] == 1))

      m_ind_traits[, "mi_event"] <-
        m_ind_traits[, "mi_event"] *
          (m_ind_traits[, "mi_hist"] == 0)

      mi2 <- weibull_event2(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        health_outcome = "mi2"
      )

      m_ind_traits[, "mi_event"] <-
        (m_ind_traits[, "mi_hist"] == 0) *
        m_ind_traits[, "mi_event"] +
        (m_ind_traits[, "mi_hist"] == 1) * mi2
    } else if (event == "stroke") {
      stroke1 <- weibull_event2(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        health_outcome = "stroke_1"
      )
      stroke2 <- weibull_event2(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        health_outcome = "stroke_2"
      )

      m_ind_traits[, "stroke_event"] <-
        (stroke1 * (m_ind_traits[, "stroke_hist"] == 0)) +
        (stroke2 * (m_ind_traits[, "stroke_hist"] == 1))

      m_ind_traits[, "stroke_event"] <-
        m_ind_traits[, "stroke_event"] *
          (m_ind_traits[, "stroke_hist"] == 0)
    } else if (event == "ulcer") {
      m_ind_traits[, "ulcer_event"] <-
        logistic_event2(
          m_ind_traits,
          m_coef_ukpds_ind_traits,
          m_coef_ukpds_other_ind_traits,
          health_outcome = "ulcer"
        )

      m_ind_traits[, "ulcer_event"] <-
        m_ind_traits[, "ulcer_event"] *
          (m_ind_traits[, "ulcer_hist"] == 0)
    } else {
      m_ind_traits[, paste0(event, "_event")] <-
        weibull_event2(
          m_ind_traits,
          m_coef_ukpds_ind_traits,
          m_coef_ukpds_other_ind_traits,
          health_outcome = event
        )

      m_ind_traits[, paste0(events, "_event")] <-
        m_ind_traits[, paste0(events, "_event")] *
          (m_ind_traits[, paste0(events, "_hist")] == 0)
    }
  }

  return(m_ind_traits)
}
