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
update_health_eventsC_ <- function(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits) {
  # Ensure m_ind_traits remains a matrix
  if (!is.matrix(m_ind_traits)) {
    stop("m_ind_traits must be a matrix")
  }

  # Define events
  events <- c("amp", "blindness", "chf", "esrd", "ihd", "mi", "stroke", "ulcer")

  # Precompute event and history column names once
  v_event_cols <- paste0(events, "_event")
  v_history_cols <- paste0(events, "_hist")

  # Update history columns in one vectorized call
  m_ind_traits[, v_history_cols] <- pmax(
    m_ind_traits[, v_history_cols],
    m_ind_traits[, v_event_cols]
  )

  # Create a lookup vector where names are the coefficient names.
  coef_names <- colnames(m_coef_ukpds_ind_traits)
  # This creates a named integer vector: e.g. "amp1_no_ulcer" = 1, etc.
  coef_index <- setNames(match(coef_names, coef_names), coef_names)

  # Precompute history columns
  ulcer_hist <- m_ind_traits[, "ulcer_hist"]
  amp_hist <- m_ind_traits[, "amp_hist"]
  mi_hist <- m_ind_traits[, "mi_hist"]
  stroke_hist <- m_ind_traits[, "stroke_hist"]

  # Randomize event order
  randomized_events <- sample(events)

  for (event in randomized_events) {
    if (event == "amp") {
      amp1_no_ulcer <- weibull_eventC(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        coef_index[["amp1_no_ulcer"]]
      )
      amp1_yes_ulcer <- weibull_eventC(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        coef_index[["amp1_yes_ulcer"]]
      )

      m_ind_traits[, "amp_event"] <-
        (amp1_no_ulcer * (ulcer_hist == 0)) +
        (amp1_yes_ulcer * (ulcer_hist == 1))

      m_ind_traits[, "amp_event"] <-
        m_ind_traits[, "amp_event"] * (amp_hist == 0)

      amp2 <- weibull_eventC(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        coef_index[["amp2"]]
      )
      m_ind_traits[, "amp_event2"] <-
        amp2 * (amp_hist == 1)
    } else if (event == "mi") {
      mi1_male <- weibull_eventC(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        coef_index[["mi1_male"]]
      )
      mi1_female <- weibull_eventC(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        coef_index[["mi1_female"]]
      )

      m_ind_traits[, "mi_event"] <-
        (mi1_male * (m_ind_traits[, "female"] == 0)) +
        (mi1_female * (m_ind_traits[, "female"] == 1))

      m_ind_traits[, "mi_event"] <-
        m_ind_traits[, "mi_event"] * (mi_hist == 0)

      mi2 <- weibull_eventC(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        coef_index[["mi2"]]
      )
      m_ind_traits[, "mi_event"] <-
        (mi_hist == 0) *
        m_ind_traits[, "mi_event"] +
        (mi_hist == 1) * mi2
    } else if (event == "stroke") {
      stroke1 <- weibull_eventC(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        coef_index[["stroke_1"]]
      )
      stroke2 <- weibull_eventC(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        coef_index[["stroke_2"]]
      )

      m_ind_traits[, "stroke_event"] <-
        (stroke1 * (stroke_hist == 0)) +
        (stroke2 * (stroke_hist == 1))

      m_ind_traits[, "stroke_event"] <-
        m_ind_traits[, "stroke_event"] * (stroke_hist == 0)
    } else if (event == "ulcer") {
      m_ind_traits[, "ulcer_event"] <-
        logistic_eventC(
          m_ind_traits,
          m_coef_ukpds_ind_traits,
          m_coef_ukpds_other_ind_traits,
          coef_index[["ulcer"]]
        )

      m_ind_traits[, "ulcer_event"] <-
        m_ind_traits[, "ulcer_event"] * (ulcer_hist == 0)
    } else {
      m_ind_traits[, paste0(event, "_event")] <-
        weibull_eventC(
          m_ind_traits,
          m_coef_ukpds_ind_traits,
          m_coef_ukpds_other_ind_traits,
          coef_index[[event]]
        )

      m_ind_traits[, paste0(event, "_event")] <-
        m_ind_traits[, paste0(event, "_event")] *
          (m_ind_traits[, paste0(event, "_hist")] == 0)
    }
  }

  return(m_ind_traits)
}
