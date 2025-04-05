#' Title
#'
#' @param m_ind_traits
#' @param health_outcome
#' @param m_coef_ukpds_ind_traits
#' @param m_coef_ukpds_other_ind_traits
#'
#' @returns
#' @export
#'
#' @examples
weibull_event2 <- function(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome) {
  # Create a vector to map health_outcome to health_event since they are not the same
  v_health_outcome_event <- c(
    amp1_no_ulcer = "amp_event",
    amp1_yes_ulcer = "amp_event",
    amp2 = "amp_event2",
    mi1_male = "mi_event",
    mi1_female = "mi_event",
    mi2 = "mi_event",
    stroke_1 = "stroke_event",
    stroke_2 = "stroke_event",
    ulcer = "ulcer_event",
    blindness = "blindness_event",
    chf = "chf_event",
    esrd = "esrd_event",
    ihd = "ihd_event",
    mi = "mi_event"
  )

  health_event <- v_health_outcome_event[health_outcome]

  # Calculate patient-specific factors using model coefficients and patient data
  patient_factors <- (m_ind_traits %*%
    m_coef_ukpds_ind_traits[, health_outcome] +
    as.vector(m_coef_ukpds_other_ind_traits["lambda", health_outcome]))

  # Compute cumulative hazard at the current time step
  cum_hazard_t <- exp(patient_factors) *
    (m_ind_traits[, "diab_dur"]^(m_coef_ukpds_other_ind_traits[
      "rho",
      health_outcome
    ]))

  # Compute cumulative hazard at the next time step (by adding 1 year to diabetes duration)
  cum_hazard_t1 <- exp(patient_factors) *
    ((m_ind_traits[, "diab_dur"] + 1)^(m_coef_ukpds_other_ind_traits[
      "rho",
      health_outcome
    ]))

  # Calculate transition probability
  trans_prob <- 1 - exp(cum_hazard_t - cum_hazard_t1)

  # Simulate whether the event occurs by comparing with a random uniform value
  event <- (trans_prob > runif(nrow(m_ind_traits))) * 1 # to convert logical to numeric
  colnames(event) <- health_event
  # Return the updated matrix
  return(
    event
  )
}
