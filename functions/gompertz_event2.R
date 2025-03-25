#' Title
#'
#' @param m_ind_traits
#' @param a_coef_ukpds_ind_traits
#' @param health_outcome
#'
#' @returns
#' @export
#'
#' @examples
gompertz_event2 <- function(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome) {
  # Calculate patient-specific factors using model coefficients and patient data
  patient_factors <- (m_ind_traits %*%
    m_coef_ukpds_ind_traits[, health_outcome] +
    as.vector(m_coef_ukpds_other_ind_traits["lambda", health_outcome]))

  # Compute cumulative hazard at the current time step
  cum_hazard_t <- (1 /
    m_coef_ukpds_other_ind_traits["rho", health_outcome]) *
    exp(patient_factors) *
    (exp(m_ind_traits[, "age"] * (m_coef_ukpds_other_ind_traits["rho", health_outcome])) - 1)

  # Compute cumulative hazard at the next time step (by adding 1 year to diabetes duration)
  cum_hazard_t1 <- (1 /
    m_coef_ukpds_other_ind_traits["rho", health_outcome]) *
    exp(patient_factors) *
    (exp(
      (m_ind_traits[, "age"] + 1) *
        (m_coef_ukpds_other_ind_traits["rho", health_outcome])
    ) -
      1)

  # Calculate transition probability
  trans_prob <- 1 - exp(cum_hazard_t - cum_hazard_t1)

  # Simulate whether the event occurs by comparing with a random uniform value
  event <- (trans_prob > runif(1)) * 1
  colnames(event) <- health_outcome
  # Return the updated matrix
  return(event)
}
