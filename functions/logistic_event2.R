#' Title
#'
#' @param health_outcome
#' @param m_ind_traits
#' @param m_coef_ukpds_ind_traits
#' @param m_coef_ukpds_other_ind_traits
#'
#' @returns
#' @export
#'
#' @examples
logistic_event2 <- function(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome) {
  # Calculate patient-specific factors using model coefficients and patient data
  patient_factors <- (m_ind_traits %*%
    m_coef_ukpds_ind_traits[, health_outcome] +
    as.vector(m_coef_ukpds_other_ind_traits["lambda", health_outcome]))

  # Calculate transition probability
  trans_prob <- 1 - (exp(-patient_factors) / (1 + exp(-patient_factors)))^1

  # Simulate whether the event occurs by comparing with a random uniform value
  event <- (trans_prob > runif(nrow(m_ind_traits))) * 1 # to convert logical to numeric

  # Return the value
  return(event)
}
