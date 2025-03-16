# Weibull distribution function for event occurrence
# Logistic regression for binary event prediction


#' Calculate Transition Probability Based on a Weibull Model and Update Patient State
#' Note: An exponential model is a special case of the Weibull model where the shape
#' parameter (ρ) is set to 1, meaning the hazard function remains constant over time,
#' resulting in a constant rate of event occurrence rather than a time-dependent rate.
#'
#' This function calculates patient-specific factors, cumulative hazards,
#' and the transition probability for a given health outcome (e.g., "ihd").
#' The function updates the provided `m_ind_traits` matrix with the event occurrence
#' at the specified time step.
#'
#' @param m_ind_traits A matrix containing patient characteristics over time.
#' @param a_coef_ukpds_ind_traits A 3D array of coefficients used for calculating risk.
#' @param health_outcome A character string specifying the health outcome equation (e.g., "ihd").
#' @param health_event A character string specifying the health outcome event in the patient trace.
#' @param time_step An integer indicating the row in `m_ind_traits` to use for calculations.
#'
#' @return Whether the event occurred.
#' @export
weibull_event <- function(a_ind_traits, a_coef_ukpds_ind_traits, health_outcome, health_event, time_step) {
  # Calculate patient-specific factors using model coefficients and patient data
  patient_factors <- (a_ind_traits[, , time_step] %*% a_coef_ukpds_ind_traits[, health_outcome, 1] +
    as.vector(a_coef_ukpds_other_ind_traits["lambda", health_outcome, 1]))

  # Compute cumulative hazard at the current time step
  cum_hazard_t <- exp(patient_factors) * (a_ind_traits[, "diab_dur", time_step]^(a_coef_ukpds_other_ind_traits["rho", health_outcome, 1]))

  # Compute cumulative hazard at the next time step (by adding 1 year to diabetes duration)
  cum_hazard_t1 <- exp(patient_factors) * ((a_ind_traits[, "diab_dur", time_step] + 1)^(a_coef_ukpds_other_ind_traits["rho", health_outcome, 1]))

  # Calculate transition probability
  trans_prob <- 1 - exp(cum_hazard_t - cum_hazard_t1)

  # Simulate whether the event occurs by comparing with a random uniform value
  event <- trans_prob > runif(1)

  # Return the updated matrix
  return(event)
}
