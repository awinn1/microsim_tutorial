#' Calculate Transition Probability Based on a Logistic Regression and Update Patient State
#'
#' This function calculates patient-specific factors, cumulative hazards,
#' and the transition probability for a given health outcome (e.g., "ihd").
#' The function updates the provided `m_ind_traits` matrix with the event occurrence
#' at the specified time step.
#'
#' @param a_ind_traits A 3d array containing patient characteristics over time (slices).
#' @param a_coef_ukpds_ind_traits A 3D array of coefficients used for calculating risk.
#' @param health_outcome A character string specifying the health outcome equation (e.g., "ihd").
#' @param health_event A character string specifying the health outcome event in the patient trace.
#' @param time_step An integer indicating the row in `m_ind_traits` to use for calculations.
#'
#' @return Whether the event occurred.
#' @export
logistic_event <- function(a_ind_traits, a_coef_ukpds_ind_traits, health_outcome, health_event, time_step) {
  # Calculate patient-specific factors using model coefficients and patient data
  patient_factors <- (a_ind_traits[, , time_step] %*% a_coef_ukpds_ind_traits[, health_outcome, 1] +
    as.vector(a_coef_ukpds_other_ind_traits["lambda", health_outcome, 1]))

  # Calculate transition probability
  trans_prob <- 1 - (exp(-patient_factors) / (1 + exp(-patient_factors)))^1

  # Simulate whether the event occurs by comparing with a random uniform value
  event <- trans_prob > runif(1)


  # Return the value
  return(event)
}
