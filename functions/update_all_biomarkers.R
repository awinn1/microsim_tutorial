# Combine risk factor functions into a single pipeline
# Update patient data over time

#' Update Multiple Biomarkers in a Transition Matrix
#'
#' This function updates multiple biomarker values in the transition matrix for a given time step.
#'
#' @param a_ind_traits The patient trace, an array containing patient data with biomarker and event columns.
#' @param a_coef_ukpds_ind_traits A coefficient matrix containing biomarker and event equations.
#' @param time_step An integer representing the current time step.
#' @param next_row An integer indicating the row in `m_ind_traits` to update with new biomarker values.
#'
#' @return The updated transition matrix `m_ind_traits` with new biomarker values in the specified row.
#'
#' @examples
#' # Example usage
#' \dontrun{
#' m_ind_traits <- update_all_biomarkers(m_ind_traits, a_coef_ukpds_ind_traits, time_step = 1, next_row = 2)
#' }
#'
#' @export

update_all_biomarkers <- function(a_ind_traits, a_coef_ukpds_ind_traits, time_step, next_row) {
  # predict the next period (and perform transformations as needed)
  # the biomarkers use real values of variables, but the event equations use transformed variables
  a_ind_traits[, "a1c", next_row] <- biomarker(a_ind_traits, a_coef_ukpds_ind_traits, biomarker_eq = "hba1c", time_step = time_step)
  a_ind_traits[, "sbp_real", next_row] <- biomarker(a_ind_traits, a_coef_ukpds_ind_traits, biomarker_eq = "sbp", time_step = time_step)
  a_ind_traits[, "sbp", next_row] <- a_ind_traits[, "sbp_real", next_row] / 10
  a_ind_traits[, "ldl_real", next_row] <- biomarker(a_ind_traits, a_coef_ukpds_ind_traits, biomarker_eq = "ldl", time_step = time_step)
  a_ind_traits[, "ldl", next_row] <- a_ind_traits[, "ldl_real", next_row] * 10
  a_ind_traits[, "hdl_real", next_row] <- biomarker(a_ind_traits, a_coef_ukpds_ind_traits, biomarker_eq = "hdl", time_step = time_step)
  a_ind_traits[, "hdl", next_row] <- a_ind_traits[, "hdl_real", next_row] * 10
  a_ind_traits[, "bmi", next_row] <- biomarker(a_ind_traits, a_coef_ukpds_ind_traits, biomarker_eq = "bmi", time_step = time_step)
  a_ind_traits[, "heart_rate_real", next_row] <- biomarker(a_ind_traits, a_coef_ukpds_ind_traits, biomarker_eq = "heart_rate", time_step = time_step)
  a_ind_traits[, "heart_rate", next_row] <- a_ind_traits[, "heart_rate_real", next_row] / 10
  a_ind_traits[, "wbc", next_row] <- biomarker(a_ind_traits, a_coef_ukpds_ind_traits, biomarker_eq = "wbc", time_step = time_step)
  a_ind_traits[, "heamo", next_row] <- biomarker(a_ind_traits, a_coef_ukpds_ind_traits, biomarker_eq = "haem", time_step = time_step)

  # Update lag and first occurrence columns
  a_ind_traits[, "a1c_lag", next_row] <- a_ind_traits[, "a1c", time_step]
  a_ind_traits[, "bmi_lag", next_row] <- a_ind_traits[, "bmi", time_step]
  a_ind_traits[, "bmi_lt18_5", next_row] <- as.integer(a_ind_traits[, "bmi", next_row] < 18.5)
  a_ind_traits[, "bmi_gte25", next_row] <- as.integer(a_ind_traits[, "bmi", next_row] >= 25)

  a_ind_traits[, "hdl_lag", next_row] <- a_ind_traits[, "hdl_real", time_step]
  a_ind_traits[, "heart_rate_lag", next_row] <- a_ind_traits[, "heart_rate_real", time_step]

  # Check if this is functioning as a spline
  a_ind_traits[, "ldl_gt35", next_row] <- as.integer(a_ind_traits[, "ldl_real", next_row] > 35) / 10
  a_ind_traits[, "ldl_lag", next_row] <- a_ind_traits[, "ldl_real", time_step]
  a_ind_traits[, "sbp_lag", next_row] <- a_ind_traits[, "sbp_real", time_step]
  a_ind_traits[, "wbc_lag", next_row] <- a_ind_traits[, "wbc", time_step]

  # Update additional values
  a_ind_traits[, "egfr", next_row] <- a_ind_traits[, "egfr", 1]
  a_ind_traits[, "egfr_real", next_row] <- a_ind_traits[, "egfr_real", 1]
  a_ind_traits[, "egfr_lt60", next_row] <- a_ind_traits[, "egfr_lt60", 1]
  a_ind_traits[, "egfr_gte60", next_row] <- a_ind_traits[, "egfr_gte60", 1]
  a_ind_traits[, "albumin_mm", next_row] <- a_ind_traits[, "albumin_mm", 1]


  # Return updated matrix
  return(a_ind_traits)
}
