#' Calculate Biomarkers
#'
#' This function calculates patient-specific factors to predict the time path of a biomarker.
#'
#' @param a_ind_traits An array containing patient characteristics over time.
#' @param a_coef_ukpds_ind_traits A 3D array of coefficients used for calculating risk.
#' @param biomarker_eq A character string specifying the health outcome equation (e.g., "ihd").
#' @param time_step An integer indicating the row in `m_ind_traits` to use for calculations.
#'
#' @return The updated biomarker is stored.
#' @export
biomarker <- function(a_ind_traits, a_coef_ukpds_ind_traits, biomarker_eq, time_step) {
  # Calculate patient-specific factors using model coefficients and patient data
  updated_biomarker <- (a_ind_traits[, , time_step] %*% a_coef_ukpds_ind_traits[, biomarker_eq, 1] +
    a_coef_ukpds_other_ind_traits["lambda", biomarker_eq, 1])

  return(updated_biomarker)
}



(a_ind_traits[, , "cycle_0"] %*% a_coef_ukpds_ind_traits[, "ihd", 1] +
  a_coef_ukpds_other_ind_traits["lambda", "ihd", 1])

a_ind_traits[, , "cycle_0"] %>% head()
(a_ind_traits[, , "cycle_0"] %*% a_coef_ukpds_ind_traits[, "ihd", 1]) %>%
  head() + a_coef_ukpds_other_ind_traits["lambda", "ihd", 1] %>%
  head()
