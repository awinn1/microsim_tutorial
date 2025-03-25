#' Predict biomarker values of the next cycle
#'
#' This is my modification of Aaron's functions `biomarker()` and `update_all_biomarkers()`
#'
#' @param m_ind_traits
#' @param m_coef_ukpds_ind_traits
#' @param v_biomarkers
#' @param m_coef_ukpds_other_ind_traits
#'
#' @returns
#' @export
#'
#' @examples
update_biomarkers2 <- \(
  m_ind_traits,
  m_coef_ukpds_ind_traits,
  m_coef_ukpds_other_ind_traits,
  v_biomarkers = c(
    # name = equation
    a1c = "hba1c",
    sbp_real = "sbp",
    ldl_real = "ldl",
    hdl_real = "hdl",
    bmi = "bmi",
    heart_rate_real = "heart_rate",
    wbc = "wbc",
    heamo = "haem"
  )
) {
  # Extract lagged values
  m_ind_traits[, c(
    "a1c_lag",
    "bmi_lag",
    "hdl_lag",
    "heart_rate_lag",
    "ldl_lag",
    "sbp_lag",
    "wbc_lag"
  )] <-
    m_ind_traits[, c(
      "a1c",
      "bmi",
      "hdl_real",
      "heart_rate_real",
      "ldl_real",
      "sbp_real",
      "wbc"
    )]

  step1 <- m_ind_traits %*% m_coef_ukpds_ind_traits[, v_biomarkers]
  step2 <- m_coef_ukpds_other_ind_traits["lambda", v_biomarkers]

  m_ind_traits[, names(v_biomarkers)] <- step1 + step2[col(step1)]

  # biomarkers transformations
  # divide by 10 transformation
  m_ind_traits[, c("sbp", "heart_rate")] <- m_ind_traits[, c(
    "sbp_real",
    "heart_rate_real"
  )] /
    10

  # multiply by 10 transformation
  m_ind_traits[, c("ldl", "hdl")] <- m_ind_traits[, c(
    "ldl_real",
    "hdl_real"
  )] *
    10

  # other transformation
  m_ind_traits[, "bmi_lt18_5"] <- as.integer(m_ind_traits[, "bmi"] < 18.5)
  m_ind_traits[, "bmi_gte25"] <- as.integer(m_ind_traits[, "bmi"] >= 25)
  m_ind_traits[, "ldl_gt35"] <- as.integer(m_ind_traits[, "ldl_real"] > 35) / 10

  return(m_ind_traits)
}
