#' Title
#'
#' @param df_UKPDS_coef
#' @param df_ukpds_pop
#' @param discount_rate
#' @param q_baseline
#' @param q_blindness
#' @param q_amp
#' @param q_chf
#' @param q_esrd
#' @param q_ihd
#' @param q_mi
#' @param q_stroke
#' @param q_ulcer
#' @param c_baseline
#' @param c_blindness_e
#' @param c_blindness_c
#' @param c_amp_e
#' @param c_amp_c
#' @param c_chf_e
#' @param c_chf_c
#' @param c_esrd_e
#' @param c_esrd_c
#' @param c_ihd_e
#' @param c_ihd_c
#' @param c_mi_e
#' @param c_mi_c
#' @param c_stroke_e
#' @param c_stroke_c
#' @param c_ulcer_e
#' @param c_ulcer_c
#' @param model
#' @param cycles
#'
#' @returns
#' @export
#'
#' @examples
run_microsim_modelC <- \(
  df_UKPDS_coef,
  df_ukpds_pop,
  # a_ind_traits,
  # a_coef_ukpds_ind_traits,
  # a_other_ind_traits,
  # a_coef_ukpds_other_ind_traits,
  discount_rate = 0.03,
  q_baseline = 0.785,
  q_blindness = -0.074,
  q_amp = -0.280,
  q_chf = -0.108,
  q_esrd = -0.204,
  q_ihd = -0.090,
  q_mi = -0.055,
  q_stroke = -0.164,
  q_ulcer = -0.170,
  c_baseline = 1990,
  c_blindness_e = 4247,
  c_blindness_c = 2206,
  c_amp_e = 15153,
  c_amp_c = 5328,
  c_chf_e = 5650,
  c_chf_c = 4277,
  c_esrd_e = 43359,
  c_esrd_c = 43359,
  c_ihd_e = 14001,
  c_ihd_c = 3550,
  c_mi_e = 9518,
  c_mi_c = 3424,
  c_stroke_e = 10755,
  c_stroke_c = 3534,
  c_ulcer_e = 7076,
  c_ulcer_c = 1072,
  model = "abdullah",
  cycles = 20
) {
  # Replace NAs with 0s to avoid missing values in calculations
  df_UKPDS_coef[is.na(df_UKPDS_coef)] <- 0
  
  # Extract parameter names (used as row names)
  v_coef_names <- df_UKPDS_coef$Parameter # Get row names from the 'Parameter' column
  
  # Determine the number of parameters (rows)
  n_coef_names <- length(v_coef_names) # Count the number of parameters
  
  # Extract factor names (used as column names), excluding the first column
  v_factors_names <- as.vector(colnames(df_UKPDS_coef[-1])) # Get column names excluding 'Parameter'
  
  # Determine the number of factors (columns)
  n_equa_names <- length(v_factors_names) # Count the number of factors
  
  # allow for bootstrapped coefficients
  boot <- 1 # Abdullah: Only one value for now?
  rep_names <- paste0("boot_rep_", 1:boot)
  
  # create an array that holds onto everything!
  a_coef_ukpds <- array(
    data = NA,
    dim = c(n_coef_names, n_equa_names, boot),
    dimnames = list(v_coef_names, v_factors_names, rep_names)
  )
  
  
  
  
  # fill in the array with coefficents from the dataset (The first depth layer is exactly the same as your UKPDS_coef. We can just import it directly)
  # Use UKPDS_coef as our first depth layer in the 3D array
  a_coef_ukpds[, , 1] <- as.matrix(df_UKPDS_coef[, -1])
  
  
  # Remove lambda, rho, and death and create an array for individual traits and another array for those three parameters
  # This separation is based on Rob's suggestion
  a_coef_ukpds_ind_traits <- a_coef_ukpds[1:62, , "boot_rep_1", drop = FALSE]
  a_coef_ukpds_other_ind_traits <- a_coef_ukpds[63:65, , , drop = FALSE]
  
  
  num_i <- nrow(df_ukpds_pop) # number of simulated individuals
  # Define the number of time points
  
  
  
  v_ids <- paste("id", 1:num_i, sep = "_")
  
  v_cycles <- paste("cycle", 0:cycles, sep = "_")
  
  nbr_coef_names <- n_coef_names
  
  # Create an array with columns for each variable, row for each person, and a
  # slice for each period
  
  a_all_ind_traits <- array(
    data = NA,
    dim = c(num_i, nbr_coef_names, cycles + 1),
    dimnames = list(v_ids, v_coef_names, v_cycles)
  )
  
  
  # need this to be the same number of columns as the coefficient table is long/rows
  # print(dim(a_all_ind_traits)) # to verify the dimensions
  # print(dimnames(a_all_ind_traits)) # to verify the dimension names
  
  # Fill patient trace
  # Create a vector for time invariant characteristics
  v_time_invariant_characteristics <- c(
    "age_diag",
    "black",
    "indian",
    "female",
    "smoke",
    "a1c_first",
    "bmi_first",
    "heart_rate_first",
    "ldl_first",
    "hdl_first",
    "sbp_first",
    "wbc_first",
    "heamo_first"
  )
  
  a_all_ind_traits[
    ,
    v_time_invariant_characteristics,
  ] <- as.matrix(df_ukpds_pop[, v_time_invariant_characteristics])
  
  # break array up into stuff individual traits
  a_ind_traits <- a_all_ind_traits[, 1:62, ]
  # remaining array that captures lambda, rho and death
  a_other_ind_traits <- a_all_ind_traits[, 63:65, ]
  
  # The Aaron's version of the model has an issue but I still don't know exactly what is it?
  if (model == "abdullah") {
    for (time_step in 2:(cycles + 1)) {
      a_ind_traits[, , time_step] <- update_biomarkers2(
        m_ind_traits = a_ind_traits[, , time_step - 1],
        m_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits[, , 1],
        m_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits[, , 1]
      )
      
      a_ind_traits[, "age", time_step] <- a_ind_traits[, "age", time_step - 1] + 1
      a_ind_traits[, "diab_dur", time_step] <- a_ind_traits[, "diab_dur", time_step - 1] + 1
      a_ind_traits[, "diab_dur_log", time_step] <- log(a_ind_traits[, "diab_dur", time_step])
      
      
      a_ind_traits[, , time_step] <- update_health_eventsC(
        a_ind_traits[, , time_step],
        m_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits[, , 1],
        m_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits[, , 1]
      )
      
      a_other_ind_traits[, , time_step] <- mortalityC(
        m_ind_traits = a_ind_traits[, , time_step],
        m_other_ind_traits = a_other_ind_traits[, , time_step],
        m_other_ind_traits_previous = a_other_ind_traits[, , time_step - 1],
        m_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits[, , 1],
        m_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits[, , 1]
      )
    }
  }
  return(
    list(
      a_ind_traits = a_ind_traits,
      a_other_ind_traits = a_other_ind_traits
    )
  )
}

#' Title
#'
#' @param df_UKPDS_coef
#' @param df_ukpds_pop
#' @param discount_rate
#' @param q_baseline
#' @param q_blindness
#' @param q_amp
#' @param q_chf
#' @param q_esrd
#' @param q_ihd
#' @param q_mi
#' @param q_stroke
#' @param q_ulcer
#' @param c_baseline
#' @param c_blindness_e
#' @param c_blindness_c
#' @param c_amp_e
#' @param c_amp_c
#' @param c_chf_e
#' @param c_chf_c
#' @param c_esrd_e
#' @param c_esrd_c
#' @param c_ihd_e
#' @param c_ihd_c
#' @param c_mi_e
#' @param c_mi_c
#' @param c_stroke_e
#' @param c_stroke_c
#' @param c_ulcer_e
#' @param c_ulcer_c
#' @param model
#' @param cycles
#'
#' @returns
#' @export
#'
#' @examples
run_microsim_modelC_ <- \(
  df_UKPDS_coef,
  df_ukpds_pop,
  # a_ind_traits,
  # a_coef_ukpds_ind_traits,
  # a_other_ind_traits,
  # a_coef_ukpds_other_ind_traits,
  discount_rate = 0.03,
  q_baseline = 0.785,
  q_blindness = -0.074,
  q_amp = -0.280,
  q_chf = -0.108,
  q_esrd = -0.204,
  q_ihd = -0.090,
  q_mi = -0.055,
  q_stroke = -0.164,
  q_ulcer = -0.170,
  c_baseline = 1990,
  c_blindness_e = 4247,
  c_blindness_c = 2206,
  c_amp_e = 15153,
  c_amp_c = 5328,
  c_chf_e = 5650,
  c_chf_c = 4277,
  c_esrd_e = 43359,
  c_esrd_c = 43359,
  c_ihd_e = 14001,
  c_ihd_c = 3550,
  c_mi_e = 9518,
  c_mi_c = 3424,
  c_stroke_e = 10755,
  c_stroke_c = 3534,
  c_ulcer_e = 7076,
  c_ulcer_c = 1072,
  model = "abdullah",
  cycles = 20
) {
  # Replace NAs with 0s to avoid missing values in calculations
  df_UKPDS_coef[is.na(df_UKPDS_coef)] <- 0
  
  # Extract parameter names (used as row names)
  v_coef_names <- df_UKPDS_coef$Parameter # Get row names from the 'Parameter' column
  
  # Determine the number of parameters (rows)
  n_coef_names <- length(v_coef_names) # Count the number of parameters
  
  # Extract factor names (used as column names), excluding the first column
  v_factors_names <- as.vector(colnames(df_UKPDS_coef[-1])) # Get column names excluding 'Parameter'
  
  # Determine the number of factors (columns)
  n_equa_names <- length(v_factors_names) # Count the number of factors
  
  # allow for bootstrapped coefficients
  boot <- 1 # Abdullah: Only one value for now?
  rep_names <- paste0("boot_rep_", 1:boot)
  
  # create an array that holds onto everything!
  a_coef_ukpds <- array(
    data = NA,
    dim = c(n_coef_names, n_equa_names, boot),
    dimnames = list(v_coef_names, v_factors_names, rep_names)
  )
  
  
  
  
  # fill in the array with coefficents from the dataset (The first depth layer is exactly the same as your UKPDS_coef. We can just import it directly)
  # Use UKPDS_coef as our first depth layer in the 3D array
  a_coef_ukpds[, , 1] <- as.matrix(df_UKPDS_coef[, -1])
  
  
  # Remove lambda, rho, and death and create an array for individual traits and another array for those three parameters
  # This separation is based on Rob's suggestion
  a_coef_ukpds_ind_traits <- a_coef_ukpds[1:62, , "boot_rep_1", drop = FALSE]
  a_coef_ukpds_other_ind_traits <- a_coef_ukpds[63:65, , , drop = FALSE]
  
  
  num_i <- nrow(df_ukpds_pop) # number of simulated individuals
  # Define the number of time points
  
  
  
  v_ids <- paste("id", 1:num_i, sep = "_")
  
  v_cycles <- paste("cycle", 0:cycles, sep = "_")
  
  nbr_coef_names <- n_coef_names
  
  # Create an array with columns for each variable, row for each person, and a
  # slice for each period
  
  a_all_ind_traits <- array(
    data = NA,
    dim = c(num_i, nbr_coef_names, cycles + 1),
    dimnames = list(v_ids, v_coef_names, v_cycles)
  )
  
  
  # need this to be the same number of columns as the coefficient table is long/rows
  # print(dim(a_all_ind_traits)) # to verify the dimensions
  # print(dimnames(a_all_ind_traits)) # to verify the dimension names
  
  # Fill patient trace
  # Create a vector for time invariant characteristics
  v_time_invariant_characteristics <- c(
    "age_diag",
    "black",
    "indian",
    "female",
    "smoke",
    "a1c_first",
    "bmi_first",
    "heart_rate_first",
    "ldl_first",
    "hdl_first",
    "sbp_first",
    "wbc_first",
    "heamo_first"
  )
  
  
  a_all_ind_traits[
    ,
    v_time_invariant_characteristics,
  ] <- as.matrix(df_ukpds_pop[, v_time_invariant_characteristics])
  
  # Fill baseline values of the first cycle
  a_all_ind_traits[, , 1] <- as.matrix(df_ukpds_pop[, -1])
  
  
  # break array up into stuff individual traits
  a_ind_traits <- a_all_ind_traits[, 1:62, ]
  # remaining array that captures lambda, rho and death
  a_other_ind_traits <- a_all_ind_traits[, 63:65, ]
  
  
  # The Aaron's version of the model has an issue but I still don't know exactly what is it?
  if (model == "abdullah") {
    for (time_step in 2:(cycles + 1)) {
      a_ind_traits[, , time_step] <- update_biomarkers2(
        m_ind_traits = a_ind_traits[, , time_step - 1],
        m_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits[, , 1],
        m_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits[, , 1]
      )
      
      a_ind_traits[, "age", time_step] <- a_ind_traits[, "age", time_step - 1] + 1
      a_ind_traits[, "diab_dur", time_step] <- a_ind_traits[, "diab_dur", time_step - 1] + 1
      a_ind_traits[, "diab_dur_log", time_step] <- log(a_ind_traits[, "diab_dur", time_step])
      
      
      a_ind_traits[, , time_step] <- update_health_eventsC_(
        a_ind_traits[, , time_step],
        m_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits[, , 1],
        m_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits[, , 1]
      )
      
      a_other_ind_traits[, , time_step] <- mortalityC(
        m_ind_traits = a_ind_traits[, , time_step],
        m_other_ind_traits = a_other_ind_traits[, , time_step],
        m_other_ind_traits_previous = a_other_ind_traits[, , time_step - 1],
        m_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits[, , 1],
        m_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits[, , 1]
      )
    }
  }
  return(
    list(
      a_ind_traits = a_ind_traits,
      a_other_ind_traits = a_other_ind_traits
    )
  )
}
