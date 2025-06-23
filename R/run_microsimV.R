#' @title Calculate Transition Probability Based on a Gompertz Model and Update Patient State (vectorised)
#'
#' @description This function calculates patient-specific factors, cumulative hazards,
#' and the transition probability for mortality.
#'
#' @param m_ind_traits A matrix to store patient data.
#' @param m_coef_ukpds_ind_traits A matrix of coefficients for individual‐level risk equations.
#' @param m_coef_ukpds_other_ind_traits A matrix of coefficients `c("lambda", "rho", "death"` used for calculating risk.
#' @param health_outcome A character string specifying the health outcome equation (e.g., "ihd").
#'
#' @return Vector of logical values for whether the event would occur.
#'
#' @export
#'
gompertz_eventV <- function(
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
  event <- (trans_prob > stats::runif(nrow(m_ind_traits))) * 1 # to convert logical to numeric
  colnames(event) <- health_outcome

  # Return the updated matrix
  return(event)
}

#' @title Calculate Transition Probability Based on a Logistic Regression and Update Patient State (vectorised)
#'
#' @description This function calculates patient-specific factors, cumulative hazards,
#' and the transition probability for a given health outcome (e.g., "ihd").
#' The function updates the provided `m_ind_traits` matrix with the event occurrence
#' at the specified time step.
#'
#' @param m_ind_traits A matrix to store patient data.
#' @param m_coef_ukpds_ind_traits A matrix of coefficients for individual‐level risk equations.
#' @param m_coef_ukpds_other_ind_traits A matrix of coefficients `c("lambda", "rho", "death"` used for calculating risk.
#' @param health_outcome A character string specifying the health outcome equation (e.g., "ihd").
#'
#' @return Vector of logical values for whether the event would occur.
#'
#' @export
#'
logistic_eventV <- function(
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
  event <- (trans_prob > stats::runif(nrow(m_ind_traits))) * 1 # to convert logical to numeric

  # Return the value
  return(event)
}

#' @title Calculate Mortality Events for a Given Time Step (vectorised)
#'
#' @description This function calculates mortality events for a given time step
#' based on new health events and medical history using Gompertz and logistic models.
#'
#' @param m_ind_traits A matrix containing patient-level data, including health event history.
#' @param m_other_ind_traits A matrix containing lmabda, rhos and death.
#' @param m_other_ind_traits_previous A matrix containing preceding cycle lmabda, rhos and death.
#' @param m_coef_ukpds_ind_traits A matrix of coefficients for individual‐level risk equations.
#' @param m_coef_ukpds_other_ind_traits A matrix of coefficients `c("lambda", "rho", "death"` used for calculating risk.
#'
#' @return The updated `m_ind_traits` matrix with the mortality status recorded for the specified time step.
#'
#' @export
#'
mortalityV <- function(
    m_ind_traits,
    m_other_ind_traits,
    m_other_ind_traits_previous,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits) {

  # Calculate new health event occurrence and prior history

  # Define events of interest
  events <- c("amp", "blindness", "chf", "esrd", "ihd", "mi", "stroke", "ulcer")

  # Create vectors containing events and event-history names:
  v_event_cols <- paste0(events, "_event")
  v_hist_cols <- paste0(events, "_hist")

  # Get the maximum across those columns, for the given time_step
  # Calculate any new health event
  new_event <- max(m_ind_traits[, v_event_cols])
  # Calculate any prior history of health events
  any_history <- max(m_ind_traits[, v_hist_cols])

  # Determine event-history combinations
  nhne <- new_event == 0 & any_history == 0 # No history, no event
  yhne <- new_event == 0 & any_history == 1 # Yes history, no event
  nhye <- new_event == 1 & any_history == 0 # No history, new event
  yhye <- new_event == 1 & any_history == 1 # Yes history, new event

  # Mortality calculations using Gompertz and logistic models
  death_nhne <- gompertz_eventV(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome = "death_nhne"
  )

  death_yhne <- gompertz_eventV(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome = "death_yhne"
  )

  death_nhye <- logistic_eventV(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome = "death_1st_event"
  )

  death_yhye <- logistic_eventV(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome = "death_yhye"
  )

  # Calculate new mortality status
  new_death <- nhne * death_nhne + yhne * death_yhne + nhye * death_nhye +
    yhye * death_yhye

  # Update the mortality status in the matrix for the given time step
  m_other_ind_traits[, "death"] <- new_death + m_other_ind_traits_previous[, "death"]

  return(m_other_ind_traits)
}

#' @title Calculate Transition Probability Based on a Weibull Model (vectorised)
#'
#' @description This function calculates patient-specific factors, cumulative hazards,
#' and the transition probability for a given health outcome (e.g., "ihd").
#' The function updates the provided `m_ind_traits` matrix with the event occurrence
#' at the specified time step.
#'
#' @param m_ind_traits A matrix to store patient data.
#' @param m_coef_ukpds_ind_traits A matrix of coefficients for individual‐level risk equations.
#' @param m_coef_ukpds_other_ind_traits A matrix of coefficients `c("lambda", "rho", "death"` used for calculating risk.
#' @param health_outcome A character string specifying the health outcome equation (e.g., "ihd").
#'
#' @return Vector of logical values for whether the event would occur.
#'
#' @export
#'
weibull_eventV <- function(
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
  patient_factors <- m_ind_traits %*%
    m_coef_ukpds_ind_traits[, health_outcome] +
    as.vector(m_coef_ukpds_other_ind_traits["lambda", health_outcome])

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
  event <- (trans_prob > stats::runif(nrow(m_ind_traits))) * 1 # to convert logical to numeric
  colnames(event) <- health_event

  # Return the updated matrix
  return(event)
}

#' @title Predict biomarker values of the next cycle (vectorised)
#'
#' @description This function calculates patient-specific factors to predict the time path of
#' a biomarker (linear progression of risk factors).
#'
#' @param m_ind_traits A matrix to store patient data.
#' @param m_coef_ukpds_ind_traits A matrix of coefficients for individual‐level risk equations.
#' @param m_coef_ukpds_other_ind_traits A matrix of coefficients `c("lambda", "rho", "death"` used for calculating risk.
#' @param v_biomarkers A character vector specifying the names of the biomarkers to be updated.
#'
#' @return The updated transition matrix `m_ind_traits` with new biomarker values in the specified row.
#'
#' @export
#'
update_biomarkersV <- function(
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
      heamo = "haem")) {
  # Extract lagged values
  m_ind_traits[, c("a1c_lag", "bmi_lag", "hdl_lag", "heart_rate_lag", "ldl_lag",
                   "sbp_lag", "wbc_lag")] <-
    m_ind_traits[, c("a1c", "bmi", "hdl_real", "heart_rate_real", "ldl_real",
                     "sbp_real", "wbc")]

  step1 <- m_ind_traits %*% m_coef_ukpds_ind_traits[, v_biomarkers]
  step2 <- m_coef_ukpds_other_ind_traits["lambda", v_biomarkers]

  m_ind_traits[, names(v_biomarkers)] <- step1 + step2[col(step1)]

  # biomarkers transformations
  # divide by 10 transformation
  m_ind_traits[, c("sbp", "heart_rate")] <- m_ind_traits[, c(
    "sbp_real",
    "heart_rate_real"
  )] / 10

  # multiply by 10 transformation
  m_ind_traits[, c("ldl", "hdl")] <- m_ind_traits[, c(
    "ldl_real",
    "hdl_real"
  )] * 10

  # other transformation
  m_ind_traits[, "bmi_lt18_5"] <- as.integer(m_ind_traits[, "bmi"] < 18.5)
  m_ind_traits[, "bmi_gte25"] <- as.integer(m_ind_traits[, "bmi"] >= 25)
  m_ind_traits[, "ldl_gt35"] <- as.integer(m_ind_traits[, "ldl_real"] > 35) / 10

  return(m_ind_traits)
}

#' @title Update Health Events Over Time Steps (vectorised)
#'
#' @param m_ind_traits A matrix to store patient data.
#' @param m_coef_ukpds_ind_traits A matrix of coefficients for individual‐level risk equations.
#' @param m_coef_ukpds_other_ind_traits A matrix of coefficients `c("lambda", "rho", "death"` used for calculating risk.
#'
#' @return The updated transition matrix `m_ind_traits` with new biomarker values in the specified row.
#'
#' @export
#'
update_health_eventsV <- function(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits) {
  # Ensure m_ind_traits remains a matrix
  if (!is.matrix(m_ind_traits)) {
    stop("m_ind_traits must be a matrix")
  }

  # Initialize event variables and update history
  events <- c("amp", "blindness", "chf", "esrd", "ihd", "mi", "stroke", "ulcer")

  # create event and history column names once and save each group in a vector
  v_event_cols <- paste0(events, "_event")
  v_history_cols <- paste0(events, "_hist")

  # Update history columns in one vectorized call
  m_ind_traits[, v_history_cols] <- pmax(
    m_ind_traits[, v_history_cols],
    m_ind_traits[, v_event_cols]
  )

  # Randomize event order
  randomized_events <- sample(events)

  for (event in randomized_events) {
    if (event == "amp") {
      amp1_no_ulcer <- weibull_eventV(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        health_outcome = "amp1_no_ulcer"
      )
      amp1_yes_ulcer <- weibull_eventV(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        health_outcome = "amp1_yes_ulcer"
      )

      m_ind_traits[, "amp_event"] <-
        (amp1_no_ulcer * (m_ind_traits[, "ulcer_hist"] == 0)) +
        (amp1_yes_ulcer * (m_ind_traits[, "ulcer_hist"] == 1))

      # Ensure that this is a new event
      m_ind_traits[, "amp_event"] <-
        m_ind_traits[, "amp_event"] *
        (m_ind_traits[, "amp_hist"] == 0)

      # Calculate amp2 event
      amp2 <- weibull_eventV(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        health_outcome = "amp2"
      )

      m_ind_traits[, "amp_event2"] <-
        amp2 * (m_ind_traits[, "amp_hist"] == 1)
    } else if (event == "mi") {
      mi1_male <- weibull_eventV(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        health_outcome = "mi1_male"
      )
      mi1_female <- weibull_eventV(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        health_outcome = "mi1_female"
      )

      m_ind_traits[, "mi_event"] <-
        (mi1_male * (m_ind_traits[, "female"] == 0)) +
        (mi1_female * (m_ind_traits[, "female"] == 1))

      m_ind_traits[, "mi_event"] <-
        m_ind_traits[, "mi_event"] *
        (m_ind_traits[, "mi_hist"] == 0)

      mi2 <- weibull_eventV(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        health_outcome = "mi2"
      )

      m_ind_traits[, "mi_event"] <-
        (m_ind_traits[, "mi_hist"] == 0) *
        m_ind_traits[, "mi_event"] +
        (m_ind_traits[, "mi_hist"] == 1) * mi2
    } else if (event == "stroke") {
      stroke1 <- weibull_eventV(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        health_outcome = "stroke_1"
      )
      stroke2 <- weibull_eventV(
        m_ind_traits,
        m_coef_ukpds_ind_traits,
        m_coef_ukpds_other_ind_traits,
        health_outcome = "stroke_2"
      )

      m_ind_traits[, "stroke_event"] <-
        (stroke1 * (m_ind_traits[, "stroke_hist"] == 0)) +
        (stroke2 * (m_ind_traits[, "stroke_hist"] == 1))

      m_ind_traits[, "stroke_event"] <-
        m_ind_traits[, "stroke_event"] *
        (m_ind_traits[, "stroke_hist"] == 0)
    } else if (event == "ulcer") {
      m_ind_traits[, "ulcer_event"] <-
        logistic_eventV(
          m_ind_traits,
          m_coef_ukpds_ind_traits,
          m_coef_ukpds_other_ind_traits,
          health_outcome = "ulcer"
        )

      m_ind_traits[, "ulcer_event"] <-
        m_ind_traits[, "ulcer_event"] *
        (m_ind_traits[, "ulcer_hist"] == 0)
    } else {
      m_ind_traits[, paste0(event, "_event")] <-
        weibull_eventV(
          m_ind_traits,
          m_coef_ukpds_ind_traits,
          m_coef_ukpds_other_ind_traits,
          health_outcome = event
        )

      m_ind_traits[, paste0(events, "_event")] <-
        m_ind_traits[, paste0(events, "_event")] *
        (m_ind_traits[, paste0(events, "_hist")] == 0)
    }
  }

  return(m_ind_traits)
}

#' @title Run a microsimulation of UKPDS patient outcomes (vectorised)
#'
#' @description
#' This version of the function calculates the per-cycle costs and QALYs inside
#' the main simulation loop for all patients simultaneously.
#'
#' @param num_i Integer. Number of individual patients to simulate.
#' @param num_cycles Integer. Number of cycles (e.g., years) to run each patient.
#' @param m_ukpds_pop Data frame or list. Baseline population characteristics used by `initialize_patients()`.
#' @param a_coef_ukpds_ind_traits Array. Coefficients for individual‐level risk equations.
#' @param a_coef_ukpds_other_ind_traits A 3D array of coefficients `c("lambda", "rho", "death"` used for calculating risk.
#' @param v_coef_names Character vector. A vector of coefficient names.
#' @param discount_rate Numeric. Annual discount rate (e.g., 0.03 for 3%).
#' @param u_baseline Numeric. Utility weight for baseline (no decrement).
#' @param u_blindness Numeric. Utility decrement for blindness (event or history).
#' @param u_amp Numeric. Utility decrement for amputation (event or history).
#' @param u_chf Numeric. Utility decrement for congestive heart failure (event or history).
#' @param u_esrd Numeric. Utility decrement for end‐stage renal disease (event or history).
#' @param u_ihd Numeric. Utility decrement for ischemic heart disease (event or history).
#' @param u_mi Numeric. Utility decrement for myocardial infarction (event or history).
#' @param u_stroke Numeric. Utility decrement for stroke (event or history).
#' @param u_ulcer Numeric. Utility decrement for foot ulcer (event or history).
#' @param c_baseline Numeric. Cost when alive and free of events.
#' @param c_blindness_e Numeric. Cost for a blindness event.
#' @param c_blindness_c Numeric. Cost for ongoing blindness history.
#' @param c_amp_e Numeric. Cost for an amputation event.
#' @param c_amp_c Numeric. Cost for ongoing amputation history.
#' @param c_chf_e Numeric. Cost for a CHF event.
#' @param c_chf_c Numeric. Cost for ongoing CHF history.
#' @param c_esrd_e Numeric. Cost for an ESRD event.
#' @param c_esrd_c Numeric. Cost for ongoing ESRD history.
#' @param c_ihd_e Numeric. Cost for an IHD event.
#' @param c_ihd_c Numeric. Cost for ongoing IHD history.
#' @param c_mi_e Numeric. Cost for an MI event.
#' @param c_mi_c Numeric. Cost for ongoing MI history.
#' @param c_stroke_e Numeric. Cost for a stroke event.
#' @param c_stroke_c Numeric. Cost for ongoing stroke history.
#' @param c_ulcer_e Numeric. Cost for an ulcer event.
#' @param c_ulcer_c Numeric. Cost for ongoing ulcer history.
#' @param seed_no Numeric. Specifying the random number generator seed number.
#'
#' @return A data.frame with `num_i` rows and columns:
#' \describe{
#'   \item{cost}{Total undiscounted cost per patient.}
#'   \item{qalys}{Total undiscounted QALYs per patient.}
#'   \item{disc_costs}{Total discounted cost per patient.}
#'   \item{disc_qalys}{Total discounted QALYs per patient.}
#' }
#'
#' @examples
#' \dontrun{
#' resultsV <- run_microsimV(
#'     num_i = 25,
#'     num_cycles = 20,
#'     m_ukpds_pop = m_ukpds_pop,
#'     a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
#'     a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
#'     v_coef_names = v_coef_names,
#'     discount_rate = 0.03,
#'     u_baseline = 0.785,
#'     u_blindness = -0.074,
#'     u_amp = -0.280,
#'     u_chf = -0.108,
#'     u_esrd = -0.204,
#'     u_ihd = -0.090,
#'     u_mi = -0.055,
#'     u_stroke = -0.164,
#'     u_ulcer = -0.170,
#'     c_baseline = 1990,
#'     c_blindness_e = 4247,
#'     c_blindness_c = 2206,
#'     c_amp_e = 15153,
#'     c_amp_c = 5328,
#'     c_chf_e = 5650,
#'     c_chf_c = 4277,
#'     c_esrd_e = 43359,
#'     c_esrd_c = 43359,
#'     c_ihd_e = 14001,
#'     c_ihd_c = 3550,
#'     c_mi_e = 9518,
#'     c_mi_c = 3424,
#'     c_stroke_e = 10755,
#'     c_stroke_c = 3534,
#'     c_ulcer_e = 7076,
#'     c_ulcer_c = 1072,
#'     seed_no = 1234)
#' }
#'
#' @export
#'
run_microsimV <- function(
    num_i,
    num_cycles,
    m_ukpds_pop = m_ukpds_pop,
    a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
    v_coef_names = v_coef_names,
    discount_rate = 0.03,
    u_baseline = 0.785,
    u_blindness = -0.074,
    u_amp = -0.280,
    u_chf = -0.108,
    u_esrd = -0.204,
    u_ihd = -0.090,
    u_mi = -0.055,
    u_stroke = -0.164,
    u_ulcer = -0.170,
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
    seed_no = 1234) {

  set.seed(seed_no)

  # --- 1. Data and Array Preparation ---
  m_ukpds_pop <- m_ukpds_pop[1:num_i, ]
  v_ids <- paste("id", 1:num_i, sep = "_")
  v_cycles <- paste("cycle", 0:num_cycles, sep = "_")

  a_all_ind_traits <- array(
    data = NA,
    dim = c(num_i, length(v_coef_names), num_cycles + 1),
    dimnames = list(v_ids, v_coef_names, v_cycles)
  )

  a_all_ind_traits[, , 1] <- m_ukpds_pop[, -1]

  a_ind_traits <- a_all_ind_traits[, 1:62, ]
  a_other_ind_traits <- a_all_ind_traits[, 63:65, ]

  # --- 2. Pre-allocate Matrices for Storing Cycle-by-Cycle Results ---
  sim_cycles <- 1:num_cycles
  m_cycle_costs <- matrix(NA, nrow = num_i, ncol = num_cycles, dimnames = list(v_ids, sim_cycles))
  m_cycle_qalys <- matrix(NA, nrow = num_i, ncol = num_cycles, dimnames = list(v_ids, sim_cycles))
  m_disc_cycle_costs <- matrix(NA, nrow = num_i, ncol = num_cycles, dimnames = list(v_ids, sim_cycles))
  m_disc_cycle_qalys <- matrix(NA, nrow = num_i, ncol = num_cycles, dimnames = list(v_ids, sim_cycles))

  # --- 3. Main Simulation Loop with In-Loop Outcome Calculation ---
  for (time_step in 2:(num_cycles + 1)) {

    # --- Health State and Biomarker Updates ---
    # Note: Vectorized update functions (update_biomarkersV, etc.) are assumed to exist
    a_ind_traits[, , time_step] <- update_biomarkersV(
      m_ind_traits = a_ind_traits[, , time_step - 1],
      m_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits[, , 1],
      m_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits[, , 1]
    )

    a_ind_traits[, "age", time_step] <- a_ind_traits[, "age", time_step - 1] + 1
    a_ind_traits[, "diab_dur", time_step] <- a_ind_traits[, "diab_dur", time_step - 1] + 1
    a_ind_traits[, "diab_dur_log", time_step] <- log(a_ind_traits[, "diab_dur", time_step])

    a_ind_traits[, , time_step] <- update_health_eventsV(
      a_ind_traits[, , time_step],
      m_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits[, , 1],
      m_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits[, , 1]
    )

    a_other_ind_traits[, , time_step] <- mortalityV(
      m_ind_traits = a_ind_traits[, , time_step],
      m_other_ind_traits = a_other_ind_traits[, , time_step],
      m_other_ind_traits_previous = a_other_ind_traits[, , time_step - 1],
      m_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits[, , 1],
      m_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits[, , 1]
    )

    # --- Vectorized Cost and QALY Calculation for the CURRENT cycle ---
    t <- time_step - 1 # Get current simulation cycle number (1, 2, ..., 20)

    # Vector (length = num_i) indicating if patient is alive in the current cycle
    v_alive <- a_other_ind_traits[, "death", time_step] == 0

    # Calculate costs for the current cycle
    v_cost_t <- v_alive * (
      c_baseline +
        a_ind_traits[, "blindness_event", time_step] * c_blindness_e +
        a_ind_traits[, "blindness_hist", time_step]  * c_blindness_c +
        (a_ind_traits[, "amp_event", time_step] | a_ind_traits[, "amp_event2", time_step]) * c_amp_e +
        a_ind_traits[, "amp_hist", time_step] * c_amp_c +
        a_ind_traits[, "chf_event", time_step] * c_chf_e +
        a_ind_traits[, "chf_hist", time_step] * c_chf_c +
        a_ind_traits[, "esrd_event", time_step] * c_esrd_e +
        a_ind_traits[, "esrd_hist", time_step] * c_esrd_c +
        a_ind_traits[, "ihd_event", time_step] * c_ihd_e +
        a_ind_traits[, "ihd_hist", time_step] * c_ihd_c +
        a_ind_traits[, "mi_event", time_step] * c_mi_e +
        a_ind_traits[, "mi_hist", time_step] * c_mi_c +
        a_ind_traits[, "stroke_event", time_step] * c_stroke_e +
        a_ind_traits[, "stroke_hist", time_step] * c_stroke_c +
        a_ind_traits[, "ulcer_event", time_step] * c_ulcer_e +
        a_ind_traits[, "ulcer_hist", time_step] * c_ulcer_c
    )

    # Calculate QALYs for the current cycle
    m_decrements_t <- cbind(
      (a_ind_traits[, "blindness_event", time_step] |
         a_ind_traits[, "blindness_hist", time_step]) * u_blindness,
      (a_ind_traits[, "amp_event", time_step] |
         a_ind_traits[, "amp_event2", time_step] |
         a_ind_traits[, "amp_hist", time_step]) * u_amp,
      (a_ind_traits[, "chf_event", time_step] |
         a_ind_traits[, "chf_hist", time_step]) * u_chf,
      (a_ind_traits[, "esrd_event", time_step] |
         a_ind_traits[, "esrd_hist", time_step]) * u_esrd,
      (a_ind_traits[, "ihd_event", time_step] |
         a_ind_traits[, "ihd_hist", time_step]) * u_ihd,
      (a_ind_traits[, "mi_event", time_step] |
         a_ind_traits[, "mi_hist", time_step]) * u_mi,
      (a_ind_traits[, "stroke_event", time_step] |
         a_ind_traits[, "stroke_hist", time_step]) * u_stroke,
      (a_ind_traits[, "ulcer_event", time_step] |
         a_ind_traits[, "ulcer_hist", time_step]) * u_ulcer
    )
    v_min_decrements_t <- apply(m_decrements_t, 1, min)
    v_qaly_t <- v_alive * (u_baseline + v_min_decrements_t)

    # Store results in their respective matrices
    m_cycle_costs[, t] <- v_cost_t
    m_cycle_qalys[, t] <- v_qaly_t

    # Calculate and store discounted results
    disc_factor_t <- 1 / (1 + discount_rate)^t
    m_disc_cycle_costs[, t] <- v_cost_t * disc_factor_t
    m_disc_cycle_qalys[, t] <- v_qaly_t * disc_factor_t
  }

  # --- 4. Summarize and Return Final Results ---

  # Use rowSums to get the total per patient over all num_cycles
  results <- data.frame(
    id = v_ids,
    total_cost = rowSums(m_cycle_costs),
    total_qaly = rowSums(m_cycle_qalys),
    total_disc_cost = rowSums(m_disc_cycle_costs),
    total_disc_qaly = rowSums(m_disc_cycle_qalys)
  )

  return(results)
}
