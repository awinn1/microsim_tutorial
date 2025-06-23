#' Initialize baseline values for multiple patients
#'
#' @param num_patients Numeric. The total number of patients to process.
#' @param m_ukpds_pop A Matrix containing patient characteristics.
#' @param m_ind_traits A matrix to store patient data.
#'
#' @return The updated matrix with initialized patient data.
#'
#' @export
#'
initialize_patients <- function(
    num_patients,
    m_ukpds_pop,
    m_ind_traits) {

  patient <- num_patients

  # 1. Create a vector of column names containing the individual characteristics you want to copy:
  v_ind_traits <- c(
    "age", "age_diag", "black", "indian", "female",
    "diab_dur", "diab_dur_log", "smoke",
    "a1c", "a1c_lag", "a1c_first",
    "bmi", "bmi_lt18_5", "bmi_gte25", "bmi_lag", "bmi_first",
    "egfr", "egfr_lt60", "egfr_gte60",
    "hdl", "hdl_lag", "hdl_first",
    "heart_rate", "heart_rate_lag", "heart_rate_first",
    "ldl", "ldl_gt35", "ldl_lag", "ldl_first",
    "albumin_mm", "sbp", "sbp_lag", "sbp_first",
    "wbc", "wbc_lag", "wbc_first",
    "amp_event", "amp_event2", "amp_hist",
    "atria_fib", "blindness_event", "blindness_hist",
    "chf_event", "chf_hist",
    "esrd_event", "esrd_hist",
    "ihd_event", "ihd_hist",
    "mi_event", "mi_hist",
    "pvd_event", "stroke_event", "stroke_hist",
    "ulcer_event", "ulcer_hist"
  )

  # 2. Assign all these columns in a single step.
  m_ind_traits[1, v_ind_traits] <- m_ukpds_pop[patient, v_ind_traits]

  # 3. Handle any variables that aren't in m_ukpds_pop.
  m_ind_traits[1, "heamo"] <- 15
  m_ind_traits[1, "heamo_first"] <- 15

  # Event history tracking
  event_vars <- c("amp_event", "amp_event2", "amp_hist", "atria_fib",
                  "blindness_event", "blindness_hist", "chf_event", "chf_hist",
                  "esrd_event", "esrd_hist", "ihd_event", "ihd_hist",
                  "mi_event", "mi_hist", "pvd_event", "stroke_event",
                  "stroke_hist", "ulcer_event", "ulcer_hist")

  for (var in event_vars) {
    m_ind_traits[1, var] <-  m_ukpds_pop[patient, var]
  }

  m_ind_traits[1,"sbp_real"]<- m_ind_traits[1,"sbp"]*10
  m_ind_traits[1,"egfr_real"]<- m_ind_traits[1,"egfr"]*10
  m_ind_traits[1,"hdl_real"]<- m_ind_traits[1,"hdl"]/10
  m_ind_traits[1,"heart_rate_real"]<- m_ind_traits[1,"heart_rate"]*10
  m_ind_traits[1,"ldl_real"]<- m_ind_traits[1,"ldl"]/10

  # Atrial Fib and PVD do not update
  m_ind_traits[, "atria_fib"] <- m_ind_traits[1, "atria_fib"]
  m_ind_traits[, "pvd_event"] <- m_ind_traits[1, "pvd_event"]

  return(m_ind_traits)
}

#' Calculate Biomarkers
#'
#' This function calculates patient-specific factors to predict the time path of
#' a biomarker (linear progression of risk factors).
#'
#' @param m_ind_traits A matrix containing patient characteristics over time.
#' @param a_coef_ukpds_ind_traits A 3D array of coefficients used for calculating risk.
#' @param a_coef_ukpds_other_ind_traits A 3D array of coefficients `c("lambda", "rho", "death"` used for calculating risk.
#' @param biomarker_eq A character string specifying the health outcome equation (e.g., "ihd").
#' @param time_step An integer indicating the row in `m_ind_traits` to use for calculations.
#'
#' @return The updated biomarker is stored.
#'
#' @export
#'
biomarker <- function(
    m_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    biomarker_eq,
    time_step) {

  # Calculate patient-specific factors using model coefficients and patient data
  updated_biomarker <- m_ind_traits[max(1,time_step-1),] %*%
    a_coef_ukpds_ind_traits[,  biomarker_eq, 1] +
    a_coef_ukpds_other_ind_traits["lambda",  biomarker_eq, 1]

  return(updated_biomarker)
}

#' Update Multiple Biomarkers in a Transition Matrix
#'
#' This function updates multiple biomarker values in the transition matrix for a given time step.
#'
#' @param m_ind_traits The patient trace, a matrix containing patient data with biomarker and event columns.
#' @param a_coef_ukpds_ind_traits A 3D array of coefficients used for calculating risk.
#' @param a_coef_ukpds_other_ind_traits A 3D array of coefficients `c("lambda", "rho", "death"` used for calculating risk.
#' @param time_step An integer representing the current time step.
#' @param next_row An integer indicating the row in `m_ind_traits` to update with new biomarker values.
#'
#' @return The updated transition matrix `m_ind_traits` with new biomarker values in the specified row.
#'
#' @examples
#' \dontrun{
#' # Example usage
#' m_ind_traits <- update_all_biomarkers(
#'   m_ind_traits,
#'   a_coef_ukpds_ind_traits,
#'   time_step = 1,
#'   next_row = 2
#' )
#' }
#'
#' @export
#'
update_all_biomarkers <- function(
    m_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    time_step,
    next_row) {
  # predict the next period (and perform transformations as needed)
  # the biomarkers use real values of variables, but the event equations use transformed variables
  m_ind_traits[next_row, "a1c"] <- biomarker(
    m_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    biomarker_eq = "hba1c",
    time_step = time_step
  )
  m_ind_traits[next_row, "sbp_real"] <- biomarker(
    m_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    biomarker_eq = "sbp",
    time_step = time_step
  )
  m_ind_traits[next_row, "sbp"] <- m_ind_traits[next_row, "sbp_real"] /10
  m_ind_traits[next_row, "ldl_real"] <- biomarker(
    m_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    biomarker_eq = "ldl",
    time_step = time_step
  )
  m_ind_traits[next_row, "ldl"] <- m_ind_traits[next_row, "ldl_real"] * 10
  m_ind_traits[next_row, "hdl_real"] <- biomarker(
    m_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    biomarker_eq = "hdl",
    time_step = time_step
  )
  m_ind_traits[next_row, "hdl"] <- m_ind_traits[next_row, "hdl_real"] * 10
  m_ind_traits[next_row, "bmi"] <- biomarker(
    m_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    biomarker_eq = "bmi",
    time_step = time_step
  )
  m_ind_traits[next_row, "heart_rate_real"] <- biomarker(
    m_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    biomarker_eq = "heart_rate",
    time_step = time_step
  )
  m_ind_traits[next_row, "heart_rate"] <- m_ind_traits[next_row, "heart_rate_real"] /10
  m_ind_traits[next_row, "wbc"] <- biomarker(
    m_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    biomarker_eq = "wbc",
    time_step = time_step
  )
  m_ind_traits[next_row, "heamo"] <- biomarker(
    m_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    biomarker_eq = "haem",
    time_step = time_step
  )

  # Update lag and first occurrence columns
  m_ind_traits[next_row, "a1c_lag"] <- m_ind_traits[time_step, "a1c"]
  m_ind_traits[next_row, "a1c_first"] <- m_ind_traits[1, "a1c"]

  m_ind_traits[next_row, "bmi_lag"] <- m_ind_traits[time_step, "bmi"]
  m_ind_traits[next_row, "bmi_lt18_5"] <- as.integer(m_ind_traits[next_row, "bmi"] < 18.5)
  m_ind_traits[next_row, "bmi_gte25"] <- as.integer(m_ind_traits[next_row, "bmi"] >= 25)
  m_ind_traits[next_row, "bmi_first"] <- m_ind_traits[1, "bmi"]

  m_ind_traits[next_row, "hdl_lag"] <- m_ind_traits[time_step, "hdl_real"]
  m_ind_traits[next_row, "hdl_first"] <- m_ind_traits[1, "hdl_real"]

  m_ind_traits[next_row, "heart_rate_lag"] <- m_ind_traits[time_step, "heart_rate_real"]
  m_ind_traits[next_row, "heart_rate_first"] <- m_ind_traits[1, "heart_rate_real"]

  m_ind_traits[next_row, "ldl_gt35"] <- as.integer(m_ind_traits[next_row, "ldl_real"] > 35) /10
  m_ind_traits[next_row, "ldl_lag"] <- m_ind_traits[time_step, "ldl_real"]
  m_ind_traits[next_row, "ldl_first"] <- m_ind_traits[1, "ldl_real"]

  m_ind_traits[next_row, "sbp_lag"] <- m_ind_traits[time_step, "sbp_real"]
  m_ind_traits[next_row, "sbp_first"] <- m_ind_traits[1, "sbp_real"]

  m_ind_traits[next_row, "wbc_lag"] <- m_ind_traits[time_step, "wbc"]
  m_ind_traits[next_row, "wbc_first"] <- m_ind_traits[1, "wbc"]

  m_ind_traits[next_row, "heamo_first"] <- m_ind_traits[1, "heamo"]

  m_ind_traits[next_row, "egfr"] <- m_ind_traits[1, "egfr"]
  m_ind_traits[next_row, "egfr_real"] <- m_ind_traits[1, "egfr_real"]
  m_ind_traits[next_row, "egfr_lt60"] <- m_ind_traits[1, "egfr_lt60"]
  m_ind_traits[next_row, "egfr_gte60"] <- m_ind_traits[1, "egfr_gte60"]
  m_ind_traits[next_row, "albumin_mm"] <- m_ind_traits[1, "albumin_mm"]

  # Return updated matrix
  return(m_ind_traits)
}

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
#' @param a_coef_ukpds_other_ind_traits A 3D array of coefficients `c("lambda", "rho", "death"` used for calculating risk.
#' @param health_outcome A character string specifying the health outcome equation (e.g., "ihd").
#' @param health_event A character string specifying the health outcome event in the patient trace.
#' @param time_step An integer indicating the row in `m_ind_traits` to use for calculations.
#'
#' @return Whether the event occurred.
#'
#' @export
#'
weibull_event <- function(
    m_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    health_outcome,
    health_event,
    time_step) {

  # Calculate patient-specific factors using model coefficients and patient data
  patient_factors <- m_ind_traits[time_step,] %*%
    a_coef_ukpds_ind_traits[, health_outcome, 1] +
    as.vector(a_coef_ukpds_other_ind_traits["lambda", health_outcome, 1])

  # Compute cumulative hazard at the current time step
  cum_hazard_t <- exp(patient_factors) *
    (m_ind_traits[time_step, "diab_dur"]^(a_coef_ukpds_other_ind_traits["rho", health_outcome, 1]))

  # Compute cumulative hazard at the next time step (by adding 1 year to diabetes duration)
  cum_hazard_t1 <- exp(patient_factors) *
    ((m_ind_traits[time_step, "diab_dur"] + 1)^(a_coef_ukpds_other_ind_traits["rho", health_outcome, 1]))

  # Calculate transition probability
  trans_prob <- 1 - exp(cum_hazard_t - cum_hazard_t1)

  # Simulate whether the event occurs by comparing with a random uniform value
  event <- trans_prob > stats::runif(1)

  # Return the updated matrix
  return(event)
}

#' Calculate Transition Probability Based on a Logistic Regression and Update Patient State
#'
#' This function calculates patient-specific factors, cumulative hazards,
#' and the transition probability for a given health outcome (e.g., "ihd").
#' The function updates the provided `m_ind_traits` matrix with the event occurrence
#' at the specified time step.
#'
#' @param m_ind_traits A matrix containing patient characteristics over time.
#' @param a_coef_ukpds_ind_traits A 3D array of coefficients used for calculating risk.
#' @param a_coef_ukpds_other_ind_traits A 3D array of coefficients `c("lambda", "rho", "death"` used for calculating risk.
#' @param health_outcome A character string specifying the health outcome equation (e.g., "ihd").
#' @param health_event A character string specifying the health outcome event in the patient trace.
#' @param time_step An integer indicating the row in `m_ind_traits` to use for calculations.
#'
#' @return Whether the event occurred.
#'
#' @export
#'
logistic_event <- function(
    m_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    health_outcome,
    health_event,
    time_step) {

  # Calculate patient-specific factors using model coefficients and patient data
  patient_factors <- m_ind_traits[time_step,] %*%
    a_coef_ukpds_ind_traits[, health_outcome, 1] +
    as.vector(a_coef_ukpds_other_ind_traits["lambda", health_outcome, 1])

  # Calculate transition probability
  trans_prob=1-(exp(-patient_factors)/(1+exp(-patient_factors)))^1

  # Simulate whether the event occurs by comparing with a random uniform value
  event <- trans_prob > stats::runif(1)

  # Return the value
  return(event)
}

#' @title Update Health Events Over Time Steps
#' @description This function updates health events in a patient data matrix (`m_ind_traits`) by applying Weibull
#' and logistic event functions in a randomized order across multiple time steps.
#'
#' @param m_ind_traits A matrix containing patient-level data, including health event history.
#' @param a_coef_ukpds_ind_traits A coefficient matrix used in Weibull and logistic event calculations.
#' @param a_coef_ukpds_other_ind_traits A 3D array of coefficients `c("lambda", "rho", "death"` used for calculating risk.
#' @param time_step An integer indicating the current time step to update events.
#'
#' @return Updated `m_ind_traits` matrix with event and history values updated for the given time step.
#'
#' @examples
#' \dontrun{
#' m_ind_traits <- update_health_events(m_ind_traits, a_coef_ukpds_ind_traits, time_step = 1)
#' }
#'
#' @export
#'
update_health_events <- function(
    m_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    time_step) {

  # Ensure m_ind_traits remains a matrix
  if (!is.matrix(m_ind_traits)) {
    stop("m_ind_traits must be a matrix.")
  }

  # Initialize event variables and update history
  events <- c("amp", "blindness", "chf", "esrd", "ihd", "mi", "stroke", "ulcer")

  # create event and history column names once and save each group in a vector
  v_event_cols      <- paste0(events, "_event")
  v_history_cols    <- paste0(events, "_hist")

  # Update history columns in one vectorized call
  m_ind_traits[time_step, v_event_cols] <- 0
  m_ind_traits[time_step, v_history_cols] <- pmax(
    m_ind_traits[max(1, time_step - 1), v_history_cols],
    m_ind_traits[max(1, time_step - 1), v_event_cols]
  )

  m_ind_traits[time_step, "amp_event2"] <- 0
  # Randomize event order
  randomized_events <- sample(events)

  for (events in randomized_events) {
    if (events == "amp") {
      amp1_no_ulcer <- weibull_event(
        m_ind_traits,
        a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits,
        health_outcome = "amp1_no_ulcer",
        health_event = "amp_event",
        time_step = time_step
      )
      amp1_yes_ulcer <- weibull_event(
        m_ind_traits,
        a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits,
        health_outcome = "amp1_yes_ulcer",
        health_event = "amp_event",
        time_step = time_step
      )

      m_ind_traits[time_step, "amp_event"] <-
        (amp1_no_ulcer * (m_ind_traits[time_step, "ulcer_hist"] == 0)) +
        (amp1_yes_ulcer * (m_ind_traits[time_step, "ulcer_hist"] == 1))

      #ensure that this is a new event
      m_ind_traits[time_step, "amp_event"] <-
        m_ind_traits[time_step, "amp_event"] *
        (m_ind_traits[time_step, "amp_hist"] == 0)

      amp2 <- weibull_event(
        m_ind_traits,
        a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits,
        health_outcome = "amp2",
        health_event = "amp_event2",
        time_step = time_step
      )
      m_ind_traits[time_step, "amp_event2"] <- 0
      m_ind_traits[time_step, "amp_event2"] <- amp2 *
        (m_ind_traits[time_step, "amp_hist"] == 1)

    } else if (events == "mi") {
      mi1_male <- weibull_event(
        m_ind_traits,
        a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits,
        health_outcome = "mi1_male",
        health_event = "mi_event",
        time_step = time_step
      )
      mi1_female <- weibull_event(
        m_ind_traits,
        a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits,
        health_outcome = "mi1_female",
        health_event = "mi_event",
        time_step = time_step
      )

      m_ind_traits[time_step, "mi_event"] <-
        (mi1_male * (m_ind_traits[time_step, "female"] == 0)) +
        (mi1_female * (m_ind_traits[time_step, "female"] == 1))
      m_ind_traits[time_step, "mi_event"] <-
        m_ind_traits[time_step, "mi_event"] *
        (m_ind_traits[time_step, "mi_hist"] == 0)

      mi2 <- weibull_event(
        m_ind_traits,
        a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits,
        health_outcome = "mi2",
        health_event = "mi_event",
        time_step = time_step
      )
      m_ind_traits[time_step, "mi_event"] <-
        (m_ind_traits[time_step, "mi_hist"] == 0) *
        m_ind_traits[time_step, "mi_event"] +
        (m_ind_traits[time_step, "mi_hist"] == 1) * mi2

    } else if (events == "stroke") {
      stroke1 <- weibull_event(
        m_ind_traits,
        a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits,
        health_outcome = "stroke_1",
        health_event = "stroke_event",
        time_step = time_step
      )
      stroke2 <- weibull_event(
        m_ind_traits,
        a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits,
        health_outcome = "stroke_2",
        health_event = "stroke_event",
        time_step = time_step
      )

      m_ind_traits[time_step, "stroke_event"] <-
        (stroke1 * (m_ind_traits[time_step, "stroke_hist"] == 0)) +
        (stroke2 * (m_ind_traits[time_step, "stroke_hist"] == 1))
      m_ind_traits[time_step, "stroke_event"] <-
        m_ind_traits[time_step, "stroke_event"] *
        (m_ind_traits[time_step, "stroke_hist"] == 0)

    } else if (events == "ulcer") {
      m_ind_traits[time_step, "ulcer_event"] <- logistic_event(
        m_ind_traits,
        a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits,
        health_outcome = "ulcer",
        health_event = "ulcer_event",
        time_step = time_step
      )
      m_ind_traits[time_step, "ulcer_event"] <-
        m_ind_traits[time_step, "ulcer_event"] *
        (m_ind_traits[time_step, "ulcer_hist"] == 0)

    } else {
      m_ind_traits[time_step, paste0(events, "_event")] <- weibull_event(
        m_ind_traits,
        a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits,
        health_outcome = events,
        health_event = paste0(events, "_event"),
        time_step = time_step
      )

      m_ind_traits[time_step,  paste0(events, "_event")] <-
        m_ind_traits[time_step,  paste0(events, "_event")] *
        (m_ind_traits[time_step,  paste0(events, "_hist")] == 0)
    }
  }

  return(m_ind_traits)
}

#' Calculate Transition Probability Based on a Gompertz Model and Update Patient State
#'
#' This function calculates patient-specific factors, cumulative hazards,
#' and the transition probability for mortality.
#' The function updates the provided `m_ind_traits` matrix with the event occurrence
#' at the specified time step.
#'
#' @param m_ind_traits A matrix containing patient characteristics over time.
#' @param a_coef_ukpds_ind_traits A 3D array of coefficients used for calculating risk.
#' @param a_coef_ukpds_other_ind_traits A 3D array of coefficients `c("lambda", "rho", "death"` used for calculating risk.
#' @param health_outcome A character string specifying the health outcome equation (e.g., "ihd").
#' @param health_event A character string specifying the health outcome event in the patient trace.
#' @param time_step An integer indicating the row in `m_ind_traits` to use for calculations.
#'
#' @return The event occurrence stored.
#'
#' @export
#'
gompertz_event <- function(
    m_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    health_outcome,
    health_event,
    time_step) {

  # Calculate patient-specific factors using model coefficients and patient data
  patient_factors <- m_ind_traits[time_step, ] %*%
    a_coef_ukpds_ind_traits[, health_outcome, 1] +
    as.vector(a_coef_ukpds_other_ind_traits["lambda", health_outcome, 1])

  # Compute cumulative hazard at the current time step
  cum_hazard_t <- (1/a_coef_ukpds_other_ind_traits["rho", health_outcome, 1]) *
    exp(patient_factors) *
    (exp(m_ind_traits[time_step, "age"] *
           (a_coef_ukpds_other_ind_traits["rho", health_outcome, 1])) -1)

  # Compute cumulative hazard at the next time step (by adding 1 year to diabetes duration)
  cum_hazard_t1 <- (1/a_coef_ukpds_other_ind_traits["rho", health_outcome, 1]) *
    exp(patient_factors) *
    (exp((m_ind_traits[time_step, "age"]+1) *
           (a_coef_ukpds_other_ind_traits["rho", health_outcome, 1])) -1)

  # Calculate transition probability
  trans_prob <- 1 - exp(cum_hazard_t - cum_hazard_t1)

  # Simulate whether the event occurs by comparing with a random uniform value
  event <- trans_prob > stats::runif(1)

  # Return the updated matrix
  return(event)
}

#' @title Calculate Mortality Events for a Given Time Step
#'
#' @description This function calculates mortality events for a given time step
#' based on new health events and medical history using Gompertz and logistic models.
#'
#' @param m_ind_traits A matrix containing patient-level data, including health event history.
#' @param m_other_ind_traits A matrix containing lmabda, rhos and death.
#' @param a_coef_ukpds_ind_traits A coefficient matrix used in Gompertz and logistic event calculations.
#' @param a_coef_ukpds_other_ind_traits A 3D array of coefficients `c("lambda", "rho", "death"` used for calculating risk.
#' @param time_step An integer specifying the time step at which mortality should be calculated.
#'
#' @return The updated `m_ind_traits` matrix with the mortality status recorded for the specified time step.
#'
#' @examples
#' \dontrun{
#' m_ind_traits <- mortality(m_ind_traits, a_coef_ukpds_ind_traits, time_step = 5)
#' }
#'
#' @export
#'
mortality <- function(
    m_ind_traits,
    m_other_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    time_step) {

  # Calculate new health event occurrence and prior history

  # Define events of interest
  events <- c("amp", "blindness", "chf", "esrd", "ihd", "mi", "stroke", "ulcer")

  # Create vectors containing events and event-history names:
  v_event_cols <- paste0(events, "_event")
  v_hist_cols  <- paste0(events, "_hist")

  # Get the maximum across those columns, for the given time_step
  # Calculate any new health event
  new_event  <- max(m_ind_traits[time_step, v_event_cols])
  # Calculate any prior history of health events
  any_history <- max(m_ind_traits[time_step, v_hist_cols])

  # Determine event-history combinations
  nhne <- new_event == 0 & any_history == 0  # No history, no event
  yhne <- new_event == 0 & any_history == 1  # Yes history, no event
  nhye <- new_event == 1 & any_history == 0  # No history, new event
  yhye <- new_event == 1 & any_history == 1  # Yes history, new event

  # Mortality calculations using Gompertz and logistic models
  death_nhne <- gompertz_event(
    m_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    health_outcome = "death_nhne",
    health_event = "death_nhne",
    time_step = time_step
  )

  death_yhne <- gompertz_event(
    m_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    health_outcome = "death_yhne",
    health_event = "death_yhne",
    time_step = time_step
  )

  death_nhye <- logistic_event(
    m_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    health_outcome = "death_1st_event",
    health_event = "death_nhye",
    time_step = time_step
  )

  death_yhye <- logistic_event(
    m_ind_traits,
    a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits,
    health_outcome = "death_yhye",
    health_event = "death_yhye",
    time_step = time_step
  )

  # Calculate new mortality status
  new_death <- nhne * death_nhne + yhne * death_yhne + nhye * death_nhye +
    yhye * death_yhye

  # Update the mortality status in the matrix for the given time step
  m_other_ind_traits[time_step, "death"] <- new_death +
    m_other_ind_traits[max(time_step - 1, 1), "death"]

  return(m_other_ind_traits)
}

#' Run a microsimulation of UKPDS patient outcomes
#'
#' `run_microsim()` simulates a cohort of patients through a series of health-state transition cycles,
#' calculating per-cycle costs, utilities (QALYs), and discounted totals for each patient.
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
#' @details
#' The function loops over each patient, initializing their individual traits and
#' then simulating health events, biomarker updates, and mortality each cycle.
#' Utilities are calculated as baseline utility plus the most severe decrement
#' per cycle, and costs accumulate based on events and histories. All outcomes
#' are discounted by the specified `discount_rate`.
#'
#' @examples
#' \dontrun{
#' results <- run_microsim(
#'   num_i = 25,
#'   num_cycles = 50,
#'   m_ukpds_pop = m_ukpds_pop,
#'   a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
#'   a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
#'   v_coef_names = v_coef_names,
#'   discount_rate = 0.03,
#'   u_baseline = 0.785,
#'   u_blindness = -0.074,
#'   u_amp = -0.280,
#'   u_chf = -0.108,
#'   u_esrd = -0.204,
#'   u_ihd = -0.090,
#'   u_mi = -0.055,
#'   u_stroke = -0.164,
#'   u_ulcer = -0.170,
#'   c_baseline = 1990,
#'   c_blindness_e = 4247,
#'   c_blindness_c = 2206,
#'   c_amp_e = 15153,
#'   c_amp_c = 5328,
#'   c_chf_e = 5650,
#'   c_chf_c = 4277,
#'   c_esrd_e = 43359,
#'   c_esrd_c = 43359,
#'   c_ihd_e = 14001,
#'   c_ihd_c = 3550,
#'   c_mi_e = 9518,
#'   c_mi_c = 3424,
#'   c_stroke_e = 10755,
#'   c_stroke_c = 3534,
#'   c_ulcer_e = 7076,
#'   c_ulcer_c = 1072
#' )
#' head(results)
#' }
#'
#' @export
#'
run_microsim <- function(
    num_i,
    num_cycles,
    m_ukpds_pop = m_ukpds_pop,
    a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
    a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
    v_coef_names = v_coef_names,
    discount_rate,
    # utilities (baseline + decrements)
    u_baseline,
    u_blindness,
    u_amp,
    u_chf,
    u_esrd,
    u_ihd,
    u_mi,
    u_stroke,
    u_ulcer,
    # costs
    c_baseline,
    c_blindness_e,
    c_blindness_c,
    c_amp_e,
    c_amp_c,
    c_chf_e,
    c_chf_c,
    c_esrd_e,
    c_esrd_c,
    c_ihd_e,
    c_ihd_c,
    c_mi_e,
    c_mi_c,
    c_stroke_e,
    c_stroke_c,
    c_ulcer_e,
    c_ulcer_c,
    seed_no = 1234) {

  set.seed(seed_no)

  # pre‐allocate output
  cols <- c("cost", "qalys", "disc_costs", "disc_qalys")
  patient_summary <- matrix(
    data = NA,
    nrow = num_i,
    ncol = length(cols),
    dimnames = list(NULL, cols)
  )

  # initialize matrices
  m_ind_traits <- matrix(
    nrow = num_cycles + 1,
    ncol = length(v_coef_names[1:62]),
    dimnames = list(
      paste("cycle", 0:num_cycles, sep ="_"),
      v_coef_names[1:62])
  )
  m_other_ind_traits <- matrix(
    nrow = num_cycles + 1,
    ncol = length(v_coef_names[63:65]),
    dimnames = list(
      paste("cycle", 0:num_cycles, sep ="_"),
      v_coef_names[63:65])
  )

  for (patient in seq_len(num_i)) {
    # initialize
    m_ind_traits <- initialize_patients(patient, m_ukpds_pop, m_ind_traits)
    # carry forward invariants
    for (var in c("age_diag","black","indian","female","smoke")) {
      m_ind_traits[, var] <- m_ind_traits[1, var]
    }

    # Set default values for lambda, rho, and death
    m_other_ind_traits[1, "death"]  <- 0
    m_other_ind_traits[1, "lambda"] <- 0
    m_other_ind_traits[1, "rho"]    <- 1

    # simulate
    for (t in seq_len(num_cycles)) {
      alive_prev <- max(t - 1, 1)

      m_other_ind_traits[t, "death"] <- m_other_ind_traits[alive_prev, "death"]
      # Set default values for lambda and rho
      m_other_ind_traits[t, c("lambda", "rho")] <- 1

      # age & duration
      m_ind_traits[t, "age"]          <- m_ind_traits[alive_prev, "age"] + 1
      m_ind_traits[t, "diab_dur"]     <- m_ind_traits[alive_prev, "diab_dur"] + 1
      m_ind_traits[t, "diab_dur_log"] <- log(m_ind_traits[t, "diab_dur"])

      # events & risk updates
      m_ind_traits       <- update_health_events(
        m_ind_traits,
        a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits,
        time_step = t
      )
      m_other_ind_traits <- mortality(
        m_ind_traits,
        m_other_ind_traits,
        a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits,
        time_step = t
      )
      m_ind_traits       <- update_all_biomarkers(
        m_ind_traits,
        a_coef_ukpds_ind_traits,
        a_coef_ukpds_other_ind_traits,
        time_step = t,
        next_row = t + 1
      )
    }

    # summarize cycles
    cycles <- 1:num_cycles
    m_summary <- matrix(
      data = NA,
      nrow = num_cycles,
      ncol = length(cols),
      dimnames = list(cycles, cols)
    )

    for (t in cycles) {
      alive <- as.double(m_other_ind_traits[t, "death"] == 0)

      # cost accumulation
      m_summary[t, "cost"] <- alive * (
        c_baseline +
          m_ind_traits[t, "blindness_event"] * c_blindness_e +
          m_ind_traits[t, "blindness_hist"]  * c_blindness_c +
          m_ind_traits[t, "amp_event"]       * c_amp_e +
          m_ind_traits[t, "amp_event2"]      * c_amp_e +
          m_ind_traits[t, "amp_hist"]        * c_amp_c +
          m_ind_traits[t, "chf_event"]       * c_chf_e +
          m_ind_traits[t, "chf_hist"]        * c_chf_c +
          m_ind_traits[t, "esrd_event"]      * c_esrd_e +
          m_ind_traits[t, "esrd_hist"]       * c_esrd_c +
          m_ind_traits[t, "ihd_event"]       * c_ihd_e +
          m_ind_traits[t, "ihd_hist"]        * c_ihd_c +
          m_ind_traits[t, "mi_event"]        * c_mi_e +
          m_ind_traits[t, "mi_hist"]         * c_mi_c +
          m_ind_traits[t, "stroke_event"]    * c_stroke_e +
          m_ind_traits[t, "stroke_hist"]     * c_stroke_c +
          m_ind_traits[t, "ulcer_event"]     * c_ulcer_e +
          m_ind_traits[t, "ulcer_hist"]      * c_ulcer_c
      )

      # utility (QALYs this cycle = baseline + worst decrement)
      decrements <- c(
        m_ind_traits[t, "blindness_event"] * u_blindness,
        m_ind_traits[t, "blindness_hist"]  * u_blindness,
        m_ind_traits[t, "amp_event"]       * u_amp,
        m_ind_traits[t, "amp_event2"]      * u_amp,
        m_ind_traits[t, "amp_hist"]        * u_amp,
        m_ind_traits[t, "chf_event"]       * u_chf,
        m_ind_traits[t, "chf_hist"]        * u_chf,
        m_ind_traits[t, "esrd_event"]      * u_esrd,
        m_ind_traits[t, "esrd_hist"]       * u_esrd,
        m_ind_traits[t, "ihd_event"]       * u_ihd,
        m_ind_traits[t, "ihd_hist"]        * u_ihd,
        m_ind_traits[t, "mi_event"]        * u_mi,
        m_ind_traits[t, "mi_hist"]         * u_mi,
        m_ind_traits[t, "stroke_event"]    * u_stroke,
        m_ind_traits[t, "stroke_hist"]     * u_stroke,
        m_ind_traits[t, "ulcer_event"]     * u_ulcer,
        m_ind_traits[t, "ulcer_hist"]      * u_ulcer
      )
      m_summary[t, "qalys"] <- alive * (u_baseline + min(decrements, na.rm = TRUE))

      # discounting
      m_summary[t, "disc_costs"] <- m_summary[t, "cost"]  / (1 + discount_rate)^t
      m_summary[t, "disc_qalys"] <- m_summary[t, "qalys"] / (1 + discount_rate)^t
    }

    # total per patient
    patient_summary[patient, "cost"]       <- sum(m_summary[, "cost"])
    patient_summary[patient, "qalys"]      <- sum(m_summary[, "qalys"])
    patient_summary[patient, "disc_costs"] <- sum(m_summary[, "disc_costs"])
    patient_summary[patient, "disc_qalys"] <- sum(m_summary[, "disc_qalys"])
  }

  return(as.data.frame(patient_summary))
}
