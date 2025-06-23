#include "run_microSim_cpp.h"
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
//' Vectorized Mortality Simulation (C++) - erroneous version
//'
//' @description
//' Simulates mortality for a cohort of individuals for a single time step.
//' The function determines which of four different mortality risk equations to use
//' based on the individuals' health event status in the current cycle.
//'
//' @note
//' **IMPORTANT**: The current implementation classifies the *entire cohort* of
//' patients into a single event/history category based on the maximum event/history
//' status found across all patients, rather than assessing each patient's status
//' individually. This means all patients in a batch are simulated using the
//' same risk equation.
//'
//' @param m_ind_traits A numeric matrix of patient characteristics for the
//'   **current** cycle. It is used to check for new health events and history.
//'   The function assumes specific column indices for event/history flags.
//' @param m_other_ind_traits A numeric matrix for ancillary data (`lambda`, `rho`,
//'   `death`) for the **current** cycle. Its 'death' column (column 3) will be
//'   overwritten by this function.
//' @param m_other_ind_traits_previous A numeric matrix for ancillary data from the
//'   **previous** cycle. This is used to carry forward the cumulative death status.
//' @param m_coef_ukpds_ind_traits The main model coefficients matrix.
//' @param m_coef_ukpds_other_ind_traits The ancillary model coefficients matrix
//'   (containing `lambda` and `rho` values).
//'
//' @details
//' The function routes the entire cohort to one of four models based on a
//' collective check for new events or any prior history:
//' - **No History, No Event (NHNE):** Calls `gompertz_eventC` with index 22.
//' - **Yes History, No Event (YHNE):** Calls `gompertz_eventC` with index 24.
//' - **No History, Yes Event (NHYE):** Calls `logistic_eventC` with index 23.
//' - **Yes History, Yes Event (YHYE):** Calls `logistic_eventC` with index 25.
//'
//' @return
//' Returns the `m_other_ind_traits` matrix with an updated 'death' column
//' (column 3), reflecting the cumulative mortality status.
//'
//' @export
//'
// [[Rcpp::export]]
arma::mat mortalityCz(
    arma::mat& m_ind_traits,
    arma::mat& m_other_ind_traits,
    arma::mat& m_other_ind_traits_previous,
    arma::mat& m_coef_ukpds_ind_traits,
    arma::mat& m_coef_ukpds_other_ind_traits) {
  // Setup
  int n_rows = m_ind_traits.n_rows; // Number of individuals/rows

  arma::uvec event_col_ids = {43, 47, 49, 51, 53, 55, 58, 60}; // indices for event columns
  arma::uvec hist_col_ids = {45, 48, 50, 52, 54, 56, 59, 61}; // indices for history columns

  // Calculate any new health event
  double new_event = m_ind_traits.cols(event_col_ids).max();

  // Calculate any prior history of health events
  double any_history = m_ind_traits.cols(hist_col_ids).max();

  // Determine event-history combinations
  bool nhne = (new_event == 0) && (any_history == 0); // No history, no event
  bool yhne = (new_event == 0) && (any_history == 1); // Yes history, no event
  bool nhye = (new_event == 1) && (any_history == 0); // No history, new event
  bool yhye = (new_event == 1) && (any_history == 1); // Yes history, new event

  arma::mat new_death;  // To store the newly computed mortality status

  // Mortality calculations using Gompertz and logistic models
  if(nhne) {
    int health_outcome_index = 22; // Use the original (R) index for this scenario
    auto death_nhne = gompertz_eventC(
      m_ind_traits,
      m_coef_ukpds_ind_traits,
      m_coef_ukpds_other_ind_traits,
      health_outcome_index);
    new_death = arma::conv_to<arma::mat>::from(death_nhne);

  } else if(yhne) {
    int health_outcome_index = 24; // Use the original (R) index for this scenario
    auto death_yhne = gompertz_eventC(
      m_ind_traits,
      m_coef_ukpds_ind_traits,
      m_coef_ukpds_other_ind_traits,
      health_outcome_index);
    new_death = arma::conv_to<arma::mat>::from(death_yhne);

  } else if(nhye) {
    int health_outcome_index = 23; // Use the original (R) index for this scenario
    auto death_nhye = logistic_eventC(
      m_ind_traits,
      m_coef_ukpds_ind_traits,
      m_coef_ukpds_other_ind_traits,
      health_outcome_index);
    new_death = arma::conv_to<arma::mat>::from(death_nhye);

  } else if(yhye) {
    int health_outcome_index = 25; // Use the original (R) index for this scenario
    auto death_yhye = logistic_eventC(
      m_ind_traits,
      m_coef_ukpds_ind_traits,
      m_coef_ukpds_other_ind_traits,
      health_outcome_index);
    new_death = arma::conv_to<arma::mat>::from(death_yhye);

  } else {
    // In case conditions do not match, default to zero mortality
    new_death = arma::zeros<arma::mat>(n_rows, 1);
  }
  // pdate the mortality status in the matrix
  m_other_ind_traits.col(2) = new_death + m_other_ind_traits_previous.col(2);

  return m_other_ind_traits;
}

// [[Rcpp::depends(RcppArmadillo)]]
//' Vectorized Mortality Simulation (C++)
//'
//' @description
//' Simulates mortality for a cohort of individuals for a single time step.
//' The function determines which of four different mortality risk equations to use
//' based on the individuals' health event status in the current cycle.
//'
//' @param m_ind_traits A numeric matrix of patient characteristics for the
//'   **current** cycle. It is used to check for new health events and history.
//'   The function assumes specific column indices for event/history flags.
//' @param m_other_ind_traits A numeric matrix for ancillary data (`lambda`, `rho`,
//'   `death`) for the **current** cycle. Its 'death' column (column 3) will be
//'   overwritten by this function.
//' @param m_other_ind_traits_previous A numeric matrix for ancillary data from the
//'   **previous** cycle. This is used to carry forward the cumulative death status.
//' @param m_coef_ukpds_ind_traits The main model coefficients matrix.
//' @param m_coef_ukpds_other_ind_traits The ancillary model coefficients matrix
//'   (containing `lambda` and `rho` values).
//'
//' @details
//' The function routes the entire cohort to one of four models based on a
//' collective check for new events or any prior history:
//' - **No History, No Event (NHNE):** Calls `gompertz_eventC` with index 22.
//' - **Yes History, No Event (YHNE):** Calls `gompertz_eventC` with index 24.
//' - **No History, Yes Event (NHYE):** Calls `logistic_eventC` with index 23.
//' - **Yes History, Yes Event (YHYE):** Calls `logistic_eventC` with index 25.
//'
//' @return
//' Returns the `m_other_ind_traits` matrix with an updated 'death' column
//' (column 3), reflecting the cumulative mortality status.
//'
//' @export
//'
// [[Rcpp::export]]
arma::mat mortalityC(
    arma::mat& m_ind_traits,
    arma::mat& m_other_ind_traits,
    arma::mat& m_other_ind_traits_previous,
    arma::mat& m_coef_ukpds_ind_traits,
    arma::mat& m_coef_ukpds_other_ind_traits) {

  // Setup
  arma::uvec event_col_ids = {43, 47, 49, 51, 53, 55, 58, 60}; // indices for event columns
  arma::uvec hist_col_ids = {45, 48, 50, 52, 54, 56, 59, 61};  // indices for history columns

  // --- Per-Patient Status Check ---
  // Get row-wise maximums to check status for each patient individually.
  // The result is a column vector with one entry per patient.
  arma::vec new_event_status = arma::max(m_ind_traits.cols(event_col_ids), 1);
  arma::vec any_history_status = arma::max(m_ind_traits.cols(hist_col_ids), 1);

  // --- Create Logical Masks for Each Patient ---
  // A '1' indicates the patient belongs to that group, '0' otherwise.
  arma::uvec nhne_mask = (new_event_status == 0) && (any_history_status == 0);
  arma::uvec yhne_mask = (new_event_status == 0) && (any_history_status == 1);
  arma::uvec nhye_mask = (new_event_status == 1) && (any_history_status == 0);
  arma::uvec yhye_mask = (new_event_status == 1) && (any_history_status == 1);

  // --- Calculate All Four Potential Mortality Outcomes ---
  // Each function is called for all patients.
  int health_outcome_index = 22; // Use the original (R) index for this scenario
  arma::umat death_nhne = gompertz_eventC(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome_index);

  health_outcome_index = 24;
  arma::umat death_yhne = gompertz_eventC(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome_index);

  health_outcome_index = 23;
  arma::umat death_nhye = logistic_eventC(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome_index);

  health_outcome_index = 25;
  arma::umat death_yhye = logistic_eventC(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome_index);

  // --- Combine Results Using Masks ---
  // Convert masks and results to matrices for element-wise multiplication.
  // Each patient's result is selected from the correct model.
  arma::mat new_death =
    arma::conv_to<arma::mat>::from(nhne_mask) % arma::conv_to<arma::mat>::from(death_nhne) +
    arma::conv_to<arma::mat>::from(yhne_mask) % arma::conv_to<arma::mat>::from(death_yhne) +
    arma::conv_to<arma::mat>::from(nhye_mask) % arma::conv_to<arma::mat>::from(death_nhye) +
    arma::conv_to<arma::mat>::from(yhye_mask) % arma::conv_to<arma::mat>::from(death_yhye);

  // Update the mortality status in the matrix (col 2 is the 3rd column, 'death')
  // by adding the new death status to the previous cumulative total.
  m_other_ind_traits.col(2) = new_death + m_other_ind_traits_previous.col(2);

  return m_other_ind_traits;
}
