#include "run_microSim_cpp.h"
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Vectorized Weibull Event Simulation (C++)
//'
//' @description
//' Calculates the probability of an event within a one-year time step for a cohort
//' of individuals using a Weibull proportional hazards model. Event occurrence is
//' then simulated based on this probability. This function is implemented in C++
//' via RcppArmadillo for high performance.
//'
//' @param m_ind_traits A numeric matrix of individual characteristics. Each row
//'   represents an individual. **CRITICAL:** The 6th column (index 5) of this
//'   matrix must be the diabetes duration (`diab_dur`).
//' @param m_coef_ukpds_ind_traits A numeric matrix of the main model
//'   coefficients. Each column corresponds to a different health outcome.
//' @param m_coef_ukpds_other_ind_traits A numeric matrix of ancillary parameters.
//'   It must have two rows: the first for `lambda` (location) and the second for
//'   `rho` (shape). Each column corresponds to a health outcome.
//' @param health_outcome_index An integer specifying the column index (1-based)
//'   for the relevant health outcome in the coefficient matrices.
//'
//' @return
//' An integer matrix of size `n x 1` (where `n` is the number of individuals).
//' A value of 1 indicates the event occurred for that individual in the cycle,
//' and 0 indicates it did not.
//'
//' @export
//'
// [[Rcpp::export]]
arma::umat weibull_eventC(
    arma::mat& m_ind_traits,
    arma::mat& m_coef_ukpds_ind_traits,
    arma::mat& m_coef_ukpds_other_ind_traits,
    int health_outcome_index) {

  // Setup
  int n_rows = m_ind_traits.n_rows; // Number of individuals/rows
  int idx = health_outcome_index - 1; // Adjust index

  // Extract the coefficient column
  arma::vec coef = m_coef_ukpds_ind_traits.col(idx);

  // Extract constants
  double lambda = m_coef_ukpds_other_ind_traits(0, idx);
  double rho = m_coef_ukpds_other_ind_traits(1, idx);

  // Extract diabetes duration
  const arma::vec& diab_dur = m_ind_traits.col(5);
  arma::vec diab_dur1 = diab_dur + 1.0; // Duration at next time step

  // Weibull model calculations
  arma::vec patient_factors = m_ind_traits * coef;
  patient_factors += lambda;
  arma::vec patient_factors_exp = arma::exp(patient_factors);

  arma::vec cum_hazard_t  = patient_factors_exp % arma::pow(diab_dur,  rho);
  arma::vec cum_hazard_t1 = patient_factors_exp % arma::pow(diab_dur1, rho);

  arma::vec trans_prob = 1.0 - arma::exp(cum_hazard_t - cum_hazard_t1);

  // Generate a matrix of random numbers from a uniform distribution
  arma::mat random_numbers = arma::randu(n_rows, 1);

  // Determine event occurrence
  arma::umat event = trans_prob > random_numbers;

  return event;
}
