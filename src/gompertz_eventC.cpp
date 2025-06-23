#include "run_microSim_cpp.h"
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Vectorized Gompertz Event Simulation (C++)
//'
//' @description
//' Calculates the probability of an event within a one-year time step for a cohort
//' of individuals using a Gompertz proportional hazards model. This function is
//' implemented in C++ via RcppArmadillo for high performance.
//'
//' @param m_ind_traits A numeric matrix of individual characteristics, where each
//'   row represents an individual and columns are variables (e.g., age, biomarkers).
//'   The first column is assumed to be `age`.
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
arma::umat gompertz_eventC(
    arma::mat& m_ind_traits,
    const arma::mat& m_coef_ukpds_ind_traits,
    const arma::mat& m_coef_ukpds_other_ind_traits,
    int health_outcome_index) {

  // Setup
  int n_rows = m_ind_traits.n_rows; // Number of individuals/rows
  int idx = health_outcome_index - 1; // Adjust index

  // Extract the coefficient column
  arma::vec coef = m_coef_ukpds_ind_traits.col(idx);

  // Extract parameters
  double lambda = m_coef_ukpds_other_ind_traits(0, idx);
  double rho = m_coef_ukpds_other_ind_traits(1, idx);
  double inv_rho = 1.0 / rho;

  // Extract age column
  const arma::vec& age = arma::vec(m_ind_traits.col(0));
  arma::vec age1 = age + 1;

  // Gompertz model calculations -----
  arma::vec patient_factors = m_ind_traits * coef;
  patient_factors += lambda;
  arma::vec patient_factors_exp = arma::exp(patient_factors);

  // Compute cumulative hazard at time t
  arma::mat p_t0 = arma::exp(age*rho) - 1.0;
  arma::mat cum_hazard_t = inv_rho*(patient_factors_exp%p_t0);

  // Compute cumulative hazard at time t+1
  arma::mat p_t1 = arma::exp(age1*rho) - 1.0;
  arma::mat cum_hazard_t1 = inv_rho*(patient_factors_exp%p_t1);

  // Calculate transition probability
  arma::mat trans_prob = 1 - arma::exp(cum_hazard_t - cum_hazard_t1);

  // Generate a matrix of random numbers from a uniform distribution
  arma::mat random_numbers = arma::randu(n_rows, 1);

  // Determine event occurrence
  arma::umat event = trans_prob > random_numbers;

  return event;
}
