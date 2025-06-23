#include "run_microSim_cpp.h"
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Vectorized Logistic Event Simulation (C++)
//'
//' @description
//' Calculates the probability of a binary event for a cohort of individuals using
//' a logistic regression model. The event occurrence is then simulated based on
//' this probability. This function is implemented in C++ via RcppArmadillo
//' for high performance.
//'
//' @param m_ind_traits A numeric matrix of individual characteristics (predictors),
//'   where each row represents an individual.
//' @param m_coef_ukpds_ind_traits A numeric matrix of the main model
//'   coefficients. Each column corresponds to a different health outcome.
//' @param m_coef_ukpds_other_ind_traits A numeric matrix of ancillary parameters.
//'   The first row must contain the intercept (`lambda`) for each health outcome.
//' @param health_outcome_index An integer specifying the column index (1-based)
//'   for the relevant health outcome in the coefficient matrices.
//'
//' @return
//' An integer matrix of size `n x 1` (where `n` is the number of individuals).
//' A value of 1 indicates the event occurred, and 0 indicates it did not.
//'
//' @export
//'
// [[Rcpp::export]]
arma::umat logistic_eventC(
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

  // Logistic model calculations

  arma::vec patient_factors = m_ind_traits * coef;
  patient_factors += lambda;
  arma::vec patient_factors_exp = arma::exp(-patient_factors);

  // Calculate transition probability
  arma::mat trans_prob = 1 - (patient_factors_exp/ (1 + patient_factors_exp));

  // Generate a matrix of random numbers from a uniform distribution
  arma::mat random_numbers = arma::randu(n_rows, 1);

  arma::umat event = trans_prob > random_numbers;

  return event;
}
