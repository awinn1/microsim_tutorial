#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
auto weibull_eventC(
    arma::mat& m_ind_traits,  
    arma::mat& m_coef_ukpds_ind_traits,
    arma::mat& m_coef_ukpds_other_ind_traits,
    int health_outcome_index
) {
  
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