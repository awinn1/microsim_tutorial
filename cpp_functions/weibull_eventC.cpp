#include <RcppArmadillo.h>
#include <unordered_map>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
// [[Rcpp::export]]
auto weibull_eventC(
    arma::mat& m_ind_traits,  
    arma::mat& m_coef_ukpds_ind_traits,
    arma::mat& m_coef_ukpds_other_ind_traits,
    int health_outcome_index
) {
  
  // Extract the coefficient column
  arma::vec coef = m_coef_ukpds_ind_traits.col(health_outcome_index - 1);
  
  // Matrix multiplication to get the first part of patient_factors
  arma::mat p1 = m_ind_traits * coef;
  
  // Extract lambda from m_coef_ukpds_other_ind_traits related to the health_outcome
  double lambda = m_coef_ukpds_other_ind_traits(0, health_outcome_index - 1);
  
  // Calculate patient_factors
  arma::mat patient_factors = p1 + lambda;
  
  // Compute cumulative hazards at times t and t+1 
  // Exponentiate patients_factors and define diab_dur
  arma::mat patient_factors_exp = arma::exp(patient_factors);
  arma::vec diab_dur = arma::vec(m_ind_traits.col(5));
  arma::vec diab_dur1 = arma::vec(m_ind_traits.col(5)+1);
  
  // Extract rho
  double rho = m_coef_ukpds_other_ind_traits(1, health_outcome_index - 1);
  
  // Combine the variables above to calculate hazard at time t and t+1
  arma::mat cum_hazard_t = patient_factors_exp % arma::pow(diab_dur, rho);
  arma::mat cum_hazard_t1 = patient_factors_exp % arma::pow(diab_dur1, rho);
  
  // Calculate transition probability
  arma::mat diff_cum_hazard = arma::exp(cum_hazard_t - cum_hazard_t1);
  arma::mat trans_prob = 1 - diff_cum_hazard;
  
  // Generate a matrix of random numbers from a uniform distribution
  arma::mat random_numbers = arma::randu(m_ind_traits.n_rows, 1);
  
  // Assign column name
  arma::umat event = trans_prob > random_numbers;
  
  return event;
  
}