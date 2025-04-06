#include <RcppArmadillo.h>
#include <unordered_map>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
auto gompertz_eventC(
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
    
    // Exponentiate patient_factors
    arma::mat patient_factors_exp = arma::exp(patient_factors);
    
    // Extract rho and inverse rho
    double rho = m_coef_ukpds_other_ind_traits(1, health_outcome_index - 1);
    double inv_rho = 1 / arma::as_scalar(rho);
  
    // Compute cumulative hazard at time t
    arma::vec age = arma::vec(m_ind_traits.col(0));
    arma::mat p3 = arma::exp(age*rho)-1;
    arma::mat cum_hazard_t = inv_rho*(patient_factors_exp%p3);
  
    // Compute cumulative hazard at time t+1
    arma::vec age1 = arma::vec(m_ind_traits.col(0)+1);
    arma::mat p3t1 = arma::exp(age1*rho)-1;
    arma::mat cum_hazard_t1 = inv_rho*(patient_factors_exp%p3t1);
  
    // Calculate transition probability
    arma::mat diff_cum_hazard = arma::exp(cum_hazard_t - cum_hazard_t1);
    arma::mat trans_prob = 1 - diff_cum_hazard;
  
  // Generate a matrix of random numbers from a uniform distribution
  arma::mat random_numbers = arma::randu(m_ind_traits.n_rows, 1);
  
  // Assign column name
  arma::umat event = trans_prob > random_numbers;
  
  return event;
}
