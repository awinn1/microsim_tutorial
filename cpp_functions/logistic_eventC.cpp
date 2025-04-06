#include <RcppArmadillo.h>
#include <unordered_map>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
auto logistic_eventC(
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
  arma::mat patient_factors_exp = arma::exp(-patient_factors);

  // Calculate transition probability
  arma::mat trans_prob = 1 - (patient_factors_exp/ (1 + patient_factors_exp));
  
  // Generate a matrix of random numbers from a uniform distribution
    arma::mat random_numbers = arma::randu(m_ind_traits.n_rows, 1);
    
    arma::umat event = trans_prob > random_numbers;
    
    return event;
}
