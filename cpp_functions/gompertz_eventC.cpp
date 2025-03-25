#include <RcppArmadillo.h>
#include <unordered_map>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
auto gompertz_eventC(arma::mat& m_ind_traits,  
                    NumericMatrix m_coef_ukpds_ind_traits,
                    arma::mat& m_coef_ukpds_other_ind_traits,
                    CharacterVector health_outcome
                    ) {

    
  // Extract health_outcome index using its name and convert it to arma::uvec with respecting 0 indexing
  Rcpp::CharacterVector m_coef_ukpds_ind_traits_cols = colnames(m_coef_ukpds_ind_traits);
  IntegerVector idx = match(health_outcome, m_coef_ukpds_ind_traits_cols);
  arma::uvec col_idx = as<arma::uvec>(idx) - 1;
  
  
  // Generate patient_factors matrix
    // Convert m_coef_ukpds_ind_traits to arma matrix and extract the coef column
    arma::mat A = as<arma::mat>(m_coef_ukpds_ind_traits);
    arma::vec g = arma::vec(A.cols(col_idx));
    
    // Matrix multiplication to get the first part of patient_factors
    arma::mat p1 = m_ind_traits * g;
    
    // Extract lambda from m_coef_ukpds_other_ind_traits related to the health_outcome
    arma::uvec lambda = {static_cast<unsigned int>(0)};
    arma::mat p2 = arma::mat(m_coef_ukpds_other_ind_traits.submat(lambda, col_idx));
    double value = arma::as_scalar(p2);
    
    // Calculate patient_factors
    arma::mat patient_factors = p1 + value;
    // Exponentiate patient_factors
    arma::mat patient_factors_exp = arma::exp(patient_factors);
  
  // Compute cumulative hazard at time t
  arma::vec age = arma::vec(m_ind_traits.col(0));
  arma::uvec rho = {static_cast<unsigned int>(1)};
  arma::mat coef = arma::mat(m_coef_ukpds_other_ind_traits.submat(rho, col_idx));
  double inv_coef = 1 / arma::as_scalar(coef);
  arma::mat p3 = arma::exp(age*coef)-1;
  
  arma::mat cum_hazard_t = inv_coef*(patient_factors_exp%p3);
  
  // Compute cumulative hazard at time t+1
  arma::vec age1 = arma::vec(m_ind_traits.col(0)+1);
  arma::mat p3t1 = arma::exp(age1*coef)-1;

  arma::mat cum_hazard_t1 = inv_coef*(patient_factors_exp%p3t1);
  
  // Calculate transition probability
  arma::mat diff_cum_hazard = arma::exp(cum_hazard_t - cum_hazard_t1);
  arma::mat trans_prob = 1 - diff_cum_hazard;
  
  double random_number = arma::randu();
  
  // Assign column name
  // NumericMatrix event = wrap(trans_prob > random_number);
  arma::umat event = trans_prob > random_number;
  
  return event;
}
