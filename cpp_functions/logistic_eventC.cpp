#include <RcppArmadillo.h>
#include <unordered_map>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
auto logistic_eventC(arma::mat& m_ind_traits,  
                    NumericMatrix m_coef_ukpds_ind_traits,
                    arma::mat& m_coef_ukpds_other_ind_traits,
                    CharacterVector health_outcome
                    ) {
  
  // Extract health_outcome index using its name and convert it to arma::uvec with respecting 0 indexing
  Rcpp::CharacterVector m_coef_ukpds_ind_traits_cols = colnames(m_coef_ukpds_ind_traits);
  IntegerVector idx = match(health_outcome, m_coef_ukpds_ind_traits_cols);
  arma::uvec col_idx = as<arma::uvec>(idx) - 1;
  
  // Convert m_coef_ukpds_ind_traits to arma matrix and extract the coefficients column
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
  arma::mat patient_factors_exp = arma::exp(-patient_factors);

  // Calculate transition probability
  arma::mat trans_prob = 1 - (patient_factors_exp/ (1 + patient_factors_exp));
  
  // Simulate whether the event occurs by comparing with a random uniform value
    double random_number = arma::randu();
    arma::umat event = trans_prob > random_number;
  
  return event;
}
