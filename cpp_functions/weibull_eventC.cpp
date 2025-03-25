#include <RcppArmadillo.h>
#include <unordered_map>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
auto weibull_eventC(arma::mat& m_ind_traits,  
                    NumericMatrix m_coef_ukpds_ind_traits,
                    arma::mat& m_coef_ukpds_other_ind_traits,
                    CharacterVector health_outcome
                    ) {
  
  // Create a vector to map health_outcome to health_event since they are not the same
  CharacterVector v_health_outcome_event = CharacterVector::create(
    _["amp1_no_ulcer"] = "amp_event",
    _["amp1_yes_ulcer"] = "amp_event",
    _["amp2"] = "amp_event2",
    _["mi1_male"] = "mi_event",
    _["mi1_female"] = "mi_event",
    _["mi2"] = "mi_event",
    _["stroke_1"] = "stroke_event",
    _["stroke_2"] = "stroke_event",
    _["ulcer"] = "ulcer_event",
    _["blindness"] = "blindness_event",
    _["chf"] = "chf_event",
    _["esrd"] = "esrd_event",
    _["ihd"] = "ihd_event",
    _["mi"] = "mi_event"
  );
  // Extract health_event from the v_health_outcome_event vector
  Rcpp::CharacterVector health_event = v_health_outcome_event[health_outcome];
  
  // Extract health_outcome index using its name and convert it to arma::uvec with respecting 0 indexing
  Rcpp::CharacterVector m_coef_ukpds_ind_traits_cols = colnames(m_coef_ukpds_ind_traits);
  IntegerVector idx = match(health_outcome, m_coef_ukpds_ind_traits_cols);
  arma::uvec col_idx = as<arma::uvec>(idx) - 1;
  
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
  
  // Compute cumulative hazard at time t and t+1
    // Exponentiate patient_factors
  arma::mat patient_factors_exp = arma::exp(patient_factors);
  arma::vec diab_dur = arma::vec(m_ind_traits.col(5));
  arma::vec diab_dur1 = arma::vec(m_ind_traits.col(5)+1);
  
    // Extract rho
  arma::uvec rho = {static_cast<unsigned int>(1)};
  arma::mat rho2 = arma::mat(m_coef_ukpds_other_ind_traits.submat(rho, col_idx));
  double value2 = arma::as_scalar(rho2);
  
  // Combine the variables above to calculate hazard at time t and t+1
  arma::mat cum_hazard_t = patient_factors_exp%arma::pow(diab_dur, value2);
  arma::mat cum_hazard_t1 = patient_factors_exp%arma::pow(diab_dur1, value2);
  
  // Calculate transition probability
  arma::mat diff_cum_hazard = arma::exp(cum_hazard_t - cum_hazard_t1);
  arma::mat trans_prob = 1 - diff_cum_hazard;
  
  double random_number = arma::randu();
  
  // Assign column name
  NumericMatrix event = wrap(trans_prob > random_number);
  colnames(event) = health_event;
  // arma::umat logical_mat = trans_prob > random_number;
  
  return event;
}
