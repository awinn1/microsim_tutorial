#include <RcppArmadillo.h>
#include <unordered_map>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


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


// [[Rcpp::export]]
auto mortalityC(
    arma::mat& m_ind_traits,  
    arma::mat& m_other_ind_traits,
    arma::mat& m_other_ind_traits_previous,
    arma::mat& m_coef_ukpds_ind_traits,
    arma::mat& m_coef_ukpds_other_ind_traits
) {
  
  arma::uvec event_col_ids = {43, 47, 49, 51, 53, 55, 58, 60}; // indices for event columns
  arma::uvec hist_col_ids = {45, 48, 50, 52, 54, 56, 59, 61}; // indices for history columns
  
  // Calculate any new health event
  double new_event = m_ind_traits.cols(event_col_ids).max();
  
  // Calculate any prior history of health events
  double any_history = m_ind_traits.cols(hist_col_ids).max();
  
  // Determine event-history combinations
  double nhne = (new_event == 0) && (any_history == 0); // No history, no event
  double yhne = (new_event == 0) && (any_history == 1); // Yes history, no event
  double nhye = (new_event == 1) && (any_history == 0); // No history, new event
  double yhye = (new_event == 1) && (any_history == 1); // Yes history, new event
  
  // Mortality calculations using Gompertz and logistic models
  int health_outcome_index = 22; // Put the original index not the cpp index
  auto death_nhne = gompertz_eventC(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome_index
  );
  
  int health_outcome_index2 = 24; // Put the original index not the cpp index
  auto death_yhne = gompertz_eventC(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome_index2
  );

  int health_outcome_index3 = 23; // Put the original index not the cpp index
  auto death_nhye = logistic_eventC(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome_index3
  );

  int health_outcome_index4 = 25; // Put the original index not the cpp index
  auto death_yhye = logistic_eventC(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome_index4
  );  
  
  // Calculate new mortality status
  arma::mat new_death =
    nhne * arma::conv_to<arma::mat>::from(death_nhne) +
    yhne * arma::conv_to<arma::mat>::from(death_yhne) +
    nhye * arma::conv_to<arma::mat>::from(death_nhye) +
    yhye * arma::conv_to<arma::mat>::from(death_yhye);
  
  // pdate the mortality status in the matrix
  m_other_ind_traits.col(2) = new_death + m_other_ind_traits_previous.col(2);
  
  return m_other_ind_traits;
  
  
}