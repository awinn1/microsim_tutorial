#include <RcppArmadillo.h>
#include <unordered_map>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


auto gompertz_eventC(
    arma::mat& m_ind_traits,  
    const arma::mat& m_coef_ukpds_ind_traits,
    const arma::mat& m_coef_ukpds_other_ind_traits,
    int health_outcome_index
) {
  
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


auto logistic_eventC(
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


// [[Rcpp::export]]
auto mortalityC(
    arma::mat& m_ind_traits,  
    arma::mat& m_other_ind_traits,
    arma::mat& m_other_ind_traits_previous,
    arma::mat& m_coef_ukpds_ind_traits,
    arma::mat& m_coef_ukpds_other_ind_traits
) {
  // Setup
  int n_rows = m_ind_traits.n_rows; // Number of individuals/rows
  
  arma::uvec event_col_ids = {43, 47, 49, 51, 53, 55, 58, 60}; // indices for event columns
  arma::uvec hist_col_ids = {45, 48, 50, 52, 54, 56, 59, 61}; // indices for history columns
  
  // Calculate any new health event
  double new_event = m_ind_traits.cols(event_col_ids).max();
  
  // Calculate any prior history of health events
  double any_history = m_ind_traits.cols(hist_col_ids).max();
  
  // Determine event-history combinations
  bool nhne = (new_event == 0) && (any_history == 0); // No history, no event
  bool yhne = (new_event == 0) && (any_history == 1); // Yes history, no event
  bool nhye = (new_event == 1) && (any_history == 0); // No history, new event
  bool yhye = (new_event == 1) && (any_history == 1); // Yes history, new event
  
  arma::mat new_death;  // To store the newly computed mortality status
  
  // Mortality calculations using Gompertz and logistic models
  if(nhne) {
    int health_outcome_index = 22; // Use the original (R) index for this scenario
    auto death_nhne = gompertz_eventC(
      m_ind_traits,
      m_coef_ukpds_ind_traits,
      m_coef_ukpds_other_ind_traits,
      health_outcome_index
    );
    new_death = arma::conv_to<arma::mat>::from(death_nhne);
    
  } else if(yhne) {
    int health_outcome_index = 24; // Use the original (R) index for this scenario
    auto death_yhne = gompertz_eventC(
      m_ind_traits,
      m_coef_ukpds_ind_traits,
      m_coef_ukpds_other_ind_traits,
      health_outcome_index
    );
    new_death = arma::conv_to<arma::mat>::from(death_yhne);
    
  } else if(nhye) {
    int health_outcome_index = 23; // Use the original (R) index for this scenario
    auto death_nhye = logistic_eventC(
      m_ind_traits,
      m_coef_ukpds_ind_traits,
      m_coef_ukpds_other_ind_traits,
      health_outcome_index
    );
    new_death = arma::conv_to<arma::mat>::from(death_nhye);
    
  } else if(yhye) {
    int health_outcome_index = 25; // Use the original (R) index for this scenario
    auto death_yhye = logistic_eventC(
      m_ind_traits,
      m_coef_ukpds_ind_traits,
      m_coef_ukpds_other_ind_traits,
      health_outcome_index
    );
    new_death = arma::conv_to<arma::mat>::from(death_yhye);
    
  } else {
    // In case conditions do not match, default to zero mortality
    new_death = arma::zeros<arma::mat>(n_rows, 1);
  }
  // pdate the mortality status in the matrix
  m_other_ind_traits.col(2) = new_death + m_other_ind_traits_previous.col(2);
  
  return m_other_ind_traits;
  
  
}