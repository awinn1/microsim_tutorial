#include <RcppArmadillo.h>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <random>

// [[Rcpp::depends(RcppArmadillo)]]

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

// [[Rcpp::export]]
Rcpp::NumericMatrix update_health_eventsC(
    Rcpp::NumericMatrix m_ind_traits,
    Rcpp::NumericMatrix m_coef_ukpds_ind_traits,
    Rcpp::NumericMatrix m_coef_ukpds_other_ind_traits
) {
  // Extract dimnames
  Rcpp::List dn_traits = m_ind_traits.attr("dimnames");
  Rcpp::CharacterVector trait_names = dn_traits[1];
  Rcpp::List dn_coefs = m_coef_ukpds_ind_traits.attr("dimnames");
  Rcpp::CharacterVector coef_names = dn_coefs[1];
  
  // Build index maps
  std::unordered_map<std::string,int> trait_idx;
  for (int i = 0; i < trait_names.size(); ++i) {
    trait_idx[std::string(trait_names[i])] = i;
  }
  std::unordered_map<std::string,int> coef_idx;
  for (int i = 0; i < coef_names.size(); ++i) {
    coef_idx[std::string(coef_names[i])] = i + 1;
  }
  
  // Convert to Armadillo matrices
  arma::mat traits = Rcpp::as<arma::mat>(m_ind_traits);
  arma::mat coef_ind = Rcpp::as<arma::mat>(m_coef_ukpds_ind_traits);
  arma::mat coef_other = Rcpp::as<arma::mat>(m_coef_ukpds_other_ind_traits);
  
  // Define events
  std::vector<std::string> events = {"amp","blindness","chf","esrd","ihd","mi","stroke","ulcer"};
  
  // 1) Update history columns
  for (const auto &ev : events) {
    traits.col(trait_idx[ev + "_hist"]) = arma::max(
      traits.col(trait_idx[ev + "_hist"]),
      traits.col(trait_idx[ev + "_event"])
    );
  }
  
  // Cache history vectors
  arma::vec ulcer_hist = traits.col(trait_idx["ulcer_hist"]);
  arma::vec amp_hist    = traits.col(trait_idx["amp_hist"]);
  arma::vec mi_hist     = traits.col(trait_idx["mi_hist"]);
  arma::vec stroke_hist = traits.col(trait_idx["stroke_hist"]);
  
  // 2) Randomize event order using std::shuffle
  std::random_device rd;
  std::mt19937 rng(rd());
  std::shuffle(events.begin(), events.end(), rng);
  
  // 3) Process each event
  for (const auto &ev : events) {
    if (ev == "amp") {
      arma::umat a1_no = weibull_eventC(traits, coef_ind, coef_other, coef_idx["amp1_no_ulcer"]);
      arma::umat a1_yes= weibull_eventC(traits, coef_ind, coef_other, coef_idx["amp1_yes_ulcer"]);
      arma::vec amp_evt = arma::conv_to<arma::vec>::from(
        a1_no % (ulcer_hist == 0) + a1_yes % (ulcer_hist == 1)
      );
      arma::vec mask_amp0 = arma::conv_to<arma::vec>::from(amp_hist == 0);
      amp_evt %= mask_amp0;
      traits.col(trait_idx["amp_event"]) = amp_evt;
      
      arma::umat a2 = weibull_eventC(traits, coef_ind, coef_other, coef_idx["amp2"]);
      arma::vec amp2 = arma::conv_to<arma::vec>::from(a2);
      arma::vec mask_amp1 = arma::conv_to<arma::vec>::from(amp_hist == 1);
      amp2 %= mask_amp1;
      traits.col(trait_idx["amp_event2"]) = amp2;
      
    } else if (ev == "mi") {
      arma::umat m1_m = weibull_eventC(traits, coef_ind, coef_other, coef_idx["mi1_male"]);
      arma::umat m1_f = weibull_eventC(traits, coef_ind, coef_other, coef_idx["mi1_female"]);
      arma::vec mi_evt = arma::conv_to<arma::vec>::from(
        m1_m % (traits.col(trait_idx["female"]) == 0) +
          m1_f % (traits.col(trait_idx["female"]) == 1)
      );
      arma::vec mask_mi0 = arma::conv_to<arma::vec>::from(mi_hist == 0);
      mi_evt %= mask_mi0;
      traits.col(trait_idx["mi_event"]) = mi_evt;
      
      arma::umat m2 = weibull_eventC(traits, coef_ind, coef_other, coef_idx["mi2"]);
      arma::vec mi2 = arma::conv_to<arma::vec>::from(m2);
      arma::vec mask_mi1 = arma::conv_to<arma::vec>::from(mi_hist == 1);
      arma::vec combined = mask_mi0 % mi_evt + mask_mi1 % mi2;
      traits.col(trait_idx["mi_event"]) = combined;
      
    } else if (ev == "stroke") {
      arma::umat s1 = weibull_eventC(traits, coef_ind, coef_other, coef_idx["stroke_1"]);
      arma::umat s2 = weibull_eventC(traits, coef_ind, coef_other, coef_idx["stroke_2"]);
      arma::vec strk = arma::conv_to<arma::vec>::from(
        s1 % (stroke_hist == 0) + s2 % (stroke_hist == 1)
      );
      arma::vec mask_strk0 = arma::conv_to<arma::vec>::from(stroke_hist == 0);
      strk %= mask_strk0;
      traits.col(trait_idx["stroke_event"]) = strk;
      
    } else if (ev == "ulcer") {
      arma::umat ul = logistic_eventC(traits, coef_ind, coef_other, coef_idx["ulcer"]);
      arma::vec ulc = arma::conv_to<arma::vec>::from(ul);
      arma::vec mask_ulc0 = arma::conv_to<arma::vec>::from(ulcer_hist == 0);
      ulc %= mask_ulc0;
      traits.col(trait_idx["ulcer_event"]) = ulc;
      
    } else {
      arma::umat evm = weibull_eventC(traits, coef_ind, coef_other, coef_idx.at(ev));
      arma::vec evv = arma::conv_to<arma::vec>::from(evm);
      arma::vec mask_evhist0 = arma::conv_to<arma::vec>::from(traits.col(trait_idx[ev + "_hist"]) == 0);
      evv %= mask_evhist0;
      traits.col(trait_idx[ev + "_event"]) = evv;
    }
  }
  
  // Return result preserving dimnames
  Rcpp::NumericMatrix out = Rcpp::wrap(traits);
  out.attr("dimnames") = m_ind_traits.attr("dimnames");
  return out;
}
