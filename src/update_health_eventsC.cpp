#include "run_microSim_cpp.h"
#include <RcppArmadillo.h>
#include <string>
#include <unordered_map>
#include <random>
#include <algorithm>

//' Vectorized Health Event Simulation with Randomization (C++)
//'
//' @description
//' Simulates the occurrence of a set of new health events for a cohort of
//' individuals over a single time step.
//'
//' This function is a core part of the simulation engine and performs three
//' key operations:
//' 1. Updates the health history status for each patient.
//' 2. Processes the list of possible health events in a **random order** to
//'    prevent any bias from a fixed sequence.
//' 3. Simulates each event using the appropriate risk engine (`weibull_eventC` or
//'    `logistic_eventC`), applying special logic for complex events.
//'
//' @details
//' This function uses the column names from the input matrices to find the
//' correct data, making it more robust than functions that rely on hardcoded
//' numeric indices.
//'
//' The following events have special conditional logic:
//' - **`amp`**: First amputation risk depends on ulcer history. Second amputation
//'   risk depends on first amputation history.
//' - **`mi`**: First MI risk is different for males and females. Second MI risk
//'   is different from the first.
//' - **`stroke`**: Second stroke risk is different from the first.
//' - **`ulcer`**: Uses a logistic model instead of a Weibull model.
//'
//' For all events, a new event can only occur if the patient does not have a
//' pre-existing history of that same event.
//'
//' @param m_ind_traits A numeric matrix with named columns containing patient
//'   characteristics. It must include columns for event status (e.g., `"mi_event"`)
//'   and event history (e.g., `"mi_hist"`) for all simulated events.
//' @param m_coef_ukpds_ind_traits A numeric matrix of the main model coefficients,
//'   with named columns corresponding to the risk equations (e.g., `"mi1_male"`).
//' @param m_coef_ukpds_other_ind_traits A numeric matrix of ancillary parameters
//'   (`lambda`, `rho`) with corresponding named columns.
//'
//' @return
//' A `NumericMatrix` with the same dimensions and dimnames as the input
//' `m_ind_traits`, but with updated values in the `_event` and `_hist` columns
//' reflecting the outcomes of the current simulation cycle.
//'
//' @export
//'
// [[Rcpp::export]]
Rcpp::NumericMatrix update_health_eventsC(
    Rcpp::NumericMatrix m_ind_traits,
    Rcpp::NumericMatrix m_coef_ukpds_ind_traits,
    Rcpp::NumericMatrix m_coef_ukpds_other_ind_traits) {

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
      arma::umat a1_no = weibull_eventC(
        traits,
        coef_ind,
        coef_other,
        coef_idx["amp1_no_ulcer"]);

      arma::umat a1_yes= weibull_eventC(
        traits,
        coef_ind,
        coef_other,
        coef_idx["amp1_yes_ulcer"]);

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
