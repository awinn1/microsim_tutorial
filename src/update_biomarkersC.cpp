#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Vectorized Biomarker and Lag Variable Update (C++)
//'
//' @description
//' Updates the biomarker values for a cohort of individuals for a single time
//' step. This function performs three main tasks in sequence:
//' 1. Updates lag variables (e.g., `a1c_lag`) with the current biomarker values.
//' 2. Predicts new values for a set of 8 core biomarkers using a linear model.
//' 3. Computes several transformations based on these new biomarker values (e.g.,
//'    `sbp_real` to `sbp`, creating BMI categories).
//'
//' @details
//' **CRITICAL:** This function is highly optimized and operates on the assumption
//' of a fixed matrix structure. It uses hardcoded column indices to read,
//' update, and transform values. The input matrix `m_ind_traits` **must**
//' conform to the specific column order that the model expects. Modifying the
//' input matrix structure without updating the internal indices in this function
//' will lead to incorrect calculations.
//'
//' @param m_ind_traits A numeric matrix of patient characteristics. This matrix
//'   serves as both input (containing values from the previous cycle) and output
//'   (it is modified in-place with updated values for the current cycle).
//' @param m_coef_ukpds_ind_traits A matrix of the main model coefficients. The
//'   function will use the first 8 columns for its predictions.
//' @param m_coef_ukpds_other_ind_traits A matrix of ancillary parameters. The
//'   function will use the first row of the first 8 columns as the `lambda`
//'   (intercept) values.
//'
//' @return
//' The updated `m_ind_traits` matrix with new values for biomarkers, lag
//' variables, and transformed variables for the current simulation cycle.
//'
// [[Rcpp::export]]
arma::mat update_biomarkersC(
    arma::mat& m_ind_traits,
    arma::mat& m_coef_ukpds_ind_traits,
    arma::mat& m_coef_ukpds_other_ind_traits) {

  // Assign lag values
  arma::uvec row_ids = arma::regspace<arma::uvec>(0, m_ind_traits.n_rows - 1);
  arma::uvec lag_col_ids = {9, 14, 21, 25, 30, 35, 39}; // indices for lag columns
  arma::uvec _col_ids = {8, 11, 23, 27, 32, 37, 38}; // indices for values

  m_ind_traits.submat(row_ids, lag_col_ids) = m_ind_traits.submat(row_ids, _col_ids);

  arma::mat biomarkers = m_ind_traits * m_coef_ukpds_ind_traits.cols(0,7);
  arma::mat lambda = m_coef_ukpds_other_ind_traits.submat(0,0,0,7);

  biomarkers.each_row() += lambda;

  // Assign the predictions back to m_ind_traits
  arma::uvec col_ids = {8, 37, 32, 23, 11, 27, 38, 41}; // indices
  m_ind_traits.submat(row_ids, col_ids) = biomarkers;

  // Biomarkers transformations
  // Divide by 10 transformations (sbp and heart_rate)
  arma::uvec col_ids2 = {37,27}; // indices of sbp_real and heart_rate_real
  arma::uvec col_ids3 = {34,24}; // indices of sbp and heart_rate
  m_ind_traits.submat(row_ids, col_ids3) = m_ind_traits.submat(row_ids, col_ids2) / 10  ;

  // Multiply by 10 transformations (ldl and hdl)
  arma::uvec col_ids4 = {32, 23}; // indices of ldl_real and hdl_real
  arma::uvec col_ids5 = {28, 20}; // indices of ldl and hdl
  m_ind_traits.submat(row_ids, col_ids5) = m_ind_traits.submat(row_ids, col_ids4) * 10  ;

  // Other transformations
  m_ind_traits.col(12) = arma::conv_to<arma::vec>::from(m_ind_traits.col(11) < 18.5);
  m_ind_traits.col(13) = arma::conv_to<arma::vec>::from(m_ind_traits.col(11) >= 25);
  m_ind_traits.col(29) = arma::conv_to<arma::vec>::from(m_ind_traits.col(32) > 35) / 10;

  return m_ind_traits;
}
