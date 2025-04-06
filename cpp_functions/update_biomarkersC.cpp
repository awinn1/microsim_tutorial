#include <RcppArmadillo.h>
#include <unordered_map>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
auto update_biomarkersC(
    arma::mat& m_ind_traits,  
    arma::mat& m_coef_ukpds_ind_traits,
    arma::mat& m_coef_ukpds_other_ind_traits
) {
  
  // Creating a copy for testing to avoid overwriting (delete it later)
  arma::mat m_ind_traits_copy = m_ind_traits;
  
  // Assign lag values
  arma::uvec row_ids = arma::regspace<arma::uvec>(0, m_ind_traits.n_rows - 1);
  arma::uvec lag_col_ids = {9, 14, 21, 25, 30, 35, 39}; // indices for lag columns
  arma::uvec _col_ids = {8, 11, 23, 27, 32, 37, 38}; // indices for values
  
  m_ind_traits_copy.submat(row_ids, lag_col_ids) = m_ind_traits_copy.submat(row_ids, _col_ids);
  
  arma::mat biomarkers = m_ind_traits * m_coef_ukpds_ind_traits.cols(0,7);
  arma::mat lambda = m_coef_ukpds_other_ind_traits.submat(0,0,0,7);
  
  biomarkers.each_row() += lambda;
  
  // Assign the predictions back to m_ind_traits
  arma::uvec col_ids = {8, 37, 32, 23, 11, 27, 38, 41}; // indices 
  m_ind_traits_copy.submat(row_ids, col_ids) = biomarkers;
  // 
  // 
  // Biomarkers transformations
  // Divide by 10 transformations (sbp and heart_rate)
  arma::uvec col_ids2 = {37,27}; // indices of sbp_real and heart_rate_real
  arma::uvec col_ids3 = {34,24}; // indices of sbp and heart_rate
  m_ind_traits_copy.submat(row_ids, col_ids3) = m_ind_traits_copy.submat(row_ids, col_ids2) / 10  ;

  // Multiply by 10 transformations (ldl and hdl)
  arma::uvec col_ids4 = {32, 23}; // indices of ldl_real and hdl_real
  arma::uvec col_ids5 = {28, 20}; // indices of ldl and hdl
  m_ind_traits_copy.submat(row_ids, col_ids5) = m_ind_traits_copy.submat(row_ids, col_ids4) * 10  ;
  
  // Other transformations
  m_ind_traits_copy.col(12) = arma::conv_to<arma::vec>::from(m_ind_traits_copy.col(11) < 18.5);
  m_ind_traits_copy.col(13) = arma::conv_to<arma::vec>::from(m_ind_traits_copy.col(11) >= 25);
  m_ind_traits_copy.col(29) = arma::conv_to<arma::vec>::from(m_ind_traits_copy.col(32) > 35) / 10;
  
  return m_ind_traits_copy;
  
}