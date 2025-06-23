#ifndef RUN_MICROSIM_CPP_H
#define RUN_MICROSIM_CPP_H

#include <RcppArmadillo.h>

arma::umat gompertz_eventC(
    arma::mat& m_ind_traits,
    const arma::mat& m_coef_ukpds_ind_traits,
    const arma::mat& m_coef_ukpds_other_ind_traits,
    int health_outcome_index);

arma::umat logistic_eventC(
    arma::mat& m_ind_traits,
    arma::mat& m_coef_ukpds_ind_traits,
    arma::mat& m_coef_ukpds_other_ind_traits,
    int health_outcome_index);

arma::mat mortalityC(
    arma::mat& m_ind_traits,
    arma::mat& m_other_ind_traits,
    arma::mat& m_other_ind_traits_previous,
    arma::mat& m_coef_ukpds_ind_traits,
    arma::mat& m_coef_ukpds_other_ind_traits);

arma::umat weibull_eventC(
    arma::mat& m_ind_traits,
    arma::mat& m_coef_ukpds_ind_traits,
    arma::mat& m_coef_ukpds_other_ind_traits,
    int health_outcome_index);

#endif
