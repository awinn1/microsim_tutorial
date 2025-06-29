% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_microsimV.R
\name{gompertz_eventV}
\alias{gompertz_eventV}
\title{Calculate Transition Probability Based on a Gompertz Model and Update Patient State (vectorised)}
\usage{
gompertz_eventV(
  m_ind_traits,
  m_coef_ukpds_ind_traits,
  m_coef_ukpds_other_ind_traits,
  health_outcome
)
}
\arguments{
\item{m_ind_traits}{A matrix to store patient data.}

\item{m_coef_ukpds_ind_traits}{A matrix of coefficients for individual‐level risk equations.}

\item{m_coef_ukpds_other_ind_traits}{A matrix of coefficients \verb{c("lambda", "rho", "death"} used for calculating risk.}

\item{health_outcome}{A character string specifying the health outcome equation (e.g., "ihd").}
}
\value{
Vector of logical values for whether the event would occur.
}
\description{
This function calculates patient-specific factors, cumulative hazards,
and the transition probability for mortality.
}
