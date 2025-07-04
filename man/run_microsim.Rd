% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_microsim.R
\name{run_microsim}
\alias{run_microsim}
\title{Run a microsimulation of UKPDS patient outcomes}
\usage{
run_microsim(
  num_i,
  num_cycles,
  m_ukpds_pop = m_ukpds_pop,
  a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
  a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
  v_coef_names = v_coef_names,
  discount_rate,
  u_baseline,
  u_blindness,
  u_amp,
  u_chf,
  u_esrd,
  u_ihd,
  u_mi,
  u_stroke,
  u_ulcer,
  c_baseline,
  c_blindness_e,
  c_blindness_c,
  c_amp_e,
  c_amp_c,
  c_chf_e,
  c_chf_c,
  c_esrd_e,
  c_esrd_c,
  c_ihd_e,
  c_ihd_c,
  c_mi_e,
  c_mi_c,
  c_stroke_e,
  c_stroke_c,
  c_ulcer_e,
  c_ulcer_c,
  seed_no = 1234
)
}
\arguments{
\item{num_i}{Integer. Number of individual patients to simulate.}

\item{num_cycles}{Integer. Number of cycles (e.g., years) to run each patient.}

\item{m_ukpds_pop}{Data frame or list. Baseline population characteristics used by \code{initialize_patients()}.}

\item{a_coef_ukpds_ind_traits}{Array. Coefficients for individual‐level risk equations.}

\item{a_coef_ukpds_other_ind_traits}{A 3D array of coefficients \verb{c("lambda", "rho", "death"} used for calculating risk.}

\item{v_coef_names}{Character vector. A vector of coefficient names.}

\item{discount_rate}{Numeric. Annual discount rate (e.g., 0.03 for 3\%).}

\item{u_baseline}{Numeric. Utility weight for baseline (no decrement).}

\item{u_blindness}{Numeric. Utility decrement for blindness (event or history).}

\item{u_amp}{Numeric. Utility decrement for amputation (event or history).}

\item{u_chf}{Numeric. Utility decrement for congestive heart failure (event or history).}

\item{u_esrd}{Numeric. Utility decrement for end‐stage renal disease (event or history).}

\item{u_ihd}{Numeric. Utility decrement for ischemic heart disease (event or history).}

\item{u_mi}{Numeric. Utility decrement for myocardial infarction (event or history).}

\item{u_stroke}{Numeric. Utility decrement for stroke (event or history).}

\item{u_ulcer}{Numeric. Utility decrement for foot ulcer (event or history).}

\item{c_baseline}{Numeric. Cost when alive and free of events.}

\item{c_blindness_e}{Numeric. Cost for a blindness event.}

\item{c_blindness_c}{Numeric. Cost for ongoing blindness history.}

\item{c_amp_e}{Numeric. Cost for an amputation event.}

\item{c_amp_c}{Numeric. Cost for ongoing amputation history.}

\item{c_chf_e}{Numeric. Cost for a CHF event.}

\item{c_chf_c}{Numeric. Cost for ongoing CHF history.}

\item{c_esrd_e}{Numeric. Cost for an ESRD event.}

\item{c_esrd_c}{Numeric. Cost for ongoing ESRD history.}

\item{c_ihd_e}{Numeric. Cost for an IHD event.}

\item{c_ihd_c}{Numeric. Cost for ongoing IHD history.}

\item{c_mi_e}{Numeric. Cost for an MI event.}

\item{c_mi_c}{Numeric. Cost for ongoing MI history.}

\item{c_stroke_e}{Numeric. Cost for a stroke event.}

\item{c_stroke_c}{Numeric. Cost for ongoing stroke history.}

\item{c_ulcer_e}{Numeric. Cost for an ulcer event.}

\item{c_ulcer_c}{Numeric. Cost for ongoing ulcer history.}

\item{seed_no}{Numeric. Specifying the random number generator seed number.}
}
\value{
A data.frame with \code{num_i} rows and columns:
\describe{
\item{cost}{Total undiscounted cost per patient.}
\item{qalys}{Total undiscounted QALYs per patient.}
\item{disc_costs}{Total discounted cost per patient.}
\item{disc_qalys}{Total discounted QALYs per patient.}
}
}
\description{
\code{run_microsim()} simulates a cohort of patients through a series of health-state transition cycles,
calculating per-cycle costs, utilities (QALYs), and discounted totals for each patient.
}
\details{
The function loops over each patient, initializing their individual traits and
then simulating health events, biomarker updates, and mortality each cycle.
Utilities are calculated as baseline utility plus the most severe decrement
per cycle, and costs accumulate based on events and histories. All outcomes
are discounted by the specified \code{discount_rate}.
}
\examples{
\dontrun{
results <- run_microsim(
  num_i = 25,
  num_cycles = 50,
  m_ukpds_pop = m_ukpds_pop,
  a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
  a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits,
  v_coef_names = v_coef_names,
  discount_rate = 0.03,
  u_baseline = 0.785,
  u_blindness = -0.074,
  u_amp = -0.280,
  u_chf = -0.108,
  u_esrd = -0.204,
  u_ihd = -0.090,
  u_mi = -0.055,
  u_stroke = -0.164,
  u_ulcer = -0.170,
  c_baseline = 1990,
  c_blindness_e = 4247,
  c_blindness_c = 2206,
  c_amp_e = 15153,
  c_amp_c = 5328,
  c_chf_e = 5650,
  c_chf_c = 4277,
  c_esrd_e = 43359,
  c_esrd_c = 43359,
  c_ihd_e = 14001,
  c_ihd_c = 3550,
  c_mi_e = 9518,
  c_mi_c = 3424,
  c_stroke_e = 10755,
  c_stroke_c = 3534,
  c_ulcer_e = 7076,
  c_ulcer_c = 1072
)
head(results)
}

}
