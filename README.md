
# microsim.tutorial

<!-- badges: start -->
<!-- badges: end -->

This repo contains code and materials for the Microsimulation Using R tutorial, covering patient-level data structures, Rcpp-based model functions, and simulation result analysis. Includes example datasets, R scripts, and documentation to guide users through building efficient microsimulation models in health economics.

Please note that the repository (web folder) name is microsim_tutorial but the package is called microsim.tutorial.

## Installation

You can install the development version of microsim.tutorial like so:

``` r
devtools::install_github(repo = "github.com/awinn1/microsim_tutorial")
library(microsim.tutorial)
```

## The un-Vectorized R version

``` r
## Micorsim R

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
```

## The Vectorized R version

``` r
## Micorsim R - vectorised

resultsV <- run_microsimV(
    num_i = 25,
    num_cycles = 20,
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
    c_ulcer_c = 1072,
    seed_no = 1234)
```

## The C++ (Rcpp) Version

```
## Micorsim Rcpp

resultsC <- run_microsimC(
    num_i = 25,
    num_cycles = 20,
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
    c_ulcer_c = 1072,
    seed_no = 1234)
```

## Parallel Execution of PSA

```
# PSA

psa_results <- microsim.tutorial::run_psa(
  df_ukpds_bootstraps = microsim.tutorial::df_ukpds_bootstraps,
  df_ukpds_coef = microsim.tutorial::df_ukpds_coef,
  m_ukpds_pop = microsim.tutorial::m_ukpds_pop,
  num_i = 25,
  num_cycles = 20,
  parallel = TRUE
)
```
