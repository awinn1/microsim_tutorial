# # Load functions ------------------------------------------------------------------------------------------------

invisible(lapply(fs::dir_ls("cpp_functions"), Rcpp::sourceCpp))
invisible(lapply(fs::dir_ls("functions"), source))


# Model code used in the test -------------------------------------------------------------------------------------

## Coefficients dataset
df_UKPDS_coef <- read.csv("data/ukpds_coef.csv")

## Population dataset
df_ukpds_pop <- read.csv("data/population.csv")

cycles <- 20


# Replace NAs with 0s to avoid missing values in calculations
df_UKPDS_coef[is.na(df_UKPDS_coef)] <- 0

# Extract parameter names (used as row names)
v_coef_names <- df_UKPDS_coef$Parameter # Get row names from the 'Parameter' column

# Determine the number of parameters (rows)
n_coef_names <- length(v_coef_names) # Count the number of parameters

# Extract factor names (used as column names), excluding the first column
v_factors_names <- as.vector(colnames(df_UKPDS_coef[-1])) # Get column names excluding 'Parameter'

# Determine the number of factors (columns)
n_equa_names <- length(v_factors_names) # Count the number of factors

# allow for bootstrapped coefficients
boot <- 1 # Abdullah: Only one value for now?
rep_names <- paste0("boot_rep_", 1:boot)

# create an array that holds onto everything!
a_coef_ukpds <- array(
  data = NA,
  dim = c(n_coef_names, n_equa_names, boot),
  dimnames = list(v_coef_names, v_factors_names, rep_names)
)




# fill in the array with coefficents from the dataset (The first depth layer is exactly the same as your UKPDS_coef. We can just import it directly)
# Use UKPDS_coef as our first depth layer in the 3D array
a_coef_ukpds[, , 1] <- as.matrix(df_UKPDS_coef[, -1])


# Remove lambda, rho, and death and create an array for individual traits and another array for those three parameters
# This separation is based on Rob's suggestion
a_coef_ukpds_ind_traits <- a_coef_ukpds[1:62, , "boot_rep_1", drop = FALSE]
a_coef_ukpds_other_ind_traits <- a_coef_ukpds[63:65, , , drop = FALSE]


num_i <- nrow(df_ukpds_pop) # number of simulated individuals
# Define the number of time points



v_ids <- paste("id", 1:num_i, sep = "_")

v_cycles <- paste("cycle", 0:cycles, sep = "_")

nbr_coef_names <- n_coef_names

# Create an array with columns for each variable, row for each person, and a
# slice for each period

a_all_ind_traits <- array(
  data = NA,
  dim = c(num_i, nbr_coef_names, cycles + 1),
  dimnames = list(v_ids, v_coef_names, v_cycles)
)


# need this to be the same number of columns as the coefficient table is long/rows
# print(dim(a_all_ind_traits)) # to verify the dimensions
# print(dimnames(a_all_ind_traits)) # to verify the dimension names

# Fill patient trace
# Create a vector for time invariant characteristics
v_time_invariant_characteristics <- c(
  "age_diag",
  "black",
  "indian",
  "female",
  "smoke",
  "a1c_first",
  "bmi_first",
  "heart_rate_first",
  "ldl_first",
  "hdl_first",
  "sbp_first",
  "wbc_first",
  "heamo_first"
)


a_all_ind_traits[
  ,
  v_time_invariant_characteristics,
] <- as.matrix(df_ukpds_pop[, v_time_invariant_characteristics])

# Fill baseline values of the first cycle
a_all_ind_traits[, , 1] <- as.matrix(df_ukpds_pop[, -1])


# break array up into stuff individual traits
a_ind_traits <- a_all_ind_traits[, 1:62, ]
# remaining array that captures lambda, rho and death
a_other_ind_traits <- a_all_ind_traits[, 63:65, ]




# Inputs ----------------------------------------------------------------------------------------------------------

m_ind_traits <- a_ind_traits[, , 1]
m_coef_ukpds_ind_traits <- a_coef_ukpds_ind_traits[, , 1]
m_coef_ukpds_other_ind_traits <- a_coef_ukpds_other_ind_traits[, , 1]
health_outcome <- "amp1_no_ulcer"


# update_biomarkersC ----------------------------------------------------------------------------------------------
# I did not notice much improvement with update_biomarkersC()
microbenchmark::microbenchmark(
  Cpp = update_biomarkersC(
    m_ind_traits = m_ind_traits,
    m_coef_ukpds_ind_traits = m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits = m_coef_ukpds_other_ind_traits
  ),
  R = update_biomarkers2(
    m_ind_traits = m_ind_traits,
    m_coef_ukpds_ind_traits = m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits = m_coef_ukpds_other_ind_traits
  ),
  times = 1000
)

# weibull_event -----------------------------------------------------------------------------------------------
microbenchmark::microbenchmark(
  weibull_Cpp_result = weibull_eventC(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    match(health_outcome, colnames(m_coef_ukpds_ind_traits))
  ),
  weibull_R_result = weibull_event2(
    m_ind_traits = m_ind_traits,
    m_coef_ukpds_ind_traits = m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits = m_coef_ukpds_other_ind_traits,
    health_outcome = health_outcome
  ),
  times = 100
)

# logistic_event --------------------------------------------------------------------------------------------------
microbenchmark::microbenchmark(
  logistic_Cpp_result = logistic_eventC(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    match("ulcer", colnames(m_coef_ukpds_ind_traits))
  ),
  logistic_R_result = logistic_event2(
    m_ind_traits = m_ind_traits,
    m_coef_ukpds_ind_traits = m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits = m_coef_ukpds_other_ind_traits,
    health_outcome = "ulcer"
  ),
  times = 100
)



# gompertz_event --------------------------------------------------------------------------------------------------------
microbenchmark::microbenchmark(
  gompertz_R_result = gompertz_event2(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    health_outcome = "death_nhne"
  ),
  gompertz_Cpp_result = gompertz_eventC(
    m_ind_traits,
    m_coef_ukpds_ind_traits,
    m_coef_ukpds_other_ind_traits,
    match("death_nhne", colnames(m_coef_ukpds_ind_traits))
  ),
  times = 100
)


# update_health_eventsC --------------------------------------------------------------------------------------------------------
profvis::profvis({
  # Replace NAs with 0s to avoid missing values in calculations
  df_UKPDS_coef[is.na(df_UKPDS_coef)] <- 0
  
  # Extract parameter names (used as row names)
  v_coef_names <- df_UKPDS_coef$Parameter # Get row names from the 'Parameter' column
  
  # Determine the number of parameters (rows)
  n_coef_names <- length(v_coef_names) # Count the number of parameters
  
  # Extract factor names (used as column names), excluding the first column
  v_factors_names <- as.vector(colnames(df_UKPDS_coef[-1])) # Get column names excluding 'Parameter'
  
  # Determine the number of factors (columns)
  n_equa_names <- length(v_factors_names) # Count the number of factors
  
  # allow for bootstrapped coefficients
  boot <- 1 # Abdullah: Only one value for now?
  rep_names <- paste0("boot_rep_", 1:boot)
  
  # create an array that holds onto everything!
  a_coef_ukpds <- array(
    data = NA,
    dim = c(n_coef_names, n_equa_names, boot),
    dimnames = list(v_coef_names, v_factors_names, rep_names)
  )
  
  
  
  
  # fill in the array with coefficents from the dataset (The first depth layer is exactly the same as your UKPDS_coef. We can just import it directly)
  # Use UKPDS_coef as our first depth layer in the 3D array
  a_coef_ukpds[, , 1] <- as.matrix(df_UKPDS_coef[, -1])
  
  
  # Remove lambda, rho, and death and create an array for individual traits and another array for those three parameters
  # This separation is based on Rob's suggestion
  a_coef_ukpds_ind_traits <- a_coef_ukpds[1:62, , "boot_rep_1", drop = FALSE]
  a_coef_ukpds_other_ind_traits <- a_coef_ukpds[63:65, , , drop = FALSE]
  
  
  num_i <- nrow(df_ukpds_pop) # number of simulated individuals
  # Define the number of time points
  
  
  
  v_ids <- paste("id", 1:num_i, sep = "_")
  
  v_cycles <- paste("cycle", 0:cycles, sep = "_")
  
  nbr_coef_names <- n_coef_names
  
  # Create an array with columns for each variable, row for each person, and a
  # slice for each period
  
  a_all_ind_traits <- array(
    data = NA,
    dim = c(num_i, nbr_coef_names, cycles + 1),
    dimnames = list(v_ids, v_coef_names, v_cycles)
  )
  
  
  # need this to be the same number of columns as the coefficient table is long/rows
  # print(dim(a_all_ind_traits)) # to verify the dimensions
  # print(dimnames(a_all_ind_traits)) # to verify the dimension names
  
  # Fill patient trace
  # Create a vector for time invariant characteristics
  v_time_invariant_characteristics <- c(
    "age_diag",
    "black",
    "indian",
    "female",
    "smoke",
    "a1c_first",
    "bmi_first",
    "heart_rate_first",
    "ldl_first",
    "hdl_first",
    "sbp_first",
    "wbc_first",
    "heamo_first"
  )
  
  
  a_all_ind_traits[
    ,
    v_time_invariant_characteristics,
  ] <- as.matrix(df_ukpds_pop[, v_time_invariant_characteristics])
  
  # Fill baseline values of the first cycle
  a_all_ind_traits[, , 1] <- as.matrix(df_ukpds_pop[, -1])
  
  
  # break array up into stuff individual traits
  a_ind_traits <- a_all_ind_traits[, 1:62, ]
  # remaining array that captures lambda, rho and death
  a_other_ind_traits <- a_all_ind_traits[, 63:65, ]
})

microbenchmark::microbenchmark(  
  tmpR = update_health_eventsC_(
    a_ind_traits[, , time_step],
    m_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits[, , 1],
    m_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits[, , 1]
  ),
  tmpC = update_health_eventsC(
    a_ind_traits[, , time_step],
    m_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits[, , 1],
    m_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits[, , 1]
  ),
  times = 2
)
# Unit: milliseconds
# expr      min       lq     mean   median       uq      max neval cld
# tmpR 245.1453 245.1453 249.5122 249.5122 253.8791 253.8791     2 a  
# tmpC 203.1318 203.1318 207.0469 207.0469 210.9621 210.9621     2  b 

profvis::profvis({			
  run_microsim_modelC(			
    df_UKPDS_coef = df_UKPDS_coef,			
    df_ukpds_pop = df_ukpds_pop,			
    cycles = 20,			
    model = "abdullah"			
  )			
})
