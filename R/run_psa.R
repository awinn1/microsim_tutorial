#' @title Update a Coefficient Data Frame with a Bootstrap Replicate
#'
#' @description
#' This function takes a main coefficient data frame and updates its values for
#' a set of biomarkers using a single row (a specific bootstrap replicate) from
#' a bootstrap results data frame.
#'
#' @param df_ukpds_coef A data frame containing the model's coefficients. It must
#' have a character column named `Parameter` and columns for each biomarker
#' (e.g., "hba1c", "sbp").
#' @param df_ukpds_bootstraps A data frame where each row is a bootstrap replicate and
#' columns correspond to the bootstrapped parameter estimates (e.g.,
#' "re_hb_cons", "re_sbp_fem").
#' @param i An integer specifying the row index of the bootstrap replicate to use
#' from `df_ukpds_bootstraps`.
#'
#' @return
#' A modified `df_ukpds_coef` data frame with updated values from the specified
#' bootstrap replicate.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assume df_coef is your main coefficient data frame and df_boots contains
#' # 100 bootstrap replicates. To update the coefficients using the 5th
#' # bootstrap replicate:
#'
#' updated_df_coef <- update_coef_with_bootstrap(df_coef, df_boots, 5)
#' }
update_coef_with_bootstrap <- function(
    df_ukpds_coef = microsim.tutorial::df_ukpds_coef,
    df_ukpds_bootstraps = microsim.tutorial::df_ukpds_bootstraps,
    i) {

  # Create a mapping structure. Each element name is the target biomarker column.
  # It contains two vectors:
  #   - target_params: The rows to update in df_ukpds_coef$Parameter
  #   - source_cols: The columns to draw values from in df_ukpds_bootstraps
  param_mapping <- list(
    hba1c = list(
      target_params = c("lambda", "indian", "female", "black", "a1c_lag", "diab_dur_log", "a1c_first"),
      source_cols = c("re_hb_cons", "re_hb_asian", "re_hb_fem", "re_hb_black", "re_hb_lag_hba1c", "re_hb_log_year", "re_hb_base_hba1c")
    ),
    sbp = list(
      target_params = c("lambda", "female", "indian", "sbp_lag", "sbp_first"),
      source_cols = c("re_sbp_sbp", "re_sbp_fem", "re_sbp_asian", "re_sbp_lag_sbp", "re_sbp_base_sbp")
    ),
    ldl = list(
      target_params = c("lambda", "female", "black", "indian", "ldl_lag", "diab_dur_log", "ldl_first"),
      source_cols = c("re_ldl_cons", "re_ldl_fem", "re_ldl_black", "re_ldl_asian", "re_ldl_lag_ldl", "re_ldl_log_year", "re_ldl_base_ldl")
    ),
    hdl = list(
      target_params = c("lambda", "female", "black", "hdl_lag", "hdl_first"),
      source_cols = c("re_hdl_cons", "re_hdl_fem", "re_hdl_black", "re_hdl_lag_hdl", "re_hdl_base_hdl")
    ),
    bmi = list(
      target_params = c("lambda", "female", "black", "indian", "bmi_lag", "diab_dur_log", "bmi_first"),
      source_cols = c("re_bmi_cons", "re_bmi_fem", "re_bmi_black", "re_bmi_asian", "re_bmi_lag_bmi", "re_bmi_log_year", "re_bmi_base_bmi")
    ),
    heart_rate = list(
      target_params = c("lambda", "female", "heart_rate_lag", "diab_dur_log", "heart_rate_first"),
      source_cols = c("re_hra_cons", "re_hra_fem", "re_hra_lag_hrate", "re_hra_log_year", "re_hra_base_hrate")
    ),
    wbc = list(
      target_params = c("lambda", "female", "black", "wbc_lag", "diab_dur_log", "wbc_first"),
      source_cols = c("re_wbc_cons", "re_wbc_fem", "re_wbc_black", "re_wbc_lag_wbc", "re_wbc_log_year", "re_wbc_base_wbc")
    ),
    haem = list(
      target_params = c("lambda", "female", "black", "diab_dur_log", "heamo_first"),
      source_cols = c("re_hae_cons", "re_hae_fem", "re_hae_black", "re_hae_log_year", "re_hae_base_haem")
    )
  )

  # Loop through the mapping for each biomarker
  for (biomarker_name in names(param_mapping)) {
    # Get the mapping for the current biomarker
    current_map <- param_mapping[[biomarker_name]]

    # Extract the target and source names
    target_params <- current_map$target_params
    source_cols <- current_map$source_cols

    # Find the row indices in the target data frame
    target_row_indices <- match(target_params, df_ukpds_coef$Parameter)

    # Extract the values from the bootstrap data frame
    bootstrap_values <- df_ukpds_bootstraps[i, source_cols] |> unlist()

    # Perform the assignment
    df_ukpds_coef[target_row_indices, biomarker_name] <- bootstrap_values
  }

  # Return the modified data frame
  return(df_ukpds_coef)
}

#' @title Create UKPDS Coefficient Arrays from Data Frames
#'
#' @description
#' This function processes a base UKPDS coefficient data frame and an optional
#' set of bootstrap replicates to produce the 3D arrays needed for simulation.
#' It handles the creation of arrays with the correct dimensions and populates
#' them with data from each bootstrap replicate.
#'
#' @details
#' The function performs the following steps:
#' 1. Determines the number of bootstrap replicates from `df_bootstraps`.
#' 2. Creates a main 3D array to hold all coefficients for all replicates.
#' 3. Loops through each bootstrap replicate, calling the
#'    `update_coef_with_bootstrap` function to get the correct parameters for
#'    that replicate, and populates a "slice" of the 3D array.
#' 4. Splits the main array into two separate arrays as required by the model:
#'    one for the main individual traits and one for other parameters
#'    (lambda, rho, death).
#'
#' @param df_ukpds_coef A data frame containing the base model coefficients, with
#' a `Parameter` column.
#' @param df_bootstraps An optional data frame where each row is a bootstrap
#' replicate. If `NULL`, only the base coefficients are used.
#'
#' @return
#' A list containing two 3D arrays:
#' \describe{
#'   \item{a_coef_ukpds_ind_traits}{An array for the main individual traits,
#'     based on the first bootstrap replicate.}
#'   \item{a_coef_ukpds_other_ind_traits}{An array for lambda, rho, and death,
#'     containing data from all bootstrap replicates.}
#' }
#'
#' @seealso `update_coef_with_bootstrap()`
#'
#' @export
#'
create_ukpds_arrays <- function(
    df_ukpds_coef = microsim.tutorial::df_ukpds_coef,
    df_bootstraps = NULL) {

  # --- 1. Initial Setup and Dimension Calculation ---

  df_ukpds_coef[is.na(df_ukpds_coef)] <- 0

  v_coef_names <- df_ukpds_coef$Parameter
  n_coef_names <- length(v_coef_names)

  v_factors_names <- as.vector(colnames(df_ukpds_coef[-1]))
  n_equa_names <- length(v_factors_names)

  # Determine the number of bootstrap replicates
  boot <- if (is.null(df_bootstraps)) 1 else nrow(df_bootstraps)
  rep_names <- paste0("boot_rep_", 1:boot)

  # --- 2. Create and Populate the Main 3D Array ---

  a_coef_ukpds <- array(
    data = NA,
    dim = c(n_coef_names, n_equa_names, boot),
    dimnames = list(v_coef_names, v_factors_names, rep_names)
  )

  # Loop through each bootstrap replicate to populate the array
  for (i in 1:boot) {
    # Use the base coefficients for the first replicate if no bootstraps are provided
    # Otherwise, always update with the specified bootstrap replicate
    temp_coef_df <- if (is.null(df_bootstraps)) {
      df_ukpds_coef
    } else {
      update_coef_with_bootstrap(df_ukpds_coef, df_bootstraps, i)
    }

    # Assign the updated matrix to the i-th slice of the 3D array
    a_coef_ukpds[, , i] <- as.matrix(temp_coef_df[, -1])
  }


  # --- 3. Split the Array and Return ---

  # Split the main array into two parts based on the model's requirements
  # Note: The `ind_traits` part only uses the first bootstrap replicate
  a_coef_ukpds_ind_traits <- a_coef_ukpds[1:62, , 1, drop = FALSE]

  # The `other_ind_traits` part keeps all bootstrap replicates
  a_coef_ukpds_other_ind_traits <- a_coef_ukpds[63:65, , , drop = FALSE]

  return(
    list(
      a_coef_ukpds_ind_traits = a_coef_ukpds_ind_traits,
      a_coef_ukpds_other_ind_traits = a_coef_ukpds_other_ind_traits
    )
  )
}

#' reduce computation time.
#'
#' @details
#' The function performs the following steps:
#' 1.  Sets up a parallel backend if `parallel = TRUE`, automatically detecting
#'     and managing the number of CPU cores to use. It provides informative
#'     messages and warnings about the parallel setup.
#' 2.  Splits the bootstrap parameter sets among the parallel workers.
#' 3.  For each parameter set (i.e., each row in `df_ukpds_bootstraps`), it:
#'     a.  Updates the base coefficients using `update_coef_with_bootstrap()`.
#'     b.  Creates the necessary 3D arrays using `create_ukpds_arrays()`.
#'     c.  Runs the core microsimulation model (`run_microsimC()`).
#' 4.  The results are aggregated into a nested list.
#'
#' @param df_ukpds_bootstraps A data frame where each row represents a complete
#'   set of bootstrapped parameters for one PSA iteration.
#' @param df_ukpds_coef A data frame containing the base model coefficients, which
#'   will be used as a template and updated in each iteration.
#' @param m_ukpds_pop A matrix containing the baseline patient population
#'   characteristics for the microsimulation.
#' @param num_i An integer specifying the number of individuals to simulate in
#'   each microsimulation run.
#' @param num_cycles An integer specifying the number of cycles (e.g., years)
#'   for each microsimulation run.
#' @param parallel A logical value. If `TRUE` (the default), the PSA will be run
#'   in parallel across multiple CPU cores.
#' @param n_cores An optional integer specifying the number of cores to use for
#'   parallel execution. If `NULL`, the function will attempt to use half of the
#'   available cores.
#'
#' @return
#' A nested list containing the results. The top-level list has one element per
#' parallel worker. Each of these elements is a list containing the results from
#' `run_microsimC` for each PSA iteration processed by that worker.
#'
#' @seealso `update_coef_with_bootstrap()`, `create_ukpds_arrays()`, `run_microsimC()`
#'
#' @examples
#' \dontrun{
#' # This is a conceptual example as it is computationally intensive and
#' # requires several large data objects to be defined.
#'
#' # Assume df_boots, df_coef, and pop_matrix are all correctly defined.
#'
#' # Run the PSA for 100 individuals across all bootstrap sets in parallel
#' psa_results <- microsim.tutorial::run_psa(
#'   df_ukpds_bootstraps = microsim.tutorial::df_ukpds_bootstraps,
#'   df_ukpds_coef = microsim.tutorial::df_ukpds_coef,
#'   m_ukpds_pop = microsim.tutorial::m_ukpds_pop,
#'   num_i = 25,
#'   num_cycles = 20,
#'   parallel = TRUE
#' )
#' }
#'
#' @export
#'
run_psa <- function(
    df_ukpds_bootstraps = microsim.tutorial::df_ukpds_bootstraps,
    df_ukpds_coef = microsim.tutorial::df_ukpds_coef,
    m_ukpds_pop = microsim.tutorial::m_ukpds_pop,
    num_i,
    num_cycles = 20,
    parallel = TRUE,
    n_cores = NULL) {

  if(isTRUE(parallel)) {
    # Detect the number of cores in the machine being used
    n_all_cores <- parallel::detectCores()

    # Use half of the existing cores if none specified by the user
    if(is.null(n_cores)) {
      if(n_all_cores > 4) {
        n_cores <- floor(n_all_cores / 2)

        message(
          paste0("Using ", n_cores, " cores to run PSA in parallel")
        )
      } else {
        n_cores <- 1

        warning(
          paste0("Device has ", n_cores, " cores. PSA will not run in parallel")
        )
      }
    } else {
      if(n_cores >= n_all_cores) {

        warning(
          paste0(
            "Requested ", n_cores, " cores to run PSA in parallel, but only ",
            n_all_cores, " cores are detectable in this device."
          )
        )

        if(n_all_cores > 4) {

          n_cores <- floor(n_all_cores / 2)
        } else {

          n_cores <- n_all_cores - 2
        }

        message(
          paste0("Using ", n_cores, " cores to run PSA in parallel")
        )
      }
    }

    n_cores <- min(n_cores, nrow(df_ukpds_bootstraps))

    # Create a cluster using the cores
    clust <- parallel::makeCluster(n_cores)
    # # Export the global environment variables to each core.
    # parallel::clusterExport(
    #   cl = clust,
    #   varlist =  ls(envir = .GlobalEnv),
    #   envir = .GlobalEnv
    # )
    # # Export the function environment variables to each core.
    # parallel::clusterExport(
    #   cl = clust,
    #   varlist =  ls(),
    #   envir = environment()
    # )

    # Split the PSA data.frame for efficient parallel computation
    l_split_df_ukpds_bootstraps <- split(
      x = as.data.frame(df_ukpds_bootstraps),
      f = rep(
        x = 1:n_cores,
        length.out = nrow(df_ukpds_bootstraps)
      )
    )

    # Run the model PSA in parallel
    start_time <- Sys.time()
    cat(paste("Running PSA at", start_time |> format("%Y-%m-%d-%H:%M") , "\n"))
    a_psa_results <- parallel::parLapply(
      X = l_split_df_ukpds_bootstraps,
      cl = clust,
      fun = function(df_ukpds_bootstraps) {

        lapply(
          X = 1:nrow(df_ukpds_bootstraps),
          FUN = function(i) {

            # Update coefficients dataset
            df_ukpds_coef <- microsim.tutorial::update_coef_with_bootstrap(
              df_ukpds_coef = microsim.tutorial::df_ukpds_coef,
              df_ukpds_bootstraps = microsim.tutorial::df_ukpds_bootstraps,
              i
            )

            # Initialize coefficients arrays
            l_coef_arrays <- microsim.tutorial::create_ukpds_arrays(
              df_ukpds_coef = microsim.tutorial::df_ukpds_coef
            )

            # Run the model
            res <- microsim.tutorial::run_microsimC(
              num_i = num_i,
              num_cycles = num_cycles,
              m_ukpds_pop = microsim.tutorial::m_ukpds_pop,
              a_coef_ukpds_ind_traits = l_coef_arrays$a_coef_ukpds_ind_traits,
              a_coef_ukpds_other_ind_traits = l_coef_arrays$a_coef_ukpds_other_ind_traits,
              v_coef_names = microsim.tutorial::v_coef_names
            )

            return(res)
          }
        )
      }
    )
  } else {
    lapply(
      X = 1:nrow(df_ukpds_bootstraps),
      FUN = function(i) {

        # Update coefficients dataset
        df_ukpds_coef <- microsim.tutorial::update_coef_with_bootstrap(
          df_ukpds_coef = microsim.tutorial::df_ukpds_coef,
          df_ukpds_bootstraps = microsim.tutorial::df_ukpds_bootstraps,
          i
        )

        # Initialize coefficients arrays
        l_coef_arrays <- microsim.tutorial::create_ukpds_arrays(
          df_ukpds_coef = microsim.tutorial::df_ukpds_coef
        )

        # Run the model
        res <- microsim.tutorial::run_microsimC(
          num_i = num_i,
          num_cycles = num_cycles,
          m_ukpds_pop = microsim.tutorial::m_ukpds_pop,
          a_coef_ukpds_ind_traits = l_coef_arrays$a_coef_ukpds_ind_traits,
          a_coef_ukpds_other_ind_traits = l_coef_arrays$a_coef_ukpds_other_ind_traits,
          v_coef_names = microsim.tutorial::v_coef_names
        )

        return(res)
      }
    )
  }
}
