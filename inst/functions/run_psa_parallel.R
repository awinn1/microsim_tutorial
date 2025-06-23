#' Title
#'
#' @param df_UKPDS_coef 
#' @param df_ukpds_pop 
#' @param df_bootstraps 
#' @param cycles 
#' @param parallel 
#'
#' @returns
#' @export
#'
#' @examples
run_psa <- function(
    df_UKPDS_coef,
    df_ukpds_pop,
    df_bootstraps,
    cycles = 20,
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
    
    n_cores <- min(n_cores, nrow(df_bootstraps))
    
    # Create a cluster using the cores
    clust <- parallel::makeCluster(n_cores)
    # Export the global environment variables to each core.
    parallel::clusterExport(
      cl = clust,
      varlist =  ls(envir = .GlobalEnv),
      envir = .GlobalEnv
    )
    # Export the function environment variables to each core.
    parallel::clusterExport(
      cl = clust,
      varlist =  ls(),
      envir = environment()
    )
    
    # Split the PSA data.frame for efficient parallel computation
    l_split_df_bootstraps <- split(
      x = as.data.frame(df_bootstraps),
      f = rep(
        x = 1:n_cores,
        length.out = nrow(df_bootstraps)
      )
    )
    
    # Run the model PSA in parallel
    start_time <- Sys.time()
    cat(paste("Running PSA at", start_time |> format("%Y-%m-%d-%H:%M") , "\n"))
    a_psa_results <- parallel::parLapply(
      X = l_split_df_bootstraps,
      cl = clust,
      fun = function(df_bootstraps) {
        
        lapply(
          X = 1:nrow(df_bootstraps),
          FUN = function(i) {
            
            # A1c
            df_UKPDS_coef[
              match(c("lambda", "indian", "female", "black", "a1c_lag", "diab_dur_log", "a1c_first"), df_UKPDS_coef$Parameter),
              "hba1c"
            ] <- df_bootstraps[
              i,
              c("re_hb_cons", "re_hb_asian", "re_hb_fem", "re_hb_black", "re_hb_lag_hba1c",
                "re_hb_log_year", "re_hb_base_hba1c")
            ] |> unlist()
            
            # SBP
            df_UKPDS_coef[
              match(c("lambda", "female", "indian", "sbp_lag", "sbp_first"), df_UKPDS_coef$Parameter),
              "sbp"
            ] <- df_bootstraps[
              i,
              c("re_sbp_sbp", "re_sbp_fem", "re_sbp_asian", "re_sbp_lag_sbp", "re_sbp_base_sbp")
            ] |> unlist() # re_sbp_sbp NOT re_sbp_cons
            
            #LDL
            df_UKPDS_coef[
              match(c("lambda", "female", "black", "indian", "ldl_lag", "diab_dur_log", "ldl_first"), df_UKPDS_coef$Parameter),
              "ldl"
            ] <- df_bootstraps[
              i,
              c("re_ldl_cons", "re_ldl_fem", "re_ldl_black", "re_ldl_asian", "re_ldl_lag_ldl",
                "re_ldl_log_year", "re_ldl_base_ldl")
            ] |> unlist()
            
            #HDL
            df_UKPDS_coef[
              match(c("lambda", "female", "black", "hdl_lag", "hdl_first"), df_UKPDS_coef$Parameter),
              "hdl"
            ] <- df_bootstraps[
              i,
              c("re_hdl_cons", "re_hdl_fem", "re_hdl_black", "re_hdl_lag_hdl", "re_hdl_base_hdl")
            ] |> unlist()
            
            #BMI
            df_UKPDS_coef[
              match(c("lambda", "female", "black", "indian", "bmi_lag", "diab_dur_log", "bmi_first"), df_UKPDS_coef$Parameter),
              "bmi"
            ] <- df_bootstraps[
              i,
              c("re_bmi_cons", "re_bmi_fem", "re_bmi_black", "re_bmi_asian", "re_bmi_lag_bmi",
                "re_bmi_log_year", "re_bmi_base_bmi")
            ] |> unlist()
            
            # <3 rate
            df_UKPDS_coef[
              match(c("lambda", "female", "heart_rate_lag", "diab_dur_log", "heart_rate_first"), df_UKPDS_coef$Parameter),
              "heart_rate"
            ]  <- df_bootstraps[
              i,
              c("re_hra_cons", "re_hra_fem", "re_hra_lag_hrate", "re_hra_log_year",
                "re_hra_base_hrate")
            ] |> unlist()
            
            # white blood cell count
            df_UKPDS_coef[
              match(c("lambda", "female", "black", "wbc_lag", "diab_dur_log", "wbc_first"), df_UKPDS_coef$Parameter),
              "wbc"
            ] <- df_bootstraps[
              i,
              c("re_wbc_cons", "re_wbc_fem", "re_wbc_black", "re_wbc_lag_wbc", "re_wbc_log_year",
                "re_wbc_base_wbc")
            ] |> unlist()
            
            # haem
            df_UKPDS_coef[
              match(c("lambda", "female", "black", "diab_dur_log", "heamo_first"), df_UKPDS_coef$Parameter),
              "haem" 
            ] <- df_bootstraps[
              i,
              c("re_hae_cons", "re_hae_fem", "re_hae_black", "re_hae_log_year",
                "re_hae_base_haem")
            ] |> unlist()
            
            # Run the model
            res <- run_microsim_modelC(
              df_UKPDS_coef = df_UKPDS_coef,
              df_ukpds_pop = df_ukpds_pop,
              cycles = 20
            )
            
            return(res)
          }
        )
      }
    )
  } else {
    lapply(
      X = 1:nrow(df_bootstraps),
      FUN = function(i) {
        
        # A1c
        df_UKPDS_coef[
          match(c("lambda", "indian", "female", "black", "a1c_lag", "diab_dur_log", "a1c_first"), df_UKPDS_coef$Parameter),
          "hba1c"
        ] <- df_bootstraps[
          i,
          c("re_hb_cons", "re_hb_asian", "re_hb_fem", "re_hb_black", "re_hb_lag_hba1c",
            "re_hb_log_year", "re_hb_base_hba1c")
        ] |> unlist()
        
        # SBP
        df_UKPDS_coef[
          match(c("lambda", "female", "indian", "sbp_lag", "sbp_first"), df_UKPDS_coef$Parameter),
          "sbp"
        ] <- df_bootstraps[
          i,
          c("re_sbp_sbp", "re_sbp_fem", "re_sbp_asian", "re_sbp_lag_sbp", "re_sbp_base_sbp")
        ] |> unlist() # re_sbp_sbp NOT re_sbp_cons
        
        #LDL
        df_UKPDS_coef[
          match(c("lambda", "female", "black", "indian", "ldl_lag", "diab_dur_log", "ldl_first"), df_UKPDS_coef$Parameter),
          "ldl"
        ] <- df_bootstraps[
          i,
          c("re_ldl_cons", "re_ldl_fem", "re_ldl_black", "re_ldl_asian", "re_ldl_lag_ldl",
            "re_ldl_log_year", "re_ldl_base_ldl")
        ] |> unlist()
        
        #HDL
        df_UKPDS_coef[
          match(c("lambda", "female", "black", "hdl_lag", "hdl_first"), df_UKPDS_coef$Parameter),
          "hdl"
        ] <- df_bootstraps[
          i,
          c("re_hdl_cons", "re_hdl_fem", "re_hdl_black", "re_hdl_lag_hdl", "re_hdl_base_hdl")
        ] |> unlist()
        
        #BMI
        df_UKPDS_coef[
          match(c("lambda", "female", "black", "indian", "bmi_lag", "diab_dur_log", "bmi_first"), df_UKPDS_coef$Parameter),
          "bmi"
        ] <- df_bootstraps[
          i,
          c("re_bmi_cons", "re_bmi_fem", "re_bmi_black", "re_bmi_asian", "re_bmi_lag_bmi",
            "re_bmi_log_year", "re_bmi_base_bmi")
        ] |> unlist()
        
        # <3 rate
        df_UKPDS_coef[
          match(c("lambda", "female", "heart_rate_lag", "diab_dur_log", "heart_rate_first"), df_UKPDS_coef$Parameter),
          "heart_rate"
        ]  <- df_bootstraps[
          i,
          c("re_hra_cons", "re_hra_fem", "re_hra_lag_hrate", "re_hra_log_year",
            "re_hra_base_hrate")
        ] |> unlist()
        
        # white blood cell count
        df_UKPDS_coef[
          match(c("lambda", "female", "black", "wbc_lag", "diab_dur_log", "wbc_first"), df_UKPDS_coef$Parameter),
          "wbc"
        ] <- df_bootstraps[
          i,
          c("re_wbc_cons", "re_wbc_fem", "re_wbc_black", "re_wbc_lag_wbc", "re_wbc_log_year",
            "re_wbc_base_wbc")
        ] |> unlist()
        
        # haem
        df_UKPDS_coef[
          match(c("lambda", "female", "black", "diab_dur_log", "heamo_first"), df_UKPDS_coef$Parameter),
          "haem" 
        ] <- df_bootstraps[
          i,
          c("re_hae_cons", "re_hae_fem", "re_hae_black", "re_hae_log_year",
            "re_hae_base_haem")
        ] |> unlist()
        
        # Run the model
        res <- run_microsim_modelC(
          df_UKPDS_coef = df_UKPDS_coef,
          df_ukpds_pop = df_ukpds_pop,
          cycles = 20
        )
        
        return(res)
      }
    )
  }
}
