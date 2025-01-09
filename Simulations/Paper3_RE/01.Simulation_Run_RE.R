# Title: Simulating Survival Data with Recurrent Events and Terminal Death
# Description: [Brief description of the purpose and objectives of the simulation]
# Author: Christina Zhou
# Date: 07.02.2023

setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper3_RE")
source("02.Simulation_Parameters_RE.R")
start_time = Sys.time()

# Generate column names based on methods for result
column_names = NULL
column_names <- lapply(all_methods, function(method) {
  if (method == "obs") {
    c(
      paste(method, "survival", sep = "_"),
      paste(method, "endpoint", sep = "_"),
      paste("time", method, sep = ".")
    )
  } else {
    c(
      paste(method, "survival", sep = "_"),
      paste(method, "endpoint", sep = "_"),
      paste(method, "n_phase2", sep = "_"),
      paste("time", method, sep = ".")
    )
  }
})
# Flatten the list of column names
column_names <- c(unlist(column_names),"train_cens")
# Define a custom sorting function
custom_sort <- function(names) {
  order(
    !grepl("_survival$", names),
    !grepl("_endpoint$", names),
    !grepl("_n_phase2$", names),
    !grepl("train_cens", names),
    !grepl("^time_", names),
    !grepl("^czmk_", names),
    !grepl("^zom_", names),
    !grepl("^obs_", names),
    names
  )
}
sorted_column_names <- column_names[custom_sort(c(column_names))]
print(sorted_column_names)
print(filename)
filename.tmp <- gsub("\\.rds", "_tmp.rds", filename)

survival_val.fn <- function(data){
  # mean (across people) truncated survival time (terminal event) or
  # mean survival probability at t0 = crit.eval.os
  mean(data$obs_time, na.rm = TRUE)
}
endpoint_val.fn <- function(data, idName, epName, txName) {
  if (endpoint == "RE") {
    # Calculate the mean frequency function at tau
    mff_tau_df <- data %>%
      group_by(!!sym(idName)) %>%
      summarize(
        Number_RE = sum(!!sym(epName)),
        Trt = mean(!!sym(txName))
      ) %>% as.data.frame()
    print(head(mff_tau_df))

    # Calculate the mean across people
    mean_value <- mean(mff_tau_df$Number_RE, na.rm = TRUE)

    # Return both the dataset and the mean value
    return(list(mff_tau_df = mff_tau_df, mean_value = mean_value))
  }
}

# need to source CHECKPackagesScripts.R where we want to check_RE = 1
init_seed = 1995
idName0 = "ID"
epName0 = "IndR"
txName0 = "Trt"
nodesize = 3
mindeath = round(sqrt(c(nodesize)), 0)
splittest1 = "mean_surv"
splittest2 = "gray_re"
endpoint = "RE"
# rm(arg_list);rm(arg.czmk.test);rm(arg.czmk.train)
arg_list = list(N=n,
                G=G,
                # lambda_0D=lambda_0D,lambda_0R=lambda_0R,beta_R=beta_R,beta_D=beta_D,
                # omega_D=omega_D,omega_R=omega_R,gamma_D=gamma_D,gamma_R=gamma_R,
                ztype = 0, # covariate distribution
                zparam = 0.3, # covariate distribution parameter
                # zseed = 2025, # don't need - delete: jan 1, 2025
                ctype=ctype,
                cparam=censor_rate,
                censor_min=censor_min,
                censor_max=censor_max,
                # gaptype=gaptype, # alpha1/2 depends on treatment, so moved gapparam1/2 to F01.DynamicsRE.R, as of Jan 3, 2025
                # gapparam1=gapparam1,gapparam2=gapparam2, #gapparams are the rhos
                ncov = ncov,
                tau0=tau0,
                predHazardFn_D = predHazardFn_D,
                predHazardFn_R = predHazardFn_R,
                predPropensityFn = predPropensityFn # list of predictor functions
)

arg.obs.train <- c(arg_list)
arg.czmk.train = list(endPoint = endpoint,
                      idName = idName0,
                      epName = epName0,
                      txName = txName0,
                      criticalValue1 = criterion_phase1,
                      criticalValue2 = criterion_phase2,
                      evalTime = crit_t0_eval,
                      splitRule1 = splittest1,
                      splitRule2 = splittest2,
                      ERT = FALSE,
                      uniformSplit = TRUE,
                      replace = FALSE,
                      randomSplit = 0.2,
                      nTree = 300,
                      pooled = FALSE,
                      tol1 = tol1,
                      stratifiedSplit = 0.1)
arg.czmk.test = c(arg_list,
                  evaluate = TRUE,
                  crit.eval.surv = criterion_phase1,
                  crit.eval.endpoint = criterion_phase2)
# update args that alr exist in arg.czmk.test
arg.czmk.test$N = n.eval
arg.czmk.test$ctype = 99
arg.obs.no.censor <- arg.zom.test <- arg.czmk.test
# sim = 1

# Initialize an empty list to collect datasets for all simulations
all_sims_data.mff <- list()

# Define the function that runs a single simulation
# run_simulation <- function(sim){

  # Create the data frame with the 'sim' column
  # result <- data.frame(sim.no = sim)
  result <- data.frame(sim = 1:n.sim)
  # Add columns to the result data frame
  result[sorted_column_names] <- NA
  attr(result, "criterion_phase1") <- list(criterion = criterion_phase1) #, crit.value = crit.value_phase1)
  attr(result, "criterion_phase2") <- list(criterion = criterion_phase2) #, crit.value = crit.value_phase2)

  # flow: obs (1) -> optimal (n.mc), obs.no.censor (n.mc)
  print(Sys.time())

  ######################################################################
  ######################################################################
  ######################################################################
  ######################################################################
  ######################################################################
  ### simulation
  for (sim in 1:n.sim){ # for-loop for sims (if NOT parallelization)
  # message("starting run_simulation for sim#", sim)

  cat("\n\n##################################################################\n")
  cat("##################################################################\n")
  cat("################## Simulation ",sim, "##################\n")
  cat("#########      Endpoint: ", endpoint, "      #########\n")
  cat("##################################################################\n")
  cat("##################################################################\n")
  train_seed = sim*10000 + init_seed*3
  test_seed = train_seed + 30306

  set.seed(init_seed*sim+505)
  u1_train = runif(n)
  u1_test = runif(n.eval)
  u2_train = runif(n)
  u2_test = runif(n.eval)
  u3_train = runif(n)
  u3_test = runif(n.eval)
  arg.czmk.test$u1 <- arg.zom.test$u1 <- arg.obs.no.censor$u1 <- u1_test
  arg.czmk.test$u2 <- arg.zom.test$u2 <- arg.obs.no.censor$u2 <- u2_test
  arg.czmk.test$u3 <- arg.zom.test$u3 <- arg.obs.no.censor$u3 <- u3_test
  arg.obs.train$u1 = u1_train
  arg.obs.train$u2 = u2_train
  arg.obs.train$u3 = u3_train

  # arg.obs.train$zseed = round(arg.obs.train$zseed*sim)
  # arg.czmk.test$zseed = round(test_seed*sim+1)
  # arg.zom.test$zseed = round(test_seed*sim+2)
  # arg.obs.no.censor$zseed = round(test_seed*sim+3)

  cat ("%%% Training Data for", sim_data_type, "Simulation:",sim,"%%%\n")
  tt(1)

  # if (sim_data_type == "RE"){
    message("Recurrent Events Survival Data Simulation")
    message("using train_seed (", train_seed, ") to generate training data")
    set.seed(train_seed)
    sim.train = do.call(gdata_RE, arg.obs.train)
    obs.times_train <<- times_act
    ph1_obs_train <<-pred.hazard1
    gap1_obs_train <<- gaptime1
    tt_obs_train <<- tt_fail
    df_recurr = sim.train$dataset_recurrent; #View(df_recurr)
    df_surv = sim.train$dataset_survival; head(df_surv)
    name = sprintf("%s_%s",sim.train$name, sim_data_type); #print(name)
    # update this to include name_surv
    assign(name, df_recurr)
    assign(sprintf("df_sim%s", sim), df_surv)
    # head_recurr = head(df_recurr); head_recurr
    # kable(head_recurr, format = "latex", caption = "Recurrent Events Dataset Example")

    # Calculate the number of censored observations
    n_censored <- sum(df_surv$indD == 0)  # Count of censored observations
    # Calculate the total number of observations
    n_total <- nrow(df_surv)  # Total number of rows in the dataset
    # Calculate the percentage of censored data
    training_censored <- n_censored / n_total
    cat("Percentage of censored data:", training_censored*100, "%\n")
    result[sim, "train_cens"] = round(training_censored,3)

    # Output --------------------------------------------------------------------
    if (savingrds == TRUE){
      folder_path <- sprintf("./2_pipeline/%s", date_folder)
      # Create the folder if it doesn't exist
      if (!dir.exists(folder_path)) {
        dir.create(folder_path, recursive = TRUE)
      }
      write_csv(df_recurr, file = sprintf("%s/01_recurr_%s.csv", folder_path, name))
      write_csv(df_surv, file = sprintf("%s/01_surv_%s.csv", folder_path, name))
    }
    # }

  data_to_use = df_recurr
  tau = max(data_to_use$R_closed)

  timePointsSurvival = data_to_use %>%
    filter(IndD == 1) %>%
    dplyr::select(R_closed) %>%
    distinct() %>%
    unlist(use.names = F); length(timePointsSurvival)

  timePointsEndpoint = data_to_use %>%
    filter(IndR == 1) %>%
    dplyr::select(R_closed) %>%
    distinct() %>%
    unlist(use.names = F); length(timePointsEndpoint)

  model1 = paste0("Surv(R_closed, IndD) ~ ",
                  paste(paste0("Z", 1:ncov, ""), collapse = " + ")) %>%
    as.formula()
  model2 = paste0("Surv(L_open, R_closed, IndR) ~ ",
                  paste(paste0("Z", 1:ncov, ""), collapse = " + ")) %>%
    as.formula()
  models_RE = list(model1, model2)

  a0 = data_to_use %>% filter(Trt == 0); #dim(a0)
  a0_ids = a0 %>% pull(ID) %>% unique(); #length(a0_ids)
  a1 = data_to_use %>% filter(Trt == 1); #dim(a1)
  a1_ids = a1 %>% pull(ID) %>% unique(); #length(a1_ids)

  mff_to_remove <- c("zom.mff.df", "obs.mff.df", "czmk.mff.df")
  for (obj in mff_to_remove) {
    if (exists(obj, envir = .GlobalEnv)) {
      rm(list = obj, envir = .GlobalEnv)
    }
  }

  # # # obs policy value (testing)
  message("using test_seed to generate obs testing data")
  set.seed(test_seed)
  obs.data.rep <- do.call(gdata_RE, arg.obs.no.censor) # no censoring for eval sets
  # obs.times_test <<- times_act
  # ph1_obs_test <<-pred.hazard1
  # gap1_obs_test <<- gaptime1
  # tt_obs_test <<- tt_fail
  # rep_obs <<- obs.data.rep
  # Define your variables in a list for cleaner handling
  variables_to_assign <- list(
    times_test = times_act,
    ph1_test = pred.hazard1,
    ph2_test = pred.hazard2,
    gap1_test = gaptime1,
    tt_test = tt_fail
  )
  # Assign each variable using sprintf
  for (name in names(variables_to_assign)) {
    assign(sprintf("obs_sim%s_%s", sim, name), variables_to_assign[[name]], envir = .GlobalEnv)
  }
  assign(sprintf("rep_obs_sim%s", sim), obs.data.rep, envir = .GlobalEnv)

  obs.test.df_recurr = obs.data.rep$dataset_recurrent;
  obs.test.df_surv = obs.data.rep$dataset_survival;
  # result["obs_survival"] = survival_val.fn(obs.test.df_surv)
  result[sim, "obs_survival"] = survival_val.fn(obs.test.df_surv)
  obs.mff.result <- endpoint_val.fn(data = obs.test.df_recurr, idName0, epName0, txName0)
  obs.mff.df <- cbind(simulation = sim, obs.mff.result$mff_tau_df,
                      survival = obs.test.df_surv$obs_time,
                      method = "observed")
  # result["obs_endpoint"] = obs.mff.result$mean_value
  result[sim, "obs_endpoint"] = obs.mff.result$mean_value
  # View(obs.mff.df); View(result)

    # test.name = sprintf("Test.%s_%s",test.sim$name, sim_data_type); print(name)
  # # update this to include name_surv
  # assign(test.name, test.df_recurr)
  # # test.head_recurr = head(test.df_recurr); test.head_recurr
  # # kable(test.head_recurr, format = "latex", caption = "Test Set Recurrent Events Dataset Example")
  #
  # # Define folder path
  # test.folder_path <- sprintf("./2_pipeline/test/%s", date_folder)
  #
  # # Create the folder if it doesn't exist
  # if (!dir.exists(test.folder_path)) {
  #   dir.create(test.folder_path, recursive = TRUE)
  # }
  # # Write CSV files
  # write_csv(test.df_recurr, file = sprintf("%s/test.01_recurr_%s.csv", test.folder_path, test.name))
  # write_csv(test.df_surv, file = sprintf("%s/test.01_surv_%s.csv", test.folder_path, test.name))
  # result["time.obs"] <- tt(2, reset = TRUE, unit = "min")["elapsed"]
  result[sim, "time.obs"] <- tt(2, reset = TRUE, unit = "min")["elapsed"]
  rm(obs.data.rep); gc()

  # czmk
  cat("\n******************************\n")
  # estimation
  cat ("1. czmk for simulation", sim, "\n")
  if (!skip.czmk) {
    cat ("  1. czmk - Policy estimation for RE Simulation",sim,"\n")
    # new package itrSurv
    set.seed(train_seed + 1)
    optimal.czmk <- do.call(itrSurv::itrSurv,
                            c(arg.czmk.train,
                              list(data = data_to_use,
                                   models = models_RE,
                                   timePointsSurvival = timePointsSurvival,
                                   timePointsEndpoint = timePointsEndpoint,
                                   tau = tau0,
                                   mTry = sqrt(ncov),
                                   minEventEnd = 3L,
                                   minEventSurv = 3L,
                                   nodeSizeEnd = 6L,
                                   nodeSizeSurv = 6L)))
    # czmk.times_train <<- times_act
    # ph1_czmk_train <<-pred.hazard1
    # gap1_czmk_train <<- gaptime1
    # tt_czmk_train <<- tt_fail
    czmk.error <- class(optimal.czmk)[1] == "try-error"
    arg.czmk.test$policy <- if (!czmk.error) optimal.czmk
    policy_czmk <<- arg.czmk.test$policy
    #training
    predd_surv_train_czmk <<- predd_surv
    predd_ep_train_czmk <<- predd_ep
    # if (!czmk.error) result["czmk_n_phase2"] <- mean(arg.czmk.test$policy@phaseResults[["SurvivalPhase1Results"]]@optimal@Ratio_Stopping_Ind == 0) #result[sim, "czmk_n_phase2"]
    if (!czmk.error) result[sim,"czmk_n_phase2"] <- mean(arg.czmk.test$policy@phaseResults[["SurvivalPhase1Results"]]@optimal@Ratio_Stopping_Ind == 0)
    rm(optimal.czmk); gc()

    cat ("  \n 1. czmk - Evaluation for RE Simulation",sim,"\n")
    if (!czmk.error) {
      set.seed(test_seed)
      # TO DO: MAKE SURE THE COVARIATES ARE THE SAME FOR CZMK AND ZOM TEST SET
      czmk.data.rep <- do.call(gdata_RE, arg.czmk.test); #head(czmk.data.rep$dataset_survival$Z1)
      czmk.test.df_recurr = czmk.data.rep$dataset_recurrent; #View(czmk.test.df_recurr)
      czmk.test.df_surv = czmk.data.rep$dataset_survival; #View(czmk.test.df_surv)
      # czmk.times_test <<- times_act
      # ph1_czmk_test <<- pred.hazard1
      # gap1_czmk_test <<- gaptime1
      # tt_czmk_test <<- tt
      # rep_czmk <<- czmk.data.rep

      # Define your variables in a list for cleaner handling
      variables_to_assign <- list(
        times_test = times_act,
        ph1_test = pred.hazard1,
        ph2_test = pred.hazard2,
        gap1_test = gaptime1,
        tt_test = tt_fail,
        predd_surv_eval = predd_surv,
        predd_ep_eval = predd_ep,
        test.df_recurr = czmk.test.df_recurr,
        test.df_surv = czmk.test.df_surv
      )
      # Assign each variable using sprintf
      for (name in names(variables_to_assign)) {
        assign(sprintf("czmk_sim%s_%s", sim, name), variables_to_assign[[name]], envir = .GlobalEnv)
      }
      assign(sprintf("rep_czmk_sim%s", sim), czmk.data.rep, envir = .GlobalEnv)

      # result["czmk_survival"] = survival_val.fn(czmk.test.df_surv)
      result[sim,"czmk_survival"] = survival_val.fn(czmk.test.df_surv)
      # czmk.test.df_recurr %>% group_by(ID) %>% summarize(Number_RE = sum(IndR), Trt = mean(Trt))
      czmk.mff.result <<- endpoint_val.fn(data = czmk.test.df_recurr, idName0, epName0, txName0)
      czmk.mff.df <<- cbind(simulation = sim, czmk.mff.result$mff_tau_df,
                            survival = czmk.test.df_surv$obs_time,
                            method = "czmk")
      # View(czmk.mff.df); View(result)
      # result["czmk_endpoint"] = czmk.mff.result$mean_value #result[sim, "czmk_endpoint"]
      result[sim, "czmk_endpoint"] = czmk.mff.result$mean_value
    } # end of if czmk.error
    # result["time.czmk"] <- tt(2, reset = TRUE, units = "mins")["elapsed"] #result[sim, "time.czmk"]
    result[sim, "time.czmk"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]
    arg.czmk.test$policy <- NULL; gc()
    rm(czmk.data.rep); gc()
  } # if skip.czmk

  cat("\n******************************\n")
  # estimation
  cat ("2. Estimation - zero-order model for RE Simulation",sim,"\n")
  if (!skip.zom) {
    cat ("  2. zom - Policy estimation for RE Simulation",sim,"\n")
    set.seed(train_seed + 2)
    optimal.zom <- do.call(itrSurv::itrSurv, c(arg.czmk.train,
                                               list(data = data_to_use,
                                                    models = models_RE,
                                                    timePointsSurvival = timePointsSurvival,
                                                    timePointsEndpoint = timePointsEndpoint,
                                                    tau = tau0,
                                                    mTry = 1,
                                                    minEventEnd = 1L,
                                                    minEventSurv = 1L,
                                                    nodeSizeEnd = 1e+9,
                                                    nodeSizeSurv = 1e+9)))
    # zom.times_train <<- times_act
    # ph1_zom_train <<-pred.hazard1
    # gap1_zom_train <<- gaptime1
    # tt_zom_train <<- tt
    zom.error <- class(optimal.zom)[1] == "try-error"
    arg.zom.test$policy <- if (!zom.error) optimal.zom
    policy_zom <- arg.zom.test$policy
    # if (!zom.error) result["zom_n_phase2"] <- mean(arg.zom.test$policy@phaseResults[["SurvivalPhase1Results"]]@optimal@Ratio_Stopping_Ind == 0) # result[sim, "zom_n_phase2"]
    if (!zom.error) result[sim, "zom_n_phase2"] <- mean(arg.zom.test$policy@phaseResults[["SurvivalPhase1Results"]]@optimal@Ratio_Stopping_Ind == 0)
    rm(optimal.zom); gc()
    if (!zom.error){
      cat ("  6. zero-order model - Evaluation for RE Simulation",sim,"\n")
      set.seed(test_seed)
      zom.data.rep <- do.call(gdata_RE, arg.zom.test); head(zom.data.rep$dataset_survival$Z1)
      # zom.times_test <<- times_act
      # ph1_zom_test <<-pred.hazard1
      # ph2_zom_test <<-pred.hazard2
      # gap1_zom_test <<- gaptime1
      # tt_zom_test <<- tt
      # Define your variables in a list for cleaner handling
      variables_to_assign <- list(
        times_test = times_act,
        ph1_test = pred.hazard1,
        ph2_test = pred.hazard2,
        gap1_test = gaptime1,
        tt_test = tt_fail
      )
      # Assign each variable using sprintf
      for (name in names(variables_to_assign)) {
        assign(sprintf("zom_sim%s_%s", sim, name), variables_to_assign[[name]], envir = .GlobalEnv)
      }
      assign(sprintf("rep_zom_sim%s", sim), zom.data.rep, envir = .GlobalEnv)
      zom.test.df_recurr = zom.data.rep$dataset_recurrent; #View(czmk.test.df_recurr)
      zom.test.df_surv = zom.data.rep$dataset_survival; #View(czmk.test.df_surv)
      # result["zom_survival"] = survival_val.fn(zom.test.df_surv) #result[sim, "zom_survival"]
      result[sim, "zom_survival"] = survival_val.fn(zom.test.df_surv)
      zom.mff.result <- endpoint_val.fn(data = zom.test.df_recurr, idName0, epName0, txName0)
      zom.mff.df <- cbind(simulation = sim, zom.mff.result$mff_tau_df,
                          survival = zom.test.df_surv$obs_time,
                          method = "zom")
      # result["zom_endpoint"] = zom.mff.result$mean_value #result[sim, "zom_endpoint"]
      result[sim, "zom_endpoint"] = zom.mff.result$mean_value
    }
    # result["time.zom"] <- tt(2, reset = TRUE, units = "mins")["elapsed"] #result[sim, "time.zom"]
    result[sim, "time.zom"] <- tt(2, reset = TRUE, units = "mins")["elapsed"]
    arg.zom.test$policy <- NULL; gc()
    rm(zom.data.rep); gc()
    } # if !skip.zom

  # Initialize an empty list to collect datasets
  mff.dataset_list <- list()
  # Check if each dataset exists and add it to the list
  if (exists("czmk.mff.df")) mff.dataset_list[["czmk"]] <- czmk.mff.df
  if (exists("zom.mff.df")) mff.dataset_list[["zom"]] <- zom.mff.df
  if (exists("obs.mff.df")) mff.dataset_list[["obs"]] <- obs.mff.df
  # Stack the datasets
  mff.stacked_data <- do.call(rbind, mff.dataset_list)
  row.names(mff.stacked_data) <- NULL
  # Optionally add a column to indicate the simulation number
  mff.stacked_data$simulation <- sim
  # Append to the list of all simulations' data
  all_sims_data.mff[[sim]] <- mff.stacked_data

  ### saving and cleaning
  if (savingrds == TRUE){
    message('saving rds tmp')
    saveRDS(result, filename.tmp) # saving the temporary results
    if (!dir.exists(paste0(dir_rds,"/mff"))) {
      dir.create(paste0(dir_rds,"/mff"), recursive = TRUE)
    }
    filename_stacked.tmp = paste0(dir_rds,"/mff/stacked.mff_sim", sim, ".rds")
    # saveRDS(mff.stacked_data, filename_stacked.tmp)
    saveRDS(all_sims_data.mff, filename_stacked.tmp)
    message('saved rds stacked')
  }
  gc()
  cat("---------------------------------------------------\n")
  cat("------------ End of Simulation", sim, "------------\n")
  print(result[sim,])
  View(result)
  cat("---------------------------------------------------\n")
  } # for (sim in 1:n.sim) but removed for parallelizing
  # return(result) # this is for parallel
  # } # this is for run_simulation for parallel

# Combine all simulations into one big dataset for MFF
mff_allsims0 <- do.call(rbind, all_sims_data.mff)
row.names(mff_allsims0) <- NULL
mff_allsims <- mff_allsims0 %>% mutate(RE_lived = Number_RE/survival)
# View(mff_allsims)
if (savingrds == TRUE){
  write.csv(mff_allsims, paste0(dir_rds,"/mff/mff_allsims.csv"), row.names = FALSE)
}

# mff_allsims %>%
#   group_by(simulation, Number_RE, method) %>%
#   summarize(count = n()) %>%
#   View()
num_re_prop = mff_allsims %>%
  group_by(simulation, Number_RE, method) %>%
  summarize(count = n(), .groups = "drop") %>%  # Step 1: Calculate `count`
  group_by(Number_RE, method) %>%
  summarize(mean_count = count/n.eval)  # Step 2: Calculate mean
if (savingrds == TRUE){
  write.csv(num_re_prop, paste0(dir_rds,"/mff/num_re.csv"), row.names = FALSE)
}



# # Number of cores to use for parallelization
# num_cores <- 5
# results_list = list()
# # Run the simulations in parallel
# results_list <- mclapply(2, run_simulation, mc.cores = num_cores)


message("\n END OF SIMS \n")
# result
end_time = Sys.time()

message("End of Script: 01.Simulation_Run_RE.R")
sprintf("Overall Time Took: %s", round(end_time - start_time,2))

View(result)

print("end of script")
# End of script -------------------------------------------------------------
#
#
mff_allsims %>%
  group_by(method) %>%
  summarise(
    # Summary for Number_RE
    mff_mean = mean(Number_RE, na.rm = TRUE),
    mff_sd = sd(Number_RE, na.rm = TRUE),
    mff_min = min(Number_RE, na.rm = TRUE),
    mff_25 = quantile(Number_RE, 0.25, na.rm = TRUE),  # 25th percentile
    mff_median = median(Number_RE, na.rm = TRUE),  # Median
    mff_75 = quantile(Number_RE, 0.75, na.rm = TRUE),  # 75th percentile
    mff_max = max(Number_RE, na.rm = TRUE))
mff_allsims %>%
  group_by(method) %>%
  summarise(
    # Summary for survival
    surv_mean = mean(survival, na.rm = TRUE),
    surv_sd = sd(survival, na.rm = TRUE),
    surv_min = min(survival, na.rm = TRUE),
    surv_25 = quantile(survival, 0.25, na.rm = TRUE),  # 25th percentile
    surv_median = median(survival, na.rm = TRUE),  # Median
    surv_75 = quantile(survival, 0.75, na.rm = TRUE),  # 75th percentile
    surv_max = max(survival, na.rm = TRUE)
  )

# mff_allsims %>%
#   filter(simulation == 20) %>%
#   group_by(method) %>%
#   summarize(meanre = mean(Number_RE),
#             meansurv = mean(survival),
#             percent_trt1 = mean(Trt))

# library(reshape2)
# # Create summary heatmap
# heatmap_data <- mff.stacked_data %>%
#   group_by(method, ID) %>%
#   summarize(mean_TotalEvents = mean(Number_RE))
#
# ggplot(heatmap_data,
#        aes(x = method, y = as.factor(ID), fill = mean_TotalEvents)) +
#   geom_tile() +
#   # facet_wrap(~method) +
#   theme_minimal() +
#   labs(
#     title = "Heatmap of Total Recurrent Events (averaged across simulations) per Treatment and Person",
#     x = "Recommended Treatment",
#     y = "Person",
#     fill = "Total Events"
#   ) +
#   scale_fill_gradient(
#     name = "Total Events",
#     low = "white",   # Light color for low values
#     high = "blue"   # Dark color for high values
#   )

mff_allsims %>% 
  group_by(method, Trt) %>% 
  summarise(per_trt = n()/(n.eval*n.sim), 
            mean_surv = mean(survival), 
            mean_RE = mean(Number_RE), 
            mean_RE_lived = mean(RE_lived))

