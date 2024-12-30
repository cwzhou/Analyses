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
column_names <- unlist(column_names)
# Define a custom sorting function
custom_sort <- function(names) {
  order(
    !grepl("_survival$", names),
    !grepl("_endpoint$", names),
    !grepl("_n_phase2$", names),
    !grepl("^time_", names),
    !grepl("^czmk_", names),
    !grepl("^zom_", names),
    !grepl("^obs_", names),
    names
  )
}
sorted_column_names <- column_names[custom_sort(column_names)]

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
    
    # Calculate the mean across people
    mean_value <- mean(mff_tau_df$Number_RE, na.rm = TRUE)
    
    # Return both the dataset and the mean value
    return(list(mff_tau_df = mff_tau_df, mean_value = mean_value))
  }
}


# need to source CHECKPackagesScripts.R where we want to check_RE = 1
init_seed = 2025
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
                ctype=ctype,
                cparam=censor_rate,
                censor_min=censor_min,
                censor_max=censor_max,
                gaptype=gaptype,gapparam1=gapparam1,gapparam2=gapparam2, #gapparams are the rhos
                ncov = ncov,
                tau0=tau0,
                predHazardFn_D = predHazardFn_D,
                predHazardFn_R = predHazardFn_R,
                predPropensityFn = predPropensityFn # list of predictor functions
)
arg.obs.train <- arg_list
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
sim = 2
sim = 10 # testing malloc

message("starting run_simulation for sim#", sim)
  
  cat("\n\n#################################\n")
  cat("######### Simulation ",sim, "#########\n")
  cat("      Endpoint: ", endpoint, "      \n")
  cat("#################################\n")
  
  train_seed = sim*10000 + init_seed
  test_seed = train_seed + 10
  
  cat ("%%% Training Data for", sim_data_type, "Simulation:",sim,"%%%\n")
  tt(1)
  
  message("using train_seed (", train_seed, ") to generate training data")
  set.seed(train_seed)
  # if (sim_data_type == "RE"){
  message("Recurrent Events Survival Data Simulation")
  sim.train = do.call(gdata_RE, arg.obs.train)
  df_recurr = sim.train$dataset_recurrent; #View(df_recurr)
  df_surv = sim.train$dataset_survival; #View(df_surv)
  name = sprintf("%s_%s",sim.train$name, sim_data_type); print(name)
  # update this to include name_surv
  assign(name, df_recurr)
  # head_recurr = head(df_recurr); head_recurr
  # kable(head_recurr, format = "latex", caption = "Recurrent Events Dataset Example")
  
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
  
  # # # obs policy value (testing)
  message("using test_seed to generate obs testing data")
  set.seed(test_seed)
  obs.data.rep <- do.call(gdata_RE, arg.obs.no.censor) # no censoring for eval sets
  rep_obs <<- obs.data.rep
  obs.test.df_recurr = obs.data.rep$dataset_recurrent;
  obs.test.df_surv = obs.data.rep$dataset_survival;
  # result["obs_survival"] = survival_val.fn(obs.test.df_surv)
  result[sim, "obs_survival"] = survival_val.fn(obs.test.df_surv)
  obs.mff.result <- endpoint_val.fn(data = obs.test.df_recurr, idName0, epName0, txName0)
  obs.mff.df <- cbind(simulation = sim, obs.mff.result$mff_tau_df, method = "observed")
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
