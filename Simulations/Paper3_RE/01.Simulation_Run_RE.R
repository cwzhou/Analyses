# Title: Simulating Survival Data with Recurrent Events and Terminal Death
# Description: [Brief description of the purpose and objectives of the simulation]
# Author: Christina Zhou
# Date: 07.02.2023

# setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Github/Paper_1/Paper1_Sims/GenerateData/1_code/0_simulations/RE")
setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper3_RE")

start_time = Sys.time()

# Specify the methods and skip.methods
all_methods <- c("czmk", "zom", "obs");
skip_method <- !c(TRUE, TRUE, TRUE);
n.methods <- length(all_methods)
column_names = NULL

# Generate column names based on methods for result
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

survival_val.fn <- function(data){
  # mean truncated survival time (terminal event) or
  # mean overall survival probability at t0 = crit.eval.os
  mean(data$St_eval, na.rm = TRUE)
}
endpoint_val.fn <- function(data) {
  # mean frequency function at tau
  mean(data$MFF_eval, na.rm = TRUE)
}

set.seed(111)

# Libraries --------------------------------------------------------------
source("02.Simulation_Libraries_RE.R")

# Functions -----------------------------------------------------------------
source("02.Simulation_Functions_RE.R")

# Hyperparameters  -----------------------------------------------------------------
source("02.Simulation_Parameters_RE.R")

# need to source CHECKPackagesScripts.R where we want to check_RE = 1
train_seed = 2025
ncov = 2
nodesize = 3
mindeath = round(sqrt(c(nodesize)), 0)
criterion_phase1 = "mean"
criterion_phase2 = "mean"
splittest1 = "mean_surv"
splittest2 = "gray_re"
endpoint = "RE"
arg.czmk = list(endPoint = endpoint,
                 idName = "ID",
                 epName = "IndR",
                 txName = "Trt",
                 criticalValue1 = criterion_phase1,
                 criticalValue2 = criterion_phase2,
                 evalTime = 1,
                 splitRule1 = splittest1,
                 splitRule2 = splittest2,
                 ERT = FALSE,
                 uniformSplit = TRUE,
                 replace = FALSE,
                 randomSplit = 0.2,
                 nTree = 300,
                 pooled = FALSE,
                 tol1 = c(0.1,0.1),
                 stratifiedSplit = 0.1)
arg.obs <- arg.obs.no.censor <- arg.zom <- arg.czmk
arg.obs.no.censor$ctype <- 99

sim = 1
# Define the function that runs a single simulation
# run_simulation <- function(sim){
  message("starting run_simulation")
  ### simulation
  r00 = data.frame(sim.no = sim) #data.frame(sim = n.sim)
  r00[column_names] <- NA
  
  # Create the data frame with the 'sim' column
  result <- data.frame(sim.no = sim) #data.frame(sim = 1:n.sim)
  # Add columns to the result data frame
  result[sorted_column_names] <- NA
  attr(result, "criterion_phase1") <- list(criterion = criterion_phase1) #, crit.value = crit.value_phase1)
  attr(result, "criterion_phase2") <- list(criterion = criterion_phase2) #, crit.value = crit.value_phase2)
  
  # trt_result <- data.frame(sim.no = rep(sim, each = n.methods), #rep(1:n.sim, each = n.methods),
  #                          method = NA,
  #                          surv_A = NA, surv_B = NA,
  #                          endpoint_A = NA, endpoint_B = NA) # each = 3 means repeating b/c theres 3 methods rn (obs, czmk, zom)
  
  ######################################################################
  ######################################################################
  ######################################################################
  ######################################################################
  ######################################################################
  
  # flow: obs (1) -> optimal (n.mc), obs.no.censor (n.mc)
  print(Sys.time())
  # for (sim in 1:n.sim) {
  
  cat("\n\n#################################")
  cat("######### Simulation ",sim, "#########")
  cat("      Endpoint: ", endpoint, "      ")
  cat("#################################\n")

  
  train_seed = sim*10000 + 2025
  
  cat ("1. Data for Simulation",sim,":", sim_data_type, "\n")
  tt(1)
  
  # obs.data
  message("using train_seed for generating training data")
  set.seed(train_seed)
  
  # Main Simulation -----------------------------------------------------------
  if (sim_data_type == "RE"){
    message("Recurrent Events Survival Data Simulation")
    sim = gdata_RE(N=N,G=G,
                   lambda_0D=lambda_0D,lambda_0R=lambda_0R,beta_R=beta_R,beta_D=beta_D,
                   omega_D=omega_D,omega_R=omega_R,
                   ztype=ztype,zparam=zparam,ctype=ctype,cparam=cparam,
                   gaptype=gaptype,gapparam1=gapparam1,gapparam2=gapparam2,
                   num_A=num_A,tau=tau,seed1=train_seed)
    df_recurr = sim$dataset_recurrent; #View(df_recurr)
    df_surv = sim$dataset_survival; #View(df_surv)
    name = sprintf("%s_%s",sim$name, sim_data_type); print(name)
    # update this to include name_surv
    assign(name, df_recurr)
    # head_recurr = head(df_recurr); head_recurr
    # kable(head_recurr, format = "latex", caption = "Recurrent Events Dataset Example")
    
    # Output --------------------------------------------------------------------
    # Save or export simulation results, plots, or any other relevant output
    # Get today's date in the desired format (e.g., YYYY-MM-DD)
    today_date <- format(Sys.Date(), "%Y-%m-%d")
    
    # Define folder path
    folder_path <- sprintf("./2_pipeline/%s", today_date)
    
    # Create the folder if it doesn't exist
    if (!dir.exists(folder_path)) {
      dir.create(folder_path, recursive = TRUE)
    }
    
    # Write CSV files
    write_csv(df_recurr, file = sprintf("%s/01_recurr_%s.csv", folder_path, name))
    write_csv(df_surv, file = sprintf("%s/01_surv_%s.csv", folder_path, name))
  }
  
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
  
  model1 = "Surv(R_closed, IndD) ~ Cov" %>% as.formula()
  model2 = "Surv(L_open, R_closed, IndR) ~ Cov" %>% as.formula()
  models_RE = list(model1, model2)
  
  a0 = data_to_use %>% filter(Trt == 0); #dim(a0)
  a0_ids = a0 %>% pull(ID) %>% unique(); #length(a0_ids)
  a1 = data_to_use %>% filter(Trt == 1); #dim(a1)
  a1_ids = a1 %>% pull(ID) %>% unique(); #length(a1_ids)
  
  # # # obs policy value
  message("using train_seed+10 to generate obs testing data")
  test_seed = train_seed + 10
  
  test.sim = gdata_RE(N=N,G=G,
                 lambda_0D=lambda_0D,lambda_0R=lambda_0R,beta_R=beta_R,beta_D=beta_D,
                 omega_D=omega_D,omega_R=omega_R,
                 ztype=ztype,zparam=zparam,ctype=99,cparam=cparam,
                 gaptype=gaptype,gapparam1=gapparam1,gapparam2=gapparam2,
                 num_A=num_A,tau=tau,seed1=test_seed)
  test.df_recurr = test.sim$dataset_recurrent; #View(test.df_recurr)
  test.df_surv = test.sim$dataset_survival; #View(test.df_surv)
  test.name = sprintf("Test.%s_%s",test.sim$name, sim_data_type); print(name)
  # update this to include name_surv
  assign(test.name, test.df_recurr)
  # test.head_recurr = head(test.df_recurr); test.head_recurr
  # kable(test.head_recurr, format = "latex", caption = "Test Set Recurrent Events Dataset Example")
  
  # Define folder path
  test.folder_path <- sprintf("./2_pipeline/test/%s", today_date)
  
  # Create the folder if it doesn't exist
  if (!dir.exists(test.folder_path)) {
    dir.create(test.folder_path, recursive = TRUE)
  }
  # Write CSV files
  write_csv(test.df_recurr, file = sprintf("%s/test.01_recurr_%s.csv", test.folder_path, test.name))
  write_csv(test.df_surv, file = sprintf("%s/test.01_surv_%s.csv", test.folder_path, test.name))
  
  # result[sim, "obs_survival"]
  # result[sim, "obs_endpoint"]
  # result["obs_survival"] <- survival_val.fn(data_to_use)
  # result["obs_endpoint"] <- endpoint_val.fn(obs.data.rep)
  
  
  
  
  
  
  
set.seed(train_seed + 1)
optimal.czmk <- do.call(itrSurv::itrSurv,
                        c(arg.czmk,
                          list(data = data_to_use,
                               models = models_RE,
                               timePointsSurvival = timePointsSurvival,
                               timePointsEndpoint = timePointsEndpoint,
                               tau = tau,
                               mTry = sqrt(ncov),
                               minEventEnd = 3L,
                               minEventSurv = 3L,
                               nodeSizeEnd = 6L,
                               nodeSizeSurv = 6L)))
# czmk.error <- class(optimal.czmk)[1] == "try-error"
# arg.czmk$policy <- if (!czmk.error) optimal.czmk
# View(arg.czmk$policy)
# policy_czmk <<- arg.czmk$policy
# #training
# predd_surv_czmk <<- predd_surv
# predd_ep_czmk <<- predd_ep
# 
#set.seed(train_seed + 2)
#optimal.zom <- do.call(itrSurv::itrSurv, c(arg.czmk,
#                                           list(mTry = 1,
#                                                data = data_to_use,
#                                                models = models_RE,
#                                                timePointsSurvival = timePointsSurvival,
#                                                timePointsEndpoint = timePointsEndpoint,
#                                                tau = tau,
#                                                minEventEnd = 1L,
#                                                minEventSurv = 1L,
#                                                nodeSizeEnd = 1e9,
#                                                nodeSizeSurv = 1e9)))
# 
# 
# zom.error <- class(optimal.zom)[1] == "try-error"
# arg.zom$policy <- if (!zom.error) optimal.zom
# policy_zom <- arg.zom$policy
# if (!zom.error) result["zom_n_phase2"] <- mean(arg.zom$policy@phaseResults[["SurvivalPhase1Results"]]@optimal@Ratio_Stopping_Ind == 0) # result[sim, "zom_n_phase2"]
# rm(optimal.zom); gc()
# }
# 



print("end of script")


# End of script -------------------------------------------------------------
