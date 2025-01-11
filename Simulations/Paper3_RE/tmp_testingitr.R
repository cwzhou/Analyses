setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper3_RE")
# Libraries --------------------------------------------------------------
source("02.Simulation_Libraries_RE.R")
# Functions -----------------------------------------------------------------
source("02.Simulation_Functions_RE.R")
local = 1
savingrds = FALSE
date_folder = "2025-02-01"
n.eval = 1000
sim = 10
n.sim = 1000
sim_data_type = "RE"
endpoint = sim_data_type
tau0 = 10
ntree1 = 300
##### Gap Time Hyperparameters #####
G = 10 #5 # total gap times
# now defined in F01.DynamicsRE.R
# gaptype = 0 # failure vs gap time indicator
# gapparam1 = 0.1 #rho for failure # alpha1/2 depends on treatment, so moved gapparam1/2 to F01.DynamicsRE.R, as of Jan 3, 2025
# gapparam2 = 0.7 #rho for gap times

# Specify the methods and skip.methods
all_methods <- c("czmk", "zom", "obs");
skip_method <- c(!TRUE, !TRUE, !TRUE);
n.methods <- length(all_methods)
# Loop to create logical SKIP objects for each method and assign skip_method
assign_skip_function(all_methods, skip_method)

### 0. Get the setting number.
arg <- commandArgs(trailingOnly = TRUE)
num.args = 14
num.params = 8
if (length(arg) < num.args) {
  arg = c(2, 1,
          rep(1,num.params),
          2, 1, 1, 1) #c(2,rep(1, num.args-1)) # by default, with RE endpoint as default
  warning(sprintf("commandArgs was not provided. Set as c(%s).",
                  toString(arg)))
}
names(arg)[1:num.args] = c("endpoint", "censor",
                           "beta_D","beta_R",
                           "gamma_D","gamma_R",
                           "omega_D","omega_R",
                           "lambda_0D","lambda_0R",
                           "propensity", "size",
                           "crit_surv", "crit_endpoint")
# print("arg:")
# print(arg)

arg1 <- as.numeric(arg[1]) # 1..3 # 1=CR,2=RE,3=MC
arg2 <- as.numeric(arg[2]) # 1,2 #censoring = 20%, censoring = 50%
arg3 <- as.numeric(arg[3]) # 1..4 4 possible beta_D combinations
arg4 <- as.numeric(arg[4]) # 1..4 4 possible beta_R combinations
arg5 <- as.numeric(arg[5]) # 1..4 4 possible gamma_D combinations
arg6 <- as.numeric(arg[6]) # 1..4 4 possible gamma_R combinations
arg7 <- as.numeric(arg[7]) # 1..4 4 possible omega_D combinations
arg8 <- as.numeric(arg[8]) # 1..4 4 possible omega_R combinations
arg9 <- as.numeric(arg[9]) # 1..4 4 possible lambda_0D combinations
arg10 <- as.numeric(arg[10]) # 1..4 4 possible lambda_0R combinations
arg11 <- as.numeric(arg[11]) # 1..2 2 possible propensity combinations
arg12 <- as.numeric(arg[12]) # 1..2 2 possible size combinations
arg13 <- as.numeric(arg[13]) # 1..3 3 possible crit for surv (mean, mean.prob.combo, prob)
arg14 <- as.numeric(arg[14]) # 1..3 3 possible crit for ep (mean, mean.prob.combo, prob)
arg.date <- if (is.na(arg[num.args+1]) | arg[num.args+1] == "") date_folder else as.character(arg[num.args+1])

### 1. setup
# default settings
if (arg1 == 1){
  endpoint = "CR"
} else if (arg1 == 2){
  endpoint = "RE"
} else if (arg1 == 3){
  endpoint = "MC"
} else{
  stop("Endpoint must be 1=CR, 2=RE, or 3=MC")
}
default <- list(n.eval = n.eval,
                n.sim = n.sim,
                tau = tau0,
                endpoint = endpoint)

# arg2 censor
censor <- list(low.censoring = list(
  # we want about 20% censoring
  ctype = 1, # uniform censoring
  censor_min = 0,
  censor_max = tau0/3, # higher is less censoring for unif
  censor_rate = 0  # not used for unif
),
high.censoring = list(
  # we want about 50% censoring
  ctype = 1, # uniform censoring
  censor_min = 0,
  censor_max = 0.7*tau0, # lower is more censoring for unif
  censor_rate = 0 # not used for unif
))


# arg3-10 betas, gammas, omegas, lambda0s
if (endpoint == "RE"){
  # no intercept b/c survival
  # Covariate effects for survival
  betasD <- list(
    betaD1 = list(
      betaD.hazard0 = c(log(2.2), log(1.2), log(1.7), log(2.5), log(1.7)),
      betaD.hazard1 = c(log(1.4), log(2.7), log(1.5), log(1.9), log(1.1))  # Covariate effects for Treatment 0 (5 parameters)
      # Covariate effects for Treatment 1 (5 parameters)
    ),
    betaD2 = list(
      betaD.hazard0 = c(log(1.5), log(1.2), log(0.7), log(1.0), log(1.3), log(1.5), log(1.0), log(1.6), log(1.3), log(0.9)),  # Covariate effects for Treatment 0 (10 parameters)
      betaD.hazard1 = c(log(1.0), log(0.8), log(1.3), log(1.1), log(0.9), log(1.2), log(1.4), log(1.3), log(1.1), log(1.2))   # Covariate effects for Treatment 1 (10 parameters)
    )
  )
  # betasD <- list(
  #   betaD1 = list(
  #     betaD.hazard0 = c(log(1.3), log(0.9), log(1.3), log(1.2), log(1.8)),  # Covariate effects for Treatment 0 (5 parameters)
  #     betaD.hazard1 = c(log(2.1), log(0.2), log(1.1), log(0.95), log(1.15))  # Covariate effects for Treatment 1 (5 parameters)
  #   ),
  #   betaD2 = list(
  #     betaD.hazard0 = c(log(1.5), log(1.2), log(0.7), log(1.0), log(1.3), log(1.5), log(1.0), log(1.6), log(1.3), log(0.9)),  # Covariate effects for Treatment 0 (10 parameters)
  #     betaD.hazard1 = c(log(1.0), log(0.8), log(1.3), log(1.1), log(0.9), log(1.2), log(1.4), log(1.3), log(1.1), log(1.2))   # Covariate effects for Treatment 1 (10 parameters)
  #   )
  # )

  # Covariate effects for recurrence
  # betasR <- list(
  #   betaR1 = list(
  #     betaR.hazard0 = c(log(2.8), log(3.9), log(0.2), log(1.05), log(1.1)), #beta.hazard0 = c(log(0.8), log(3.9), log(0.2), log(0.05), log(1.1)), #c(log(1.8), log(1.9), log(1.2), log(1.5), log(1.1)),  # Covariate effects for Treatment 0 (5 parameters)
  #     betaR.hazard1 = c(log(2.1), log(2.6), log(1.1), log(1.25), log(0.31)) #beta.hazard1 = c(log(3.1), log(5.6), log(2.1), log(1.25), log(0.11)) #c(log(1.1), log(0.6), log(2.1), log(1.25), log(1.5))   # Covariate effects for Treatment 1 (5 parameters)
  #   ),
  #   betaR2 = list(
  #     betaR.hazard0 = c(log(1.5), log(1.8), log(1.3), log(1.4), log(1.2), log(1.0), log(1.4), log(1.6), log(1.1), log(0.9)),  # Covariate effects for Treatment 0 (10 parameters)
  #     betaR.hazard1 = c(log(1.0), log(1.2), log(1.7), log(1.4), log(1.3), log(1.1), log(1.5), log(1.6), log(1.3), log(1.2))   # Covariate effects for Treatment 1 (10 parameters)
  #   )
  # )

  betasR <- list(
    betaR1 = list(
      betaR.hazard0 = c(log(7.4), log(2.8), log(5.3), log(4.4), log(3.2)), #beta.hazard0 = c(log(0.8), log(3.9), log(0.2), log(0.05), log(1.1)), #c(log(1.8), log(1.9), log(1.2), log(1.5), log(1.1)),  # Covariate effects for Treatment 0 (5 parameters)
      betaR.hazard1 = c(log(3.8), log(20.5), log(12.1), log(15.7), log(16)) #beta.hazard1 = c(log(3.1), log(5.6), log(2.1), log(1.25), log(0.11)) #c(log(1.1), log(0.6), log(2.1), log(1.25), log(1.5))   # Covariate effects for Treatment 1 (5 parameters)
    ),
    betaR2 = list(
      betaR.hazard0 = c(log(1.5), log(1.8), log(1.3), log(1.4), log(1.2), log(1.0), log(1.4), log(1.6), log(1.1), log(0.9)),  # Covariate effects for Treatment 0 (10 parameters)
      betaR.hazard1 = c(log(1.0), log(1.2), log(1.7), log(1.4), log(1.3), log(1.1), log(1.5), log(1.6), log(1.3), log(1.2))   # Covariate effects for Treatment 1 (10 parameters)
    )
  )

  # Interaction between treatment and covariates for survival
  # gammaD: Represents the interaction effect between treatment and covariates on the survival hazard.
  gammaD <- list(
    # gammaD1 = list(
    #   gammaD.hazard0 = c(log(0.5), log(2.5), log(1.4), log(1.3), log(0.7)),  # Interaction effects for Treatment 0 (5 parameters)
    #   gammaD.hazard1 = c(log(1.1), log(1.05), log(0.5), log(1.2), log(1.1))   # Interaction effects for Treatment 1 (5 parameters)
    # ),
    gammaD1 = list(
      gammaD.hazard0 = c(log(1), log(1.0), log(1.0), log(1.0), log(1.0)),  # Interaction effects for Treatment 0 (5 parameters)
      gammaD.hazard1 = c(log(0.1), log(1.4), log(0.05), log(0.01), log(2.01))   # Interaction effects for Treatment 1 (5 parameters)
    ),
    gammaD2 = list(
      gammaD.hazard0 = c(log(1.2), log(0.9), log(1.1), log(1.5), log(1.3), log(1.0), log(1.2), log(1.1), log(1.4), log(1.2)),  # Interaction effects for Treatment 0 (10 parameters)
      gammaD.hazard1 = c(log(0.8), log(1.0), log(1.3), log(1.1), log(1.0), log(1.2), log(1.3), log(1.5), log(1.4), log(1.0))   # Interaction effects for Treatment 1 (10 parameters)
    )
  )

  # Interaction between treatment and covariates for recurrence
  # gammaR: Represents the interaction effect between treatment and covariates on the recurrence hazard.
  gammaR <- list(
    gammaR1 = list(
      gammaR.hazard0 = c(log(1), log(1), log(1), log(1), log(1)),  # Interaction for Treatment 0 recurrence (5 parameters)
      gammaR.hazard1 = c(log(12), log(0.1), log(6), log(10), log(3.3))  # Interaction for Treatment 1 recurrence (5 parameters)
    ),
    gammaR2 = list(
      gammaR.hazard0 = c(log(1.3), log(1.2), log(1.6), log(1.1), log(1.4), log(1.0), log(1.3), log(1.5), log(1.2), log(1.1)),  # Interaction for Treatment 0 recurrence (10 parameters)
      gammaR.hazard1 = c(log(1.0), log(1.1), log(1.4), log(1.2), log(1.3), log(1.2), log(1.5), log(1.6), log(1.3), log(1.4))   # Interaction for Treatment 1 recurrence (10 parameters)
    )
  )

  # Parameters for interaction: disease for treatment
  # omegaD: Represents the direct treatment effect on the survival hazard, independent of covariates.
  omegaD <- list(
    omegaD1 = list(
      omegaD.hazard0 = log(1), #log(1.3),  # Treatment 0 survival effect
      omegaD.hazard1 = log(3) #log(3.5)   # Treatment 1 survival effect
    ),
    omegaD2 = list(  # Second list (identical to omegaD1 for now)
      omegaD.hazard0 = log(3),  # Treatment 0 survival effect
      omegaD.hazard1 = log(3)   # Treatment 1 survival effect
    )
  )

  # Parameters for treatment+covariate interaction
  # Parameters for recurrence hazard (treatment effects)
  # omegaR: Represents the direct treatment effect on the recurrence hazard, independent of covariates.
  omegaR <- list(
    omegaR1 = list(
      omegaR.hazard0 = log(1), #log(3.3), #log(2.3), #log(3.3) # Treatment 0 recurrence effect
      omegaR.hazard1 = log(0.1) #log(2.7) #log(1.5) #log(1.3)  # Treatment 1 recurrence effect
    ),
    omegaR2 = list(  # Second list (identical to omegaR1 for now)
      omegaR.hazard0 = log(3),  # Treatment 0 recurrence effect
      omegaR.hazard1 = log(3)   # Treatment 1 recurrence effect
    )
  )

  # Baseline hazards for survival (disease-related, same for both lists)
  lambda0D <- list(
    lambda0D1 = list(
      lambda0D.hazard0 = c(1),  # Consistent baseline hazard
      lambda0D.hazard1 = c(1)   # Consistent baseline hazard
    ),
    lambda0D2 = list(  # Second list (identical to lambda0D1 for now)
      lambda0D.hazard0 = c(1),  # Consistent baseline hazard
      lambda0D.hazard1 = c(1)   # Consistent baseline hazard
    )
  )

  # Baseline hazards for recurrence (same for both lists)
  lambda0R <- list(
    lambda0R1 = list(
      lambda0R.hazard0 = c(1),  # Consistent baseline hazard
      lambda0R.hazard1 = c(1)   # Consistent baseline hazard
    ),
    lambda0R2 = list(  # Second list (identical to lambda0R1 for now)
      lambda0R.hazard0 = c(1),  # Consistent baseline hazard
      lambda0R.hazard1 = c(1)   # Consistent baseline hazard
    )
  )

  ncovD.list <- lapply(betasD, function(x) length(x$betaD.hazard1) )
  ncovR.list <- lapply(betasR, function(x) length(x$betaR.hazard1) )
  if (ncovD.list$betaD1 == ncovR.list$betaR1){
    ncov.list <- ncovD.list
  } else{
    stop("02.Simulation_Parameters_RE.R: number of covariate for betas D and R don't align.")
  }
}

# arg5 propensity
propensity <-   # (int), covariate (1~5)
  list(obs = list(beta.propensity = function(p) c(-rep(0.5, p))), #no unmeasured confounder
       rct  = list(beta.propensity = function(p) c(rep(0, p))))  # RCT

# arg6 size
size <- list(small.sample.size = list(n = 300),
             large.sample.size = list(n = 1000))

# arg7 and 8 crit
criterion_phase1 = "mean"
criterion_phase2 = "mean"
crit_t0_eval = 1
mean_tol1 = c(0.1,0) # this is for differences in years so we don't want it to be too big
crit_surv <- list(crit1 = list(criterion_phase1 = criterion_phase1,
                               crit.value_phase1 = crit_t0_eval,
                               value_phase1 = "truncated survival mean E[T]",
                               tol1 = mean_tol1)
)
if (endpoint == "RE"){
  crit3 = list(criterion_phase2 = criterion_phase2,
               crit.value_phase2 = crit_t0_eval,
               value_phase2 = "mean frequency function at time tau"
  )
  crit_endpoint <- list(crit1 = crit3#,
                        # crit2 = crit4
  )
}

setting = c(all_methods = all_methods,
            n.methods = n.methods,
            arg = list(arg),
            default,
            censor[[arg2]],
            # beta
            betasD[[arg3]], ncov = ncov.list[[arg3]],
            betasR[[arg4]],
            # gamma
            gammaD[[arg5]],
            gammaR[[arg6]],
            # omega
            omegaD[[arg7]],
            omegaR[[arg8]],
            # lambda0 (baseline)
            lambda0D[[arg9]],
            lambda0R[[arg10]],
            # rest of args
            propensity[[arg11]],
            size[[arg12]],
            crit_surv[[arg13]],
            crit_endpoint[[arg14]]
)

### 2. put a selected setting into the global environment
cat("setting (endpoint, censor, beta_D, beta_R,
    gamma_D, gamma_R, omega_D, omega_R, lambda_0D, lambda_0R,
    propensity, size, crit_phase1, crit_phase2) ", arg, "\n")
list2env(setting, envir = globalenv())
beta.propensity <- beta.propensity(ncov)

# complete the scores
predPropensityFn <- function(covariate) {
  # need to make sure this is right for RE setting..
  # message("predPropensityFn")
  # removed intercept - 12/5/24
  plogis(cbind(covariate) %*% beta.propensity)
}

# HAZARD FUNCTIONS
if (endpoint == "RE"){
  # lambda_D = lambda_0D*exp(t(beta_D)%*%z + omega_D * A + A * t(gamma_D) %*% z) #no intercept for z b/c survival
  # lambda_R = lambda_0R*exp(t(beta_R)%*%z + omega_R * A + A * t(gamma_R) %*% z)
  predHazardFn_D <- function(action, covariate) {
    ifelse(action == 1,
           lambda0D.hazard1 * exp(cbind(covariate) %*% betaD.hazard1 + omegaD.hazard1 * action + action * cbind(covariate) %*% gammaD.hazard1),
           lambda0D.hazard0 * exp(cbind(covariate) %*% betaD.hazard0 + omegaD.hazard0 * action + action * cbind(covariate) %*% gammaD.hazard0)
           # lambda0D.hazard1 * cbind(covariate) %*% exp(beta.hazard1) + omegaD.hazard1 * action + action * cbind(covariate) %*% gammaD.hazard1,
           # lambda0D.hazard0 * cbind(covariate) %*% exp(beta.hazard0)  + omegaD.hazard0 * action + action * cbind(covariate) %*% gammaD.hazard0
    )
  }
  predHazardFn_R <- function(action, covariate) {
    ifelse(action == 1,
           lambda0R.hazard1 * exp(cbind(covariate) %*% betaR.hazard1 + omegaR.hazard1 * action + action * cbind(covariate) %*% gammaR.hazard1),
           lambda0R.hazard0 * exp(cbind(covariate) %*% betaR.hazard0 + omegaR.hazard0 * action + action * cbind(covariate) %*% gammaR.hazard0)
           # lambda0R.hazard1 * cbind(covariate) %*% exp(beta.hazard1) + omegaR.hazard1 * action + action * cbind(covariate) %*% gammaR.hazard1,
           # lambda0R.hazard0 * cbind(covariate) %*% exp(beta.hazard0) + omegaR.hazard0 * action + action * cbind(covariate) %*% gammaR.hazard0
    )
  }
}
if (local == 1) {
  base_dir <- "./output/"
  figure_dir <- "./figure/"
  dir_rds = sprintf("./output/%s", date_folder)
} else {
  base_dir <- "/work/users/c/w/cwzhou/Proj3RE/output/"
  figure_dir <- "/work/users/c/w/cwzhou/Proj3RE/figure/"
  dir_rds = sprintf("%s/%s", base_dir, date_folder)
}
dir_fig = dir_rds %>% gsub("output/", "figure/", .)
if (!dir.exists(base_dir)) {
  dir.create(base_dir, recursive = TRUE)
}
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}
if (!dir.exists(dir_rds)) {
  dir.create(dir_rds, recursive = TRUE)
}
if (!dir.exists(dir_fig)) {
  dir.create(dir_fig, recursive = TRUE)
}
# print(dir_rds)
# print(dir_fig)
if (local == 0){
  if (endpoint == "CR"){
    dir_rds_tmp = sprintf("/users/c/w/cwzhou/Dissertation/Paper_1/output/%s/%s",
                          generate_failure_method,
                          date_folder)
  } else if (endpoint == "RE"){
    dir_rds_tmp = sprintf("/work/users/c/w/cwzhou/Proj3RE/output/%s",
                          date_folder)
  } else{
    message("endpoint DNE")
  }
  if (savingrds == TRUE){
    if (!dir.exists(dir_rds_tmp)) dir.create(dir_rds_tmp)
  }
}
if (savingrds == TRUE){
  if (!dir.exists(dir_rds)) dir.create(dir_rds)
  if (!dir.exists(dir_fig)) dir.create(dir_fig)
}
if (endpoint == "CR"){
  if (generate_failure_method == "fine_gray"){
    filename = paste0(dir_rds,"/simResult_", generate_failure_method,
                      "_censor", arg2,
                      "_nCauses", arg3, "_cause1prob", arg9,
                      "_beta", arg4, "_prop", arg5,
                      "_n", arg6, "_critS", arg7, "_critE", arg8, ".rds")
  } else if (generate_failure_method == "simple_exp"){
    filename = paste0(dir_rds,"/simResult_", generate_failure_method,
                      "_censor", arg2,
                      "_nCauses", arg3,
                      "_beta", arg4, "_prop", arg5,
                      "_n", arg6, "_critS", arg7, "_critE", arg8, ".rds")
  } else{
    stop("we dont have this generate_failure_method coded yet.")
  }

} else if (endpoint == "RE"){
  filename = paste0(dir_rds, "/simResult_RE",
                    "_censor", arg2, "_prop", arg11,
                    "_n", arg12,
                    "_betaD.", arg3,
                    "_gammaD.", arg5,
                    "_omegaD.", arg7,
                    "_lambda0D.", arg9,
                    ".rds")

} else{
  message("endpoint DNE")
}

### 3. Run the simulation
cv.nodesize = FALSE
message("End of Simulation_Parameters_RE.R")

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
                      nTree = ntree1,
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
# for (sim in 1:n.sim){ # for-loop for sims (if NOT parallelization)
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
  df_recurr = sim.train$dataset_recurrent; ##View(df_recurr)
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
    if (local == 1){
      folder_path <- sprintf("./2_pipeline/%s", date_folder)
    } else{
      folder_path <- sprintf("/work/users/c/w/cwzhou/Proj3RE/2_pipeline/%s",
                             date_folder)
    }
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


  # tau0 = max(timePointsSurvival)

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
  # #View(obs.mff.df); #View(result)

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
      czmk.test.df_recurr = czmk.data.rep$dataset_recurrent; ##View(czmk.test.df_recurr)
      czmk.test.df_surv = czmk.data.rep$dataset_survival; ##View(czmk.test.df_surv)
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
      # #View(czmk.mff.df); #View(result)
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
      zom.test.df_recurr = zom.data.rep$dataset_recurrent; ##View(czmk.test.df_recurr)
      zom.test.df_surv = zom.data.rep$dataset_survival; ##View(czmk.test.df_surv)
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
  mff.stacked_data0 <- do.call(rbind, mff.dataset_list)
  row.names(mff.stacked_data0) <- NULL
  mff.stacked_data = mff.stacked_data0 %>%
    mutate(simulation = sim,
           RE_lived = Number_RE/survival)
  # Append to the list of all simulations' data
  all_sims_data.mff[[sim]] <- mff.stacked_data

  message("Simulation #", sim)
  mff.stacked_data %>%
    group_by(method, Trt) %>%
    summarise(per_trt = round(n()/(n.eval)*100,2),
              mean_surv = mean(survival),
              mean_RE = mean(Number_RE),
              mean_RE_lived = mean(RE_lived)) %>%
    print()
  message("Simulation #", sim)

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
  # #View(result)
  cat("---------------------------------------------------\n")
# } # for (sim in 1:n.sim) but removed for parallelizing
# return(result) # this is for parallel
# } # this is for run_simulation for parallel

# Combine all simulations into one big dataset for MFF
mff_allsims <- do.call(rbind, all_sims_data.mff)
row.names(mff_allsims) <- NULL
# #View(mff_allsims)
if (savingrds == TRUE){
  write.csv(mff_allsims, paste0(dir_rds,"/mff/mff_allsims.csv"), row.names = FALSE)
}

# mff_allsims %>%
#   group_by(simulation, Number_RE, method) %>%
#   summarize(count = n()) %>%
#   #View()
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

# #View(result)

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
  group_by(simulation,method, Trt) %>%
  summarise(per_trt = round(n()/(n.eval)*100,2),
            mean_surv = mean(survival),
            mean_RE = mean(Number_RE),
            # mean_RE_yr = mean(Number_RE/survival),
            mean_RE_lived = mean(RE_lived))

mff_allsims %>%
  group_by(method, Trt) %>%
  summarise(per_trt = round(n()/(n.eval)*100,2),
            #per_trt = round(n()/(n.eval*n.sim)*100,2),
            mean_surv = mean(survival),
            mean_RE = mean(Number_RE),
            # mean_RE_yr = mean(Number_RE/survival),
            mean_RE_lived = mean(RE_lived))

mff_allsims %>%
  group_by(simulation, method) %>%
  summarise(per_trt1 = mean(Trt)*100,
            mean_RE = mean(Number_RE),
            mean_surv = mean(survival),
            # mean_RE_yr = mean(Number_RE/survival),
            mean_RE_lived = mean(RE_lived))

mff_allsims %>%
  group_by(method) %>%
  summarise(mean_surv = mean(survival),
            mean_RE = mean(Number_RE),
            mean_RE_lived = mean(RE_lived))

print("end of script")

