# sbatch -p general -N 1 --mem=100GB -n 2 -t 10-07:00:00 --mail-type=end --mail-user=cwzhou@email.unc.edu --wrap="Rscript 01.Simulation_Run_RE.R > REoutput_20250111_500nsims_300ntree.txt"

# Title: Setting hyperparameters for Simulation Scripts
# Description: [Brief description of the purpose and objectives of the simulation]
# Author: Christina Zhou
# Date: 11.12.2024

### DATA TYPE ###
local = 0
if (local == 1){
  setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper3_RE")
} else{
  setwd("/nas/longleaf/home/cwzhou/Dissertation/Analyses/Simulations/Paper3_RE")
}
# Libraries --------------------------------------------------------------
source("02.Simulation_Libraries_RE.R")
# Functions -----------------------------------------------------------------
source("02.Simulation_Functions_RE.R")

savingrds = TRUE
date_folder = "2025-01-11"
n.eval = 1000
n.sim = 500
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

  # arg1 <- as.numeric(arg[1]) # 1..3 # 1=CR,2=RE,3=MC
  # arg2 <- as.numeric(arg[2]) # 1,2 #censoring = 20%, censoring = 50%
  # arg3 <- as.numeric(arg[3]) # 1..4 4 possible beta_D combinations
  # arg4 <- as.numeric(arg[4]) # 1..4 4 possible beta_R combinations
  # arg5 <- as.numeric(arg[5]) # 1..4 4 possible gamma_D combinations
  # arg6 <- as.numeric(arg[6]) # 1..4 4 possible gamma_R combinations
  # arg7 <- as.numeric(arg[7]) # 1..4 4 possible omega_D combinations
  # arg8 <- as.numeric(arg[8]) # 1..4 4 possible omega_R combinations
  # arg9 <- as.numeric(arg[9]) # 1..4 4 possible lambda_0D combinations
  # arg10 <- as.numeric(arg[10]) # 1..4 4 possible lambda_0R combinations
  # arg11 <- as.numeric(arg[11]) # 1..2 2 possible propensity combinations
  # arg12 <- as.numeric(arg[12]) # 1..2 2 possible size combinations
  # arg13 <- as.numeric(arg[13]) # 1..3 3 possible crit for surv (mean, mean.prob.combo, prob)
  # arg14 <- as.numeric(arg[14]) # 1..3 3 possible crit for ep (mean, mean.prob.combo, prob)
  # arg.date <- if (is.na(arg[num.args+1]) | arg[num.args+1] == "") date_folder else as.character(arg[num.args+1])


} else{
  message("endpoint DNE")
}

# print(filename)


### 3. Run the simulation
cv.nodesize = FALSE
# if (skip.zom == "FALSE"){
#   skip.czmk <- FALSE
# }

# ### More description
# n       # training sample size
# n.eval  # evaluation sample size
# n.sim   # number of simulation replicates
# ncovD   # number of covariates for terminal event (should agree with betaD coefficients' dimension)
# ncovR   # number of covariates for recurrent event (should agree with betaR coefficients' dimension)
# tau  # total study length
# beta1.hazard1

message("End of 02.Simulation_Parameters_RE.R")


# if (sim_data_type == "RE"){
#   ##### Study Design Hyperparameters #####
#   seed1 = 2024
#   N = 80 # Number of subjects
#   tau0 = 1 # stop time/end of study
#   G = 3 # total gap times
#   ##### Covariate-Related Hyperparameters #####
#   num_A = 2 # Number of treatments
#   lambda_0D = 1 # baseline hazard for disease
#   lambda_0R = 1 # baseline hazard for recurrence
#   beta_D = log(2) # parameters for disease for covariates
#   beta_R = log(2) # parameters for recurrence for covariates
#   omega_D = log(3) #parameters for disease for treatment
#   omega_R = log(3) #parameters for treatment+covariate interaction
#   ztype = 0 # covariate distribution
#   zparam = 0.3 # covariate distribution parameter
#
#   ## TREATMENT COVARIATE INTERACTION ##
#
#   ##### Censoring Hyperparameters #####
#   ctype = 1 # censoring distribution
#   cparam = 2 # censoring distribution parameter
#   ##### Gap Time Hyperparameters #####
#   gaptype = 0 # failure vs gap time indicator
#   gapparam1 = 0.1 #rho for failure
#   gapparam2 = 0.7 #rho for gap times
#
#   print(sprintf("Parameters: N=%s, tau0=%s, G=%s, num_A=%s, lambda_0D=%s, lambda_0R=%s, beta_D=%s, beta_R=%s, omega_D=%s, omega_R=%s, ztype=%s, zparam=%s, ctype=%s, cparam=%s, gaptype=%s, gapparam1=%s, gapparam2=%s.",
#                 N, tau0, G, num_A,
#                 lambda_0D, lambda_0R, beta_D, beta_R,
#                 omega_D, omega_R,
#                 ztype, zparam, ctype, cparam, gaptype, gapparam1, gapparam2))
# }
