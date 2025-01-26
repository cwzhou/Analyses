# things to note:
# - censoring for training; no censoring for testing; truncated at tau for both (differs by sim generate failure setting)
# right now, in cr01 t0_pmcr is set to 0.2 regardless of other parameters.
 # local = 1 # uncomment for CR02.Simulation_Summary.R plots # comment for running sims
local = 0 # local = 0 for cluster
#
# # below is only needed if running this script directly. comment out if running CR01.Simulation_Run.R script.
# # uncomment if running alone (aka for CR02.Simulation_Summary.R
# if (local == 1){
#    setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper1_CR")
#  } else{
#    setwd("/nas/longleaf/home/cwzhou/Dissertation/Analyses/Simulations/Paper1_CR")
#  }

#### libraries and functions
source("F01.Simulation_Functions.R") # calls libraries

savingrds = TRUE
# date_folder = "2024-09-09" # "2024-08-31" #Sys.Date() # "2024-08-20/24" #"2024-02-27" # this is the most recent date with results; # very old date: "2024-02-18"
date_folder = "2025-01-26" #"2024-09-13" #Sys.Date()
n.eval = 1000 #n.eval = 10000
n.sim = 10
mean_tol1 = c(0.07,0) # this is for differences in years so we don't want it to be too big
prob_tol1 = c(0.15, 0.01)
combo_tol1 = c(mean_tol1[1], prob_tol1[1], mean_tol1[2], prob_tol1[2])
generate_failure_method = c("simple_exp","fine_gray") #"simple_exp" # "fine_gray"
generate_failure_method = generate_failure_method[2]

if (generate_failure_method == "simple_exp"){
  crit_t0_eval = 1 #1 year (we dont use days bc its calculated using the rates which was for years)
} else if (generate_failure_method == "fine_gray"){
  crit_t0_eval = 1 # 1 year (we dont use days bc its calculated using the rates which was for years)
} else{
  stop("crit_t0_eval for generate_failure_method not defined in CR00.Simulation_Parameters.R.")
}

# Specify the methods and skip.methods
all_methods <- c("czmk", "csk", "pmcr", "aipwe", "zom", "obs");
skip_method <- !c(TRUE, TRUE, !TRUE, TRUE, TRUE, TRUE);
# skip_method <- c(!TRUE, TRUE, TRUE, TRUE, TRUE, TRUE);

#### Run this Script FOR CR. Change name later.

# Specify the methods and hyperparameters
## ADD THIS IN SOMEWHERE (fix so that these are not needed for CR02.Simulation_Run.R)
n.methods <- length(all_methods)
n.phases = 2
n.causes = 2
# CR setting so we want to include a second model for priority cause! (1 here)
# Specify the priority cause
priority_cause = 1
# we assume priority_cause is always called cause 1 # required

# Loop to create logical SKIP objects for each method and assign skip_method
assign_skip_function(all_methods, skip_method)

### 0. Get the setting number.
arg <- commandArgs(trailingOnly = TRUE)
if (length(arg) < 9) {
  arg = c(1, 1, 1, 1, 1, 1, 1, 1, 1) # by default
  warning(sprintf("commandArgs was not provided. Set as c(%s).",
                  toString(arg)))
}
names(arg)[1:9] = c("endpoint", "censor", "ncauses", "beta",
                    "propensity", "size", "crit_surv",
                    "crit_endpoint", "cause1prob")
print("arg:")
print(arg)

arg1 <- as.numeric(arg[1]) # 1..3 # 1=CR,2=RE,3=MC
arg2 <- as.numeric(arg[2]) # 1,2 #censoring = 20%, censoring = 50%
arg3 <- as.numeric(arg[3]) # 1..inf, number of causes
arg4 <- as.numeric(arg[4]) # 1..4 4 possible beta combinations
arg5 <- as.numeric(arg[5]) # 1..2 2 possible propensity combinations
arg6 <- as.numeric(arg[6]) # 1..2 2 possible size combinations
arg7 <- as.numeric(arg[7]) # 1..3 3 possible crit for surv (mean, mean.prob.combo, prob)
arg8 <- as.numeric(arg[8]) # 1..3 3 possible crit for ep (mean, mean.prob.combo, prob)
arg9 <- as.numeric(arg[9]) # 1..2 2 possible probabilities for cause1prob
arg.date <- if (is.na(arg[10]) | arg[10] == "") date_folder else as.character(arg[10])

### 1. setup

# default settings
if (arg1 == 1){
  endpoint = "CR"
} else if (arg1 == 2){
  endpoint = "RE"
  # ncauses = 1
} else if (arg1 == 3){
  endpoint = "MC"
  # ncauses = 1
} else{
  stop("Endpoint must be 1=CR, 2=RE, or 3=MC")
}

if (generate_failure_method == "simple_exp"){
  tau.1 = 2
  # message("IN F01.DynamicsCR.R: we have f2_constant = 0.3 to make cause2 prevalance less and f1_constant.")
  # message("simple_exp with exp censoring")
  default <- list(n.eval = n.eval,
                  n.sim = n.sim,
                  tau = tau.1,#365,#1 year - we use days here b/c censoring is for days #2.5,
                  generate_failure_method = generate_failure_method,
                  endpoint = endpoint)
} else if (generate_failure_method == "fine_gray"){
  tau.1 = 3
  # message("fine_gray with uniform censoring")
  default <- list(n.eval = n.eval,
                  n.sim = n.sim,
                  tau = tau.1,#365,
                  generate_failure_method = generate_failure_method,
                  endpoint = endpoint)
} else{
  stop("generated failure time setting not specified")
}


# arg7 and 8 crit
crit_surv <- list(crit1 = list(criterion_phase1 = "mean",
                               crit.value_phase1 = crit_t0_eval,#NULL,
                               value_phase1 = "truncated overall survival mean E[T]",
                               tol1 = mean_tol1),
                  crit2 = list(criterion_phase1 = "mean.prob.combo",
                               crit.value_phase1 = crit_t0_eval, #* 2200,
                               value_phase1 = "S_0(t_0)",
                               tol1 = combo_tol1)
                  )
if (endpoint == "CR"){
  crit3 = list(criterion_phase2 = "mean",
               crit.value_phase2 = crit_t0_eval, #NULL,
               value_phase2 = "truncated PC-CIF mean PC-CIF[T] \n or Years Lost due to PC"
  )
  crit4 = list(criterion_phase2 = "mean.prob.combo",
               crit.value_phase2 = crit_t0_eval, #* 2200,
               value_phase2 = "PC-CIF(t_0)"
  )
  crit_endpoint <- list(crit1 = crit3,
                        crit2 = crit4)
}

# arg9
if (endpoint == "CR" & generate_failure_method == "fine_gray"){
  cause1_prob <- list(small.cause1prob = list(cause1prob = 0.2),
                      large.cause1prob = list(cause1prob = 0.8))
} else{
  cause1_prob <- list(cause1.prob = list(cause1prob = 1))
}

# arg3 ncauses
if (endpoint == "CR"){
  ncauses <- list(small.ncauses = list(M = 2))#,
  # large.ncauses = list(M = 4))
} else{
  ncauses = list(ncauses = list(M = 1L))
}

# arg6 size
size <- list(small.sample.size = list(n = 300),
             large.sample.size = list(n = 700))

# arg2 size
if (generate_failure_method == "fine_gray"){
  censor <- list(low.censoring = list(
    # we want about 20% censoring
    ctype = 1, # uniform censoring
    censor_min = 0,
    censor_max = tau.1+tau.1/2, # higher is less censoring for unif
    censor_rate = 0  # not used for unif
    ),
  high.censoring = list(
    # we want about 50% censoring
    ctype = 1, # uniform censoring
    censor_min = 0,
    censor_max = 0.7*tau.1, # lower is more censoring for unif
    censor_rate = 0 # not used for unif
    ))
  } else if (generate_failure_method == "simple_exp"){
    censor <- list(low.censoring = list(ctype = 0, # exp censoring
                                      censor_min = 0, # not used for exp
                                      censor_max = 0,# not used for exp
                                      censor_rate = 0.5 # lower censoring: range
                                  ),
                   high.censoring = list(ctype = 0, # exp censoring
                                         censor_min = 0, # not used for exp
                                         censor_max = 0, # not used for exp
                                         censor_rate = 3.3 # higher is more censoring for exp
                                   ))
  } else{
    stop("generate_failure_method is only coded up for simple exponential and fine-gray setting right now.")
  }

# arg5 propensity
propensity <-   # (int), covariate (1~5)
  list(obs = list(beta.propensity = function(p) c(0, -rep(0.5, p))), #no unmeasured confounder
       rct  = list(beta.propensity = function(p) c(0, rep(0, p))))  # RCT

# arg4 beta
if (endpoint == "CR"){

  if (default$generate_failure_method == "simple_exp"){
    # print("simple exp coeff")
    # currently using runif covariates

    betas <- list(
      # currently using runif covariates

      # we want difference between cause1 and cause2 (csk);
      # also want difference between 0 and 1 within cause1 (zom)

      # this works when tau = 2
      # N = 300 works better than N = 1000
      beta1 = list(
        beta1.hazard0 = c(0, -1.5,-1.9,-0.3), # (int), covariate (1~3)
        beta1.hazard1 = c(0, -0.6,-1.9,1.5),
        beta2.hazard0 = c(0, 1,-0.1,0.4),#c(0, 1.1,-1.3,0.3), #c(0,-0.1,-0.2),
        beta2.hazard1 = c(0, 1.2, 0.2,-1)),#c(0, -0.6,0.5,0.3)),

      # 6 covariates
      beta2 = list(
        beta1.hazard0 = c(0,-0.5,-1.2,1.5,0.1,0.6), #c(0, -0.5,-0.5,-0.5,-0.3,0.3), # (int), covariate (1~6)
        beta1.hazard1 = c(0,0.3,-0.4,-0.2,0.9,1.6), #c(0, -0.1,0.6,-0.9,0.5,-0.3),
        beta2.hazard0 = c(0,-0.6,1.4,1.5,2,1), #c(0, 0.1,0.3,-0.6,0.3,0.8),
        beta2.hazard1 = c(0,1.5,-2.1,-0.8,0.1,1)) #c(0, -0.3,-0.3,-0.2,-0.2,0.2))
    )

  } else if (default$generate_failure_method == "fine_gray"){
    # print("fine-gray coeff")
    # currently using runif covariates

    # temporary while working on testing dataset
    betas <- list(
      beta1 = list(
        beta1.hazard0 = c(0,0.1,0.3),#c(0,-1,-1.4),
        beta1.hazard1 = c(0,-1.8,-1.5),#c(0,0.8,0.7),
        beta2.hazard0 = c(0,-1.1,-0.3),#c(0,-0.2,1.2), #c(0,0.2,0.8),
        beta2.hazard1 = c(0,-0.2,1.2)),#c(0,-0.3,-2)),
      beta2 = list(
        beta1.hazard0 = c(0,0.9,-0.7,0.6,-0.1,0.5),
        beta1.hazard1 = c(0,-1.2,0.1,-0.5,-0.2,-1.2),
        beta2.hazard0 = c(0,-0.5,-0.4,-1.1,-0.8,-0.4), #c(0,-0.1,-0.2),
        beta2.hazard1 = c(0,0.4,-0.3,0.6,0.2,1.1))
    )
  } else{
    stop("generate failure method not specified")
  }

  alpha1 <- 0 #log(lambda0)
  alpha2 <- 0 #log(lambda0)

  ncov.list <- lapply(betas, function(x) length(x$beta1.hazard1) -1 ) # removing intercept
} else{
  stop("betas not coded up yet for endpoint: ", endpoint)
}

# setting1 n.boot 50 / n 300
# this setting is for the first of each thing, since arg = all 1s
# args 1- 9: "endpoint", "censor", "ncauses", "beta", "propensity",
# "size", "crit_surv", "crit_endpoint", "cause1prob" (for CR)
if (endpoint == "CR"){
  if (generate_failure_method == "fine_gray"){
    setting = c(all_methods = list(all_methods),
                n.methods = n.methods,
                arg = list(arg), default, ncauses[[arg3]],
                censor[[arg2]],
                cause1_prob[[arg9]],
                betas[[arg4]],
                ncov = ncov.list[[arg4]],
                propensity[[arg5]],
                size[[arg6]],
                crit_surv[[arg7]],
                crit_endpoint[[arg8]])
  } else if (generate_failure_method == "simple_exp"){
    setting = c(all_methods = list(all_methods),
                n.methods = n.methods,
                arg = list(arg), default, ncauses[[arg3]],
                censor[[arg2]],
                # cause1_prob[[arg9]],
                betas[[arg4]],
                ncov = ncov.list[[arg4]],
                propensity[[arg5]],
                size[[arg6]],
                crit_surv[[arg7]],
                crit_endpoint[[arg8]])
  } else{
    stop ("generate failure method DNE")
  }
} else{
  setting = c(all_methods = all_methods,
              n.methods = n.methods,
              arg = list(arg), default, ncauses[[arg3]],
              censor[[arg2]],
              betas[[arg4]], ncov = ncov.list[[arg4]],
              propensity[[arg5]], size[[arg6]],
              crit_surv[[arg7]], crit_endpoint[[arg8]])
}


### 2. put a selected setting into the global environment
cat("setting (endpoint, censor, ncauses, beta, propensity,
    size, crit_phase1, crit_phase2, cause1_prob) ", arg, "\n")
list2env(setting, envir = globalenv())
beta.propensity <- beta.propensity(ncov)


dir_rds = sprintf("./output/%s/%s", generate_failure_method, date_folder)
dir_fig = dir_rds %>% gsub("output/", "figure/", .)
print(dir_rds)
print(dir_fig)
if (local == 0){
  dir_rds_tmp = sprintf("/users/c/w/cwzhou/Dissertation/Paper_1/output/%s/%s",
                        generate_failure_method,
                        date_folder)
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

} else{
  filename = paste0(dir_rds, "/simResult_", generate_failure_method,
                    "_censor", arg2,
                    "_beta", arg4, "_prop", arg5,
                    "_n", arg6, "_critS", arg7, "_critE", arg8, ".rds")
}

print(filename)

# complete the scores
predPropensityFn <- function(covariate) {
  # message("predPropensityFn")
  plogis(cbind(1, covariate) %*% beta.propensity)
}
if (endpoint == "CR"){
  # alpha1/2 correspond to baseline hazard (to make failure times higher)
  # alpha1=alpha2=0 for now; if i change that then need to include in "settings"
  predHazardFn <- function(action, covariate, cause) {
    # message("predHazardFn")
    if (cause == 1){
      # message("cause1")
      ifelse(action == 1,
             alpha1 + cbind(1, covariate) %*% beta1.hazard1, #cause1 treatment1
             alpha1 + cbind(1, covariate) %*% beta1.hazard0) #cause1 treatment 0
      # mat_beta = c(beta11,beta12,beta13)
      # # print(dim(mat_cov))
      # # print(mat_beta)
      # # print(dim(mat_beta))
      # alpha1 + mat_cov %*% mat_beta
    } else if (cause == 2){
      # message("cause2")
      ifelse(action == 1,
             alpha2 + cbind(1, covariate) %*% beta2.hazard1, #cause2 treatment1
             alpha2 + cbind(1, covariate) %*% beta2.hazard0) #cause2 treatment0
      # mat_beta = c(beta21,beta22,beta23)
      # alpha2 + mat_cov %*% mat_beta
    } else if (cause == 3){
      # message("cause3")
      ifelse(action == 1,
             alpha2 + cbind(1, covariate) %*% beta3.hazard1, #cause3 treatment1
             alpha2 + cbind(1, covariate) %*% beta3.hazard0) #cause3 treatment0
      # mat_beta = c(beta21,beta22,beta23)
      # alpha2 + mat_cov %*% mat_beta
    }
    else{
      stop("CURRENTLY ONLY CODED FOR 2 (or 3?) CAUSES")
    }

  }
} else{
  predHazardFn <- function(action, covariate, cause = NULL) {
    ifelse(action == 1,
           cbind(covariate) %*% beta.hazard1,
           cbind(covariate) %*% beta.hazard0)
  }
}

### 3. Run the simulation
cv.nodesize = FALSE
if (skip.zom == "FALSE"){
  skip.czmk <- FALSE
}

# ### More description
# n       # training sample size
# n.eval  # evaluation sample size
# n.sim   # number of simulation replicates
# ncov   # number of covariates (should agree with beta coefficients' dimension)
# tau  # total study length
# beta1.hazard1

message("End of CR00.Simulation_Parameters.R")

