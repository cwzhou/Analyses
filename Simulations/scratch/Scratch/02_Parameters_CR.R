#### C21.simulation_body.R is to be run by source("C21.simulation_run.R")

### 0. Get the setting number.
arg <- commandArgs(trailingOnly = TRUE)
if (length(arg) < 7) {
  warning("commandArgs was not provided. Being set as c(1,2,1,1,1,1,1).")
  arg = c(1, 2, 1, 1, 1, 1, 1) # by default
}
names(arg)[1:7] = c("endpoint", "ncauses", "beta", "propensity", "size", "crit_surv", "crit_endpoint")
arg1 <- as.numeric(arg[1]) # 1..3 # 1=CR,2=RE,3=MC
if (arg1 == 1){
  arg2 <- as.numeric(arg[2]) # 1..inf # number of causes. 
} else{
  arg2 <- as.numeric(1) #arg2=1 when endpoint != CR
}
arg3 <- as.numeric(arg[3]) # 1..4 4 possible beta combinations
arg4 <- as.numeric(arg[4]) # 1..2 2 possible propensity combinations
arg5 <- as.numeric(arg[5]) # 1..2 2 possible size combinations
arg6 <- as.numeric(arg[6]) # 1..3 3 possible crit for survival (prob, mean, prob.mean.combo)
arg7 <- as.numeric(arg[7]) # 1..3 3 possible crit for endpoint (prob, mean, prob.mean.combo)
arg.date <- if (is.na(arg[8]) | arg[8] == "") Sys.Date() else as.character(arg[8])


### 1. setup

# default settings
if (arg1 == 1){
  endpoint = "CR"
  ncauses = arg2 
  cause1_prob <- list(small.cause1_prob = list(cause1_prob = 0.5),
               large.cause1_prob = list(cause1_prob = 0.7))
  arg = c(arg,1)
  arg8 = as.numeric(arg[8])
  names(arg)[8] = "cause1_prob"
  name = sprintf("%s%sp%s", endpoint,ncauses,cause1_prob[[arg8]])
} else if (arg1 == 2){
  endpoint = "RE"
  ncauses = 1
  name = endpoint
} else if (arg1 == 3){
  endpoint = "MC"
  ncauses = 1
  name = endpoint
} else{
  stop("Endpoint must be 1=CR,2=RE, or 3=MC")
}
default <- list(n.eval = 1, n.sim = 1, tau = 10, 
                censor_rate = 0.2, 
                endpoint = endpoint, ncauses = ncauses)

# arg4 crit
if (endpoint == "CR"){
  crit3 = list(criterion_phase2 = "mean", crit.value_phase2 = NULL, value_phase2 = "truncated CIF mean CIF[T]")
  crit4 = list(criterion_phase2 = "surv.mean", crit.value_phase2 = 5, value_phase2 = "CIF(5)")
}
crit_surv <- list(crit1 = list(criterion_phase1 = "mean", crit.value_phase1 = NULL, value_phase1 = "truncated survival mean E[T]"),
                  crit2 = list(criterion_phase1 = "surv.mean", crit.value_phase1 = 5, value_phase1 = "S(5)"))
crit_endpoint <- list(crit1 = crit3,
                      crit2 = crit4)

# arg3 size
size <- list(small.sample.size = list(n = 10),
             large.sample.size = list(n = 50))

# arg2 propensityx
propensity <-   # (int), covariate (1~5)
  list(obs = list(beta.propensity = function(p) c(0, -rep(0.5, p))),  # no unmeasured confounder
       rct  = list(beta.propensity = function(p) c(0, rep(0, p))))    # RCT

# arg1 beta
if (endpoint == "CR"){
  betas <- list(
    # ncov = 5 and 2 CRs
    beta1 =
      list(beta1.hazard0 = c(c(1,1,1, 1, 1)),              # covariate (1~5)
           beta1.hazard1 = c(c(2,2,1,-1,-1)),              # covariate (1~5)
           beta2.hazard0 = c(0.2 * c(1,1,1, 1, 1)),      # covariate (1~5)
           beta2.hazard1 = c(0.2 * c(1,1,1,-1,-1)))     # covariate (1~5)
  )
  ncov.list <- lapply(betas, function(x) length(x$beta1.hazard0)) # no intercept included
} else{
  stop("betas not coded up yet for endpoint: ", endpoint)
}


# setting1 n.boot 50 / n 300
# this setting is for the first of each thing, since arg = all 1s
# args 1- 7: "endpoint", "ncauses", "beta", "propensity", "size", "crit_surv", "crit_endpoint"
if (endpoint == "CR"){
  setting = c(arg = list(arg), default, betas[[arg3]], ncov = ncov.list[[arg3]],
              propensity[[arg4]], size[[arg5]], crit_surv[[arg6]], crit_endpoint[[arg7]], cause1_prob[[arg8]]) 
} else{
  setting = c(arg = list(arg), default, betas[[arg3]], ncov = ncov.list[[arg3]],
              propensity[[arg4]], size[[arg5]], crit_surv[[arg6]], crit_endpoint[[arg7]]) 
}


### 2. put a selected setting into the global environment    
cat("setting (endpoint, beta, propensity, size, crit_phase1, crit_phase2, opt cause1_prob) ", arg, "\n")
list2env(setting, envir = globalenv())
beta.propensity <- beta.propensity(ncov)
filename = paste0("output/simResult_", arg.date, "_", name, "_beta", arg3, "_prop", arg4, 
                  "_n", arg5, "_critS", arg6, "_critE", arg7, ".rds")

# complete the scores
predPropensityFn <- function(covariate) {
  plogis(cbind(1, covariate) %*% beta.propensity)
}
if (endpoint == "CR"){
  predHazardFn <- function(action, covariate, cause) {
    if (cause == 1){
      ifelse(action == 1,
             covariate %*% beta1.hazard1,
             covariate %*% beta1.hazard0)
    } else if (cause == 2){
      ifelse(action == 1, 
             covariate %*% beta2.hazard1,
             covariate %*% beta2.hazard0)
    } else{
      stop("CURRENTLY ONLY CODED FOR 2 CAUSES")
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
skip.csk <- skip.gk <- skip.dw <- skip.zom <- TRUE;
cv.nodesize = FALSE

### More description
# n       # training sample size
# n.eval  # evaluation sample size
# n.sim   # number of simulation replicates
# ncov   # number of covariates (should agree with beta coefficients' dimension)
# tau  # total study length

# source("Testi_Simulation_Body.R")