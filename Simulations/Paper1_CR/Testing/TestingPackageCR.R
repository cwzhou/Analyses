samplesize.tmp = 300

# obs.data = read.csv("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper1_CR/MISC/example_datadf.csv")

# things to note:
# - censoring for training; no censoring for testing; truncated at tau for both (differs by sim generate failure setting)
# right now, in cr01 t0_pmcr is set to 0.2 regardless of other parameters.
local = 1 #0 # local = 0 for cluster
#### libraries and functions
source("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper1_CR/F01.Simulation_Functions.R") # calls libraries

date_folder = Sys.Date()
n.eval = 1000 #n.eval = 10000
n.sim = 1 #n.sim = 200
mean_tol1 = c(0.10,0)
prob_tol1 = c(0.3, 0.01)
combo_tol1 = c(mean_tol1[1], prob_tol1[1], mean_tol1[2], prob_tol1[2])
generate_failure_method = c("simple_exp","fine_gray") #"simple_exp" # "fine_gray"
generate_failure_method = generate_failure_method[2]

if (generate_failure_method == "simple_exp"){
  crit_t0_eval = 1 # one year?
} else if (generate_failure_method == "fine_gray"){
  crit_t0_eval = 0.5 # half year?
} else{
  stop("crit_t0_eval for generate_failure_method not defined in CR00.Simulation_Parameters.R.")
}

# Specify the methods and skip.methods
all_methods <- c("czmk", "csk", "pmcr", "aipwe", "zom", "obs");
# skip_method <- !c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE);
skip_method <- c(!TRUE, TRUE, TRUE, TRUE, TRUE, TRUE);
savingrds = FALSE

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
if (length(arg) < 8) {
  arg = c(1, 1, 1, 1, 1, 1, 1, 1) # by default
  warning(sprintf("commandArgs was not provided. Set as c(%s).",
                  toString(arg)))
}
names(arg)[1:8] = c("endpoint", "ncauses", "beta",
                    "propensity", "size", "crit_surv",
                    "crit_endpoint", "cause1prob")
# print("arg:")
# print(arg)

arg1 <- as.numeric(arg[1]) # 1..3 # 1=CR,2=RE,3=MC
arg2 <- as.numeric(arg[2]) # 1..inf, number of causes
arg3 <- as.numeric(arg[3]) # 1..4 4 possible beta combinations
arg4 <- as.numeric(arg[4]) # 1..2 2 possible propensity combinations
arg5 <- as.numeric(arg[5]) # 1..2 2 possible size combinations
arg6 <- as.numeric(arg[6]) # 1..3 3 possible crit for surv (mean, mean.prob.combo, prob)
arg7 <- as.numeric(arg[7]) # 1..3 3 possible crit for ep (mean, mean.prob.combo, prob)
arg8 <- as.numeric(arg[8]) # 1..2 2 possible probabilities for cause1prob
arg.date <- if (is.na(arg[9]) | arg[9] == "") date_folder else as.character(arg[9])

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
  stop("Endpoint must be 1=CR,2=RE, or 3=MC")
}

if (generate_failure_method == "simple_exp"){
  # message("IN F01.multiPhaseDynamicsCR.R: we have f2_constant = 0.3 to make cause2 prevalance less and f1_constant.")
  message("simple_exp with exp censoring")
  default <- list(n.eval = n.eval,
                  n.sim = n.sim,
                  tau = 2,#2.5,
                  ctype = 0, # exp censoring
                  censor_min = 0,
                  censor_max = 0.5,
                  censor_rate = 0.8,#0.6, # lower is less censoring for exp
                  generate_failure_method = generate_failure_method,
                  endpoint = endpoint)
} else if (generate_failure_method == "fine_gray"){
  message("fine_gray with uniform censoring")
  default <- list(n.eval = n.eval,
                  n.sim = n.sim,
                  tau = 4,
                  ctype = 1, # uniform censoring
                  censor_min = 0,
                  censor_max = 3, # higher is less censoring for unif
                  censor_rate = 0.1,
                  generate_failure_method = generate_failure_method,
                  endpoint = endpoint)
} else{
  stop("generated failure time setting not specified")
}


# arg6 and 7 crit
crit_surv <- list(crit1 = list(criterion_phase1 = "area",
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

# arg8
if (endpoint == "CR" & generate_failure_method == "fine_gray"){
  cause1_prob <- list(small.cause1prob = list(cause1prob = 0.3),
                      large.cause1prob = list(cause1prob = 0.7))
} else{
  cause1_prob <- list(cause1.prob = list(cause1prob = 1))
}

# arg2 ncauses
if (endpoint == "CR"){
  ncauses <- list(small.ncauses = list(M = 2))#,
  # large.ncauses = list(M = 4))
} else{
  ncauses = list(ncauses = list(M = 1L))
}

# arg5 size
size <- list(small.sample.size = list(n = samplesize.tmp),
             large.sample.size = list(n = 200))

# arg4 propensity
propensity <-   # (int), covariate (1~5)
  list(obs = list(beta.propensity = function(p) c(0, -rep(0.5, p))), #no unmeasured confounder
       rct  = list(beta.propensity = function(p) c(0, rep(0, p))))  # RCT

# arg3 beta
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
        beta1.hazard0 = c(0, -1.6,-1.2,1.5), # (int), covariate (1~3)
        beta1.hazard1 = c(0, 0.3,-0.4,-0.2),
        beta2.hazard0 = c(0, 1.1,-1.3,0.3), #c(0,-0.1,-0.2),
        beta2.hazard1 = c(0, -0.6,0.5,0.3)),

      # 6 covariates
      beta2 = list(
        beta1.hazard0 = c(0, -0.5,-0.5,-0.5,-0.3,0.3), # (int), covariate (1~6)
        beta1.hazard1 = c(0, -0.1,0.6,-0.9,0.5,-0.3),
        beta2.hazard0 = c(0, 0.1,0.3,-0.6,0.3,0.8),
        beta2.hazard1 = c(0, -0.3,-0.3,-0.2,-0.2,0.2))
    )

  } else if (default$generate_failure_method == "fine_gray"){
    # print("fine-gray coeff")
    # currently using runif covariates

    # temporary while working on testing dataset
    betas <- list(
      beta1 = list(
        beta1.hazard0 = c(0,-1,-1.4),
        beta1.hazard1 = c(0,0.8,0.7),
        beta2.hazard0 = c(0,-0.2,1.2), #c(0,0.2,0.8),
        beta2.hazard1 = c(0,-0.3,-2)),
      beta2 = list(
        beta1.hazard0 = c(0,-0.9,-0.7,-0.6,0.1,-0.5),
        beta1.hazard1 = c(0,0.5,-0.1,0.5,0.2,-0.2),
        beta2.hazard0 = c(0,0.5,0.4,0.2,-0.2,-0.2), #c(0,-0.1,-0.2),
        beta2.hazard1 = c(0,-0.4,-0.3,-0.6,-0.2,0.5))
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
# args 1- 7: "endpoint", "ncauses", "beta", "propensity",
# "size", "crit_surv", "crit_endpoint"
if (endpoint == "CR"){
  if (generate_failure_method == "fine_gray"){
    setting = c(all_methods = list(all_methods),
                n.methods = n.methods,
                arg = list(arg), default, ncauses[[arg2]],
                cause1_prob[[arg8]],
                betas[[arg3]], ncov = ncov.list[[arg3]],
                propensity[[arg4]], size[[arg5]],
                crit_surv[[arg6]], crit_endpoint[[arg7]])
  } else if (generate_failure_method == "simple_exp"){
    setting = c(all_methods = list(all_methods),
                n.methods = n.methods,
                arg = list(arg), default, ncauses[[arg2]],
                # cause1_prob[[arg8]],
                betas[[arg3]], ncov = ncov.list[[arg3]],
                propensity[[arg4]], size[[arg5]],
                crit_surv[[arg6]], crit_endpoint[[arg7]])
  } else{
    stop ("generate failure method DNE")
  }

} else{
  setting = c(all_methods = all_methods,
              n.methods = n.methods,
              arg = list(arg), default, ncauses[[arg2]],
              betas[[arg3]], ncov = ncov.list[[arg3]],
              propensity[[arg4]], size[[arg5]],
              crit_surv[[arg6]], crit_endpoint[[arg7]])
}


### 2. put a selected setting into the global environment
cat("setting (endpoint, ncauses, beta, propensity,
    size, crit_phase1, crit_phase2, cause1_prob) ", arg, "\n")
list2env(setting, envir = globalenv())
beta.propensity <- beta.propensity(ncov)


dir_rds = sprintf("./output/%s/%s", generate_failure_method, date_folder)
dir_fig = dir_rds %>% gsub("output/", "figure/", .)
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
                      "_nCauses", arg2, "_cause1prob", arg8,
                      "_beta", arg3, "_prop", arg4,
                      "_n", arg5, "_critS", arg6, "_critE", arg7, ".rds")
  } else if (generate_failure_method == "simple_exp"){
    filename = paste0(dir_rds,"/simResult_", generate_failure_method,
                      "_nCauses", arg2,
                      "_beta", arg3, "_prop", arg4,
                      "_n", arg5, "_critS", arg6, "_critE", arg7, ".rds")
  } else{
    stop("we dont have this generate_failure_method coded yet.")
  }

} else{
  filename = paste0(dir_rds, "/simResult_", generate_failure_method,
                    "_beta", arg3, "_prop", arg4,
                    "_n", arg5, "_critS", arg6, "_critE", arg7, ".rds")
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
message("End of CR00.Simulation_Parameters.R")

#### This script is run automatically from 02_Simulation_Parameters_CR.R. DO not run alone.
setting_seed = 0

start_time = Sys.time()
column_names = NULL

if (criterion_phase1 %in% c("mean", "area", "surv.area")) {
  crit.eval.os = "mean"
  # overall_survival_val.fn <- function(data){
  #   mean(subset(data$event.time, data$status > 0),
  #        na.rm = TRUE)
  #   }
} else if (criterion_phase1 %in% c("prob", "mean.prob.combo")) {
  crit.eval.os = "prob"
  # overall_survival_val.fn <- function(data) {
  #   mean(subset(data$event.time,
  #               data$status > 0) >= crit.value_phase1,
  #        na.rm = TRUE)}
}
# NOTE: below is assuming endpoint = CR aka we want mean LESS THAN critvalue
if (criterion_phase2 %in% c("mean", "area", "surv.area")) {
  crit.eval.cif = "mean"
} else if (criterion_phase2 %in% c("prob", "mean.prob.combo")) {
  crit.eval.cif = "prob"
  # endpoint_val.fn <- function(data, priority_cause) {
  #   mean(subset(data$event.time,
  #               # MAKE SURE THIS IS RIGHT since survival is >=
  #               data$status == priority_cause) <= crit.value_phase2,
  #        na.rm = TRUE)}
}

print(filename)
filename.tmp <- gsub("\\.rds", "_tmp.rds", filename)
if (local == 0){
  extracted_part <- gsub(".*/([[:digit:]]{4}-[[:digit:]]{2}-[[:digit:]]{2})/(.*).rds", "\\2", filename)
  filename.true2 <- sprintf("%s/%s_true2.rds", dir_rds_tmp, extracted_part)
} else{
  filename.true2 <- gsub("\\.rds", "_true2.rds", filename)
}
### simulation
r00 = data.frame(sim = n.sim)
colnames1 = lapply(all_methods, function(method) {
  c(
    paste(method, "training_os_trtprop", sep = "_"),
    paste(method, "training_cause1_trtprop", sep = "_"),
    paste(method, "testing_os_trtprop", sep = "_")
  )
})
# Flatten the list of column names
column_names1 <- unlist(colnames1)
r00[column_names1] <- NA

# Create the data frame with the 'sim' column
result <- data.frame(sim = 1:n.sim)
# Generate column names based on methods
column_names <- lapply(all_methods, function(method) {
  if (method == "obs" | method == "csk" | method == "pmcr" | method == "aipwe") {
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
    !grepl("^csk_", names),
    !grepl("^pmcr_", names),
    !grepl("^aipwe_", names),
    !grepl("^zom_", names),
    !grepl("^obs_", names),
    names
  )
}
sorted_column_names <- column_names[custom_sort(column_names)]
# print(sorted_column_names)
# Add columns to the result data frame
result[sorted_column_names] <- NA
attr(result, "criterion_phase1") <- list(criterion = criterion_phase1, crit.value = crit.value_phase1)
attr(result, "criterion_phase2") <- list(criterion = criterion_phase2, crit.value = crit.value_phase2)

trt_result <- data.frame(sim = rep(1:n.sim, each = n.methods),
                         method = NA,
                         surv_A = NA, surv_B = NA,
                         endpoint_A = NA, endpoint_B = NA) # each = 3 means repeating b/c theres 3 methods rn (obs, czmk, zom)

if (endpoint != "CR" | generate_failure_method != "fine_gray"){
  message("setting cause1prob to NULL b/c not in fine-gray setting")
  cause1prob = NULL
}
if (endpoint == "CR"){
  arg_list =   list(
    N = n, tau = tau, # structural parameters
    ztype = 2, #uniform covariates
    ctype=ctype,cparam=censor_rate,censor_min=censor_min,censor_max=censor_max,
    ncov = ncov,
    # M = M, #number of causes # NOTE WE CURRENTLY ONLY HAVE GDATA_CR SUPPORT 2 CAUSES b/c of Fine-Gray simulation setting
    mass_p = cause1prob,
    generate_failure_method = generate_failure_method,
    predHazardFn = predHazardFn,
    predPropensityFn = predPropensityFn, # list of predictor functions,
    priority_cause = priority_cause
  )
} else{
  list(
    N = n, tau = tau, tick = 0.01,  # structural parameters
    censor_rate = censor_rate,
    ncov = ncov, #M = M,
    predHazardFn = predHazardFn,
    predPropensityFn = predPropensityFn, # list of predictor functions
    hidden_data = TRUE, printFlag = FALSE
  )
}

nodesize = 5
mindeath = round(sqrt(c(nodesize)), 0)
tau = tau

arg.obs <- arg.czmk <- arg.csk <-
  arg.pmcr <- arg.aipwe <-
  arg.obs.no.censor <- arg.zom <- arg_list

# for NOW only: no censoring

# evaluation datasets dont have censoring
arg.czmk$ctype <- arg.csk$ctype <-
  arg.pmcr$ctype <- arg.aipwe$ctype <-
  arg.zom$ctype <- arg.obs.no.censor$ctype <- 99 # no censoring

arg.czmk$N <- arg.csk$N <-
  arg.pmcr$N <- arg.aipwe$N <-
  arg.obs.no.censor$N <- arg.zom$N <- n.eval

arg.czmk$crit.eval.cif <- arg.csk$crit.eval.cif <-
  arg.pmcr$crit.eval.cif <- arg.aipwe$crit.eval.cif <-
  arg.obs.no.censor$crit.eval.cif <- arg.zom$crit.eval.cif <- crit.eval.cif
arg.czmk$crit.eval.os <- arg.csk$crit.eval.os <-
  arg.pmcr$crit.eval.os <- arg.aipwe$crit.eval.os <-
  arg.obs.no.censor$crit.eval.os <- arg.zom$crit.eval.os <- crit.eval.os

#tau/2 is because of itrSurv::VerifySurvivalTime.R (line 53 etc) used in criticalValue.R
if (is.null(crit.value_phase2) || is.na(crit.value_phase2) || length(crit.value_phase2) == 0){
  t0.cif = tau/2;
  t0_pmcr = 0.1 #tau/2
  t0_aipwe = 0.1 # these need to be tuned...
} else{
  t0.cif = crit.value_phase2;
  t0_pmcr = 0.1 #crit.value_phase2
  t0_aipwe = 0.1 # this needs to be tuned..
}
if (is.null(crit.value_phase1) || is.na(crit.value_phase1) || length(crit.value_phase2) == 0){
  t0.os = tau/2;
} else{
  t0.os = crit.value_phase1;
}
arg.czmk$t0.cif <- arg.csk$t0.cif <-
  arg.pmcr$t0.cif <- arg.aipwe$t0.cif <-
  arg.obs.no.censor$t0.cif <- arg.zom$t0.cif <- t0.cif;
arg.czmk$t0.os <- arg.csk$t0.os <-
  arg.pmcr$t0.os <-arg.aipwe$t0.os <-
  arg.obs.no.censor$t0.os <- arg.zom$t0.os <- t0.os;

arg.czmk$evaluate <- arg.csk$evaluate <-
  arg.pmcr$evaluate <- arg.aipwe$evaluate <-
  arg.obs.no.censor$evaluate <- arg.zom$evaluate <- TRUE

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################

sim = 1

  cat("\n\n#################################")
  cat("\n######### Simulation ",sim, "#########")
  cat("\nEndpoint: ",endpoint)
  if (endpoint == "CR"){
    cat("\nNumber of CR: ",M)
  }
  cat("\n#################################\n")

  if (sim == 1){
    sim_count = sim
    prev_count = 0
  } else{
    sim_count = sim*2+prev_count
    prev_count = prev_count + 1
  }

  train_seed = sim*10000 + 123

  cat ("1. Data for Simulation",sim, ":",generate_failure_method,"\n")
  tt(1)

  # obs.data
  set.seed(train_seed + 0)
  u1_train = runif(n)
  u1_test = runif(n.eval)
  u2_train = runif(n)
  u2_test = runif(n.eval)
  u3_train = runif(n)
  u3_test = runif(n.eval)

  arg.obs$u1 = u1_train; arg.obs$u2 = u2_train; arg.obs$u3 = u3_train;
  arg.czmk$u1 <- arg.csk$u1 <- arg.pmcr$u1 <- arg.aipwe$u1 <- arg.zom$u1 <- arg.obs.no.censor$u1 <- u1_test
  arg.czmk$u2 <- arg.csk$u2 <- arg.pmcr$u2 <- arg.aipwe$u2 <- arg.zom$u2 <- arg.obs.no.censor$u2 <- u2_test
  arg.czmk$u3 <- arg.csk$u3 <- arg.pmcr$u3 <- arg.aipwe$u3 <- arg.zom$u3 <- arg.obs.no.censor$u3 <- u3_test
  # print(head(u1_train,1));print(head(u1_test,1));print(head(u2_train,1));print(head(u2_train,1))

  message("using train_seed for generating training data")
  set.seed(train_seed)
  obs.data <- do.call(gdata_CR, arg.obs);
  obs_1 <<- obs.data

  # View(obs.data)
  # print(head(obs_1$event.time))
  data.df <<- obs.data %>%
    mutate(D.0 = ifelse(status> 0,1,0), # event from any cause indicator
           D.1 = ifelse(status==1,1,0), # event from cause 1 indicator
           D.2 = ifelse(status==2,1,0), # event from cause 2 indicator
           obs_time = event.time,
           Trt = action) %>%
    dplyr::select(obs_time, Trt, status, D.0, D.1, D.2, contains("Z")) %>%
    arrange(obs_time)
  # dplyr::select(-c(Time_Censor, Time_Failure1, Time_Failure2, Time_Tau, obs_time_failureCR,indD))

  # create unique observed failure times for survival phase 1
  timePointsSurvival = data.df %>%
    filter(D.0 == 1) %>%
    dplyr::select(obs_time) %>%
    unlist(use.names = FALSE) %>%
    unique()
  timePointsEndpoint = timePointsSurvival # we can do this bc CR is subset of OS failure times
  # create unique observed endpoint times for endpoint phase 2

  cat("\n******************************\n")
  # estimation
  cat ("2. czmk for simulation", sim, "\n")

    cat ("  2. czmk - Policy estimation for Simulation",sim, ":",generate_failure_method,"\n")
    # new package itrSurv
    priority_vector <- c(0, priority_cause)
    # Use lapply to create a list of sprintf statements
    models_itr <- lapply(priority_vector, function(x) {
      paste0(sprintf("Surv(obs_time, D.%s) ~ ", ifelse(x == 0, "0", as.character(x))),
             paste(paste0("Z", 1:ncov, ""), collapse = " + ")) %>% as.formula
    })
    # print(crit.value_phase1)
    arg.czmk2 = list(data = data.df,
                     endPoint = "CR",
                     txName = paste("Trt"),
                     epName = "status",
                     models = models_itr,
                     timePointsSurvival = timePointsSurvival,
                     timePointsEndpoint = timePointsEndpoint,
                     tau = tau,
                     # timePoints = "uni",
                     # nTimes = 100,
                     criticalValue1 = criterion_phase1,
                     criticalValue2 = criterion_phase2,
                     evalTime = crit.value_phase1,
                     splitRule1 = "mean_surv",
                     splitRule2 = "gray_cr",
                     ERT = TRUE, uniformSplit = TRUE, replace = FALSE,
                     randomSplit = 0.2, nTree = 1,#300,
                     pooled = FALSE,
                     tol1 = tol1,
                     stratifiedSplit = 0.1)
    set.seed(train_seed + 1)
optimal.czmk <- do.call(itrSurv, c(arg.czmk2, list(mTry = sqrt(ncov),
                                                   nodeSize = nodesize,
                                                   minEvent = mindeath )))
# print(data.df)
