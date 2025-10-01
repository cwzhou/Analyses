# setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper1_CR")

# things to note:
# - censoring for training; no censoring for testing; truncated at tau for both (differs by sim generate failure setting)
# right now, in cr01 t0_pmcr is set to 0.2 regardless of other parameters.
# local = 1 # uncomment for CR02.Simulation_Summary.R plots # comment for running sims
local = 1 # local = 0 for cluster
parallel = 0 # parallel = 0 for NOT parallel code
revision = 0 #jasa revision round 1
# revision = 1 uses CR02.Simulation_Summary_Revision.R
# revision = 0 uses CR02.Simulation_Summary.R

init_seed = 2025 #init_seed = 353 was used for revision = 1

# # below is only needed if running this script directly. comment out if running CR01.Simulation_Run.R script.
# # uncomment if running alone (aka for CR02.Simulation_Summary.R
# if (local == 1){
#    setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper1_CR")
#  } else{
#    setwd("/nas/longleaf/home/cwzhou/Dissertation/Analyses/Simulations/Paper1_CR")
#  }

#### libraries and functions
source("F01.Simulation_Functions.R") # calls libraries

# sbatch -p general -N 1 --mem 15G -n 1 -t 6-11:00:00 --mail-type=end --mail-user=cwzhou@email.unc.edu --wrap="Rscript CR01.Simulation_Run.R"

savingrds = TRUE
# date_folder = "2024-09-09" # "2024-08-31" #Sys.Date() 
#10 and 20 are local for revision = 1 and 100 sims; 2025-07-21 is cluster for revision = 1 1000 sims, 10000 neval 
#"2025-02-10" this is the original submission 
#"2024-09-13" this is an old one

if (revision == 1){
  # only thing you change
  multiplier <- 25  # e.g. 1 → first 50 sims, 2 → next 50 sims, etc.
  
  # fixed settings
  n.eval <- 10000
  n.sim  <- 50     # number of sims per block
  step   <- 50     # step size for n.sim_start calculation
  
  # auto-calculated
  start_date <- as.Date("2025-09-01")  # baseline date for multiplier = 1
  date_folder <- as.character(start_date + (multiplier - 1)) 
  
  n.sim_start <- (multiplier - 1) * step + 1
  n.sim_end   <- n.sim_start + n.sim - 1
  
  # check
  sim_deets = list(
    multiplier   = multiplier,
    date_folder  = date_folder,
    n.sim_start  = n.sim_start,
    n.sim_end    = n.sim_end,
    n.eval = n.eval
  ); sim_deets
  
  # 
  
  # calculate number of multipliers needed to reach 1000 sims
  max_multiplier <- ceiling(1000 / n.sim)
  
  # create table
  sim_schedule <- data.frame(
    multiplier   = 1:max_multiplier
  ) %>%
    mutate(
      date_folder = as.character(start_date + (multiplier - 1)),
      n.sim_start = (multiplier - 1) * step + 1,
      n.sim_end   = n.sim_start + n.sim - 1
    )
  sim_schedule
  
  # n.eval = 10000 #10000 #n.eval = 10000
  # date_folder = "2025-08-09"; 
  # # n.sim = 100 #500
  # n.sim_start = 701
  # # "2025-08-02" is 1-100 sims;
  # # "2025-08-03" is 101-200 sims;
  # # "2025-08-04" is 201-300 sims;
  # # "2025-08-05" is 301-400 sims;
  # # "2025-08-06" is 401-500 sims;
  # # "2025-08-07" is 501-600 sims;
  # # "2025-08-08" is 601-700 sims;
  # # "2025-08-09" is 701-800 sims;
  # # "2025-08-10" is 801-900 sims;
  # # "2025-08-11" is 901-1000 sims;
  # n.sim = 100
  # n.sim_end = n.sim_start - 1 + n.sim # n.sim
  
} else{
  n.eval = 1000
  n.sim = 100
  date_folder = "2025-09-22"
  n.sim_start = 1
  n.sim_end = n.sim_start - 1 + n.sim
  
  # check
  sim_deets = list(
    date_folder  = date_folder,
    n.sim_start  = n.sim_start,
    n.sim_end    = n.sim_end,
    n.eval = n.eval
  ); sim_deets
  
}

mean_tol1 = c(0.1,0) #c(0.07,0) # this is for differences in years so we don't want it to be too big
prob_tol1 = c(0.15, 0.01) # IGNORE THIS, WE DONT USE, but keep in code since fortran isn't updated to ignore
combo_tol1 = c(mean_tol1[1], prob_tol1[1], mean_tol1[2], prob_tol1[2])
generate_failure_method = c("simple_exp","fine_gray") 
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
skip_method <- c(!TRUE, !TRUE, !TRUE, !TRUE, !TRUE, !TRUE);
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
  if (revision == 1){
    arg = c(1, 2, 1, 1, 1, 1, 1, 1, 1) # for revision like rda
  }
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
             large.sample.size = list(n = 1000))

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
    censor_max = 0.1*tau.1, #0.7*tau.1, # lower is more censoring for unif
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
                        censor_rate = 5.3 #3.3 # higher is more censoring for exp
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
    
    if (revision == 1){
      # currently using runif covariates
      # we want difference between cause1 and cause2 (csk);
      # also want difference between 0 and 1 within cause1 (zom)
      
      sprintf("Adding a simulation study reflecting RDA for %s setting", default$generate_failure_method)
      # temporary while working on testing dataset
      betas <- list(
        beta1 = list(
          beta1.hazard0 = 
            c(0,
              2.1, 3.8, 1.3, 1.3, 2.1, 0.2,
              2.1, 2.1, 1, 1.3, 1.9, 2.1,
              1.2, 3.1, 2.1, 1, 2.3, 0.4,
              2.3, 0.1, 1.3, 2.4, 1.1, 2.7,
              0.1, 2.3, 3.4, 3.1, 0.1, 0.3),
          # c(0,
          #                 -1.30, -3.50, -2.3, -1.4, -2.1, -1.41,
          #                 -0.41, 2.50, -3.1, -2.21, -1.4, -3.1,
          #                 0.1, -1.90, 1.1, -2.1, -1.3, -1.4,
          #                 -1.1, -0.70, -1.60, -1.3, -1.90, -1.50,
          #                 -2.1, -1.1, -0.51, -0.31, -0.21, -0.21),
          
          beta1.hazard1 = 
            c(0,
              2.5, 3.8, 1.3, 1.3, 2.1, 0.2,
              2.3, 2.1, 1.3, 1.3, 1.9, 2.3,
              1.5, 3.1, 2.1, 1.1, 2.3, 0.4,
              2.7, 0.1, 1.3, 2.4, 1.1, 3.2,
              0.1, 2.3, 3.4, 3.1, 0.4, 0.3),
          # 3.3, 3.8, 3.3, 5.31, 1.1, 2.21,
          # 4.2, 2.1, 4.3, 2, 3.1, 4.3,
          # 2, 5.1, 2.1, 3.1, 3.2, 4.2,
          # 2, 4.8, 2, 5.1, 4.1, 2.2,
          # 1.51, 1.3, 6.4, 7.3, 1.41, 3.1),
          
          # c(0,
          #                 2.71, 0.61, 1.1, -2.61, 0.81, 1,
          #                 -0.81, 2, 1.5, -1.1, 3.21, 1.71,
          #                 3.1, 3.3, 2.1, -0.91, -2.21, 1.81,
          #                 2.21, -1.1, 0.5, 1.1, 2.1, -0.21,
          #                 3.1, 0.91, 0.71, 2.1, 1.1, 0.41),
          beta2.hazard0 = 
            c(0,
              2.3, 1.1, 0.31, 2.1, 0.31, 1,
              0.1, 1, 1, 1, 1, 1,
              1, 1, 1, 1, 1, 1,
              1, 1, 1, 1, 1, 1,
              1, 1, 1, 1, 1, 1),
          # c(0,
          #                 -1.1, 1.1, 0.41, 1.41, 0.51, -1.31,
          #                 0.31, 1.1, -0.51, 0.71, 0.31, -2.1,
          #                 -0.51, -0.81, 0.1, 0.1, -3.1, -1.7,
          #                 -0.41, 0.1, -3.1, -0.81, 1.6, 1.4,
          # 1.51, -0.31, -2.1, 1.1, -0.31, -2.1), #c(0,-0.1,-0.2),
          # higher beta2.hazard1 means more cause 2 because quicker failure_t2 in obs.data
          beta2.hazard1 = 
            c(0,
              1.8, 1.4, 1.1, 1.1, 2.1, 3.1,
              1.2, 1, 1, 1, 1, 1,
              1, 1, 1, 1, 1, 1,
              1, 1, 1, 1, 1, 1,
              1, 1, 1, 1, 1, 1)
          # c(0,
          #                 -2.1, 0.1, -3.8, -0.61, -3.1, -2.1,
          #                 -1.41, -2.1, -1.1, -4.1, -0.41, -0.31,
          #                 -1.1, -3.1, 0.41, -1.61, -2.1, -1.1,
          #                 1.1, -0.21, 0.61, -3.1, -2.41, -3.1,
          #                 0.41, 2.81, 0.91, 0.81, -0.31, 0.41)
          # c(0,
          # 1, 1, 1, 1, 1, 1,
          # 1, 1, 1, 1, 1, 1,
          # 1, 1, 1, 1, 1, 1,
          # 1, 1, 1, 1, 1, 1,
          # 1, 1, 1, 1, 1, 1)
        )
      )
      
      #     betas <- list(
      #       beta1 = list(
      #         beta1.hazard0 = c(0,-0.5,-1.2,1.5,0.1,0.6), #c(0, -0.5,-0.5,-0.5,-0.3,0.3), # (int), covariate (1~6)
      #         beta1.hazard1 = c(0,0.3,-0.4,-0.2,0.9,1.6), #c(0, -0.1,0.6,-0.9,0.5,-0.3),
      #         beta2.hazard0 = c(0,-0.6,1.4,1.5,2,1), #c(0, 0.1,0.3,-0.6,0.3,0.8),
      #         beta2.hazard1 = c(0,1.5,-2.1,-0.8,0.1,1)) #c(0, -0.3,-0.3,-0.2,-0.2,0.2))
      # )
    }else{
      # currently using runif covariates
      
      # we want difference between cause1 and cause2 (csk);
      # also want difference between 0 and 1 within cause1 (zom)
      
      # this works when tau = 2
      # N = 300 works better than N = 1000
      betas <- list(
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
    }
  } else if (default$generate_failure_method == "fine_gray"){
    # print("fine-gray coeff")
    # currently using runif covariates
    if (revision == 1){
      sprintf("Adding a simulation study reflecting RDA for %s setting", default$generate_failure_method)
      
      # betas <- list(
      #   beta1 = list(
      #     beta1.hazard0 = c(0,0.1,0.3),#c(0,-1,-1.4),
      #     beta1.hazard1 = c(0,-1.8,-1.5),#c(0,0.8,0.7),
      #     beta2.hazard0 = c(0,-1.1,-0.3),#c(0,-0.2,1.2), #c(0,0.2,0.8),
      #     beta2.hazard1 = c(0,-0.2,1.2)),#c(0,-0.3,-2)),
      #   beta2 = list(
      #     beta1.hazard0 = c(0,1,1.3,-0.2,0.6,0.3,2.1,0.5,-0.3,1.4,3),
      #     beta1.hazard1 = c(0,-1.2,-0.1,-0.5,-1.2,-2.2,1.2,-0.1,0.5,-0.2,-1.2),
      #     beta2.hazard0 = c(0,-0.5,-0.5,-1.4,-1.1,0.8,1.4,-2.4,-1.1,-0.8,-0.4), #c(0,-0.1,-0.2),
      #     beta2.hazard1 = c(0,0.4,1.3,1.6,1.2,1.1,-0.4,-0.3,0.6,0.2,1.1))
      # )
      
      betas <- list(
        beta1 = list(
          beta1.hazard0 = c(0,1,1.3,-0.2,0.6,0.3,
                            0.1,0.5,-0.3,1.4,3, 
                            rep(0, 31-11)),
          beta1.hazard1 = c(0,1.2,-0.1,-0.5,1.2,-2.2,
                            1.2,1.1,0.5,-0.2,-1.2, 
                            rep(0, 31-11)),
          beta2.hazard0 = c(0,-0.5,-0.5,-1.4,-1.1,0.8,
                            1.4,-2.4,-1.1,-0.8,-0.4, 
                            rep(0, 31-11)),
          beta2.hazard1 = c(0,0.4,1.3,1.6,1.2,1.1,
                            -0.4,-0.3,0.6,0.2,1.1, 
                            rep(0, 31-11))
          
          # beta1.hazard0 = c(0,
          #                   0.1, 0.3, 0.01, -0.15, -0.13, -0.17,
          #                   0.16, -1.15, -3.27, 1.16, 0.01, 0.01,
          #                   0.14, -2.13, 0.3, 1.18, 0.01, 0.01,
          #                   0.18, 0.01, -2.2, 0.54, 0, 0,
          #                   0.01, 0.01, 0.01, 0.38, 0.15, 0.12),
          # beta1.hazard1 = c(0,
          #                   -1.8, -1.5, 1.31, -0.13, 0.01, -0.17,
          #                   0.1, 0.13, 1.71, -0.3, 0, 2.1,
          #                   0, -0.16, -0.4, 0.1, 0.1, 0.11,
          #                   0.19, 0, 0.1, 1.33, 0.4, 0.2,
          #                   0.18, 0.12, 0.01, 0, 0, -0.1),
          # beta2.hazard0 = c(0,
          #                   0.15, 0.32, 2.19, 0.13, 0.14, 0.18,
          #                   -1.18, -0.13, 0.16, 0.01, 0.25, 0.19,
          #                   -0.36, 0.34, 0.32, -1.73, 0.04, 0.15,
          #                   -0.27, -0.14, -0.13, 1.21, 0.02, 0.84,
          #                   -0.14, 1.33, 0.12, 0.25, 0.01, 0.71), #c(0,-0.1,-0.2),
          # beta2.hazard1 = c(0,
          #                   -0.23,  0.35, 0.77, -1.38, 0.36, -1.98,
          #                   -0.14, -0.81, -0.52, 0.13, 0.17, 0.15,
          #                    0.31,  0.37, -0.17, 0.18, 0.12, -0.14,
          #                   -0.22,  0.83, 0.54, -1.41, 0.16, 0.01,
          #                   -0.33, -0.55, 0.62, 0.1, -0.32, 0.01)
          # #c(0,
          #                   # 1, 1, 1, 1, 1, 1,
          #                   # 1, 1, 1, 1, 1, 1,
          #                   # 1, 1, 1, 1, 1, 1,
          #                   # 1, 1, 1, 1, 1, 1,
          #                   # 1, 1, 1, 1, 1, 1)
        )
      )
      
      # betas <- list(
      #   beta1 = list(
      #     beta1.hazard0 = c(0,-0.5,-1.2,1.5,0.1,0.6), #c(0, -0.5,-0.5,-0.5,-0.3,0.3), # (int), covariate (1~6)
      #     beta1.hazard1 = c(0,0.3,-0.4,-0.2,0.9,1.6), #c(0, -0.1,0.6,-0.9,0.5,-0.3),
      #     beta2.hazard0 = c(0,-0.6,1.4,1.5,2,1), #c(0, 0.1,0.3,-0.6,0.3,0.8),
      #     beta2.hazard1 = c(0,1.5,-2.1,-0.8,0.1,1)) #c(0, -0.3,-0.3,-0.2,-0.2,0.2))
      # )
      
    } else{
      # temporary while working on testing dataset
      betas <- list(
        beta1 = list(
          beta1.hazard0 = c(0,0.1,0.3),#c(0,-1,-1.4),
          beta1.hazard1 = c(0,-1.8,-1.5),#c(0,0.8,0.7),
          beta2.hazard0 = c(0,-1.1,-0.3),#c(0,-0.2,1.2), #c(0,0.2,0.8),
          beta2.hazard1 = c(0,-0.2,1.2)),#c(0,-0.3,-2)),
        beta2 = list(
          beta1.hazard0 = c(0,1,1.3,-0.2,0.6,0.3,2.1,0.5,-0.3,1.4,3),
          beta1.hazard1 = c(0,-1.2,-0.1,-0.5,-1.2,-2.2,1.2,-0.1,0.5,-0.2,-1.2),
          beta2.hazard0 = c(0,-0.5,-0.5,-1.4,-1.1,0.8,1.4,-2.4,-1.1,-0.8,-0.4), #c(0,-0.1,-0.2),
          beta2.hazard1 = c(0,0.4,1.3,1.6,1.2,1.1,-0.4,-0.3,0.6,0.2,1.1))
      )
    }
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


dir_rds <- file.path("./output", generate_failure_method, date_folder)
dir_fig <- sub("^output", "figure", dir_rds)
print(dir_rds)
print(dir_fig)

if (savingrds) {
  # Create both output/ and figure/ directories with subfolders
  dirs_to_create <- c(
    file.path("output", generate_failure_method),
    file.path("figure", generate_failure_method)
  )
  
  for (d in dirs_to_create) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
}

if (local == 0){
  dir_rds_tmp = sprintf("/work/users/c/w/cwzhou/Proj1/output/%s/%s",
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

if (revision == 1){
  # check
  list(
    multiplier   = multiplier,
    date_folder  = date_folder,
    n.sim_start  = n.sim_start,
    n.sim_end    = n.sim_end
  )
}
