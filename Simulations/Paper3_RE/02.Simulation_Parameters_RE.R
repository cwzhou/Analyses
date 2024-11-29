# Title: Setting hyperparameters for Simulation Scripts
# Description: [Brief description of the purpose and objectives of the simulation]
# Author: Christina Zhou
# Date: 11.12.2024

### DATA TYPE ###
sim_data_type = "RE"

if (sim_data_type == "RE"){
  ##### Study Design Hyperparameters #####
  seed1 = 2024
  N = 80 # Number of subjects
  tau = 1 # stop time/end of study
  G = 3 # total gap times
  ##### Covariate-Related Hyperparameters #####
  num_A = 2 # Number of treatments
  lambda_0D = 1 # baseline hazard for disease
  lambda_0R = 1 # baseline hazard for recurrence
  beta_D = log(2) # parameters for disease for covariates
  beta_R = log(2) # parameters for recurrence for covariates
  omega_D = log(3) #parameters for disease for treatment
  omega_R = log(3) #parameters for treatment+covariate interaction
  ztype = 0 # covariate distribution
  zparam = 0.3 # covariate distribution parameter

  ## TREATMENT COVARIATE INTERACTION ##

  ##### Censoring Hyperparameters #####
  ctype = 1 # censoring distribution
  cparam = 2 # censoring distribution parameter
  ##### Gap Time Hyperparameters #####
  gaptype = 0 # failure vs gap time indicator
  gapparam1 = 0.1 #rho for failure
  gapparam2 = 0.7 #rho for gap times

  print(sprintf("Seed: %s", seed1))
  print(sprintf("Parameters: N=%s, tau=%s, G=%s, num_A=%s, lambda_0D=%s, lambda_0R=%s, beta_D=%s, beta_R=%s, omega_D=%s, omega_R=%s, ztype=%s, zparam=%s, ctype=%s, cparam=%s, gaptype=%s, gapparam1=%s, gapparam2=%s.",
                N, tau, G, num_A,
                lambda_0D, lambda_0R, beta_D, beta_R,
                omega_D, omega_R,
                ztype, zparam, ctype, cparam, gaptype, gapparam1, gapparam2))
}
if (sim_data_type == "CR"){
  ##### Study Design Hyperparameters #####
  seed1 = 2023
  N = 100 # Number of subjects
  tau = 5 # stop time/end of study
  ##### Covariate-Related Hyperparameters #####
  num_A = 2 # Number of treatments
  lambda_0D = 1 # baseline hazard for disease
  lambda_0R = 1 # baseline hazard for recurrence
  beta_D = log(2) # parameters for disease for covariates
  beta_R = log(2) # parameters for recurrence for covariates
  omega_D = log(3) #parameters for disease for treatment
  omega_R = log(3) #parameters for treatment+covariate interaction
  ztype = 0 # covariate distribution
  zparam = 0.3 # covariate distribution parameter

  ## TREATMENT COVARIATE INTERACTION ##
  ##### Censoring Hyperparameters #####
  ctype = 1 # censoring distribution
  cparam = 2 # censoring distribution parameter

  print(sprintf("Seed: %s", seed1))
  print(sprintf("Parameters: N=%s, tau=%s, G=%s, num_A=%s, lambda_0D=%s, lambda_0R=%s, beta_D=%s, beta_R=%s, omega_D=%s, omega_R=%s, ztype=%s, zparam=%s, ctype=%s, cparam=%s, gaptype=%s, gapparam1=%s, gapparam2=%s.",
                N, tau, G, num_A,
                lambda_0D, lambda_0R, beta_D, beta_R,
                omega_D, omega_R,
                ztype, zparam, ctype, cparam, gaptype, gapparam1, gapparam2))

}
