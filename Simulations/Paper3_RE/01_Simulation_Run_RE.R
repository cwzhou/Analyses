# Title: Simulating Survival Data with Recurrent Events and Terminal Death
# Description: [Brief description of the purpose and objectives of the simulation]
# Author: Christina Zhou
# Date: 07.02.2023

setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Github/Paper_1/Paper1_Sims/GenerateData/1_code/0_simulations/RE")

# Libraries --------------------------------------------------------------
source("02_Simulation_Libraries.R")

# Functions -----------------------------------------------------------------
# Define any necessary functions used in the simulation
source("02_Simulation_Functions.R")

# Hyperparameters  -----------------------------------------------------------------
# Define constants used in the simulation (e.g., parameters, initial values)
source("02_Simulation_Parameters.R")

# Main Simulation -----------------------------------------------------------
if (sim_data_type == "CR"){
  message("Competing Risks Survival Data Simulation")

  sim = gdata_CR(N=N,G=G,
              lambda_0D=lambda_0D,lambda_0R=lambda_0R,beta_R=beta_R,beta_D=beta_D,
              omega_D=omega_D,omega_R=omega_R,
              ztype=ztype,zparam=zparam,ctype=ctype,cparam=cparam,
              gaptype=gaptype,gapparam1=gapparam1,gapparam2=gapparam2,
              num_A=num_A,tau=tau,seed1=seed1)
  df = sim$dataset; View(df)

  name = sprintf("%s_%s",sim$name, sim_data_type); print(name)
  assign(name, df)
  kable(head(df), format = "latex", caption = "Competing Risks Dataset Example")
  # Output --------------------------------------------------------------------
  # Save or export simulation results, plots, or any other relevant output
  write_csv(df_recurr, file = sprintf("../../2_pipeline/CR/01_surv_%s.csv", name))
  write_csv(df_surv, file = sprintf("../../2_pipeline/CR/01_cr_%s.csv", name))

}


if (sim_data_type == "RE"){
  message("Recurrent Events Survival Data Simulation")
  sim = gdata_RE(N=N,G=G,
              lambda_0D=lambda_0D,lambda_0R=lambda_0R,beta_R=beta_R,beta_D=beta_D,
              omega_D=omega_D,omega_R=omega_R,
              ztype=ztype,zparam=zparam,ctype=ctype,cparam=cparam,
              gaptype=gaptype,gapparam1=gapparam1,gapparam2=gapparam2,
              num_A=num_A,tau=tau,seed1=seed1)
  df_recurr = sim$dataset_recurrent; View(df_recurr)
  df_surv = sim$dataset_survival; View(df_surv)
  name = sprintf("%s_%s",sim$name, sim_data_type); print(name)
  # update this to include name_surv
  assign(name, df_recurr)
  head_recurr = head(df_recurr);head_recurr
  kable(head_recurr, format = "latex", caption = "Recurrent Events Dataset Example")

  # Output --------------------------------------------------------------------
  # Save or export simulation results, plots, or any other relevant output
  write_csv(df_recurr, file = sprintf("../../2_pipeline/RE/01_surv_%s.csv", name))
  write_csv(df_surv, file = sprintf("../../2_pipeline/RE/01_recurr_%s.csv", name))

}


# End of script -------------------------------------------------------------
