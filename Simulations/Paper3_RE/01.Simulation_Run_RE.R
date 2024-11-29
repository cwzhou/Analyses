# Title: Simulating Survival Data with Recurrent Events and Terminal Death
# Description: [Brief description of the purpose and objectives of the simulation]
# Author: Christina Zhou
# Date: 07.02.2023

# setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Github/Paper_1/Paper1_Sims/GenerateData/1_code/0_simulations/RE")
setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper3_RE")

set.seed(111)

# Libraries --------------------------------------------------------------
source("02.Simulation_Libraries_RE.R")

# Functions -----------------------------------------------------------------
source("02.Simulation_Functions_RE.R")

# Hyperparameters  -----------------------------------------------------------------
source("02.Simulation_Parameters_RE.R")

print(sim_data_type)
# Main Simulation -----------------------------------------------------------
if (sim_data_type == "RE"){
  message("Recurrent Events Survival Data Simulation")
  sim = gdata_RE(N=N,G=G,
              lambda_0D=lambda_0D,lambda_0R=lambda_0R,beta_R=beta_R,beta_D=beta_D,
              omega_D=omega_D,omega_R=omega_R,
              ztype=ztype,zparam=zparam,ctype=ctype,cparam=cparam,
              gaptype=gaptype,gapparam1=gapparam1,gapparam2=gapparam2,
              num_A=num_A,tau=tau,seed1=seed1)
  df_recurr = sim$dataset_recurrent; #View(df_recurr)
  df_surv = sim$dataset_survival; #View(df_surv)
  name = sprintf("%s_%s",sim$name, sim_data_type); print(name)
  # update this to include name_surv
  assign(name, df_recurr)
  head_recurr = head(df_recurr);head_recurr
  kable(head_recurr, format = "latex", caption = "Recurrent Events Dataset Example")

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

a0 = data_to_use %>% filter(Trt == 0); dim(a0)
a0_ids = a0 %>% pull(ID) %>% unique(); length(a0_ids)
a1 = data_to_use %>% filter(Trt == 1); dim(a1)
a1_ids = a1 %>% pull(ID) %>% unique(); length(a1_ids)

# need to source CHECKPackagesScripts.R where we want to check_RE = 1
train_seed = 2025
ncov = 2
nodesize = 3
mindeath = round(sqrt(c(nodesize)), 0)
arg.czmk2 = list(data = data_to_use,
                 endPoint = "RE",
                 idName = "ID",
                 epName = "IndR",
                 txName = "Trt",
                 models = models_RE,
                 timePointsSurvival = timePointsSurvival,
                 timePointsEndpoint = timePointsEndpoint,
                 tau = tau,
                 criticalValue1 = "mean",
                 criticalValue2 = "mean",
                 evalTime = 1,
                 splitRule1 = "mean_surv",
                 splitRule2 = "gray_re",
                 ERT = FALSE,
                 uniformSplit = TRUE,
                 replace = FALSE,
                 randomSplit = 0.2,
                 nTree = 2,
                 pooled = FALSE,
                 tol1 = c(0.1,0.1),
                 stratifiedSplit = 0.1)
set.seed(train_seed + 1)
optimal.czmk <- do.call(itrSurv, c(arg.czmk2, list(mTry = sqrt(ncov),
                                                   nodeSize = nodesize,
                                                   minEvent = mindeath )))





# End of script -------------------------------------------------------------
