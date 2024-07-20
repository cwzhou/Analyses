## This is the script used to run the simulations ##
## Edit parameters in CR00.Simulation_Parameters.R script
## Plots in CR02.Simulation_Summary.R script

# call on parameters which calls on library
source("CR00.Simulation_Parameters.R")

# call on functions then runs simulations
source("CR00.Simulation_Body.R")

# source("plotting value functions.R")
# plot_values(combined_data, OS_eval)
# plot_values(combined_data, CIF_eval)

if (local == 1){
  beep()
  # data.df%>%group_by(action) %>% summarize(mean = mean(event.time))
  # data.df%>%group_by(status,action) %>% summarize(mean = mean(event.time))
  print(result)
}

message("End of CR01.Simulation_Run.R -- proceed to CR02.Simulation_Summary.R for creating plots.")
