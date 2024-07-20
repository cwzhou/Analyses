library(ggplot2)
# combined_data <- rbind(
#   transform(rep_czmk, method = "CZMK"),
#   transform(rep_csk, method = "CSK"),
#   transform(rep_pmcr, method = "PMCR"),
#   transform(rep_zom, method = "ZOM"),
#   transform(rep_obs, method = "OBS")
# )
combined_data <- NULL  # Initialize combined_data to be NULL

# Check if each dataset exists before including it in combined_data
if (exists("rep_czmk")) {
  combined_data <- rbind(combined_data, transform(rep_czmk, method = "CZMK"))
}
if (exists("rep_csk")) {
  combined_data <- rbind(combined_data, transform(rep_csk, method = "CSK"))
}
if (exists("rep_pmcr")) {
  combined_data <- rbind(combined_data, transform(rep_pmcr, method = "PMCR"))
}
if (exists("rep_zom")) {
  combined_data <- rbind(combined_data, transform(rep_zom, method = "ZOM"))
}
if (exists("rep_obs")) {
  combined_data <- rbind(combined_data, transform(rep_obs, method = "OBS"))
}
combined_data$method = factor(combined_data$method, 
                              levels = c("CZMK", "CSK", "PMCR", "ZOM", "OBS"))

plot_values <- function(data, eval) {
  eval_col <- enquo(eval)  # Convert the input column name to a quosure
  col_name <- quo_name(eval_col)  # Extract the column name
  if (grepl("cif", col_name, ignore.case = TRUE)){
    title = "Value Function: Priority Cause CIF"
  } else{
    title = "Value Function: Overall Survival"
  }
  ggplot(data, aes(x = factor(action), y = !!eval_col, color = factor(subj.id))) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.5) +
    theme_minimal() +
    theme(legend.position = "top") + 
    facet_grid(~ method) + 
    labs(x = "Treatment", 
         y = "Value Function", 
         title = title,
         color = "Simulation")
  
}


plot_values(combined_data, OS_eval)
plot_values(combined_data, CIF_eval)
print(result)
