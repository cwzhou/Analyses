library(ggplot2)

# Define the plotting function with data, id, and timePointsSurvival as arguments
plot_func <- function(list_name, id, timePointsSurvival) {
  # Create a data frame for the first set of list_name using the specified id
  data1 <- data.frame(
    x = timePointsSurvival,  # Use timePointsSurvival for x-axis
    y = list_name[["Func"]][[2]][,id],
    group = "Func 2"
  )
  
  # Create a data frame for the second set of data using the specified id
  data2 <- data.frame(
    x = timePointsSurvival,  # Use timePointsSurvival for x-axis
    y = list_name[["Func"]][[1]][,id],
    group = "Func 1"
  )
  
  # Combine the data frames
  combined_data <- rbind(data1, data2)
  
  # Plot using ggplot2 with the updated linewidth aesthetic
  ggplot(combined_data, aes(x = x, y = y, color = group)) +
    geom_line(linewidth = 1) +
    labs(x = "Time", y = "Y-axis Label", 
         title = paste("Overlay of Two Plots for id =", id)) +
    theme_minimal() +
    scale_color_manual(values = c("Func 1" = "red", "Func 2" = "blue")) +
    scale_x_continuous(breaks = timePointsSurvival)  # Set x-axis labels at timePointsSurvival
}

# Example of how to call the function with a specific data list, id, and timePointsSurvival
plot_func(predd_surv_czmk_eval, 1, timePointsSurvival)


# Example of how to call the function with id = 1
plot_func(predd_surv_czmk_eval, 1, timePointsSurvival)
plot_func(predd_surv, 1)
