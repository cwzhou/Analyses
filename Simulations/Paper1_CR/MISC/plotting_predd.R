library(ggplot2)

# Define the plotting function
plot_func <- function(list_name, id) {
  # Create a list_name frame for the first set of list_name using the specified id
  data1 <- data.frame(
    x = seq_along(list_name[["Func"]][[2]][,id]),
    y = list_name[["Func"]][[2]][,id],
    group = "Trt 1"
  )
  
  # Create a list_name frame for the second set of list_name using the specified id
  data2 <- data.frame(
    x = seq_along(list_name[["Func"]][[1]][,id]),
    y = list_name[["Func"]][[1]][,id],
    group = "Trt -1"
  )
  
  # Combine the list_name frames
  combined_data <- rbind(data1, data2)
  
  # Plot using ggplot2
  ggplot(combined_data, aes(x = x, y = y, color = group)) +
    geom_line(linewidth = 1) +
    labs(x = "X-axis Label", y = "Y-axis Label", 
         title = paste("Overlay of Two Plots for id =", id)) +
    theme_minimal() +
    scale_color_manual(values = c("Trt 1" = "red", "Trt -1" = "blue"))
}

# Example of how to call the function with id = 1
plot_func(predd_surv_czmk_eval, 1)
plot_func(predd_surv, 1)
