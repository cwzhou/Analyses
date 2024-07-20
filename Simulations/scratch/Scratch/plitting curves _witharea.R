# Load required libraries
library(ggplot2)
library(dplyr)

# Assuming your dataset is a data frame named 'df' with rows as treatments and columns as time points
# Replace this with your actual dataset

# Example data
set.seed(123)
df <- data.frame(
  treatment = factor(rep(1:2, each = 4)),
  time_point = rep(0:3, times = 2),
  my_values = c(10,3,2,1,15,11,7,6)  # Replace 'my_values' with the actual column name in your dataset
)

ggplot(df, aes(x = time_point,
               y = my_values,
               group = treatment,
               fill = treatment
               )) +
  geom_point(size = 0.5) +
  geom_ribbon(data=df,
              aes(x=time_point,
                  ymax=my_values),
              ymin=0,
              alpha=0.3) +
  scale_fill_manual(name='',
                    values=c("1" = "purple", "2" = "red")) +
  labs(x = "Time Points",
       y = "St",
       title = "Survival Curve with Area Under the Curve") +
  theme_minimal()
