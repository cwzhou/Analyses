# ------------------------------------------------------------
## another example illustrating outcome = "test"
## unique way to check reproducibility of the forest
## ------------------------------------------------------------
library(tidyverse)
library(dplyr)
library(randomForestSRC)
library(survival)
## training step
N1 <- 10  # Replace with the desired number of subjects
set.seed(542899)
data(pbc, package = "randomForestSRC")
train <- sample(1:nrow(pbc), round(nrow(pbc) * 0.50))
pbc.out <- rfsrc(Surv(days, status) ~ ., data=pbc[train, ])
## standard prediction call
pbc.train <- predict(pbc.out, pbc[train, ], outcome = "train")
##non-standard predict call: overlays the test data on the grow forest
pbc.test <- predict(pbc.out, pbc[-train, ], outcome = "test")
## check forest reproducibilility by comparing "test" predicted survival
## curves to "train" predicted survival curves for the first 3 individuals
Time <- pbc.out$time.interest

# ONE PLOT WITH ALL INDIVIDUALS
# Open a larger graphics device
options(repr.plot.width = 10, repr.plot.height = 6)  # Adjust dimensions as needed
par(mfrow = c(1, 1))  # To create two plots side by side
# Plot the survival curves for the first three subjects in the training set
matplot(Time, t(pbc.train$survival[1:N1,]), ylab = "Survival", col = 1, type = "l")
# Plot the survival curves for the first three subjects in the test set
matlines(Time, t(pbc.test$survival[1:N1,]), col = 2)
legend("topright", legend = c("Training", "Testing"), col = c(1, 2), lty = 1)
# Reset the graphics device settings
par(mfrow = c(1, 1))
options(repr.plot.width = 6, repr.plot.height = 6)  # Reset dimensions to default






### ONE PLOT PER INDIVIDUAL

# Open a larger graphics device
options(repr.plot.width = 20, repr.plot.height = 50)  # Adjust dimensions as needed
n_subjects <- N1  # Number of subjects
n_cols <- 5  # Number of columns in the grid
n_rows <- ceiling(n_subjects / n_cols)  # Calculate the number of rows

par(mfrow = c(n_rows, n_cols))  # Create a grid of plots

# Iterate over each subject and create separate plots
for (subject in 1:n_subjects) {
  # Create a new plot for the current subject
  plot(Time, pbc.train$survival[subject,], ylim = c(0, 1), type = "l",
       xlab = "Time", ylab = "Predicted Survival Probability", col = 1, main = paste("Subject", subject))
  
  # Add the survival curve for the current subject in the test set
  lines(Time, pbc.test$survival[subject,], col = 2)
}

# Add a single legend for the entire page
legend_text <- c("Training", "Testing")
legend("topright", legend = legend_text, col = c(1, 2), lty = 1)

# Reset the graphics device settings
par(mfrow = c(1, 1))
options(repr.plot.width = 6, repr.plot.height = 6)  # Reset dimensions to default






# Open a larger graphics device
options(repr.plot.width = 10, repr.plot.height = 50)  # Adjust dimensions as needed
n_subjects <- N1  # Number of subjects
n_cols <- 5  # Number of columns in the grid
n_rows <- ceiling(n_subjects / n_cols)  # Calculate the number of rows

par(mfrow = c(n_rows, n_cols))  # Create a grid of plots

# Iterate over each subject and create separate plots
for (subject in 1:n_subjects) {
  # Create a new plot for the current subject
  plot(Time, pbc.train$survival[subject,], ylim = c(0, 1), type = "l",
       xlab = "Time", ylab = "Survival", col = 1, main = paste("Subject", subject))
  
  # Add the survival curve for the current subject in the test set
  lines(Time, pbc.test$survival[subject,], col = 2)
  
  # Add a legend to the current plot
  legend("topright", legend = c("Training", "Testing"), col = c(1, 2), lty = 1)
}

# Reset the graphics device settings
par(mfrow = c(1, 1))
options(repr.plot.width = 6, repr.plot.height = 6)  # Reset dimensions to default
