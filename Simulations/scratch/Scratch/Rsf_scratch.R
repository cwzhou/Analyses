# install.packages("randomForestSRC")  # Install RSF package
library(randomForestSRC)            # Load RSF package
library(survival)

rsf_model <- rfsrc(Surv(time, status) ~ ., data = train_data)
individual_survival <- predict(rsf_model, test_data, individual = TRUE)
subject_id <- 1  # Replace with the desired subject ID
plot(individual_survival$surv[, subject_id], type = "l", xlab = "Time", ylab = "Survival Probability")


# Load example dataset
lung = survival::lung


# Fit a random survival forest model
rsf_model <- rfsrc(Surv(time, status) ~ ., data = lung)

# Get the predicted survival curves for each individual
individual_survival <- predict(rsf_model, lung, estimate.times = lung$time)

# Choose a subject ID (e.g., 1) to plot the individual survival curve
subject_id <- 1

# Plot the individual survival curve
plot(individual_survival$surv[, subject_id], type = "l", xlab = "Time", ylab = "Survival Probability")








# Split the data into training and test sets
set.seed(123)
train_index <- sample(1:nrow(lung), round(0.7*nrow(lung)), replace = FALSE)
train_data <- lung[train_index, ]
test_data <- lung[-train_index, ]


# Train the Random Survival Forest model
rsf_model <- rfsrc(Surv(time, status) ~ ., data = train_data)

# Predict individual survival curves
individual_survival <- predict(rsf_model, test_data, individual = TRUE)

# Plot individual survival curve for a specific subject
subject_id <- 3  # Replace with the desired subject ID
plot(individual_survival$surv[, subject_id], type = "l", xlab = "Time", ylab = "Survival Probability")

# Check if subject ID exists in test dataset
if (subject_id %in% test_data$index) {
  # Get the individual survival curve for the specific subject
  subject_survival <- individual_survival$surv[, subject_id]
  
  # Check if the survival curve has non-missing values
  if (any(!is.na(subject_survival))) {
    # Plot the individual survival curve
    plot(subject_survival, type = "l", xlab = "Time", ylab = "Survival Probability")
  } else {
    cat("Survival curve for subject", subject_id, "has missing values.")
  }
} else {
  cat("Subject", subject_id, "does not exist in the test dataset.")
}
