
library(tidyverse)
library(dplyr)
library(randomForestSRC)
library(survival)

# Set a random seed for reproducibility
set.seed(123) 
N1 <- 5  # Replace with the desired number of subjects
tau = tau

## training sampling
unique_subjects <- unique(df_surv$ID)
train_prop = 0.6
# Select a proportion of subjects for training
train_subjects <- sample(unique_subjects, size = train_prop * length(unique_subjects))

# Split the df_surv based on the selected subjects
train_data <- df_surv %>%
  filter(ID %in% train_subjects) %>%
  dplyr::select(ID, Cov, Trt, obs_time, indD)

test_data <- df_surv %>%
  filter(!ID %in% train_subjects) %>%
  #create counterfactual treatment group
  mutate(cf_Trt = 1-Trt) %>%
  dplyr::select(ID, Cov, Trt, cf_Trt, obs_time, indD)

# test dataset but with og treatment group
test_data_og <- test_data %>%
  dplyr::select(ID, Cov, Trt, obs_time, indD)
# same dataset but now with couterfactual treatment group
test_data_cf <- test_data %>%
  dplyr::select(ID, Cov, cf_Trt, obs_time, indD) %>%
  mutate(Trt = cf_Trt) %>%
  dplyr::select(ID, Cov, Trt, obs_time, indD)

# two ways to do this: 1) do prediction using OG trt first, then do cf trt second, then plot combined results; 2) append the two datasets together (First part is OG treatment then repeat it but with alternative trt with ONE variable TRT.) (easier for prediciton b/c u do it once; create indicator for which is OG which is CF)
# do first way
# run random survival forest on training dataset (og treatment) for terminal event (death)
df_surv.out <- rfsrc(Surv(obs_time, indD) ~ Cov+Trt, data=train_data)
## standard prediction call for test data w/ og treatment
df_surv.ogtrt <- predict(df_surv.out, newdata = test_data_og, outcome = "train")
## standard predict call for test data w/ cf treatment
df_surv.cftrt <- predict(df_surv.out, newdata = test_data_cf, outcome = "train")

## check forest reproducibilility by comparing "test-cf" predicted survival
## curves to "test-og" predicted survival curves for the first 3 individuals
Time <- df_surv.out$time.interest

# AUC per treatment per person
# length(df_surv.ogtrt$survival[1,]) # this is a step function

AreaDiffAll = numeric(0)
mat = matrix(nrow = nrow(test_data), ncol = 9)
colnames(mat) = c("ID", "OG_AUS", "CF_AUS", "Trt", "CF_Trt", "Opt_AUS", "Opt_Trt", "Opt_Trt_Ind", "NonOpt_AUS")
for (subj in 1:nrow(test_data)){
  surv_og = df_surv.ogtrt$survival
  surv_cf = df_surv.cftrt$survival
  print(subj)
  Base = c(0, Time, tau)
  # Calculate the difference between adjacent elements of Base
  Base2 = diff(Base)
  # calculate AUC by getting area under first step (rectangle) + area under curve of next step
  Height_OG = c(1, surv_og[subj,])
  # print(head(Height_OG))
  Height_CF = c(1, surv_cf[subj,])
  # print(head(Height_CF))
  # Multiply the Height_OG values with the calculated differences of Time
  AreaOfRectangles_OG = Height_OG * Base2
  AreaOfRectangles_CF = Height_CF * Base2
  AUS_OG = sum(AreaOfRectangles_OG)
  AUS_CF = sum(AreaOfRectangles_CF)
  AUS_opt = max(AUS_OG, AUS_CF)
  AUS_nopt = min(AUS_OG, AUS_CF)
  # NOTE : FOR NOW THIS IS CODED FOR BINARY TREATMENT ONLY!!!!
  opt_trt = ifelse(AUS_opt == AUS_OG, test_data$Trt[subj], test_data$cf_Trt[subj])
  opt_trt_ind = ifelse(AUS_opt == AUS_OG, 1, 0)
  print(opt_trt)
  line = c(subj, AUS_OG, AUS_CF, test_data$Trt[subj], test_data$cf_Trt[subj], AUS_opt, opt_trt, opt_trt_ind, AUS_nopt)
  print(line)
  # AreaDiff = AUS_OG - AUS_CF
  # AreaDiffAll[subj] = AreaDiff # OG - CF
  # print(AreaDiff)
  mat[subj,] = line
}

# CALCULATE area for each person and for each treatment
# TAKE DIFFERENCE IN THE TWO AREAS # compare the two and see if it fits the criteria
eps0 = 0.1
(1-eps0)
mat = as.data.frame(mat)
mat_R = mat %>%
  mutate(Ratio = NonOpt_AUS/Opt_AUS) %>%
  #if ratio >= 0.9 then ratio_stop_ind = 1 (aka we stop b/c we know which treatment to pick; otherwise ind = 0 aka we keep going)
  mutate(Ratio_Stop_Ind = ifelse(Ratio >= (1-eps0), 1, 0)) #%>% filter(Ratio_Stop_Ind == 1)



# ONE PLOT WITH ALL INDIVIDUALS
# Open a larger graphics device
options(repr.plot.width = 10, repr.plot.height = 6)  # Adjust dimensions as needed
par(mfrow = c(1, 1))  # To create two plots side by side
# Plot the survival curves for the first three subjects in the training set
matplot(Time, t(df_surv.ogtrt$survival[1:N1,]), ylab = "Survival", col = 1, type = "l")
# Plot the survival curves for the first three subjects in the test set
matlines(Time, t(df_surv.cftrt$survival[1:N1,]), col = 2)
legend("topright", legend = c("OG", "CF"), col = c(1, 2), lty = 1)
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
  plot(Time, df_surv.ogtrt$survival[subject,], ylim = c(0, 1), type = "l",
       xlab = "Time", ylab = "Predicted Survival Probability", col = 1, main = paste("Subject", subject))
  
  # Add the survival curve for the current subject in the test set
  lines(Time, df_surv.cftrt$survival[subject,], col = 2)
  # Add a single legend for the entire page
  legend_text <- c("OG", "CF")
  legend("topright", legend = legend_text, col = c(1, 2), lty = 1)
  

  }


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
  plot(Time, df_surv.ogtrt$survival[subject,], ylim = c(0, 1), type = "l",
       xlab = "Time", ylab = "Predicted Survival Probability", col = 1, main = paste("Subject", subject))
  
  # Add the survival curve for the current subject in the test set
  lines(Time, df_surv.cftrt$survival[subject,], col = 2)
  
  # Add a legend to the current plot
  legend("topright", legend = c("OG", "CF"), col = c(1, 2), lty = 1)
}

# Reset the graphics device settings
par(mfrow = c(1, 1))
options(repr.plot.width = 6, repr.plot.height = 6)  # Reset dimensions to default
