library(ggplot2)
library(survival)

survival_times <- c(0,8,10,12,15,20,25,
                    0,3,14,15,16,17,18)
status <- c(1, 1, 1, 1, 1, 1, 1,1,1,1,1,1,1,1)  # 1 for an event (death), 0 for censored
trt <- c(1,1,1,1,1,1,1,0,0,0,0,0,0,0)
# Create a data frame
data <- data.frame(time = survival_times, status = status, trt = trt)

# Create a survival object
surv_obj0 <- Surv(data$time[trt==0], data$status[trt==0])
surv_obj1 <- Surv(data$time[trt==1], data$status[trt==1])

# Fit the survival curve
surv_fit0 <- survfit(surv_obj0 ~ 1)
surv_fit1 <- survfit(surv_obj1 ~ 1)

upper0 = c(surv_fit0$surv[1], surv_fit0$surv[1:length(surv_fit0$surv)-1])
upper1 = c(surv_fit1$surv[1], surv_fit1$surv[1:length(surv_fit1$surv)-1])

# Convert survival data to a data frame for ggplot
surv_data0 <- data.frame(trt = as.character(0), time = surv_fit0$time, survival = surv_fit0$surv, upper = upper0, lower = 0)#surv_fit$lower)
surv_data1 <- data.frame(trt = as.character(1), time = surv_fit1$time, survival = surv_fit1$surv, upper = upper1, lower = 1)#surv_fit$lower)

suv_data = rbind(surv_data0, surv_data1, by = trt)
# Plot stepwise survival curve with shaded area
ggplot(suv_data, aes(x = time, y = survival, group = trt, color = trt)) +
  geom_step() +
  # geom_ribbon(aes(ymin = lower, ymax = upper0),
  #             alpha = 0.2) +
  labs(x = "Time",
       y = "Survival Probability",
       title = "Stepwise Survival Curve with Shaded Area") +
  theme_bw()
