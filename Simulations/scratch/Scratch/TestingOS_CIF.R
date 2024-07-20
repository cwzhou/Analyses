person = 30
time = tp_surv
beta1.hazard0 = c(0, -1.6,-1.2,0.5)
beta1.hazard1 = c(0, 0.3,-0.4,-0.2)
beta2.hazard0 = c(0, 1.1,-1.3,0.3)
beta2.hazard1 = c(0, -0.6,0.5,0.3)
z1 = rep_czmk[person,]$Z1
z2 = rep_czmk[person,]$Z2
z3 = rep_czmk[person,]$Z3

lamb_cause1_trt0 = as.numeric(exp(beta1.hazard0 %*% c(1,z1,z2,z3)))
lamb_cause2_trt0 = as.numeric(exp(beta2.hazard0 %*% c(1,z1,z2,z3)))
lamb_cause1_trt1 = as.numeric(exp(beta1.hazard1 %*% c(1,z1,z2,z3)))
lamb_cause2_trt1 = as.numeric(exp(beta2.hazard1 %*% c(1,z1,z2,z3)))

if (rep_czmk[person,]$action == 1){
} else{
}

# exp(-as.numeric((lamb_cause1_trt1 + lamb_cause2_trt1))*time)
# exp(-as.numeric((lamb_cause1_trt0 + lamb_cause2_trt0))*time)
os_simple_exp = function(lam_cause1, lam_cause2, time){
  exp(-as.numeric((lam_cause1 + lam_cause2))*time)
}

# taking the integration of OS curve from 0 to tau
overs_t1 = integrate(os_simple_exp, lower = 0 , upper = tau,
                     lam_cause1 = lamb_cause1_trt1, lam_cause2 = lamb_cause2_trt1)$value
overs_t0 = integrate(os_simple_exp, lower = 0 , upper = tau,
                     lam_cause1 = lamb_cause1_trt0, lam_cause2 = lamb_cause2_trt0)$value

# calculating the exact OS curve for a population with the same covariates as ID = person
# for time points time
os_trt1 = os_simple_exp(lamb_cause1_trt1, lamb_cause2_trt1, time)
os_trt0 = os_simple_exp(lamb_cause1_trt0, lamb_cause2_trt0, time)
data <- data.frame(x = time, y0 = os_trt0, y1 = os_trt1)

ggplot(data, aes(x = x)) +
  geom_line(aes(y = y0, color = "Trt-1"), linetype = "solid") +
  geom_line(aes(y = y1, color = "Trt1"), linetype = "dashed") +
  labs(x = "Time", y = "St",
       title = sprintf("True OS Curve for a population the same as ID=%s", person),
       color = "Treatment") +
  scale_color_manual(values = c("Trt-1" = "purple", "Trt1" = "light blue"),
                     labels = c("Trt -1", "Trt 1")) +
  theme_minimal()
