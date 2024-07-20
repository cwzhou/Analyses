time = data.df$event.time
event1 = data.df$D.1
event0 = data.df$D.0
event2 = data.df$D.2
Z1 = data.df$Z1
Z2 = data.df$Z2
fit1 <- coxph(Surv(time, event1) ~ Z1 + Z2)
summary(fit1)
fit0 <- coxph(Surv(time, event0) ~ Z1 + Z2)
summary(fit0)
fit2 <- coxph(Surv(time, event2) ~ Z1 + Z2)
summary(fit2)

sfit0 = survfit(fit0, 
               newdata = list(Z1 = 0 , Z2 = 0))
                 # Z1=rep_obs$Z1, 
                 #              Z2=rep_obs$Z2))
sfit1 = survfit(fit1, 
                newdata = list(Z1 = 0, Z2 = 0))
                  #             Z1=rep_obs$Z1, 
                  #              Z2=rep_obs$Z2))
sfit2 = survfit(fit2, 
                newdata = list(Z1=0, Z2 = 0))
                               #    rep_obs$Z1, 
                               # Z2=rep_obs$Z2))

tibble(
  time=sort(time),
  # estimated0 = sfit0$surv,
  # estimated1 = sfit1$surv,
  # estimated2 = sfit2$surv,
  truth1 = 1 - pexp(time, rate1[1]),
  truth0 = 1 - pexp(time, rate0[1]),
  truth2 = 1 - pexp(time, rate2[1])
) %>% 
  gather(which, surv, -time) %>% 
  ggplot(aes(time, surv, color = which)) + 
  geom_step(size=0.8, alpha = 0.5)
