# hist(obs_1$Z1)

for (act in obs_1$action){
  obs_1$action[act] = 
    ifelse(obs_1$action[act]==-1, 
           0, 
           obs_1$action[act])
  }

cox_res = coxph(Surv(event.time, 
           rep(1,nrow(obs_1))) ~ 
        action + Z1 + Z1*action, 
      data = obs_1)
print(cox_res)

tmp1 = cox_res$coefficients[1] + cox_res$coefficients[3]*obs_1$Z1
# hist(tmp1)

