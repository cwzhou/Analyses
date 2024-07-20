cov = obs_1[,c("Z1","Z2","Z3")] %>% as.matrix()
beta1 = c(0,0,1) #c(1,1,1)
beta2 = c(0,0,-1)#c(1,-1,-1)

beta1star = cov %*% beta1
beta2star = cov %*% beta2


star = cbind(b1 = beta1star, b2 = beta2star, diff = ifelse(beta1star>beta2star,1,0),
             rbinom(nrow(beta1star),1,0.5))

head(star)
star

mean(star[,3])

