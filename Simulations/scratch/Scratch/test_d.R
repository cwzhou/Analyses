source("F01.Simulation_Functions.R")

if (endpoint == "CR"){
  arg_list =   list(
    N = n, tau = tau, # structural parameters
    # censor_rate = censor_rate,
    ncov = ncov,
    # M = M, #number of causes # NOTE WE CURRENTLY ONLY HAVE GDATA_CR SUPPORT 2 CAUSES b/c of Fine-Gray simulation setting
    mass_p = cause1prob,
    censor_max = censor_max,
    predHazardFn = predHazardFn,
    predPropensityFn = predPropensityFn # list of predictor functions
  )
} else{
  list(
    N = n, tau = tau, tick = 0.01,  # structural parameters
    censor_rate = censor_rate,
    ncov = ncov, #M = M,
    predHazardFn = predHazardFn,
    predPropensityFn = predPropensityFn, # list of predictor functions
    hidden_data = TRUE, printFlag = FALSE
  )
}

arg.obs <- arg.czmk <- arg.csk <- arg.obs.no.censor <- 
  arg.gk.lm <- arg.gk.rf <-  arg.dw <-  
  arg.zom <- arg_list

# for NOW only: no censoring
arg.obs$ctype <- 99 # comment this out when we want training data to have censoring

arg.czmk$ctype <- arg.zom$ctype <- 
  arg.csk$ctype <- arg.obs.no.censor$ctype <- 99 # evaluation datasets dont have censoring

arg.czmk$N <- arg.csk$N <- arg.obs.no.censor$N <- 
  arg.gk.lm$N <- arg.gk.rf$N <- arg.dw$N <- 
  arg.zom$N <- n.eval

nodesize = 5
mindeath = round(sqrt(c(nodesize)), 0)
tau = tau

# arg.obs$seed1 <- sim*10000
# arg.czmk$seed1 <- arg.obs.no.censor$seed1 <- 
#   arg.zom$seed1 <- arg.csk$seed1 <- sim*10000 + 10 #eval seed (generate same covariates)

# obs.data
set.seed(sim*10000)
obs.data <- do.call(gdata_CR, arg.obs)
obs_1 <<- obs.data

# # # # observed policy value
# set.seed(sim*10000 + 10)
# obs.data.rep <- do.call(gdata_CR, arg.obs.no.censor)
# rep_obs <<- obs.data.rep

data = obs_1

### starting test_d.R
d = cbind(data, Xcause1, Xcause2) %>% 
  dplyr::select(subj.id, event.time, status, action, Xcause1,Xcause2) %>%
  mutate(Xcause1 = round(Xcause1,2),
         Xcause2 = round(Xcause2,2))

# treatment A (-1 or 0) vs treatment B (1)
dA = d %>% dplyr::select(subj.id, event.time, action) %>%
  filter(action == -1) %>% 
  dplyr::select(subj.id, event.time)
dB = d %>% dplyr::select(subj.id, event.time, action) %>%
  filter(action == 1) %>% 
  dplyr::select(subj.id, event.time)
message("treatment A: ",mean(dA$event.time),"\n") #-1 (or A or 0)
message("treatment B: ",mean(dB$event.time),"\n") # 1 or B
# we want to edit cause1 treatment coefficients until treatment -1 is close to treatment 1 but slightly lower
# then, in this setting, HC will always pick trt 1 even though
library(beepr)
beep(2)

# we want cause 2 to have higher survival since we want to minimize cause 1
for (cause in 1:2){
  if (cause == 1){
    xvar = "Xcause1"
  } else{
    xvar = "Xcause2"
  }
  for (trt in c(-1,1)){
    dd = d %>%
      filter(status == cause & action == trt) %>%
      dplyr::select(subj.id, xvar)
    if (trt != 1){
      trt1 = "A"
    } else{
      trt1 = "B"
    }
    message("cause", cause, " treatment", trt1,": ",mean(dd[,2]))
  }
  message("\n")
}