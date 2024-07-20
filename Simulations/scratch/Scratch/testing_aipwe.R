arg.aipwe2 = list(data=data.df.aipwe, # b/c trt is 0/1
                  pp.v = c(2,1,2),#c(1,1,1,2),
                  tau1 = 0.15, #t0_aipwe, # t0 needs to be tuned
                  tune = c(0.001,0.01,0.5,1,
                           seq(0.1,400,
                               length.out=16))
) 
# set.seed(train_seed + 4)
aipwe.fit1 <- do.call(aipwe.fit, arg.aipwe2)
optimal.aipwe <- aipwe.fit1 #eta0-eta_{ncov}
arg.aipwe$policy = list()
arg.aipwe$policy[[1]] = "aipwe"
arg.aipwe$policy[[2]] <- if (!aipwe.error) optimal.aipwe
arg.aipwe$policy[[3]] = t0_aipwe
rm(optimal.aipwe); gc()

cat ("  \n 4. aipwe - Evaluation for Simulation",sim, ":",generate_failure_method,"\n")
set.seed(train_seed + 10)
if (!aipwe.error) {
  aipwe.data.rep <- do.call(gdata_CR, arg.aipwe)
  rep_aipwe <<- aipwe.data.rep
}
print(mean(rep_aipwe$OS_eval))
print(mean(rep_aipwe$CIF_eval))
table(rep_aipwe$action)