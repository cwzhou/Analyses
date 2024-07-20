if (skip.AIPWE != TRUE){
  ### A4. AIPWE - 2022
  set.seed(cv)
  aif = aipwe.fit(data_list = train_aipwe,
                  pp.v = ncol(train)-6, #not including obs_time, trt, and the deltas
                  tau1 = as.numeric(t0_aipwe),
                  tune = c(0.001,0.01,0.5,1,
                           seq(0.1,400,
                               length.out=16))) #eta0-eta_{ncov}

  AIPWE.i <- aif
  err.AIPWE = class(AIPWE.i)[1] == "try-error"

} else{
  err.AIPWE = TRUE
}
