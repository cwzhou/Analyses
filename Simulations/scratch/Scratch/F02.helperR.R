St2 <- function(surv, time, t.eval, tau = Inf, exponential.tail = TRUE) {
  # tau is a placeholder for compatibility
  time = c(0, time)
  surv =  c(1, surv)
  max.t = max(time)
  min.s = min(surv)
  if (t.eval > max.t) {
    message("t.eval > max(time)")
    if (exponential.tail) {
      print("exponential.tail = TRUE")
      S.t = 1 - (1-min.s)^(t.eval / max.t)
    } else {
      print("S.t = min(surv) b.c exponential.tail = FALSE")
      S.t = min.s
    }
  } else {
    message("t.eval <= max(time)")
    S.t = approxfun(time, surv)(t.eval)
  }
  return(S.t)
}

weights_aric <- function(data, 
                         weight.formula.list, 
                         weight.formula.bin.list) {
  require(dplyr) 
  require(randomForest)
  require(randomForestSRC)
  
  n = dim(data)[1]
  Tx.nm = sapply(weight.formula.list, function(x) x[2] %>% 
                   as.character %>% gsub("(factor\\()(.*)(\\))", "\\2", .) )
  prop1.model =
    randomForest(weight.formula.list[[1]], data = data, )
  propensity1 = predict(prop1.model, data, type = "prob")                   # n x n.levels matrix
  if (is.null(dim(propensity1))) {
    # in case it is binary
    propensity1 = data.frame(1 - propensity1, propensity1)
    names(propensity1) = 0:1
  }
  missing1 = is.na(propensity1[, 1]) & !is.na(data$D.0) #data$d.1             # Bookkeeping covariate missing subjects (neither censored nor failed)
  propensity1 = propensity1[cbind(rownames(data), data[, Tx.nm[1]] %>% as.character)]  # A vector of predicted pi(A|X) with NA's
  propensity1 = ifelse(is.na(propensity1), 1, propensity1)             # Replace NA's with 1's.
  propensity  = propensity1
  
  # IPCW
  Sc.hat1 <- rfsrc(Surv(time, 1 - D.0) ~ #is D.0 correct here and in line 37?
                     Age + Gender + Black + GC, 
                   data = data)
  Sc.hat1 <- sapply(1:dim(data)[1], #for each person
                    function(i) St2(Sc.hat1$survival[i, ], 
                                    Sc.hat1$time.interest, #what is time.interest
                                    t.eval = data$time[i])) #data$V.1
  
  Sc.hat <- Sc.hat1
  # clipping
  Sc.hat <- Sc.hat %>% pmax(0.05) %>% pmin(0.95)
  weight.censor = data$delta/Sc.hat # what is delta
  
  p1 <<- propensity1
  c1 <<- Sc.hat1
  
  return(list(propensity = propensity, 
              weight.censor = weight.censor))
}

