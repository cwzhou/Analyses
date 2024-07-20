library(caret); library(purrr); library(dplyr); library(survival)
library(DTRreg); library(randomForestSRC); library(dtrSurv)
# source("LF.helper.R")    # All the library and source files
# source("F00.generic.R")
# source("F21.GK1.R")
# source("F21.GK2.R")

evaluate <- function(testY,
                     weight, 
                     criterion,
                     tau) {
  if (criterion[1] == "mean") {
    mean(pmin(tau, testY) * weight)/ mean(weight)
  } else {
    mean(as.numeric(testY >= as.numeric(criterion[2])) * weight)/ mean(weight)
  }
}

getValue <- function(test, 
                     actual, 
                     estimated, 
                     propensity,
                     weight.censor,
                     criterion, 
                     tau) {
  weight = 
    apply(actual == estimated, 1, all2) %>%   
    # When there is at least one NA,
    # 1. all() does not return TRUE                     all(c(NA, NA, T)) = NA; all(c(T, T, T)) = TRUE
    # 2. all() returns NA if there is no FALSE          all(c(NA, NA, F, T)) = FALSE;  all(c(NA, NA)) = NA
    # 3. all() returns FALSE if there is at least one FALSE
    # When the second stage is not available in test set, the NA-match cases should still be counted. NA => 1.
    {ifelse(is.na(.), 1, as.numeric(.))} %>% 
    "/" (propensity)
  weight = weight * weight.censor
  evaluate(test[, "obs_time"], #test[, "T.0"] 
           weight = weight, 
           criterion = criterion, 
           tau = tau)
}

### 0. parameters
Tx.nm = "Trt"
# Tx.bin.nm = "AC" # For DW which does not admit more than two Tx arms.
Tx.nm.list = Tx.nm #paste(Tx.nm, 1:2, sep = ".")
imp = 1

## time points
tau = 40 #2700
timepoints = seq(0, sqrt(tau), length.out = 10)^2

## other parameters
nodesize = 10 #50
mindeath = round(sqrt(c(nodesize)), 0)
Ntree = 5 #300
ert = TRUE; rs = 0.2 # randomSplit = 0.2


### 0.2 all_cause / criterion
for (all_cause in c(FALSE, TRUE)) {
  all_cause_nm = if (all_cause) "_allcause" else ""
  for (value.criterion in c("mean", "surv.mean")) {
    if (value.criterion[1] != "mean") {
      value.criterion[2] = 10 # 2200 = six-year survival
      rule = "logrank"
    } else {
      rule = "mean"
    }
    
    cat("all_cause = ", all_cause,
        " value.criterion =",value.criterion[1],
        " value.s =", value.criterion[2], "tau =", tau, "\n")
    
    ## PM models
    modelPM = "Surv(obs_time, D.0) ~ age + hgb"
    form.CSK <- modelPM %>% as.formula()
    # form.GK <- modelPM %>% gsub("\\.%d", "", .) %>% gsub("\\(V, d\\)", "(V, delta)", .) %>% as.formula
    print(form.CSK)
    
    ## propensity models
    modelPr = "factor(Trt) ~ age + hgb"
    form.weight <- list(modelPr %>% gsub("%d", 1, .) %>% 
                          as.formula, modelPr %>% gsub("%d", 2, .) %>% 
                          as.formula)

    modelPrbin = "factor(AC.%d) ~ AA.%d + AS.%d + PREVHF01 + PRVCHD05 + hf + CENTER + gender + race + age + cig_yrs + bmi.%d + wth.%d + drink.%d + hypert.%d + glucose.%d + smoke.%d + hdl.%d"
    form.weight.bin <- list(modelPrbin %>% gsub("%d", 1, .) %>% as.formula, modelPrbin %>% gsub("%d", 2, .) %>% as.formula)
    print(form.weight)

    # 
    # ### 1. data preprocessing
    # dat.aric.list <- readRDS(sprintf("../data_ARIC/03.aric.comp%s.rds", all_cause_nm)) ###TBD
    # dat.aric = dat.aric.list[[imp]] %>%   # First imputed dataset.
    #   mutate(T.0 = pmin(T.0, tau), V.1 = pmin(V.1, tau), V.2 = T.0 - V.1, # Treat those administratively censored as failures for RMST.
    #          last.fu = ifelse(V.1 == tau, 1, last.fu),
    #          d.1 = ifelse(T.0 == tau & last.fu == 1, 1, d.1), d.2 = ifelse(T.0 == tau, ifelse(last.fu == 1, NA, 1), d.2),
    #          delta = ifelse(T.0 == tau, 1, delta)) %>% 
    #   mutate(ACAS.1 = ifelse(is.na(AC.1), NA, paste0(AC.1, AS.1)) %>% factor,
    #          ACAS.2 = ifelse(is.na(AC.2), NA, paste0(AC.2, AS.2)) %>% factor)
    # 
    # mean(dat.aric$d.2==0, na.rm = TRUE)
    # 
    lvls = lapply(Tx.nm.list, function(x) levels(data.df[, x]))
    for (x in 1:length(lvls)) {if (is.null(lvls[[x]])) lvls[[x]] = 0:1}

    ## cross-validation
    K = 2
    set.seed(100, kind="Mersenne-Twister")
    cv.insample  = createDataPartition(data.df[, "Trt"], p = .8, list = FALSE, times = K)
    cv.outsample = map_dfc(1:K, ~which(!(1:dim(data.df)[1] %in% cv.insample[, .]))) %>% as.matrix
    
    # skeleton
    values <- matrix(NA, K, 7, 
                     dimnames = list(1:K, c("CSK", "ZOM", "observed", 
                                            "cens", "cens.cause1", "cens.cause2",
                                            "ns.CSK")))
    
    # common parameters
    stg = 1  # stages
    
    # common file name
    nm = paste0(value.criterion[1], "_", 
                round(as.numeric(value.criterion[2]), 1), "_", 
                rule, "Split_", "tau_", tau, all_cause_nm, "_imp", imp)
    fnm = paste0("output/", Sys.Date(), "/") # folder name
    if (!dir.exists("output")) dir.create("output")
    if (!dir.exists(fnm)) dir.create(fnm)
    rds = paste0("_", Sys.Date(), ".rds")
    
    # tab = data.frame(s1 = rep(NA, K), 
    #                  s2 = NA, 
    #                  total = NA)
    
    for (cv in 1:K) {
      cat(cv, "th cv.\n")
      set.seed(cv)
      in.cv = cv.insample[, cv]   # insample index
      out.cv = cv.outsample[, cv] # outsample index
      train = data.df[in.cv, ]
      test  = data.df[out.cv,]
      values[cv, "cens"] = mean(train$D.0==0)
      values[cv, "cens.cause1"] = mean(train$D.1==0)
      values[cv, "cens.cause2"] = mean(train$D.2==0)
      
      ############################################################################################################
      #### A. rule estimation ####################################################################################
      ############################################################################################################
      
      ### A1. Cho et al (dtrSurv)
      args.CSK <- list(data = train, 
                       txName = "Trt",#Tx.nm.list,
                       models = form.CSK,
                       usePrevTime = FALSE, 
                       tau = tau, 
                       timePoints = timepoints,
                       criticalValue = value.criterion[1], 
                       evalTime = as.numeric(value.criterion[2]), 
                       splitRule = ifelse(value.criterion[1] == "mean", "mean", "logrank"),
                       ERT = ert, uniformSplit = ert, replace = !ert,
                       randomSplit = rs, nTree = Ntree, mTry = 6,#c(6, 6),
                       pooled = FALSE, # pooled = FALSE,
                       stratifiedSplit = FALSE)
      
      # actual fitting
      values[cv, "ns.CSK"] = nodesize
      set.seed(cv)
      
      CSK.aric.i <- 
        try(do.call(dtrSurv, c(args.CSK, list(nodeSize = nodesize, 
                                              minEvent = mindeath))))
      err.CSK = class(CSK.aric.i)[1] == "try-error"
      
      ### A5. zero-order model
      zom.aric.i <-
        try(do.call(dtrSurv, c(args.CSK, list(nodeSize = 1e+4, minEvent = 1e+4))))
      err.zom = class(zom.aric.i)[1] == "try-error"
      
      # if (!err.zom) {
      #   zom.pred =
      #     data.frame(
      #       zom.aric.i@stageResults[[1]]@optimal@optimalTx[1],
      #       zom.aric.i@stageResults[[2]]@optimal@optimalTx[1])
      #   tab[cv, 1] = zom.aric.i@stageResults[[1]]@optimal@optimalTx[1]
      #   tab[cv, 2] = zom.aric.i@stageResults[[2]]@optimal@optimalTx[1]
      #   tab[cv, 3] = paste0(tab[cv, 1], "-", tab[cv, 2])
      #   cat("ZOM freq table\n"); tab[, 3] %>% table %>% print
      # }
      
      ############################################################################################################
      #### B. Rules applied to a test set ########################################################################
      ############################################################################################################
      
      ### B0. skeletons for the predicted A's
      # opt.rule.pred = data.frame(A.1 = rep(NA, dim(cv.outsample)[1]), 
      #                            A.2 = NA)
      opt.rule.pred = data.frame(A.1 = rep(NA, dim(cv.outsample)[1]))
      for (q in 1) {
        opt.rule.pred[, q] = factor(NA, levels = lvls[[q]])
      }
      opt.rule.CSK <- opt.rule.zom <- opt.rule.pred

      for (q in seq_along(stg)) {
        elig = !test[, Tx.nm.list[q]] %>% is.na
        
        ### B1. the proposed method
        
        if (!err.CSK) {
          opt.CSK =
            predict(CSK.aric.i, 
                    newdata = test[elig, ], 
                    stage = q)
          opt.rule.CSK[elig, q] = factor(opt.CSK$optimal@optimalTx, levels = lvls[[q]]) # %>% as.numeric()
        }
       
        ### B5. zero-order model
        if (!err.zom) {
          opt.rule.zom[elig, q] = zom.pred[1, q]
          opt.rule.zom[elig, q] = factor(zom.pred[1, q], levels = lvls[[q]]) # %>% as.numeric()
        }
      }
      print(table(opt.rule.CSK, useNA = "always"))
      print(table(opt.rule.zom, useNA = "always"))
      
      
      ############################################################################################################
      #### C. Value estimation ###################################################################################
      ############################################################################################################
      ## weights
      # weight = weights_aric(data = test, 
      #                       weight.formula.list = form.weight, 
      #                        weight.formula.bin.list = form.weight.bin)

      test.tmp = test %>% transmute(A.1 = test[, Tx.nm.list[[1]]])
      # test.tmp.DW = test %>% transmute(A.1 = test[, paste0(Tx.bin.nm, ".1")])
      
      arg.val = list(test = test, 
                     actual = test.tmp, 
                     propensity = 1, # propensity = weight$propensity,
                     weight.censor = 1, # weight.censor = weight$weight.censor, 
                     criterion = value.criterion,
                     tau = tau)
      if (!err.CSK)  values[cv, "CSK"] = do.call(getValue, c(arg.val, list(estimated = opt.rule.CSK)))
      values[cv, "observed"] = do.call(getValue, c(arg.val[-which(names(arg.val) == "propensity")], list(estimated = test.tmp, propensity = 1)))
      if (!err.zom)  values[cv, "ZOM"] = do.call(getValue, c(arg.val, list(estimated = opt.rule.zom)))
      
      print(values[cv, ])
      print(apply(values, 2, mean, na.rm = TRUE))
      
      ############################################################################################################
      #### C. saving the results #################################################################################
      ############################################################################################################
      
      
      attr(values, "spec") = data.frame(criterion = value.criterion[1], 
                                        criterion.s = as.numeric(value.criterion[2]), 
                                        rule = rule, tau = tau,
                                        ert = ert, rs = rs, ntree = Ntree,
                                        nodesize = nodesize, mindeath = mindeath)
      # attributes(values)$spec
      if (cv==1) saveRDS(CSK.aric.i, paste0(fnm, "dtr.csk.aric_", nm, "_", cv, rds))
      saveRDS(values, paste0(fnm, "values_", nm, rds))
      
    }
    saveRDS(values, paste0(fnm, "values_", nm, rds))
    # saveRDS(tab, paste0(fnm, "tab_ZOM_", nm, rds))
  }
}

