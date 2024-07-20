#' @param p number of covariates
#' @param Sig covariate error covariance
#' @param corr cor of two error processes (rho and lambda)
#' 
#' @return event.time observed time
#' @return delta 1(min(death1, death2) <= censor)
#' @covariates(Z1-Zp) covariate values at the beginning of the stage
#' 
policy_gendata <- 
  function(n.sample = 100, tau = 10,        # structural parameters
           covariate = NULL,            
           p = 5, Sig = diag(p) + 0.2 - diag(p) * 0.2,     # covariate structure
           predHazardFn, predPropensityFn,                 # list of predictor functions
           policy = NULL,                                  # optimal rule (if !is.null, propensity scores are ignored.) for value calculation
           printFlag = TRUE                        # policy is an object of rsf.obj list.
  ) {
    
    # at.risk = arg.obs$at.risk; p = arg.obs$p; corr = arg.obs$corr;
    # n.sample = 3; n.stages = arg.obs$n.stages; tau = arg.obs$tau;
    # tick = arg.obs$tick; hidden_data = arg.obs$hidden_data; printFlag = arg.obs$printFlag;
    # predHazardFn = arg.obs$predHazardFn; predPropensityFn = arg.obs$predPropensityFn;
    # predCensorFn = arg.obs$predCensorFn; surv.previous = rep(0, n.sample)
    # rho = NULL; omega = NULL; covariate = NULL;
    # policy = NULL
    # Sig = diag(p) + 0.2 - diag(p) * 0.2
    
    require(dplyr)
    at.risk = rep(1, n.sample)
    if (length(at.risk) != n.sample)  stop ("length of at.risk should match.")
    if (!is.null(covariate)) if (dim(covariate)[1] != n.sample) stop ("length of at.risk should match.")
    
    # initial state
    if (is.null(covariate))    stop("no covariates?")
    if (is.null(colnames(covariate))) colnames(covariate) = paste0("Z", 1:p)
    # # skeleton
    output <-
      array(NA, dim = c(n.sample, 2 + 2 + 2 + p),
            dimnames = list(1:n.sample, c("subj.id", "rep.id", "event.time", "delta", "action", "at.risk",
                                         colnames(covariate))))
    output[, "subj.id"] <- 1:n.sample
    output[, "rep.id"] <- 0           # rep.id is reserved for later use (repeated random draws)
    action = rep(NA, n.sample)           # initialize action vector

    output[, "at.risk"]           <- at.risk
    output[, colnames(covariate)] <- covariate
    
    ## action                (t=0)
    if (!is.null(policy)) {
      if (attr(policy, "class") %in% c("dwSurv")) {
        # tmp.out <<- output      
        dat.wide <- dwTrans(output, n.stages = n.stages, p = p)
        # vars.tf   <- names(get_all_vars(policy$tf.mod[[stage]], dat.wide))
        # vars.blip <- names(get_all_vars(policy$blip.mod[[stage]], dat.wide))
        # vars <- unique(c(vars.tf, vars.blip))
        args <- c(list(policy,
                       newdata = dat.wide[as.logical(at.risk), ],
                       # output[as.logical(at.risk), c("subj.id", vars), stage] %>% as.data.frame,
                       stage = stage,
                       return.optimal.Q = FALSE))
        action[at.risk != 0] <- do.call(decision.rule.dw, args)$rule.fix
      } else if (attr(policy, "class") %in% c("ITRSurv")) { #### my method
        
        print("POLICY GIVEN FOR ITRSURV")
        x = output[as.logical(at.risk),] %>% output2observable()
        x[, paste0(c("T_", "delta_"), stage)] = 0  # dummy values in order for the package not to drop the NA rows.
        x = get_all_vars(policy@stageResults[[stage]]@model, x)
        args <- list(policy, newdata = x, stage = stage)
        action[at.risk != 0] <- do.call(predict, args)$optimal@optimalTx
        
        
        
      } else if (attr(policy, "class") %in% c("DTRSurv")) { #### cho's package
        x = output[as.logical(at.risk),,] %>% output2observable()
        x[, paste0(c("T_", "delta_"), stage)] = 0  # dummy values in order for the package not to drop the NA rows.
        x = get_all_vars(policy@stageResults[[stage]]@model, x)
        args <- list(policy, newdata = x, stage = stage)
        action[at.risk != 0] <- do.call(predict, args)$optimal@optimalTx
      } else if (attr(policy, "class") %in% c("GKLM")) { #### new function!
        x = output[as.logical(at.risk),,] %>% output2observable()
        action[at.risk != 0] <- 
          predict.opt.sep(policy[[stage]],  newdata = x, Tx.label = paste0("A_", stage))$optimal.Tx %>% 
          {ifelse(.==1, -1, 1)}
      } else if (attr(policy, "class") %in% c("GKRF")) { #### new function!
        x = output[as.logical(at.risk),,] %>% output2observable()
        x = x[, c(paste0("A_", stage), policy[[stage]][[1]]$xvar.names)]
        x[, paste0("A_", stage)] = 0 # dummy values
        
        action[at.risk != 0] <- 
          predict.opt.sep(policy[[stage]],  newdata = x, Tx.label = paste0("A_", stage))$optimal.Tx %>% 
          {ifelse(.==1, -1, 1)}
      }
      if (length(action[at.risk == 0]) > 0){
        action[at.risk == 0] <- NA
        stop("Why are there people not at risk at data generation??")
      }
    } else {
      propensity = predPropensityFn(covariate = covariate)
      action = suppressWarnings(rbinom(n.sample, 1, propensity) * 2 - 1) # 1 for aggressive and 0 for gentle
      # for NAs in propensity, the action is NA. Thus, the warnings are suppressed.
      if (length(action[at.risk == 0]) > 0){
        action[at.risk == 0] <- NA
        stop("Why are there people not at risk at data generation??")
      }
    }
    message("End of policy_gendata in F02.singleStage.R")
  }

