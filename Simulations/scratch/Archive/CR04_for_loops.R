####### OS #######
# weights
weight.OS = weights_rda(data = test,
                        weight.formula.list = form.weight,
                        covariates = covariates)
arg.val.OS = list(test = test, 
                  actual = test.tmp, 
                  propensity = weight.OS$propensity,
                  weight.censor = weight.OS$weight.censor, 
                  criterion = criterion_phase1,
                  tau = tau)
for (method in loop_methods) {
  err_method <- paste0("err.", method)
  opt_method <- paste0("opt.rule.", method)
  if (!get(err_method)) {
    values[cv, paste0(method, ".OS")] <- 
      do.call(getValue, 
              c(arg.val.OS, 
                list(estimated = 
                       get(opt_method))))
  }
}
values[cv, "observed.OS"] = do.call(getValue, 
                                    c(arg.val.OS[-which(names(arg.val.OS) == "propensity")],
                                      list(estimated = test.tmp, 
                                           propensity = 1)))



### B4. Goldberg-Kosorok RF
# if (!err.GKRF) {
#   opt.GKRF = predict.opt.sep(GKRF.i$survRF[[1]],
#                              newdata = test[elig, ] %>% dplyr::select(-V.2, -d.2) %>%
#                                mutate(prevTime.2 = V.1),
#                              Tx.label = paste0(Tx.nm, ".", q))
#   opt.rule.GKRF[elig, q] = lvls[[1]][opt.GKRF$optimal.Tx]
# }
# ### B5. Goldberg-Kosorok linear
# if (!err.GKLM) {
#   opt.GKLM = predict.opt.sep(GKLM.i$survRF[[1]],
#                              newdata = test[elig, ] %>% dplyr::select(-V.2, -d.2) %>%
#                                mutate(prevTime.2 = V.1),
#                              Tx.label = paste0(Tx.nm, ".", q))
#   opt.rule.GKLM[elig, q] = lvls[[1]][opt.GKLM$optimal.Tx]
# }

# ####### OS #######
# weight.OS = weights_rda(data = testing_dataset,
#                         weight.formula.list = form.weight,
#                         covariates = covariates,
#                         event_indicator_string = "D.0")
# ####### Priority Cause #######
# weight.PC = weights_rda(data = testing_dataset,
#                         weight.formula.list = form.weight,
#                         covariates = covariates,
#                         event_indicator_string = paste0("D.", priority_cause))
# arg.val.OS = list(test = test,
#                   actual = test.tmp,
#                   propensity = weight.OS$propensity,
#                   weight.censor = weight.OS$weight.censor,
#                   criterion = criterion_phase1,
#                   tau = tau)
# arg.val.PC = list(test = test,
#                   actual = test.tmp,
#                   propensity = weight.PC$propensity,
#                   weight.censor = weight.PC$weight.censor,
#                   criterion = criterion_phase2,
#                   tau = tau)

# for (method in loop_methods) {
#   err_method <- paste0("err.", method)
#   opt_method <- paste0("opt.rule.", method)
#   if (!get(err_method)) {
#     # Truncated OS Mean/Prob
#     values[cv, paste0(method, ".OS")] <-
#       do.call(getValue,
#               c(arg.val.OS,
#                 list(estimated =
#                        get(opt_method))))
#     # Truncated PC Mean/Prob (anything else is counted as censored)
#     values[cv, paste0(method, ".PC")] <-
#       do.call(getValue,
#               c(arg.val.PC,
#                 list(estimated =
#                        get(opt_method))))
#   }
# }

# if (!err.CZMK)  values[cv, "CZMK.OS"] = do.call(getValue,
#                                              c(arg.val.OS,
#                                                list(estimated = opt.rule.CZMK)))
# if (!err.CSK)  values[cv, "CSK.OS"] = do.call(getValue,
#                                            c(arg.val.OS,
#                                              list(estimated = opt.rule.CSK)))
# if (!err.PMCR)  values[cv, "PMCR.OS"] = do.call(getValue,
#                                              c(arg.val.OS,
#                                                list(estimated = opt.rule.PMCR)))
# if (!err.ZOM)  values[cv, "ZOM.OS"] = do.call(getValue,
#                                            c(arg.val.OS,
#                                              list(estimated = opt.rule.ZOM)))
# if (!err.CSKzom)  values[cv, "CSKzom.OS"] = do.call(getValue,
#                                                  c(arg.val.OS,
#                                                    list(estimated = opt.rule.CSKzom)))
# values[cv, "observed.OS"] = do.call(getValue,
#                                     c(arg.val.OS[-which(names(arg.val.OS) == "propensity")],
#                                       list(estimated = test.tmp,
#                                            propensity = 1)))