# person = which(czmk_training_recc != csk_training_recc & czmk_training_recc == 1)
which(rep_czmk$action != rep_zom$action & rep_czmk$CIF_eval > rep_zom$CIF_eval)
which(rep_czmk$action != rep_zom$action & rep_czmk$OS_eval < rep_zom$OS_eval & rep_czmk$CIF_eval > rep_zom$CIF_eval)
# which(P2_czmk$action != P2_zom$action & P2_czmk$CIF_eval > P2_zom$CIF_eval)
# P2_czmk[244,]
# person = person[3]
# person = which(rep_czmk$action!=rep_csk$action)
time_points_full = tp_surv
person = 123

trt0_timepoints = length(predd_surv_czmk_eval[["AUS_tp"]][[1]][[person]])
trt1_timepoints = length(predd_surv_czmk_eval[["AUS_tp"]][[2]][[person]])
tp_ind = max(trt0_timepoints, trt1_timepoints)
if (tp_ind == trt0_timepoints){
  time_points = predd_surv_czmk_eval[["AUS_tp"]][[1]][[person]]
} else if (tp_ind == trt1_timepoints){
  time_points = predd_surv_czmk_eval[["AUS_tp"]][[2]][[person]]
} else{
  stop("error")
}
id = plotting_indiv_predicted_curves(person, time_points, time_points_full,
                                     indiv_St_opt, indiv_St_trt0, indiv_St_trt1,
                                     indiv_CIF_opt, indiv_CIF_trt0, indiv_CIF_trt1)

# czmk_training_recc = policy_czmk@phaseResults$FinalOptimalTx_Recc
# csk_training_recc = policy_csk@stageResults[[1]]@optimal@optimalTx
# zom_training_recc = policy_zom@phaseResults$FinalOptimalTx_Recc
# obs_training = obs_1$action
# czmk_testing_recc = rep_czmk$action
# csk_testing_recc = rep_csk$action
# zom_testing_recc = rep_zom$action
# obs_testing = rep_obs$action
#
# table(test_truth$best.action, test_truth$best.status)
# table(rep_czmk$action, rep_czmk$status)
# table(rep_csk$action, rep_csk$status)
# table(rep_zom$action, rep_zom$status)
#
# table(train_truth$best.action, train_truth$best.status)
# table(czmk_training_recc, policy_czmk@call$data[,"status"])
# table(csk_training_recc, policy_csk@call$data[,"status"])
# # table(zom_training_recc, policy_zom@call$data[,"status"])
#
#
# #training dataset
# print(sprintf("trt-1 mean OS: %s", round(predd_surv_czmk[["mean"]][[1]][person],4)))
# print(sprintf("trt1 mean OS: %s", round(predd_surv_czmk[["mean"]][[2]][person],4)))
# print(sprintf("nonopt/opt ratio: %s", round(policy_czmk@phaseResults$SurvivalPhase1Results@optimal@NonOpt_Opt_Ratio[person],2)))
# print(sprintf("stop at P1? (0 = no; 1 = yes): %s", policy_czmk@phaseResults$SurvivalPhase1Results@optimal@Ratio_Stopping_Ind[person]))
# print(sprintf("trt-1 mean cause-1 CIF: %s", round(predd_ep_czmk[["mean"]][[1]][person],2)))
# print(sprintf("trt1 mean cause-1 CIF: %s", round(predd_ep_czmk[["mean"]][[2]][person],2)))
# print(czmk_training_recc[person])
# print(csk_training_recc[person])
# print(zom_training_recc[person])
# print(train_truth$best.action[person])
#
# mean(czmk_training_recc == csk_training_recc)
# mean(czmk_training_recc == train_truth$best.action)
# mean(czmk_training_recc == zom_training_recc)
# mean(csk_training_recc == train_truth$best.action)
#
# mean(czmk_testing_recc == csk_testing_recc)
# mean(czmk_testing_recc == test_truth$best.action)
# mean(csk_testing_recc == test_truth$best.action)
# mean(czmk_testing_recc == zom_testing_recc)

#testing dataset
id$St_trt_plot
if (criterion_phase1 == "mean.prob.combo"){
  print(predd_surv_czmk_eval[["Prob"]][[1]][person])
  print(predd_surv_czmk_eval[["Prob"]][[2]][person])
  print(predd_surv_czmk_eval[["Prob"]][[1]][person]/predd_surv[["mean"]][[2]][person])
  print(predd_surv_czmk_eval[["Prob"]][[2]][person]/predd_surv[["mean"]][[1]][person])
  id$CIF_trt_plot
  print(predd_ep_czmk_eval[["Prob"]][[1]][person])
  print(predd_ep_czmk_eval[["Prob"]][[2]][person])
} else{
  print(predd_surv_czmk_eval[["mean"]][[1]][person])
  print(predd_surv_czmk_eval[["mean"]][[2]][person])
  print(predd_surv_czmk_eval[["mean"]][[1]][person]/predd_surv_czmk_eval[["mean"]][[2]][person])
  print(predd_surv_czmk_eval[["mean"]][[2]][person]/predd_surv_czmk_eval[["mean"]][[1]][person])
  print(predd_surv_czmk_eval[["Stopping_Ind"]][person])
  id$CIF_trt_plot
  print(predd_ep_czmk_eval[["mean"]][[1]][person])
  print(predd_ep_czmk_eval[["mean"]][[2]][person])
}
rep_czmk$action[person]; rep_czmk$CIF_eval[person]#rep_czmk$event.time[person]; rep_czmk$status[person]
rep_csk$action[person]; rep_csk$CIF_eval[person]#rep_csk$event.time[person]; rep_csk$status[person]
rep_zom$action[person]; rep_zom$CIF_eval[person]#rep_zom$event.time[person]; rep_zom$status[person]
# test_truth$best.action[person]; test_truth$best.time[person]; test_truth$best.status[person]
#
# # indnottrue = which(policy_czmk@phaseResults[["FinalOptimalTx_Recc"]] != test_truth$best.action)
# # policy_czmk@phaseResults[["SurvivalPhase1Results"]]@optimal@NonOpt_Opt_Ratio[indnottrue]
#
# merge(test0.data.rep, test1.data.rep, by = "subj.id")[person,]
#
#
# ## training
# p2_indices = which(policy_czmk@phaseResults[["SurvivalPhase1Results"]]@optimal@Ratio_Stopping_Ind==0)
# for (p2_index in p2_indices){
#   message(p2_index)
#   # amongst those continuing to Phase2, which predicted survival AUC was greater for trt 1 over trt -1
#   ind1_train = ifelse(
#     policy_czmk@phaseResults[["SurvivalPhase1Results"]]@valueAllTx[["mean"]][[1]][p2_index] >
#       policy_czmk@phaseResults[["SurvivalPhase1Results"]]@valueAllTx[["mean"]][[2]][p2_index],
#     p2_index, NA)
#   # ind1_test = ifelse(predd_surv[["mean"]][[1]][p2_index]<predd_surv[["mean"]][[2]][p2_index], p2_index, NA)
#   if (!is.na(ind1_train)){
#     print(ind1_train)
#     # THUS, CSK WANTS TO PICK TRT 1
#     # We want to pick trt -1 so we want to see when this differs
#     # when cause1-CIF trt-1 < cause1-CIF trt1 (want minimum CIF)
#     ind2_train = ifelse(
#       policy_czmk@phaseResults[["EndPointPhase2Results"]]@valueAllTx[["mean"]][[1]][ind1_train] >
#         policy_czmk@phaseResults[["EndPointPhase2Results"]]@valueAllTx[["mean"]][[2]][ind1_train],
#       ind1_train, NA)
#     # print(ind2_train)
#     # message("train(policy):csk picks: ",policy_csk@stageResults[[1]]@optimal@optimalTx[p2_index]) # this should always give 1
#     # ind2_test = ifelse(predd_ep[["mean"]][[1]][ind1_test]<predd_ep[["mean"]][[2]][ind1_test], ind1_test, NA)
#     if (!is.na(ind2_train)){
#       print(ind2_train)
#       message("train(policy):csk picks: ",policy_csk@stageResults[[1]]@optimal@optimalTx[p2_index]) # this should always give 1
#       message("train(policy):CZMK (ME) picks: ",policy_czmk@phaseResults$FinalOptimalTx_Recc[p2_index]) # this should always give -1
#     }
#   }
# }
#
# test_truth$best.action[person];rep_czmk$action[person];rep_csk$action[person];rep_zom$action[person]
#
#
# # person = 1
# # time_points = 1:20
# # opt_data = t(indiv_St_opt)[time_points,person] %>% as.data.frame()
# # ggplot(opt_data, aes(x = time_points, y = .)) +
# #   geom_line(color = 'blue', size = 1, linetype = 'solid') +
# #   geom_point(color = 'blue', size = 3) +
# #   labs(x = 'Time Points',
# #        y = 'Values',
# #        title = 'Optimal Y for One Person Over Time') +
# #   theme_minimal()
# #
# #
# # trt0_data = as.data.frame(indiv_St_trt0)[time_points, person]
# # trt1_data = as.data.frame(indiv_St_trt1)[time_points, person]
# # trt_data_St = cbind(trt0 = trt0_data, trt1 = trt1_data) %>% as.data.frame()
# #
# # ggplot(trt_data_St, aes(x = 1:nrow(trt_data_St))) +
# #   geom_line(aes(y = trt0), color = 'red', size = 0.6, linetype = 'solid') +
# #   geom_line(aes(y = trt1), color = 'blue', size = 0.6, linetype = 'solid') +
# #   labs(x = 'Time',
# #        y = 'S(t)',
# #        title = sprintf('Survival Function for ID=%s comparing Trt0 vs Trt1', person))+
# #   scale_x_continuous(breaks = time_points) +
# #   theme_minimal()
#
