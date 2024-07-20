# supplementary table for table 2
s2a_dataset <-
  lapply(1:dim(fn)[1], function(i) {
    fn.i = fn$fn[i]
    if (file.exists(fn.i)) {
      print(sprintf("file %s exists", fn.i))
      full <- readRDS(fn.i)
      a <- full$statistics
      n.sim = full$settings$n.sim
    } else if (file.exists(gsub("\\.rds", "_tmp.rds", fn.i))) {
      print(sprintf("_TMP.RDS exists: %s", fn.i))
      full <- NULL
      a <- readRDS(gsub("\\.rds", "_tmp.rds", fn.i))
    } else {
      full <- NULL
      a <- NULL
    }
    if (is.null(a)) {
      print("a is null")
      NULL
    } else {
      print("a is not null")
      if (priority_cause == 1){
        cause_string = "cause.1"
      } else{
        cause_string = paste0("cause.", priority_cause)
      }
      # print(all_methods)
      # print(priority_cause)
      # print(cause_string)
      # print(paste0(all_methods,"_",cause_string))
      propNames = paste0(all_methods,"_",cause_string)
      # print(propNames)

      as.data.frame(a) %>%
        dplyr::select(rep = sim,
                      contains(paste0("_",cause_string))) %>%
        mutate(ncauses = fn$ncauses[i],
               cause1prob = fn$cause1prob[i], # this only exists if endpoint = CR
               beta = fn$beta[i],
               prop = fn$prop[i],
               n    = fn$n[i],
               crit = fn$critS.no[i]) %>%
        pivot_longer(names_to = "method", values_to = "priority_cause_proportion",
                     c(propNames))
    }
  }) %>%
  do.call(rbind, .) %>%
  mutate(
    method = factor(method, levels = paste0(all_methods,"_",cause_string), labels = method.nm.abc),
    ncauses = factor(ncauses, levels = 1, labels = c("1" = "number_causes=2")),
    cause1prob = factor(cause1prob, levels = 1:2, labels = c("1" = "Prob(Cause) = 0.3", "2" = "Prob(Cause) = 0.6")),
    setting = factor(beta, levels = c(1,3,2,4), # reordering so that high censoring comes last.
                     labels = c("1" = "low ncov\n same beta",
                                "3" = "low ncov\n diff beta",
                                "2" = "high ncov\n same beta",
                                "4" = "high ncov\n diff beta")),
    n = factor (n, levels = 1:2, labels = c("1" = "n=400", "2" = "n = 1000")),
    design = factor(prop, levels = 1:2, labels = c("1" = "observational", "2" = "RCT")),
    crit.label = factor(crit, levels = 1:3, labels = c("1" = sprintf("Truncated %s mean, E[T]", Phase_lab),
                                                       "2" = sprintf("%s probability at t=3, Curve(3)", Phase_lab),
                                                       "3" = sprintf("%s probability at t=??, Curve(??)", Phase_lab))))
# Calculate the proportion of people dying from cause 1 for each simulation
# Then take mean and sd across all simulations
# ^s2a and s2b plots
plot_s2 = ggplot(s2a_dataset, aes(x = method, y = priority_cause_proportion, color = method)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.1,
    color = "purple"
  ) +
  facet_grid(setting ~ design + n + cause1prob) + 
  labs(title = sprintf("S2: Means (SDs) of Proportions of Cause %s Event Across %s simulations", priority_cause, n.sim), 
       x = "Method",
       y = "Proportion of Priority Cause") +
  scale_color_discrete(labels = paste0(method.nm.abc, ": ", method.nm.formal)) +
  theme_bw()
save_plot(sprintf("%s/figures/FigS2_endpoint_%s.png", 
                  endpoint, gsub("-", "", lab.date)), 
          plot_s2, base_height = 10, base_width = 20)
message("End of CR03.Simulation_Summary_Supplementary.R")