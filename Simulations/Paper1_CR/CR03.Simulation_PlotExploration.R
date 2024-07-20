
# EXPLORING TOP PERFORMANCE (crit = 1)
if (Phase.no == 1){
  fun_label = "overall survival"
  # fun = max()
} else{
  fun_label = "cause 1 survival"
  # fun = min()
}

supp_explore = function(interest){
  df = result.comb %>%
    # filter(rep < 3) %>%
    dplyr::filter(design %in% design.filter, crit == 1) %>%
    dplyr::select(cause1prob, design, n, setting, rep, method, interest) %>%
    arrange(cause1prob, design, n, setting, rep) %>%
    as.data.frame()
  return(df)
}

# when we aren’t the top performer in terms of survival,
# what is the error / how much worse is it than the top performer in survival
# (shouldn’t be that much worse).
# ^ for each simulation, pick the largest summary statistic and then for each of the methods, take the largest minus each method.
# see if its within alpha.;
df = supp_explore("value")
df1 = df %>%
  group_by(cause1prob, design, n, setting, rep) %>%
  mutate(ind = ifelse(method == 'A' & value == max(value, Phase.no), "A",
                      ifelse(method == 'B' & value == max(value, Phase.no), "B",
                             ifelse(method == 'C' & value == max(value, Phase.no), "C",
                                    ifelse(method == 'D' & value == max(value, Phase.no), "D",
                                           NA))))) %>%
  mutate(ind_val = ifelse(!is.na(ind), value, NA)) %>%
  mutate(across(everything(), ~ifelse(is.na(.), na.omit(.), .))) %>%
  ungroup() %>%
  mutate(method = factor(method, levels = method.levels, labels = method.nm.abc)) %>%
  mutate(ind_diff = ind_val - value)

# Then do mean and variance across simulations.
# mean is the regret for survival and regret for primary cause.
# ^ each method would want a smaller mean b/c closer to best.
# Do the same for primary CR (primary cause).
df2 = df1 %>%
  # filter(!is.na(ind)) %>%
  dplyr::select(-c(ind_val,ind, value)) %>%
  arrange(cause1prob, design, n, setting, rep, method) %>%
  group_by(cause1prob, design, n, setting, method) %>%
  summarize(Mean = mean(ind_diff),
            SD = sd(ind_diff)) %>%
  ungroup()

plot_s1a = ggplot(df1, aes(x = method, y = ind_diff, color = method)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    width = 0.1,
    color = "purple"
  ) +
  # facet_grid(setting ~ design + n + cause1prob) +
  facet_grid(n) +
  labs(title = sprintf("S1: Means (SDs) of %s Differences Across %s simulations", fun_label, n.sim),
       x = "Method",
       y = "Difference in Best Method Mean and Currnet Method Mean") +
  scale_color_discrete(labels = paste0(method.nm.abc, ": ", method.nm.formal)) +
  theme_bw()

plot_s1b = ggplot(df2, aes(x = method, y = Mean, color = method)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, color = "purple") +
  # facet_grid(setting ~ design + n + cause1prob) +
  facet_grid(n) +
  labs(title = sprintf("S1b: Simplied Mean (SD) Plot of %s Differences", fun_label),
       x = "Method",
       y = "Mean") +
  scale_color_discrete(labels = paste0(method.nm.abc, ": ", method.nm.formal)) +
  theme_bw()

save_plot(sprintf("%s/figures/FigS1_%s_%s.png",
                  endpoint, crit_lab, gsub("-", "", lab.date)),
          plot_s1a, base_height = 10, base_width = 20)
# save_plot(sprintf("%s/figures/Fig1b_%s_%s.png",
#                   endpoint, crit_lab, gsub("-", "", lab.date)),
#           plot_s1b, base_height = 10, base_width = 20)