solo.plot = 1
local = 1
saving_eps = TRUE#TRUE
crit.tot = 1 # total number of critical values (for now - just mean!!)
testing_out = 1

if (local == 1){
  setwd("~/Desktop/UNC_BIOS_PhD/DissertationPhD/Thesis/Code/Analyses/Simulations/Paper1_CR")
} else{
  setwd("/nas/longleaf/home/cwzhou/Dissertation/Analyses/Simulations/Paper1_CR")
}
source("CR00.Simulation_Parameters.R") # change local in this script to 0 for cluster

library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(purrr)

lab.date = "2025-09-29" #"2025-01-02" #"2024-09-13" #"2024-09-07" #"2024-08-31" #Sys.Date()#"2024-08-20"#"2024-02-27" #"2024-02-18" #Sys.Date()  # change this for the date of RDS data you want
dir_fig = sprintf("./figure/%s/%s", generate_failure_method, lab.date)

# Sequence of dates to include
date_seq <- seq(as.Date("2025-09-01"), as.Date("2025-09-22"), by = "day")
date_seq <- date_seq[!date_seq %in% as.Date(c("2025-09-07", "2025-09-12"))]
date_seq

# date_seq = date_folder
file_name <- "simResult_fine_gray_censor2_nCauses1_cause1prob1_beta1_prop1_n1_critS1_critE1.rds"
# Read and combine all files
final_tbl <- map_dfr(date_seq, function(d) {
  if (length(date_seq) == 1){
    date_str = date_folder
  } else{
    date_str <- format(d, "%Y-%m-%d")
  }
  file_path <- file.path("./output", "fine_gray/revision", date_str, file_name)
  
  dat_list <- readRDS(file_path)
  dat = dat_list$statistics
  
  # If it's a list, combine into one data frame
  dat1 <- bind_rows(dat)
  
  dat1 %>% mutate(date = date_str)
  final_tbl = dat1
}); dim(final_tbl)
View(final_tbl)

(n.sim = nrow(final_tbl))

# # crit.no = 1
# # critS.no = 1
# # critE.no = 1
# # # below needs to edited based on what the values are - or needs to be automated
# # y_limits_prob <- c(min(0.12, 0.13), max(0.17, 0.20))
# # y_limits_mean <- c(0,1.1)#c(min(0.45, 0.50), max(0.70, 0.80))

method.levels = c(1,2,#3,4,
                  5,6)
method.nm.abc =
  c("A", "B", #"C", "D", 
    "C", "D")
method.nm.simple =
  c("czmk", "csk", #"pmcr", "aipwe", 
    "zom", "obs") # this is important for selecting variables below.
method.nm.formal =
  c("itrSurv", "dtrSurv (2023)",
    # "PMCR (2021)", "AIPWE (2021)",
    "zero-order model", "observed policy")
# "Goldberg & Kosorok (2012), RF", "Goldberg & Kosorok (2012), linear", "Simoneau et al. (2019)",

if (generate_failure_method == "fine_gray"){
  cause1prob.levels = c(1,2)
  cause1prob.labels = c(sprintf("Fine-Gray probability mass: %s",cause1_prob$small.cause1prob$cause1prob),
                        sprintf("Fine-Gray probability mass: %s", cause1_prob$large.cause1prob$cause1prob))
} else{
  cause1prob.levels = cause1_prob$cause1.prob$cause1prob
  cause1prob.labels = " "
}

censor.levels = c(1,2)
censor.labels = c("Average 56% Censoring","High Censoring (50%)")
n.levels = c(1,2)
n.labels = c(sprintf("N=%s",size$small.sample.size$n),
             sprintf("N=%s",size$large.sample.size$n))
design.levels = c(1,2)
design.labels = c("Trt: Covariate Dependent","Trt: Covariate Independent")
beta.levels = c(1)#,2)
beta.labels = c(sprintf("%s Covariates",ncov.list$beta1))#,
                # sprintf("%s Covariates",ncov.list$beta2))
# beta.levels = c(1,2)
# beta.labels = c("setting1",
#                 "setting2")
# beta.levels = c(1)
# beta.labels = c("setting1")

file_naming = function(lab.date, file_lab, crit.no){
  paste0(dir_fig,"/CR02.",file_lab,"_", gsub("-", "", lab.date), "_crit", crit.no, ".eps")
}

if (endpoint == "CR"){
  if (generate_failure_method == "fine_gray"){
    fn <-
      expand.grid(ncauses = 1,
                  censor = 1:2,
                  cause1prob = 1:length(cause1prob.levels),
                  beta = 1:length(beta.levels),
                  prop = 1:2,
                  n = 1:length(n.levels),
                  critS.no = 1:2,
                  critE.no = 1:2) %>%
      mutate(nm = paste0(beta, "-", prop, "-", n),
             fn = paste0(dir_rds,"/simResult_", generate_failure_method,
                         "_censor", censor,
                         "_nCauses", ncauses,
                         "_cause1prob", cause1prob,
                         "_beta", beta, "_prop", prop,
                         "_n", n, "_critS", critS.no,
                         "_critE", critE.no, ".rds"))
  } else if (generate_failure_method == "simple_exp"){
    fn <-
      expand.grid(ncauses = 1,
                  censor = 1:2,
                  # cause1prob = 1:length(cause1prob.levels),
                  beta = 1:length(beta.levels),
                  prop = 1:2,
                  n = 1:length(n.levels),
                  critS.no = 1:2,
                  critE.no = 1:2) %>%
      mutate(nm = paste0(beta, "-", prop, "-", n),
             fn = paste0(dir_rds,"/simResult_", generate_failure_method,
                         "_censor", censor,
                         "_nCauses", ncauses,
                         # "_cause1prob", cause1prob,
                         "_beta", beta, "_prop", prop,
                         "_n", n, "_critS", critS.no,
                         "_critE", critE.no, ".rds"))
  }
  
} else{
  fn <-
    expand.grid(ncauses = 1, censor = 1:2, beta = 1:6, prop = 1:2,
                n = 1:2, critS.no = 1:2, critE.no = 1:2) %>%
    mutate(nm = paste0(beta, "-", prop, "-", n),
           fn = paste0(endpoint,"/output/simResult_", lab.date,
                       "_censor", censor,
                       "_nCauses", ncauses,
                       "_beta", beta, "_prop", prop,
                       "_n", n, "_critS", critS.no,
                       "_critE", critE.no, ".rds"))
}
print(fn)

p.list <- p.list.solo <- list()
for (crit.no in 1:crit.tot){
  message("%%%%%%%%%%%% crit.no", crit.no, " %%%%%%%%%%%%%%%%%")
  if (crit.no == 1){
    crit_lab = "Truncated mean"
    crit_lab_1 = "a"
    # y_limits = y_limits_mean
  } else{
    crit_lab = "Probability"
    crit_lab_1 = "b"
    # y_limits = y_limits_prob
  }
  file.name.saved = file_naming(lab.date, "SimulationResults", crit.no)
  file.name.saved.solo = file_naming(lab.date, "solo.SimulationResults", crit.no)
  message(file.name.saved)
  for (Phase.no in 1:2){
    if (Phase.no == 1){
      Phase_lab = "survival"
      Phase_lab_1 = "Overall Event-Free S(t)"
    } else{
      Phase_lab = "endpoint"
      Phase_lab_1 = "Cause 1 CIF"
    }
    if (Phase.no == 1){
      crit = crit_surv
      critS.no = crit.no
      t0 = crit[[2]]$crit.value_phase1 # right now, we always want it to be set to the t0 time for when it exists aka mean.prob.combo as of now
    } else{
      crit = crit_endpoint
      critE.no = crit.no
      t0 = crit[[2]]$crit.value_phase2
    }
    print(crit_lab); print(Phase_lab)
    method.nm.simple1 = paste0(method.nm.simple, "_", Phase_lab)
    var_method = select_method_endpoints(method.nm.simple, Phase_lab)
    sim_name = "sim"
    beginning <- final_tbl %>%
      as.data.frame() %>%
      dplyr::select(rep = !!sym(sim_name),
                    all_of(var_method),
                    czmk_n_phase2, zom_n_phase2,
                    training_percent.censor, training_cause.1, training_cause.2)
    i = 1
    if (generate_failure_method == "fine_gray"){
      middle = beginning %>%
              mutate(ncauses = fn$ncauses[i],
                     censor = fn$censor[i],
                     cause1prob = fn$cause1prob[i], # this only exists if endpoint = CR
                     beta = fn$beta[i],
                     prop = fn$prop[i],
                     n    = fn$n[i],
                     crit = fn$critS.no[i])
          } else{ # simple_exp has no cause1prob
            middle = beginning %>%
              mutate(ncauses = fn$ncauses[i],
                     censor = fn$censor[i],
                     # cause1prob = fn$cause1prob[i], # this only exists if endpoint = CR
                     beta = fn$beta[i],
                     prop = fn$prop[i],
                     n    = fn$n[i],
                     crit = fn$critS.no[i])
          }
    result.comb = middle %>% pivot_longer(cols = all_of(var_method),
                                  names_to = "method",
                                  values_to = "value") %>%
      # below is ONLY for days,not years
      mutate(value = ifelse(value<5, value*365.25,value)) %>%
      mutate(
        method = factor(method,
                        levels = method.nm.simple1,
                        labels = method.nm.abc),
        censor = factor(prop,
                        levels = censor.levels,
                        labels = censor.labels),
        ncauses = factor(ncauses,
                         levels = 1:max(ncauses),
                         labels = c("1" = "2 causes")),
        setting = factor(beta,
                         levels = beta.levels,
                         labels = beta.labels),
        n = factor(n,
                   levels = n.levels,
                   labels = n.labels),
        design = factor(prop,
                        levels = design.levels,
                        labels = design.labels),
        crit.label = factor(crit,
                            levels = 1:2,
                            labels = c("1" = sprintf("Mean Truncated %s", crit_lab),
                                       "2" = sprintf("%s probability at t=%s, Curve(%s)", crit_lab, t0, t0)
                            )))
    #end of result.comb
    print("end of result.comb")
    
    if (generate_failure_method == "fine_gray"){
      result.comb1 = result.comb %>%
        mutate(cause1prob = factor(cause1prob,
                                   levels = cause1prob.levels,
                                   labels = cause1prob.labels))
      result.stat <-
        result.comb1 %>%
        aggregate(cbind(value, czmk_n_phase2, zom_n_phase2, training_percent.censor, training_cause.1, training_cause.2) ~ method +
                    ncauses + censor + cause1prob +
                    setting + n + design + crit + crit.label , data = .,
                  FUN = function(x) round(mean(x, na.rm = TRUE), 2)) %>%
        mutate(progress = paste(training_cause.1, training_cause.2, sep = " / "))
      result.stat.sd <-
        result.comb1 %>%
        aggregate(value ~ method + ncauses + censor + cause1prob + setting + n + design +
                    crit + crit.label , data = .,
                  FUN = function(x) round(sd(x, na.rm = TRUE), 3)) %>%
        rename(sd = value)
    } else{
      result.comb1 <- result.comb
      result.stat <-
        result.comb1 %>%
        aggregate(cbind(value, czmk_n_phase2, zom_n_phase2, training_percent.censor, training_cause.1, training_cause.2) ~ method +
                    ncauses + censor +
                    setting + n + design + crit + crit.label , data = .,
                  FUN = function(x) round(mean(x, na.rm = TRUE), 2)) %>%
        mutate(progress = paste(training_cause.1, training_cause.2, sep = " / "))
      result.stat.sd <-
        result.comb1 %>%
        aggregate(value ~ method + ncauses + censor + setting + n + design +
                    crit + crit.label , data = .,
                  FUN = function(x) round(sd(x, na.rm = TRUE), 3)) %>%
        rename(sd = value)
    }
    result.stat <- left_join(result.stat, result.stat.sd)
    rm(result.stat.sd)
    file.name.phase = file_naming(lab.date, Phase_lab, crit.no)
    file.name.phase.solo = file_naming(lab.date, sprintf("solo.%s",Phase_lab), crit.no)
    design.filter = c("Trt: Covariate Dependent","Trt: Covariate Independent")
    if (crit.no == 1){
      ylabs = sprintf("Mean Truncated %s", Phase_lab_1)
    } else {
      ylabs = paste0(Phase_lab_1, " ",crit_lab, " at t = ", t0)
    }
    print(sprintf("ylabs: %s", ylabs))
    result.stat.i =
      result.stat %>% filter(design %in% design.filter, crit %in% crit.no, crit.label %in% crit.label[crit.no])
    # message('LINE 213: facet_grid is most common cause for: Error in `combine_vars()`')
    if (generate_failure_method == "fine_gray"){
      p0 <-
        result.comb1 %>%
        dplyr::filter(design %in% design.filter, crit == crit.no) %>%
        ggplot(aes(x = method, y = value, group = method, color = method)) +
        facet_grid(cause1prob + setting ~ censor + n + design, scales = "free_y")
    } else{
      p0 <-
        result.comb1 %>%
        dplyr::filter(design %in% design.filter, crit == crit.no) %>%
        ggplot(aes(x = method, y = value, group = method, color = method)) +
        facet_grid(setting ~ censor + n + design, scales = "free_y")
    }
    p = p0 +
      # geom_point() +
      geom_boxplot() +
      geom_jitter(width = 0.1, height = 0) + # EPS does not support alpha.
      scale_color_discrete(labels = paste0(method.nm.abc, ": ", method.nm.formal)) +
      ylab(ylabs) +
      theme_bw() #+
    # theme(legend.position = "bottom")
    rng = suppressWarnings(layer_scales(p)$y$range$range)
    rng[3] = rng[2] - rng[1]
    rng[4] = rng[3] * 0.4 + rng[2] # y coordinate for censoring %
    rng[5] = rng[3] * 0.2 + rng[2] # y coordinate for flowchart
    
    p.list[[Phase.no]] <-
      p +
      stat_summary(aes(x = as.numeric(method),
                       y = value),
                   fun = mean,
                   geom = 'point',
                   col = "black",
                   shape = "square",
                   size = 3) +
      stat_summary(aes(x = as.numeric(method),
                       y = value,
                       label = round(..y.., 2),
                       #below is code for alternating. if you use this, comment out vjust=-1.1 below
                       # vjust = ifelse(as.numeric(method) %% 2 == 0,
                       #                2,
                       #                -1)
      ),
      fun = mean,
      geom = 'text',
      col = "black",
      vjust = -1.1,
      size = 5) +
      # Expand y-axis limits to prevent labels from being cut off
      scale_y_continuous(expand = expansion(mult = c(0, 0.1)))  # Add 10% padding above
    
    if (saving_eps == TRUE){
      ggsave(file.name.phase, p.list[[Phase.no]], device="eps", width = 12, height = 10)
      ggsave(file.name.saved %>% gsub(".eps", sprintf("_Phase%s.png", Phase.no), .) , #save as png too
             p.list[[Phase.no]],
             width = 12, height = 10)
    }
    
    
    if (solo.plot == 1){
      solo.result.comb1 = result.comb1 %>%
        filter(design %in% design.filter[1]) %>%
        filter(setting %in% "30 Covariates",
               n %in% "N=300",
               censor %in% "Average 56% Censoring")
      
      if (generate_failure_method == "fine_gray"){
        solo.result.comb1 = solo.result.comb1 %>%
          filter(cause1prob %in% "Fine-Gray probability mass: 0.2")
        if (revision == 1){
          p0.solo <-
            solo.result.comb1 %>%
            ggplot(aes(x = method, y = value, group = method, color = method))
        } else{
          p0.solo <-
            solo.result.comb1 %>%
            ggplot(aes(x = method, y = value, group = method, color = method)) +
            facet_grid(cause1prob + setting ~ censor + n + design)# , scales = "free_y")
        }
      } else{
        p0.solo <-
          solo.result.comb1 %>%
          ggplot(aes(x = method, y = value, group = method, color = method)) +
          facet_grid(setting ~ censor + n + design, scales = "free_y")
      }
      p.solo = p0.solo +
        geom_boxplot() +
        geom_jitter(width = 0.1, height = 0) + # EPS does not support alpha.
        scale_color_discrete(labels = paste0(method.nm.abc, ": ", method.nm.formal)) +
        ylab(ylabs) +
        theme_bw() #+
      # theme(legend.position = "bottom")
      rng = suppressWarnings(layer_scales(p.solo)$y$range$range)
      rng[3] = rng[2] - rng[1]
      rng[4] = rng[3] * 0.4 + rng[2] # y coordinate for censoring %
      rng[5] = rng[3] * 0.2 + rng[2] # y coordinate for flowchart
      p.list.solo[[Phase.no]] <-
        p.solo +
        stat_summary(aes(x = as.numeric(method),
                         y = value),
                     fun = mean,
                     geom = 'point',
                     col = "black",
                     shape = "square",
                     size = 2) +
        stat_summary(aes(x = as.numeric(method),
                         y = value,
                         label = round(..y.., 3),
        ),
        fun = mean,
        geom = 'text',
        col = "black",
        vjust = -1.1,
        size = 4) +
        scale_y_continuous(breaks = function(x) seq(floor(min(x)), ceiling(max(x)), by = 100),
                           # Expand y-axis limits to prevent labels from being cut off
                           expand = expansion(mult = c(0, 0.1)))  # Add 10% padding above
      
      if (saving_eps == TRUE){
        ggsave(file.name.phase.solo, p.list.solo[[Phase.no]], device="eps", width = 12, height = 10)
        ggsave(file.name.saved.solo %>% gsub(".eps", sprintf("_Phase%s.png", Phase.no), .) , #save as png too
               p.list.solo[[Phase.no]],
               width = 12, height = 10)
      }
    }
    
  } # end of Phase.no for-loop
  
  
  p.1 = p.list[[1]] +
    theme(legend.position = "none",
          axis.title.x=element_blank()) #+
  # ylim(y_limits)
  p.2 = p.list[[2]] +
    theme(legend.position = "none",
          axis.title.x=element_blank()) #+
  # ylim(y_limits)
  
  p.grid <- plot_grid(p.1,
                      p.2,
                      align = "vh",
                      axis = "tblr",
                      nrow = 1, ncol = 2)#,
  # common.legend = FALSE)
  p.grid1 <- plot_grid(p.grid,
                       get_legend(p.list[[2]] +
                                    theme(legend.direction = "horizontal",
                                          legend.key.size = unit(2, "cm"),    # Adjust the size of the legend keys
                                          legend.text = element_text(size = 12), # Adjust the size of the legend text
                                          legend.spacing.x = unit(0.1, "cm")) +
                                    guides(color = guide_legend(nrow = 1, title = "Methods"))),
                       align = "vh",
                       axis = "tblr",
                       ncol = 1,
                       nrow = 2,
                       rel_heights = c(1, 0.1))
  #we can ignore the warning message about return_all
  
  # # Extract the legend from one plot
  # legend <- get_legend(p.list[[1]] + theme(legend.position = "bottom"))
  #
  # # Combine plots side by side
  # combined <- plot_grid(p.1, p.2, align = "v", ncol = 2)
  #
  # # Add the legend underneath
  # final_plot <- plot_grid(combined,
  #                         legend,
  #                         ncol = 1,
  #                         rel_heights = c(1, 0.1))
  #
  # # Display the plot
  # print(final_plot)
  
  if (saving_eps == TRUE){
    save_plot(file.name.saved, p.grid1, base_height = 10, base_width = 20)
  }
  ggsave(file.name.saved %>% gsub(".eps", ".png", .), #save as png too
         p.grid1,
         width = 20, height = 10)
  
  if (solo.plot == 1){
    # y_limits = c(0.3,2.5)
    y_limits = c(0,1000)
    p.1.solo = p.list.solo[[1]] +
      theme(legend.position = "none",
            axis.title.x=element_blank()) +
      ylim(y_limits)
    p.2.solo = p.list.solo[[2]] +
      theme(legend.position = "none",
            axis.title.x=element_blank()) +
      ylim(y_limits)
    p.grid.solo <- plot_grid(p.1.solo,
                             p.2.solo,
                             align = "vh",
                             axis = "tblr",
                             nrow = 1, ncol = 2)#,
    # common.legend = FALSE)
    p.grid1.solo <- plot_grid(p.grid.solo,
                              get_legend(p.list.solo[[2]] +
                                           theme(legend.direction = "horizontal",
                                                 legend.key.size = unit(1, "cm"),    # Adjust the size of the legend keys
                                                 legend.text = element_text(size = 20), # Adjust the size of the legend text
                                                 legend.spacing.x = unit(0.1, "cm")) +
                                           guides(color = guide_legend(nrow = 1, title = "Methods", size=1))),
                              align = "vh",
                              axis = "tblr",
                              ncol = 1,
                              nrow = 2,
                              rel_heights = c(1, 0.1))
    save_plot(file.name.saved.solo, p.grid1.solo, base_height = 10, base_width = 20)
    ggsave(file.name.saved.solo %>% gsub(".eps", ".png", .), #save as png too
           p.grid1.solo,
           width = 20, height = 10)
  }
  
} # end of crit.tot for-loop

# EXPLORING PROPOPRTION OF PEOPLE DYING FROM CAUSE 1 (crit = 1: truncated mean)
# source("CR03.Simulation_Summary_Supplementary.R")


if (local == 1){
  # beep()
}

message("End of CR02.Simulation_Summary.R")


#############
# Columns that need conversion from days -> years
cols_years <- c(
  "czmk_survival", "csk_survival", "zom_survival", "obs_survival",
  "czmk_endpoint", "csk_endpoint", "zom_endpoint", "obs_endpoint"
)

final_tbl_means <- final_tbl %>%
  dplyr::select(-sim) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%   # mean across sims
  mutate(across(all_of(cols_years), ~ .x * 365.25)) %>%          # convert years -> days
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%          # round all numeric columns
  as_tibble()                                                     # make a nice tibble table

# View the table
View(final_tbl_means)

library(knitr)

# Create a LaTeX table
latex_table <- kable(final_tbl_means, format = "latex", booktabs = TRUE, 
                     caption = "Means Across Sims (Selected Columns in Days)")

# Print the LaTeX code
cat(latex_table)




library(dplyr)
library(knitr)
library(kableExtra)

# Columns to convert from years -> days
cols_years <- c(
  "czmk_survival", "csk_survival", "zom_survival", "obs_survival",
  "czmk_endpoint", "csk_endpoint", "zom_endpoint", "obs_endpoint"
)

# Compute means across sims, convert, round
final_tbl_means <- final_tbl %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  mutate(across(all_of(cols_years), ~ .x * 365.25)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2))) %>%
  as_tibble()

# Create LaTeX table with grouped headers
final_tbl_means %>%
  kable(format = "latex", booktabs = TRUE, longtable = FALSE, 
        caption = "Means Across Sims (Days Conversion)") %>%
  kable_styling(latex_options = c("hold_position", "scale_down")) %>%  # scale_down fits page width
  add_header_above(c(" " = 0, "Survival" = 4, "Endpoint" = 4))
