library(ggplot2)
library(tidyr)

testingg = 0

if (testingg == 1){
  dat = trtdat %>% filter(sim < 100)
} else{
  dat = trtdat
}


trtsurv = dat %>% dplyr::select(sim, method, surv_A, surv_B) %>%
  mutate(A = as.numeric(surv_A),
         B = as.numeric(surv_B)) %>%
  mutate(diff = abs(A-B)) %>%
  dplyr::select(sim, method, diff)

mean(trtsurv$diff); sd(trtsurv$diff)

trtendpoint = dat %>% dplyr::select(sim, method, endpoint_A, endpoint_B) %>%
  mutate(A = as.numeric(endpoint_A),
         B = as.numeric(endpoint_B)) %>%
  mutate(diff = abs(A-B)) %>%
  dplyr::select(sim, method, diff)

# Plotting survival outcomes
p1 = ggplot(trtsurv, aes(x = method, y = diff, color = as.factor(sim))) +
  geom_point(position = position_dodge(width = 0.8)) +
  ggtitle("Comparison of Survival Outcomes") +
  labs(x = "Method", y = "Value") +
  theme_minimal() 

p2 = ggplot(trtendpoint, aes(x = method, y = diff, color = as.factor(sim))) +
  geom_point(position = position_dodge(width = 0.8)) +
  ggtitle("Comparison of Endpoint (CIF) Outcomes") +
  labs(x = "Method", y = "Value") +
  theme_minimal() 

ggsave("CR/figure/ComparingTrtDifferences_Survival.eps", p1, device="eps", width = 10, height = 10)
ggsave("CR/figure/ComparingTrtDifferences_CIF.eps", p2, device="eps", width = 10, height = 10)

diff_plot = function(dat){
  trtsurv = dat %>% dplyr::select(sim, method, surv_A, surv_B) %>%
    mutate(A = as.numeric(surv_A),
           B = as.numeric(surv_B)) %>%
    mutate(diff = abs(A-B)) %>%
    dplyr::select(sim, method, diff)
  
  # Plotting survival outcomes
  ggplot(trtsurv, aes(x = method, y = diff, color = as.factor(sim))) +
    geom_point(position = position_dodge(width = 0.8)) +
    ggtitle("Comparison of Survival Outcomes") +
    labs(x = "Method", y = "Value") +
    theme_minimal() 


trtep = dat %>% dplyr::select(sim, method, endpoint_A, endpoint_B)
trtsurv = dat %>% dplyr::select(sim, method, surv_A, surv_B) 
trtsurv_long <- gather(trtsurv, key = "treatment", value = "value", -c(sim, method))
trtsurv_long <- trtsurv_long %>% mutate(treatment = sub(".*_", "", treatment))
trtsurv_long$value = as.numeric(trtsurv_long$value)

# trtsurv_long = trtsurv_long %>% mutate(value = ifelse(is.nan(value),999,value))


# Plotting survival outcomes
ggplot(trtsurv_long, aes(x = method, y = value, color = treatment)) +
  geom_point(position = position_dodge(width = 0.8)) +
  ggtitle("Comparison of Survival Outcomes") +
  labs(x = "Method", y = "Value") +
  facet_wrap(~ sim, scales = "free_y") +
  theme_minimal() 


# Define breaks for the y-axis ticks (including NaN if it exists)
breaks <- unique(c(trtsurv_long$value, NaN))

# Remove the y-axis tick marks for non-NaN values
your_plot <- your_plot 

# Plotting endpoint outcomes
ggplot(trtsurv_long, aes(x = method, y = value, fill = treatment)) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) +
  ggtitle("Comparison of Endpoint Outcomes") +
  labs(x = "Method", y = "Value") +
  facet_grid(treatment ~ sim, scales = "free_y") +
  theme_minimal()
