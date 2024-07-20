bad_ind = which(rep_czmk$CIF_eval > rep_zom$CIF_eval | rep_czmk$OS_eval < rep_zom$OS_eval)

# View(rep_czmk[bad_ind,])
# View(rep_zom[bad_ind,])

df = result %>% dplyr::select(sim,czmk_survival, zom_survival, czmk_endpoint, zom_endpoint)
# Reshape the data into long format
library(tidyr)
df_long <- gather(df, key = "group", value = "value", -sim) %>%
  mutate(method = ifelse(grepl("czmk", group), "czmk","zom")) %>%
  mutate(group = factor(group, levels = c("czmk_survival", "zom_survival", "czmk_endpoint", "zom_endpoint")))

# Plot comparing group1 to group2 with different colors for each group
ggplot(df_long, aes(x = group, y = value, color = method)) +
  geom_point() +
  geom_boxplot() +
  facet_wrap(~sim) +   # only if sim is low!!
  labs(x = "Group", y = "Value", title = "Comparison of Group 1 and Group 2") + 
  theme(axis.text.x = element_text(angle = 10, hjust = 1))

table(rep_zom$action)
table(rep_czmk$action)
table(rep_csk$action)
