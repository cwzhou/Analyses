# View(data.df)
# 
# data.df %>% group_by(status) %>% summarise(n = n())
# table(data.df$status, data.df$action)
# 
# # rep_czmk %>% group_by(status) %>% summarise(n = n())
# # rep_csk %>% group_by(status) %>% summarise(n = n())
# View(rep_csk)
# View(rep_czmk)
check = function(rep_obs){
  # rep_obs %>% group_by(status) %>% summarise(n = n()) %>% print()
  rep_obs %>% group_by(action) %>% summarise(n = n(), mean = mean(event.time)) %>% print()
  # print(table(rep_obs$status, rep_obs$action))
  # 
  # # test_obs = rep_obs %>% mutate(os = ifelse(status>0, 1,0))
  # # table(test_obs$os, test_obs$action)
  # 
  # rep_obs %>% dplyr::select(action, event.time) %>%
  #   arrange(action) %>%
  #   group_by(action) %>%
  #   summarise(mean = mean(event.time)) %>% print()
  # 
  # rep_obs %>% dplyr::select(status, event.time) %>%
  #   arrange(status) %>%
  #   group_by(status) %>%
  #   summarise(mean = mean(event.time)) %>% print()
  # 
  # 
  rep_obs %>% dplyr::select(action, status, event.time) %>%
    arrange(action,status) %>%
    group_by(status,action) %>%
    summarise(n = n(),mean = mean(event.time)) %>% print()

  os_mean = round(mean(rep_obs$event.time[rep_obs$status>0]))
  message("overall survival:", round(mean(rep_obs$event.time),3))
  cause1mean = rep_obs %>% filter(status == 1) %>%
    dplyr::select(event.time) %>% 
    summarise(mean = mean(event.time))
  message("cause1 survival:", round(cause1mean,3))
  
  # message("**")
  
}
# check(data.df)
check(rep_obs)
check(rep_czmk)
check(rep_csk)


