zom1 = rep_zom %>% dplyr::select(subj.id, event.time, status, action)
czmk1 = rep_czmk %>% dplyr::select(subj.id, event.time, status, action)

ds1 = merge(czmk1, zom1, by = 'subj.id')

ds2 = ds1 %>% 
mutate(indCZMK = ifelse(action.x != action.y,
                     ifelse(event.time.x > event.time.y, 1, 0),
                     1
                     ))
View(ds2)

# mean survival time for CZMK when indCZMk = 1
mean(ds2$event.time.x[ds2$indCZMK==1])
# mean survival time for ZOM when indCZMk = 1
mean(ds2$event.time.y[ds2$indCZMK==1])

# mean survival time for CZMK when indCZMk = 0
mean(ds2$event.time.x[ds2$indCZMK==0])
# mean survival time for ZOM when indCZMk = 0
mean(ds2$event.time.y[ds2$indCZMK==0])