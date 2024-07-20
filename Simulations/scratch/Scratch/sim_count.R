for (sim in 1:33){
  if (sim == 1){
      sim_count = sim
      prev_count = 0
    } else{
      sim_count = sim*2+prev_count
      prev_count = prev_count + 1
    }
  if(sim == 33){
    print(sim)
    message("last one")
    print(sim_count)
    print(prev_count)}}
