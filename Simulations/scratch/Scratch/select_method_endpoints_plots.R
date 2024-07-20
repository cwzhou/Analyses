# Generate the column names for each method
select_method_endpoints = function(method_vec, endpoint){
  method_variables <- lapply(method_vec, 
                             function(method) {
                               paste0(method, "_", endpoint)
                               }) %>% unlist()
  return(method_variables)
}

method_var = select_method_endpoints(method.nm.simple, Phase_lab)

# Pivot the data frame
beginning <- as.data.frame(a) %>%
  dplyr::select(rep = sim,
                all_of(var_method)) %>%
  pivot_longer(cols = all_of(method_var),
               names_to = "method",
               values_to = "value")
