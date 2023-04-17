# subgroups.R
# Function that counts number of rows (from a data frame) that fall into one of 
# 13 possible subgroups for three drug class possibilities. 

#------------------------------------------------------------------------------#

library(dplyr)

#------------------------------------------------------------------------------#

# mask drugs where drug order does not matter (p == 0) - makes grouping easier
# ex. tnfi  > drug2 = drug3 
# ex. drug1 = drug2 = drug3
mask.results <- function(data) {
  masked <- data %>% 
    # case when there is no drug preference (p12_ohe = 0, p23_ohe = 0)
    mutate(drug1 = ifelse(p12_ohe == 0 & p23_ohe == 0, 'drug1', drug1),
           drug2 = ifelse(p12_ohe == 0 & p23_ohe == 0, 'drug2', drug2), 
           drug3 = ifelse(p12_ohe == 0 & p23_ohe == 0, 'drug3', drug3)) %>% 
    
    # case when drug1 > drug2 = drug3 (p12_ohe = 1, p23_ohe = 0)
    mutate(drug2 = ifelse(p12_ohe == 1 & p23_ohe == 0, 'drug2', drug2),
           drug3 = ifelse(p12_ohe == 1 & p23_ohe == 0, 'drug3', drug3)) %>% 
    
    # case when drug1 = drug2 > drug3 (p12_ohe = 0, p23_ohe = 1)
    mutate(drug1 = ifelse(p12_ohe == 0 & p23_ohe == 1, 'drug1', drug1),
           drug2 = ifelse(p12_ohe == 0 & p23_ohe == 1, 'drug2', drug2))
  
  return(masked)
}

find.subgroups <- function(data){
  # overview: find number of cases of 13 possible subgroups
  # 1.  tnfi > il12 > intg
  # 2.  il12 > > 
  # 3.  intg > > 
  # 4.  tnfi > drug2 = drug3
  # ...
  # 13. drug1 = drug2 = drug3
  
  subgroups <- mask.results(data) %>% 
    # collapse to 13 possible subgroups + add count
    group_by(drug1, drug2, drug3, p12_ohe, p23_ohe) %>% 
    summarise(n = n())
  
  return(subgroups)
}

#------------------------------------------------------------------------------#