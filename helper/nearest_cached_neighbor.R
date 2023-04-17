# helper functions: find nearest cached result
nearest.value.neighbor <- function(data, var_full, size){
  # defined in Code-RShiny/cache_data/cache_data.R 
  steps = list('sm' = list('CDAI'=50, 'Age'=10, 'BMI'=5, 'CRP'=10),
               'med'= list('CDAI'=25, 'Age'= 5, 'BMI'=2, 'CRP'= 5),
               'lg' = list('CDAI'=10, 'Age'= 2, 'BMI'=1, 'CRP'= 1))
  
  var_data = list('CDAI' = list('offset'=300, 'min'=200, 'max'=450), 
                  'Age'  = list('offset'=35,  'min'=20 , 'max'=80), 
                  'BMI'  = list('offset'=20,  'min'=18 , 'max'=28), 
                  'CRP'  = list('offset'=10,  'min'=0  , 'max'=30))
  
  # get first part of var name (CDAI, Age, BMI, CRP)
  var_alias = matrix(unlist(strsplit(var_full, '_')))[1]
  
  # extract necessary helper variables
  step   = steps[[size]][[var_alias]]
  offset = var_data[[var_alias]]$offset
  min    = var_data[[var_alias]]$min
  max    = var_data[[var_alias]]$max
  
  X <- data[var_full]
  
  # real value = value + offset
  # nearest cached value = round(real/step)*step
  # re-center = -offset
  nearest <- round((X+offset)/step, 0)*step
  nearest <- ifelse(nearest < min, min, nearest)
  nearest <- ifelse(nearest > max, max, nearest)
  return(nearest-offset)
}

# main function
find.nearest.neighbor <- function(data, size){
  data$CDAI_baseline_Cent <- apply(X = data, FUN = nearest.value.neighbor,
                                   MARGIN = 1, var_full  = 'CDAI_baseline_Cent', size = size)
  
  data$Age_Cent <- apply(X = data, FUN = nearest.value.neighbor,
                         MARGIN = 1, var_full  = 'Age_Cent', size = size)
  
  data$BMI_Cent <- apply(X = data, FUN = nearest.value.neighbor,
                         MARGIN = 1, var_full  = 'BMI_Cent', size = size)
  
  data$CRP_Cent <- apply(X = data, FUN = nearest.value.neighbor,
                         MARGIN = 1, var_full  = 'CRP_Cent', size = size)
  
  if (size == 'sm') {
    # step size for BMI sm is weird 
    data <- data %>% 
      mutate(BMI_Cent = ifelse(BMI_Cent %in% c(0, 5), BMI_Cent-2, BMI_Cent))
  } 
  
  return(data)
}