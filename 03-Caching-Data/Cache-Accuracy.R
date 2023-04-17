library(dplyr)
library(lme4)
library(parallel)

#------------------------------------------------------------------------------#

setwd('UCSF/ipd-ma-cd2/publication')

## data
crohns_data1 <- read.csv("data/crohns_data.csv") %>% 
  # rename
  rename(
    Year_Cent = Year_Norm, 
    CDAI_baseline_Cent = CDAI_baseline_Norm,
    Age_Cent = Age_Norm, 
    BMI_Cent = BMI_Norm,
    CRP_Cent = CRP_Norm
  ) %>% 
  # remove columns
  dplyr::select(CDAI_baseline_Cent:Ileal, CDAI_reduction)

crohns_data1 %>% glimpse()

#------------------------------------------------------------------------------#

# some files are too large. to replicate results, please run
# 03-Cached-Data::cache_data.R to recreate these datasets locally. 
cache_sm_n100    <- read.csv('data/cache_sm_100.csv')
cache_sm_n1000   <- read.csv('data/cache_sm_1000.csv')
cache_sm_n10000  <- read.csv('data-cache/cache_sm_10000.csv')

cache_med_n100   <- read.csv('data/cache_med_100.csv')
cache_med_n1000  <- read.csv('data/cache_med_1000.csv')
cache_med_n10000 <- read.csv('data/cache_med_10000.csv')

cache_lg_n100    <- read.csv('data/cache_lg_100.csv') # 2 GB

#------------------------------------------------------------------------------#

m2_il12 <- readRDS('models/m2_il12.rds')
m2_intg <- readRDS('models/m2_intg.rds')
m2_tnfi <- readRDS('models/m2_tnfi.rds')

df_list <- list("il12" = nrow(m2_il12@frame) - 10,
                "intg" = nrow(m2_intg@frame) - 10,
                "tnfi" = nrow(m2_tnfi@frame) - 10)

source('helper/benchmark.R')               # benchmark, time
source('helper/custom_predict.R')          # lme4.predict.all
source('helper/drug_ranking.R')            # drug.ranking
source('helper/nearest_cached_neighbor.R') # find.nearest.neighbor
source('helper/subgroups.R')               # mask.results

#------------------------------------------------------------------------------#

# test out find.nearest.neighbor

# crohns_data1$CDAI_baseline_Adjusted <- apply(X = crohns_data1, 
#                                                  FUN = find.nearest.neighbor, 
#                                                  MARGIN = 1, 
#                                                  var_full  = 'CDAI_baseline_Cent',
#                                                  size = 'sm')
# crohns_data1 %>% dplyr::select(CDAI_baseline_Cent, CDAI_baseline_Adjusted) %>% head()
#
##   CDAI_baseline_Cent CDAI_baseline_Adjusted
## 1           67.28924                     50
## 2          -43.18978                    -50
## 3          -76.50760                   -100
## 4           10.90406                      0
## 5          -72.00000                    -50
## 6           60.46832                     50

#------------------------------------------------------------------------------#

## calculate nearest neighbor for different step sizes
crohns_data1_sm  <- find.nearest.neighbor(crohns_data1, 'sm')
crohns_data1_med <- find.nearest.neighbor(crohns_data1, 'med')
crohns_data1_lg  <- find.nearest.neighbor(crohns_data1, 'lg')

## find ground truth
system.time( 
  truth_n100   <- lme4.predict.all(data = crohns_data1, 
                                   m2_il12=m2_il12,
                                   m2_intg=m2_intg, 
                                   m2_tnfi=m2_tnfi,
                                   nsim = 100)
)   

system.time( 
  truth_n1000  <- lme4.predict.all(data = crohns_data1, 
                                   m2_il12=m2_il12,
                                   m2_intg=m2_intg, 
                                   m2_tnfi=m2_tnfi,
                                   nsim = 1000)
)

system.time( 
  truth_n10000 <- lme4.predict.all(data = crohns_data1, 
                                   m2_il12=m2_il12,
                                   m2_intg=m2_intg, 
                                   m2_tnfi=m2_tnfi,
                                   nsim = 10000) 
) 

## rank and mask results
truth_n100   <- mask.results( drug.ranking( truth_n100 , df_list ) )
truth_n1000  <- mask.results( drug.ranking( truth_n1000 , df_list ) )
truth_n10000 <- mask.results( drug.ranking( truth_n10000 , df_list ) )

cache_sm_n100    <- mask.results(cache_sm_n100)
cache_sm_n1000   <- mask.results(cache_sm_n1000)
cache_sm_n10000  <- mask.results(cache_sm_n10000)
cache_med_n100   <- mask.results(cache_med_n100)
cache_med_n1000  <- mask.results(cache_med_n1000)
cache_med_n10000 <- mask.results(cache_med_n10000)
cache_lg_n100    <- mask.results(cache_lg_n100)

## create truth dataset with nearest neighbor 
truth_sm_n100    <- cbind(crohns_data1_sm[,1:9], truth_n100[,11:26])
truth_sm_n1000   <- cbind(crohns_data1_sm[,1:9], truth_n1000[,11:26])
truth_sm_n10000  <- cbind(crohns_data1_sm[,1:9], truth_n10000[,11:26])

truth_med_n100   <- cbind(crohns_data1_med[,1:9], truth_n100[,11:26])
truth_med_n1000  <- cbind(crohns_data1_med[,1:9], truth_n1000[,11:26])
truth_med_n10000 <- cbind(crohns_data1_med[,1:9], truth_n10000[,11:26])

truth_lg_n100    <- cbind(crohns_data1_lg[,1:9], truth_n100[,11:26])

#------------------------------------------------------------------------------#

## compare truth vs cache results 
covariates <- c('CDAI_baseline_Cent', 'Age_Cent', 'BMI_Cent', 'CRP_Cent', 
                'HxOfTNFi', 'Sex_Male', 'SteroidUse', 'ImmUse', 'Ileal')

## SM

# inner join : covariates                       : 5703
# inner join : covariates, drug1:drug3, p12, p23: 4753 (83%)
truth_sm_n100 %>% 
  select(all_of(covariates), drug1:drug3, p12_ohe, p23_ohe) %>%
  inner_join(., 
             cache_sm_n100 %>% select(all_of(covariates), drug1:drug3, p12_ohe, p23_ohe),
             by = c(covariates, 'drug1', 'drug2', 'drug3', 'p12_ohe', 'p23_ohe')) %>% nrow()

# inner join : covariates                       : 5703
# inner join : covariates, drug1:drug3, p12, p23: 4791 (84%)
truth_sm_n1000 %>% 
  select(all_of(covariates), drug1:drug3, p12_ohe, p23_ohe) %>%
  inner_join(., 
             cache_sm_n1000 %>% select(all_of(covariates), drug1:drug3, p12_ohe, p23_ohe),
             by = c(covariates, 'drug1', 'drug2', 'drug3', 'p12_ohe', 'p23_ohe')) %>% nrow()

# inner join : covariates                       : 5703
# inner join : covariates, drug1:drug3, p12, p23: 4816 (84%)
truth_sm_n10000 %>% 
  select(all_of(covariates), drug1:drug3, p12_ohe, p23_ohe) %>%
  inner_join(., 
             cache_sm_n10000 %>% select(all_of(covariates), drug1:drug3, p12_ohe, p23_ohe),
             by = c(covariates, 'drug1', 'drug2', 'drug3', 'p12_ohe', 'p23_ohe')) %>% nrow()

#------------------------------------------------------------------------------#

# MED

# inner join : covariates                       : 5703
# inner join : covariates, drug1:drug3, p12, p23: 5028 (88%)
truth_med_n100 %>% 
  select(all_of(covariates), drug1:drug3, p12_ohe, p23_ohe) %>%
  inner_join(., 
             cache_med_n100 %>% select(all_of(covariates), drug1:drug3, p12_ohe, p23_ohe),
             by = c(covariates, 'drug1', 'drug2', 'drug3', 'p12_ohe', 'p23_ohe')) %>% nrow()

# inner join : covariates                       : 5703
# inner join : covariates, drug1:drug3, p12, p23: 5045 (88%)
truth_med_n1000 %>% 
  select(all_of(covariates), drug1:drug3, p12_ohe, p23_ohe) %>%
  inner_join(., 
             cache_med_n1000 %>% select(all_of(covariates), drug1:drug3, p12_ohe, p23_ohe),
             by = c(covariates, 'drug1', 'drug2', 'drug3', 'p12_ohe', 'p23_ohe')) %>% nrow()

# inner join : covariates                       : 5703
# inner join : covariates, drug1:drug3, p12, p23: 5082 (89%)
truth_med_n10000 %>% 
  select(all_of(covariates), drug1:drug3, p12_ohe, p23_ohe) %>%
  inner_join(., 
             cache_med_n10000 %>% select(all_of(covariates), drug1:drug3, p12_ohe, p23_ohe),
             by = c(covariates, 'drug1', 'drug2', 'drug3', 'p12_ohe', 'p23_ohe')) %>% nrow()


#------------------------------------------------------------------------------#

## LG

# inner join : covariates                       : 5703
# inner join : covariates, drug1:drug3, p12, p23: 5168 (91%)
truth_lg_n100 %>% 
  select(all_of(covariates), drug1:drug3, p12_ohe, p23_ohe) %>%
  inner_join(., 
             cache_lg_n100 %>% select(all_of(covariates), drug1:drug3, p12_ohe, p23_ohe),
             by = c(covariates, 'drug1', 'drug2', 'drug3', 'p12_ohe', 'p23_ohe')) %>% nrow()

#------------------------------------------------------------------------------#
