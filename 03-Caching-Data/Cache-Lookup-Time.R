library(dplyr)
library(lme4)
library(parallel)

#------------------------------------------------------------------------------#

setwd('UCSF/ipd-ma-cd2')

m1_plac <- readRDS('models/m1_plac.rds') 
m2_il12 <- readRDS('models/m2_il12.rds')
m2_intg <- readRDS('models/m2_intg.rds')
m2_tnfi <- readRDS('models/m2_tnfi.rds')

source('helper/benchmark.R')                # benchmark, time
source('helper/custom_predict.R')           # lme4.predict.all
source('helper/nearest_cached_neighbor.R')  # findNearestNeighbor

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
  dplyr::select(Year_Cent, CDAI_baseline_Cent:Ileal)

#------------------------------------------------------------------------------#

## BOOTSTRAP TIME

# benchmark (from benchmark.R) replicates run n = 10 times 

benchmark({ time(lme4.predict.all(data = crohns_data1, nsim = 100, ncpus = 1)) })
#    min   mean    max 
# 4.4170 4.4844 4.6110 

benchmark({ time(lme4.predict.all(data = crohns_data1, nsim = 1000, ncpus = 1)) })
#     min    mean     max 
# 27.4330 27.6731 28.1790 

benchmark({ time(lme4.predict.all(data = crohns_data1, nsim = 10000, ncpus = 1)) })
# 
#

#------------------------------------------------------------------------------#

## CACHE LOOKUP TIME

# load cached data 
cache_sm_n100    <- read.csv('data-cache/cache_sm_100.csv')
cache_med_n100   <- read.csv('data-cache/cache_med_100.csv')
cache_lg_n100    <- read.csv('data-cache/cache_lg_100.csv')  # 2 GB; takes long time to load

# sample patient
test_pat_raw <- crohns_data1[1,]
test_pat_sm  <- find.nearest.neighbor(test_pat_raw, 'sm')
test_pat_med <- find.nearest.neighbor(test_pat_raw, 'med')
test_pat_lg  <- find.nearest.neighbor(test_pat_raw, 'lg')

rbind(test_pat_raw, test_pat_sm, test_pat_med, test_pat_lg)

## SM CACHE LOOKUP
benchmark({ time(cache_sm_n100 %>% right_join(., test_pat_sm)) })
#    min   mean    max 
# 0.0160 0.0166 0.0170

## MED CACHE LOOKUP
benchmark({ time(cache_med_n100 %>% right_join(., test_pat_med)) })
#    min   mean    max 
# 0.0440 0.0452 0.0460 

## LG CACHE LOOKUP
benchmark({ time(cache_lg_n100 %>% right_join(., test_pat_lg)) })
#    min   mean    max 
# 1.8740 1.8841 1.8950
