library(dplyr)
library(tidyr)
library(lme4)     # bootMer() for bootstrapping  
library(parallel) # parallel processing

#------------------------------------------------------------------------------#

setwd('UCSF/ipd-ma-cd2/publication')

m1_plac <- readRDS('models/m1_plac.rds')
m2_il12 <- readRDS('models/m2_il12.rds')
m2_intg <- readRDS('models/m2_intg.rds')
m2_tnfi <- readRDS('models/m2_tnfi.rds')

# see helper files for more details
source('helper/custom_predict.R') # lme4.predict.all()
source('helper/drug_ranking.R')   # drug.ranking()
source('helper/subgroups.R')      # find.subgroups()

################################################################################
# Subgroup-Analysis
#   1. Bootstrap placebo and drug class predictions and prediction standard errors
#      (95% CI) for all trial-based participants (N = 5703). 
# 
#      Use custom_predict::lme4.predict.all() helper function. This function will output: 
#      plac.attrib, plac.se, il12.attrib, il12.se, intg.attrib, intg.se, tnfi.attrib, tnfi.se 
#      corresponding to each models' predicted disease reduction or CDAI reduction 
#      (fit) and standard error (95% CI). 
#
#   2. For each participant rank order drug classes by most to least effective 
#      (pairwise t-test). 
# 
#      Use drug_ranking::drug.ranking() helper function. This function will output: 
#      drug1, drug2, drug3, p12, p23 corresponding to the first, second, and third
#      most effective drug class and the pairwise t-test results between drug1-drug2 (p12)
#      and drug2-drug3 (p23). Pairwise t-test determines if one drug class is 
#      significantly more effective than another (ex. drug1 > drug2 when p12 < 0.05). 
#
#   3. Assign subgroup membership to all participants based on similar treatment
#      preferences. 
#
#      Use subgroups::find.subgroups() helper function. This function will output
#      the subgroup representation and the number of participants in that subgroup. 
################################################################################

## Data
crohns_data1 <- read.csv("data/crohns_data.csv") %>% 
  # rename
  rename(
    Year_Cent = Year_Norm, 
    CDAI_baseline_Cent = CDAI_baseline_Norm,
    Age_Cent = Age_Norm, 
    BMI_Cent = BMI_Norm,
    CRP_Cent = CRP_Norm
  ) %>% 
  # mutate
  mutate(Trial = NA) %>% 
  # remove columns
  dplyr::select(-c(X, Cohort, Visit, Endpoint, Race:CurrentOrPriorStricture))

crohns_data1 %>% glimpse()

#------------------------------------------------------------------------------#

## Bootstrap

### lme4 predict using parametric bootstrapping (use.u = TRUE) 
### if on windows, parallel process does not work; change parallel = 'no' and ncpus = 1

crohns_data1 <- lme4.predict.all(data     = crohns_data1, 
                                 m1_plac  = m1_plac,
                                 m2_il12  = m2_il12, 
                                 m2_intg  = m2_intg, 
                                 m2_tnfi  = m2_tnfi, 
                                 nsim     = 10000, 
                                 parallel = 'multicore', # = 'no' # (Windows)
                                 ncpus    = 8,           # = 1 # (Windows)
                                 seed     = 1234)

crohns_data1 %>% glimpse()

#------------------------------------------------------------------------------#

## Drug class ranking

### model degrees of freedom (for two-sample t-test evaluation : check if two 
### drug classes are significantly different from one another for a patient)
df_list <- list("il12" = nrow(m2_il12@frame) - 10, #  577 - 10
                "intg" = nrow(m2_intg@frame) - 10, # 1818 - 10
                "tnfi" = nrow(m2_tnfi@frame) - 10) # 1677 - 10

## rank drug class by patient (row)
drug_ordering <- drug.ranking(crohns_data1, df_list)

drug_ordering %>% glimpse()

#------------------------------------------------------------------------------#

## Aggregate patients into groups

### ex. il12 drug2 drug3 == Anti-IL-12/23 > (Anti-TNF = Anti-Integrin)
### In other words, 138 patients are significant responders to Anti-IL-12/23.

### ohe == one hot encoding; if p12 < 0.05 then p12_ohe = 1 else 0 

find.subgroups(drug_ordering)
#   drug1 drug2 drug3 p12_ohe p23_ohe     n
# 1 drug1 drug2 drug3       0       0  3144
# 2 drug1 drug2 il12        0       1     4
# 3 drug1 drug2 intg        0       1   355
# 4 il12  drug2 drug3       1       0   138
# 5 il12  tnfi  intg        1       1     1
# 6 tnfi  drug2 drug3       1       0  2017
# 7 tnfi  intg  il12        1       1    44

#------------------------------------------------------------------------------#

# Save
write.csv(drug_ordering, 'data/drug_ordering-n10000-useuT.csv')
