library(dplyr)    # data manipulation
library(tidyr)    # data manipulation
library(lme4)     # linear mixed effect modeling
library(lmerTest) # provides p-values on coefficients in lme4 model summaries

#------------------------------------------------------------------------------#

setwd('UCSF/ipd-ma-cd2/publication')

################################################################################
# Drug-Class-Modeling
#   1. Train placebo model
#   2. Subtract placebo attributable effect from data (isolate cdai reduction due to drug only)
#   3. Train drug models
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
  # remove columns
  dplyr::select(-c(X, Cohort, Visit, Endpoint, Race:CurrentOrPriorStricture))

crohns_data1 %>% glimpse()

#------------------------------------------------------------------------------#

## 1. Placebo model - train placebo attributable effect model

### isolate placebo data (N = 1621)
placebo_df <- crohns_data1 %>%
  filter(Group == 'Placebo') %>% 
  dplyr::select(Trial, Year_Cent:Ileal, CDAI_reduction)
  
placebo_df %>% glimpse()
  
### build formula objects
dep_var = 'CDAI_reduction'
covariate_list <- names(placebo_df)[-which(names(placebo_df) %in% c('Trial',dep_var))]
f2 <- paste0( paste(dep_var, paste(covariate_list, collapse = '+'), sep='~'), '+(1|Trial)' )
f2 

### compute model
m1_plac <- lme4::lmer(f2, data = placebo_df)

### model information
summary( lmerTest::lmer(f2, data = placebo_df) )
ranef(m1_plac)

#------------------------------------------------------------------------------#

## 2. Find 'drug reduction' = 'CDAI_reduction' - 'Placebo_attributable' 

crohns_data1 <- crohns_data1 %>%
  
  # impute (predicted) placebo reduction using placebo model
  mutate(Placebo_attributable = predict(m1_plac, newdata=crohns_data1, re.form=~0), 
         .after='CDAI_reduction') %>%
  
  # calculate drug reduction (attributable effect)
  mutate(Drug_reduction = CDAI_reduction - Placebo_attributable, 
         .after='Placebo_attributable')

crohns_data1 %>% glimpse()

#------------------------------------------------------------------------------#

## 3. Drug class models

### isolate active data (N = 4082)
active_df <- crohns_data1 %>%
  filter(Group == 'Active')

### Split into 3 training data sets by drug class type (anti-TNF, anti-IL-12/23, anti-Integrin)
h2h_TNFi <- active_df %>% filter(TNFi_Active == 1)
h2h_Il12 <- active_df %>% filter(Il12_Active == 1) 
h2h_Intg <- active_df %>% filter(Integrin_Active == 1)

### build formula objects
covariate_list_2 <- covariate_list[-1] # remove Year_Cent
dep_var = 'Drug_reduction'
f3 <- paste(dep_var, paste(covariate_list_2, collapse = '+'), sep='~')
f4 <- paste0( paste(dep_var, paste(covariate_list_2, collapse = '+'), sep='~'), '+(1|Trial)' )
f4

### compute models 
m2_tnfi <- lme4::lmer(f4, h2h_TNFi)
m2_il12 <- lme4::lmer(f4, h2h_Il12)
m2_intg <- lme4::lmer(f4, h2h_Intg)

### model information
summary( lmerTest::lmer(f4, data=h2h_TNFi) )
ranef(m2_tnfi)
summary( lmerTest::lmer(f4, data=h2h_Il12) )
ranef(m2_il12)
summary( lmerTest::lmer(f4, data=h2h_Intg) )
ranef(m2_intg)

#------------------------------------------------------------------------------#

## Save
saveRDS( m1_plac ,  file = 'models/m1_plac.rds')
saveRDS( m2_il12 ,  file = 'models/m2_il12.rds')
saveRDS( m2_intg ,  file = 'models/m2_intg.rds')
saveRDS( m2_tnfi ,  file = 'models/m2_tnfi.rds')
write.csv(crohns_data1, 'data/crohns_data1.csv')
