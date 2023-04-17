library(dplyr) # data manipulation
library(tidyr) # data manipulation
library(lme4)  # linear mixed effect modeling

#------------------------------------------------------------------------------#

setwd('UCSF/ipd-ma-cd2/publication')

m1_plac <- readRDS('models/m1_plac.rds')
m2_il12 <- readRDS('models/m2_il12.rds')
m2_intg <- readRDS('models/m2_intg.rds')
m2_tnfi <- readRDS('models/m2_tnfi.rds')

################################################################################
# Global p-values
#
# Likelihood ratio test (goodness of fit) of fit linear mixed effect 
# models (placebo, 3 drug class) and intercept-only models.   
################################################################################

## Data (view 01-Analysis/02-Drug-Class-Modeling.R for info on crohns_data1.csv)
crohns_data1 <- read.csv("data/crohns_data1.csv")

### isolate placebo data (N = 1621)
placebo_df <- crohns_data1 %>%
  filter(Group == 'Placebo') %>% 
  dplyr::select(Trial, Year_Cent:Ileal, CDAI_reduction)

### isolate active data (N = 4082)
active_df <- crohns_data1 %>%
  filter(Group == 'Active')

### Split into 3 training data sets by drug class type (anti-TNF, anti-IL-12/23, anti-Integrin)
h2h_TNFi <- active_df %>% filter(TNFi_Active == 1)
h2h_Il12 <- active_df %>% filter(Il12_Active == 1) 
h2h_Intg <- active_df %>% filter(Integrin_Active == 1)

#------------------------------------------------------------------------------#

## Intercept-only models
m0_plac <- lme4::lmer('CDAI_reduction ~ (1|Trial)', data=placebo_df)
m0_il12 <- lme4::lmer('Drug_reduction ~ (1|Trial)', data=h2h_Il12)
m0_intg <- lme4::lmer('Drug_reduction ~ (1|Trial)', data=h2h_Intg)
m0_tnfi <- lme4::lmer('Drug_reduction ~ (1|Trial)', data=h2h_TNFi)

### Likelihood ratio test (goodness of fit)
anova(m1_plac, m0_plac) # p < 2.2e-16 ***
anova(m2_tnfi, m0_tnfi) # p < 2.2e-16 ***
anova(m2_il12, m0_il12) # p < 0.001 ***
anova(m2_intg, m0_intg) # p < 0.01 **
