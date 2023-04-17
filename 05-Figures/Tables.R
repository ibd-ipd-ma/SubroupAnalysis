library(dplyr)
library(tidyr)
library(lme4)      # linear mixed effect modeling
library(sjPlot)    # Publication ready tables

#------------------------------------------------------------------------------#

setwd('UCSF/ipd-ma-cd2/publication')

m1_plac <- readRDS('models/m1_plac.rds')
m2_il12 <- readRDS('models/m2_il12.rds')
m2_intg <- readRDS('models/m2_intg.rds')
m2_tnfi <- readRDS('models/m2_tnfi.rds')

#------------------------------------------------------------------------------#

# Linear Mixed Effect Models (sjPlot)
tab_model(m1_plac, m2_il12, m2_intg, m2_tnfi,
          show.ci=F,
          show.se=T,
          pred.labels = c('Intercept','Year (Centered)','Baseline CDAI (Centered)','Age (Centered)','BMI (Centered)',
                          'CRP (mg/L) (Centered)','Sex: Male','HxOfTNFi','Steroid Use','Immunomodulator Use','Ileal Disease'),
          dv.labels = c('Placebo','Anti-Il-12/23','Anti-Integrin','Anti-TNF'),
          file='05-Figures/T2-Lienar-Mixed-Effect-Models.html')