library(dplyr)    # data manipulation
library(lme4)     # linear mixed effect models
library(merTools) # predictInterval() for mixed effect prediction interval calculations

#------------------------------------------------------------------------------#

setwd('UCSF/ipd-ma-cd2/publication')

################################################################################
# Power Analysis
#
# We found a subgroup of patients who prefer anti-IL-12/23 drugs. These patients 
# primarily consist of older (>50 years old) women (~75%). We conduct a power 
# analysis on a sample of older patients (>50 years old) in a simulated head-to-head
# clinical trial with two treatment arms: ARM1 = receives anti-IL-12/23 and ARM2 = 
# receives anti-TNF. 
################################################################################

m1_plac <- readRDS('models/m1_plac.rds')
m2_il12 <- readRDS('models/m2_il12.rds')
m2_tnfi <- readRDS('models/m2_tnfi.rds')

## data frame
crohns_data1 <- read.csv("data/crohns_data.csv") %>% 
  # rename
  rename(
    CDAI_baseline_Cent = CDAI_baseline_Norm,
    Age_Cent = Age_Norm, 
    BMI_Cent = BMI_Norm,
    CRP_Cent = CRP_Norm
  ) %>% 
  mutate(Year_Cent = mean(Year_Norm)) %>% 
  dplyr::select(Year_Cent, CDAI_baseline_Cent:Ileal)

crohns_data1 %>% glimpse()

#------------------------------------------------------------------------------#

## Helper functions

# a. randomly sample clinical arm of size 'size'
sample_cohort <- function(x, size, replace=TRUE) {
  return( x[sample(nrow(x), size, replace=replace), ] )
}

# b. after generating two treatment arms, run drug class models to see observed 
#    improvement (m2_il12 = anti-IL-12/23, m2_tnfi = anti-TNF). 
PI.helper <- function(lme4.model, data, model.name, nsim=100) {
  
  # set RE to 0 (avg), ensures equal number of columns
  # SOURCE: https://github.com/jknowles/merTools/issues/60
  data <- data %>% mutate(Trial = merTools::averageObs(lme4.model)$Trial)
  
  PI <- predictInterval(merMod = lme4.model, newdata = data,
                        level = 0.95, n.sims = nsim, which = 'fixed',
                        stat = "median", type="linear.prediction",
                        include.resid.var = TRUE)
  
  fit.name <- paste(model.name, 'pi.fit', sep='.')
  sd.name <- paste(model.name, 'pi.sd', sep='.')
  
  data[fit.name] <- PI$fit
  data[sd.name]  <- (PI$upr - PI$lwr) / 3.92
  
  data <- data %>% dplyr::select(-Trial)
  
  return(data)
}

# c. simulate a clinical trial with two treatment arms: ARM1 = receives anti-IL-12/23 and 
#    ARM2 = receives anti_TNF. For each trial arm, calculate response rate (>= 100 reduction). 
sim_trial <- function(x, arm.size=100) {
  ARM1 <- sample_cohort(x, size=arm.size)       # random sample w replacement
  ARM1 <- PI.helper(m2_il12, ARM1, 'il12', 100) # drug class fit, PI
  ARM1 <- ARM1 %>% 
    mutate(
      plac.fit = predict(m1_plac, newdata=ARM1, re.form=~0), 
      Response = 1 - pnorm(100, # clinical response
                           mean = plac.fit + il12.pi.fit,
                           sd = il12.pi.sd))
  
  ARM2 <- sample_cohort(x, size=arm.size)
  ARM2 <- PI.helper(m2_tnfi, ARM2, 'tnfi', 100)
  ARM2 <- ARM2 %>% 
    mutate(
      plac.fit = predict(m1_plac, newdata=ARM2, re.form=~0), 
      Response = 1 - pnorm(100, # clinical response
                           mean = plac.fit + tnfi.pi.fit,
                           sd = tnfi.pi.sd))
  
  return(run_test(ARM1$Response, ARM2$Response, arm.size))
}

# d. determines if two treatment arms are significantly different from one another
#    with respect to treatment outcomes (clinical response).
run_test <- function(ARM1, ARM2, arm.size) {
  
  # Confidence interval = (x1–x2) +/- t*√((sp2/n1)+(sp2/n2))
  
  n = arm.size
  mu1 = mean(ARM1) 
  sd1 = sd(ARM1)
  mu2 = mean(ARM2)
  sd2 = sd(ARM2)
  
  sp = ((n-1)*sd1^2+(n-1)*sd2^2)/(n+n-2) # pooled sd
  margin <- qt(0.975, df=n+n-1)*sqrt(sp/n + sp/n)
  
  lwr = (mu1-mu2) - margin
  upr = (mu1-mu2) + margin
  
  print(c('m1'=mu1,'s1'=sd1,'m2'=mu2,'s2'=sd2,'lwr'=lwr,'upr'=upr))
  
  # confidence interval
  return( ifelse((sign(lwr)==sign(upr)) & (mu1 > mu2), 1, 0) )
}

set.seed(1234)
replicate(5, sim_trial(x=df1, arm.size=100))

#------------------------------------------------------------------------------#

## DATASET1: AGE >= 50
df1 <- crohns_data1 %>% filter((Age_Cent+35) >= 50) 
dim(df1) # 1040 / 5703
quantile(df1$Age_Cent, probs=c(0.25,0.50,0.75))

set.seed(1234)
t1.n100 <- replicate(1000, sim_trial(x=df1, arm.size = 100))
t1.n250 <- replicate(1000, sim_trial(x=df1, arm.size = 250))
t1.n500 <- replicate(1000, sim_trial(x=df1, arm.size = 500))

mean(t1.n100) # 0.585
mean(t1.n250) # 0.873
mean(t1.n500) # 0.974

#------------------------------------------------------------------------------#

## DATASET2: AGE >=50 & SEX_MALE == 0 (FEMALE)
df3 <- crohns_data1 %>% filter((Age_Cent+35) >= 50 & Sex_Male == 0)
dim(df3) # 577 / 5703
quantile(df3$Age_Cent, probs=c(0.25,0.50,0.75))

set.seed(1234)
t3.n100 <- replicate(1000, sim_trial(x=df3, arm.size = 100))
t3.n250 <- replicate(1000, sim_trial(x=df3, arm.size = 250))
t3.n500 <- replicate(1000, sim_trial(x=df3, arm.size = 500))

mean(t3.n100) # 0.755
mean(t3.n250) # 0.969
mean(t3.n500) # 0.999

#------------------------------------------------------------------------------#

## DATASET3: AGE >=50 | SEX_MALE == 1 (MALE)
df6 <- crohns_data1 %>% filter((Age_Cent+35) >= 50 & Sex_Male == 1)
dim(df6) # 463 / 5703
quantile(df6$Age_Cent, probs=c(0.25,0.50,0.75))

set.seed(1234)
t6.n100 <- replicate(1000, sim_trial(x=df6, arm.size = 100))
t6.n250 <- replicate(1000, sim_trial(x=df6, arm.size = 250))
t6.n500 <- replicate(1000, sim_trial(x=df6, arm.size = 500))

mean(t6.n100) # 0.346
mean(t6.n250) # 0.634
mean(t6.n500) # 0.807
