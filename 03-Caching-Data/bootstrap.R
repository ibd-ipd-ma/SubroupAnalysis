# bootstrap.R 

# Helper functions to perform mixed effect (lme4) bootstrapping on variable size
# input data (in this case pre-defined grid space of possible patients).

#------------------------------------------------------------------------------#

# load models
m2_il12 <- readRDS('../models/m2_il12.rds')
m2_intg <- readRDS('../models/m2_intg.rds')
m2_tnfi <- readRDS('../models/m2_tnfi.rds')

print('Random-effect models loaded successfully.')

#------------------------------------------------------------------------------#

RunCustomBootstrap <- function(data, 
                               nsim=100, 
                               parallel='no', 
                               ncpus=1, 
                               seed=1234) {
  
  # 3. Bootstrap fit, se for each model
  data <- custom.lme4.predict(lme4.model = m2_il12,
                              data = data,
                              model.name = 'il12',
                              nsim = nsim,
                              parallel = parallel,
                              ncpus = ncpus,
                              seed = seed)
  
  data <- custom.lme4.predict(lme4.model = m2_intg,
                              data = data,
                              model.name = 'intg',
                              nsim = nsim,
                              parallel = parallel,
                              ncpus = ncpus,
                              seed = seed)
  
  data <- custom.lme4.predict(lme4.model = m2_tnfi,
                              data = data,
                              model.name = 'tnfi',
                              nsim = nsim,
                              parallel = parallel,
                              ncpus = ncpus,
                              seed = seed)
  
  print('Bootstrapped input data.')
  
  # 4. Determine drug class ranking
  data.ranked <- sort.by.drug.effectiveness(data)
  
  # 5. Calculate p-values between drug1:drug2 and drug2:drug3
  df_list <- list("il12" = nrow(m2_il12@frame) - 10, #  577 - 10
                  "intg" = nrow(m2_intg@frame) - 10, # 1818 - 10
                  "tnfi" = nrow(m2_tnfi@frame) - 10) # 1677 - 10
  
  result <- calculate.drug.p.values(data.ranked, df_list)
  
  return(result)
}

#------------------------------------------------------------------------------#

## bootstrap

# bootstrap lme4 model to obtain fit and se based on pre-specified conditions:
# use.u    = TRUE                
# parallel = 'multicore'         (Mac/Linux only)
# lme4::predict(., re.form = NA) (set random effects of prediction to 0)
#
# input: lme4.model = lme4 object
#        data       = newdata 
#        model.name = c('tnfi','il12','intg')
#        nsim       = number of simulations
#        parallel   = c('no','multicore','snow')
#        ncpus      = number of CPUs 
#        seed       = 1234 (default)
custom.lme4.predict <- function(lme4.model, data, model.name, 
                                nsim=100, 
                                parallel = 'no', 
                                ncpus = 1, 
                                seed = 1234){
  
  # return predicted values from bootstrap (random effects set to 0)
  myPred <- function(.){ predict(., newdata = data, re.form = NA) }
  
  # collapse bootstrap into median, 95% confidence interval
  sumBoot <- function(merBoot){
    return(
      data.frame(fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
                 lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
                 upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
      )
    )
  }
  
  # bootstrap 
  boot.pred <- lme4::bootMer(lme4.model, 
                             myPred, 
                             nsim     = nsim, 
                             use.u    = TRUE,         # pre-defined
                             type     = "parametric", # pre-defined
                             seed     = seed,
                             parallel = parallel, 
                             ncpus    = ncpus)
  
  # collapse bootstrap into 95% CI (fit, lwr, upr)
  boot.fit <- sumBoot( boot.pred )
  
  # append prediction to data
  data[[ paste0(model.name, '.attrib') ]] = boot.fit$fit 
  data[[ paste0(model.name, '.se') ]]     = (boot.fit$upr - boot.fit$lwr) / 3.92
  
  return(data)
}


#------------------------------------------------------------------------------#

## determine drug class ranking

# sort drugs by predicted effectiveness (drug1 > drug2 > drug3)
sort.by.drug.effectiveness <- function(data){
  
  # isolate drug class attributable effects (3)
  data.attrib <- data %>% dplyr::select(il12.attrib, intg.attrib, tnfi.attrib)
  
  result <- cbind.data.frame(data, 
                             # for each row, sort drugs by desc magnitude of effectiveness 
                             t(apply(data.attrib, 1, function(row_i){
                               sort(row_i, decreasing = TRUE)})))
  
  # rename new columns
  names(result)[(ncol(result)-2):ncol(result)] <- c("drug1.attrib", "drug2.attrib", "drug3.attrib")
  
  # map predicted effectiveness (attrib or fit) to drug class name (tnfi, il12, or intg)
  result <- result %>% 
    mutate(drug1 = ifelse(drug1.attrib == tnfi.attrib, 'tnfi', 
                          ifelse(drug1.attrib == il12.attrib, 'il12', 'intg'))) %>%
    
    mutate(drug2 = ifelse(drug2.attrib == tnfi.attrib, 'tnfi', 
                          ifelse(drug2.attrib == il12.attrib, 'il12', 'intg'))) %>% 
    
    mutate(drug3 = ifelse(drug3.attrib == tnfi.attrib, 'tnfi', 
                          ifelse(drug3.attrib == il12.attrib, 'il12', 'intg')))
  
  return(result)
}

#------------------------------------------------------------------------------#

## p-value drug class comparisons

# data = data.frame of model covariates
# args = c('drug1','drug2') or c('drug2','drug3')
# function calculates p-value of a t-score between two drug classes (args)
two.sample.t.tests <- function(data, args, df_list){
  
  # initialize 
  X  <- c() # drug class fit
  SE <- c() # drug class standard error
  df <- 0   # sum of drug class degrees of freedom
  
  # for args (drug1, drug2, and/or drug3)
  for(drug in args){
    X  <- c(X,  data[[paste0(data[[drug]], '.attrib')]]) # ex. tnfi.attrib
    SE <- c(SE, data[[paste0(data[[drug]], '.se')]])     # ex. tnfi.se
    df <- df + df_list[[data[[drug]]]]
  }
  
  # convert vectors to numeric
  X <- as.numeric(X)
  SE <- as.numeric(SE)
  
  # Calculate p-value (pt)
  # p = 2 * pt( abs(X1 - X2) / sqrt(SE1^2 + SE2^2) )
  # df = df_X1 + df_X2
  p.value <- 2 * pt(abs(X[1] - X[2]) / sqrt(SE[1]^2 + SE[2]^2),
                    df = df,
                    lower.tail = F)
  
  return(p.value)
}

# calculate p-value of t-scores between different drug comparisons
calculate.drug.p.values <- function(data, df_list){
  
  # calculate p-value of t-score between drug1 and drug2
  data['p12'] <- apply(X = data, 
                       FUN = two.sample.t.tests, 
                       MARGIN = 1,                # row-wise
                       args = c('drug1','drug2'),
                       df_list = df_list)
  
  # calculate p-value of t-score between drug2 and drug3
  data['p23'] <- apply(X = data, 
                       FUN = two.sample.t.tests, 
                       MARGIN = 1,                # row-wise
                       args = c('drug2','drug3'),
                       df_list = df_list)
  
  # one-hot encode (ohe) p12 and p23
  result <- data %>% 
    mutate(
      p12_ohe = ifelse(p12 < 0.05, 1, 0),
      p23_ohe = ifelse(p23 < 0.05, 1, 0)
    )
  
  return(result)
}

