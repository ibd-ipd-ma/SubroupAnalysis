# custom_predict.R
# Due to complexities of mixed-effect models, a custom bootstrap was created to
# estimate 95% CIs for predictions, which uses the lme4::bootMer() function. 
#
# Models, such as m2_il12, m2_intg, and m2_tnfi, were fit with random-intercepts
# (Trial of participant). However, for prediction, this random-intercept is not
# needed, so `re.form` in lme4::predict() is set to `NA` (or 0) to set the 
# random effect of prediction to 0. 
#
# The bootMer() argument `use.u` indicating whether the spherical random effects 
# should be simulated / bootstrapped as well. If TRUE, they are not changed, and 
# all inference is conditional on these values. If FALSE, new normal deviates are 
# drawn. Since, we want to accurately prognosticate the counterfactuals for any 
# given patient in any given clinical setting, we set use.u = TRUE in our 
# bootstrap (spherical random effects are not changed).
# SOURCE: https://www.rdocumentation.org/packages/lme4/versions/1.1-31/topics/bootMer
#
# Lastly, we use multiple cores (parallel) for parallel processing when executing
# the bootstrap (especially with a large nsim, ex. 1_000 or 10_000). 

#------------------------------------------------------------------------------#

library(dplyr)
library(lme4)     # bootMer() for bootstrapping  
library(parallel) # parallel processing

# m1_plac <- readRDS('models/m1_plac.rds') 
# m2_il12 <- readRDS('models/m2_il12.rds')
# m2_intg <- readRDS('models/m2_intg.rds')
# m2_tnfi <- readRDS('models/m2_tnfi.rds')

#------------------------------------------------------------------------------#

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
                                nsim  = 10000, 
                                parallel = 'multicore', 
                                ncpus = 1, 
                                use.u = TRUE, 
                                seed  = 1234){
  
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
                             use.u    = use.u,        # pre-defined
                             type     = "parametric", # pre-defined
                             seed     = seed,
                             parallel = parallel, 
                             ncpus    = ncpus)
  
  # collapse bootstrap into 95% CI (fit, lwr, upr)
  boot.fit <- sumBoot( boot.pred )
  
  # append prediction to data=
  data[[ paste0(model.name, '.attrib') ]] = boot.fit$fit 
  data[[ paste0(model.name, '.se') ]]     = (boot.fit$upr - boot.fit$lwr) / 3.92
  
  return(data)
}

#------------------------------------------------------------------------------#

# Helper function to bootstrap all models
# Ensure models are loaded prior to run
lme4.predict.all <- function(data, 
                             m1_plac, m2_il12, m2_intg, m2_tnfi,
                             nsim=100, 
                             parallel='no', 
                             ncpus=1, 
                             seed=1234) {
  
  if( !missing(m1_plac) ) {
    # find plac.fit, plac.se
    data <- custom.lme4.predict(lme4.model = m1_plac,
                                data = data,
                                model.name = 'plac',
                                nsim = nsim,
                                parallel = parallel,
                                ncpus = ncpus,
                                seed = seed)
    print("Placebo fit and SE calculated.")
  } else {
    print("Placebo model not provided (m1_plac). Fit and SE not calculated.")
  }
  
  if( !missing(m2_il12) ) {
    # find il12.fit, il12.se
    data <- custom.lme4.predict(lme4.model = m2_il12,
                                data = data,
                                model.name = 'il12',
                                nsim = nsim,
                                parallel = parallel,
                                ncpus = ncpus,
                                seed = seed)
    print("Anti-Il-12/23 fit and SE calculated.")
  } else {
    print("Anti-Il-12/23 model not loaded (m2_il12). Fit and SE not calculated.")
  }
  
  if( !missing(m2_intg) ) {
    # find intg.fit,  intg.se
    data <- custom.lme4.predict(lme4.model = m2_intg,
                                data = data,
                                model.name = 'intg',
                                nsim = nsim,
                                parallel = parallel,
                                ncpus = ncpus,
                                seed = seed)
    print("Anti-Integrin fit and SE calculated.")
  } else {
    print("Anti-Integrin model not loaded (m2_intg). Fit and SE not calculated.")
  }
  
  if( !missing(m2_tnfi) ) {
    # find tnfi.fit, tnfi.se
    data <- custom.lme4.predict(lme4.model = m2_tnfi,
                                data = data,
                                model.name = 'tnfi',
                                nsim = nsim,
                                parallel = parallel,
                                ncpus = ncpus,
                                seed = seed)
    print("Anti-TNF fit and SE calculated.")
  } else {
    print("Anti-TNF model not loaded (m2_tnfi). Fit and SE not calculated.")
  }
  
  return(data)
}
