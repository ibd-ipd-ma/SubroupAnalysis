# cache_data.R

# Function to 1) create a grid.expand of possible Crohn's Disease patients and 
# 2) predict drug class preferences (interleukin-12/23, integrin, TNFi) using
# generated models (m2_il12, m2_intg, m2_tnfi). The prediction uses a manually
# defined bootstrap, since the models are mixed-effect models and the lme4 
# package does not have a built-in predict() function. 

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(parallel))
suppressMessages(library(lme4))
library(optparse)
source('bootstrap.R')

# Run on command line to generate cached data sets:
# If using Mac/Linux, change parallel = 'multicore' and specify more ncpus based on computer specs

# Rscript cache_data.R --size sm -nsim 100 -parallel no -ncpus 1 --out cache_sm_100.csv
# Rscript cache_data.R --size sm -nsim 1000 -parallel no -ncpus 1 --out cache_sm_1000.csv
# Rscript cache_data.R --size sm -nsim 10000 -parallel no -ncpus 1 --out cache_sm_10000.csv

# Rscript cache_data.R --size med -nsim 100 -parallel no -ncpus 1 --out cache_med_100.csv
# Rscript cache_data.R --size med -nsim 1000 -parallel no -ncpus 1 --out cache_med_1000.csv
# Rscript cache_data.R --size med -nsim 10000 -parallel no -ncpus 1 --out cache_med_10000.csv

# Rscript cache_data.R --size lg -nsim 100 -parallel no -ncpus 1 --out cache_lg_100.csv
# lg nsim 1000 and 10000 broke my computers, so they're not included. 

#------------------------------------------------------------------------------#

# get input arguments
option_list = list(
  make_option(c("-s", "--size"), 
              type="character", 
              default='sm', 
              help="grid size c(sm, med, lg)", 
              metavar="character"),
  
  make_option(c("-n", "--nsim"), 
              type="numeric", 
              default=100, 
              help="bootstrap simulations [default= %default]", 
              metavar="numeric"),
  
  make_option(c("-p", "--parallel"), 
              type="character", 
              default='no', 
              help="parallel processing c(no, multicore)", 
              metavar="character"),
  
  make_option(c("-c", "--ncpus"), 
              type="numeric", 
              default=1, 
              help="number of CPUs [default= %default]", 
              metavar="numeric"),
  
  make_option(c("-o", "--out"), 
              type="character", 
              default="result.csv", 
              help="output file name [default= %default]", 
              metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#------------------------------------------------------------------------------#

# generate grid
steps = list('sm' = list('CDAI'=50, 'Age'=10, 'BMI'=5, 'CRP'=10),
             'med'= list('CDAI'=25, 'Age'= 5, 'BMI'=2, 'CRP'= 5),
             'lg' = list('CDAI'=10, 'Age'= 2, 'BMI'=1, 'CRP'= 1))

myGrid <- function(size, steps) {

  # sm:     16,128
  # med:   192,192
  # lg:  8,795,072

  data.grid <- data.frame( expand.grid(
    CDAI_baseline_Cent = seq(200, 450, steps[size][[1]]$CDAI) - 300,
    Age_Cent           = seq( 20,  80, steps[size][[1]]$Age)  - 35,
    BMI_Cent           = seq( 18,  28, steps[size][[1]]$BMI)  - 20,
    CRP_Cent           = seq(  0,  30, steps[size][[1]]$CRP)  - 10,
    Sex_Male           = c(0,1),
    HxOfTNFi           = c(0,1),
    ImmUse             = c(0,1),
    SteroidUse         = c(0,1),
    Ileal              = c(0,1)
  ) )

  return(data.grid)
}

#------------------------------------------------------------------------------#

# run bootstrap
data <- myGrid(opt$size, steps)
print(paste('Created grid of size', nrow(data)))
system.time( result <- RunCustomBootstrap(data     = data,
                                          nsim     = opt$nsim,
                                          parallel = opt$parallel,
                                          ncpus    = opt$ncpus,
                                          seed     = 1234) )

# save output csv
# move files to ./data folder
write.csv(result, opt$out)
print('Saved results.')
