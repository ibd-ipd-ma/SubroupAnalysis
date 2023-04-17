# drug_ranking.R
# Two tasks are performed:
# 1) The drug class order (from best performing:drug1 to worst performing:drug3)
#    for the three modeled drug classes - interleukin-12/23 (il12), integrin (intg), 
#    and tumor necrosis factor-alpha (tnfi) - are determined for each participant
#    (row). 
# 2) Two-sample t-tests are applied to drug pairs (drug1 vs drug2 and drug2 vs drug3)
#    to determine if a participants significantly prefers one or more drug classes
#    over another.
#
# Functions: 
# - drug.ranking: main function that calls helper functions for calculating 1, 2
#
# - sort.by.drug.effectiveness: determines drug class order (1) for participant.
#   Requires columns il12.attrib, il12.se, intg.attrib, intg.se, tnfi.attrib, and
#   tnfi.se, which can be calculated from custom_predict::custom.lme4.predict(). 
#
# - calculate.drug.p.values: calculates two-sample t-tests for drug pairs, 
#   drug1 vs drug2 and drug2 and drug3.


#------------------------------------------------------------------------------#

library(dplyr)

#------------------------------------------------------------------------------#

drug.ranking <- function(data, df_list) {
  # ranks drugs (il12, intg, tnfi) for each row by magnitude and two-sample t-test
  # 
  # requires columns: il12.attrib, il12.se, intg.attrib, il12.se, tnfi.attrib, tnfi.se
  # calculated by custom_predict.R
  
  # sort by magnitude of predicted drug class effectiveness
  # -> drug1, drug2, drug3
  result <- sort.by.drug.effectiveness(data)
  
  # calculate p-values between drug1:drug2 and drug2:drug3
  result <- calculate.drug.p.values(result, df_list)
  
  return(result)
}

#------------------------------------------------------------------------------#

sort.by.drug.effectiveness <- function(data){
  # overview: sort drugs by predicted effectiveness (drug1 > drug2 > drug3)
  
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

two.sample.t.tests <- function(data, args, df_list){
  # data = data.frame of model covariates
  # args = c('drug1','drug2') or c('drug2','drug3')
  # function calculates p-value of a t-score between two drug classes (args)
  
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

#------------------------------------------------------------------------------#