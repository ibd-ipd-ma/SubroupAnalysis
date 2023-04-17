# SubroupAnalysis

### BACKGROUND
Crohn’s disease (CD) is characterized by diverse manifestations that reflect intrinsic differences in disease biology and treatment responsiveness. Prior meta-analyses have found anti-tumor necrosis factor (TNF) drugs to be more efficacious than other drug classes. However, these findings are based on cohort-averaged effects and have ignored the role of patient-level variation in determining treatment outcomes.

### OBJECTIVE
To model personalized treatment options for patients with CD and characterize patient subgroups using participant-level data from randomized trials. 

### METHODS
We obtained participant-level data from 15 trials (N=5703) involving adults with moderate-to-severe CD (Crohn’s Disease Activity Index (CDAI) 220 to 450) treated with a biologic (anti-TNFs, anti-IL-12/23s, or anti-integrins). After normalizing the data using sequential regression and simulation1, we used linear mixed-effects regression to model CDAI reduction as a function of drug class, demographics, and disease-related features. We used these models to simulate patient-level outcomes following treatment with different drug classes. We performed hypothesis tests to compare these potential outcomes and to quantify the evidence favoring personalized treatment selection over several null models, 1) that each drug class would have the same efficacy irrespective of patient-level variation (one-size-fits-all), and 2) that pairs of drug classes would have the same efficacy in each individual patient. We used a p-value threshold of 0.05 to categorize patients into subgroups defined by treatment preferences. We performed queries on the University of California Health (UCH) Data Warehouse to measure the prevalence of these subgroups and quantify the potential real-world impact of these findings. Lastly, we prototyped a decision support tool that uses manual inputs and OMOP-formatted data to recommend treatments for individual patients.

### RESULTS
Models using patient-level features to predict drug class efficacy were superior fits to the data compared to a null model that assumes that these drugs produce the same efficacy in all patients (p<0.001). 45% (N=2559) of the cohort had statistical evidence for assignment into one of six treatment subgroups. While anti-TNFs were found to be most efficacious (or tied for first) in many patients (N=2420), the majority did not have statistical evidence favoring any one drug class (55%, N=3144). 

We identified a subgroup of older women who responded best to anti-IL-12/23, with 50% achieving a CDAI reduction  100, compared to only 3% with anti-TNF treatment. Although this subgroup corresponds to only 2% (N=139) of the trial-based cohort, 20% (N=4531) of Crohn’s patients treated at UCH (2012-2022) are women over 55. During the timeframe where all drug classes were FDA-approved, 72% of biologic-treated older women did not receive anti-IL-12/23 first line, suggesting significant room to optimize treatment allocation. 

### CONCLUSION
Our findings revise the prior literature and underscore the importance of patient-specific characteristics in treatment selection. Most patients do not appear to have superior efficacy with anti-TNFs, allowing for patients and providers to prioritize factors like safety and cost when selecting treatments. Elderly women may derive superior efficacy from anti-IL-12/23s. We have prototyped a decision support tool to advance precision medicine in CD.

## Data Access

The raw data are owned by the trial sponsors. The data may be accessed for reproduction and extension of this work following an application on the YODA and Vivli platforms and execution of a data use agreement.

## Requirements
Programming was performed in the R language (4.2.2), using the packages dplyr, tidyr, lme4, lmerTest, merTools, data.table, ggplot2, reshape2, RColorBrewer, grid, gridExtra, cowplot, sjplot, shiny, and shinydashboard.
