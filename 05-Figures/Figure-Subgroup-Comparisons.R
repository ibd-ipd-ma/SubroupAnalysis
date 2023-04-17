library(tidyverse)

library(ggplot2)
library(reshape2)     # melt()
library(RColorBrewer) # color palettes

library(grid)         # arrange plots
library(gridExtra)    # arrange plots
library(cowplot)      # common legend

#------------------------------------------------------------------------------#

setwd('UCSF/ipd-ma-cd2/publication')

# data (view 01-Analysis/03-Subgroup-Analysis.R for info on csv)
drug_ordering <- read.csv('data/drug_ordering-n10000-useuT.csv')

# see helper files for more details
source('helper/subgroups.R')      # determine patient subgroup (D1 > D2 = D3, etc.) 

################################################################################
# Subgroup Populations
################################################################################

subgroups <- find.subgroups(drug_ordering)
subgroups

# TNFI > (IL12, INTG)
prefer_tnfi <- drug_ordering %>% 
  filter(p12_ohe == 1 & drug1 %in% c('tnfi')) 

# (TNFI, IL12) > INTG
prefer_tnfi_il12 <- drug_ordering %>% 
  filter(p12_ohe == 0 & p23_ohe == 1 
         & drug1 %in% c('tnfi','il12') 
         & drug2 %in% c('tnfi','il12'))

# IL12 > (TNFI, INTG)
prefer_il12 <- drug_ordering %>% 
  filter(p12_ohe == 1 & drug1 %in% c('il12'))

#------------------------------------------------------------------------------#

drugClassSummaries <- function(subpopulation) {
  # returns N, placebo reduction average, and drug reduction average 
  # for given subpopulation - results used for bar plots
  
  # Anti-TNF summary
  tnfi_summary <- subpopulation %>% summarise(
    N = n(), 
    subgroup = 'Anti-TNF',
    placebo.avg = mean(plac.attrib), 
    drug.avg = mean(tnfi.attrib))    # change here
  
  # Anti-IL12/23
  il12_summary <- subpopulation %>% summarise(
    N = n(), 
    subgroup = 'Anti-IL-12-23',
    placebo.avg = mean(plac.attrib), 
    drug.avg = mean(il12.attrib))    # change here
  
  # Anti-Integrin
  intg_summary <- subpopulation %>% summarise(
    N = n(), 
    subgroup = 'Anti-Integrin',
    placebo.avg = mean(plac.attrib), 
    drug.avg = mean(intg.attrib))    # change here
  
  results <- rbind(data.frame(), tnfi_summary, il12_summary, intg_summary) %>% 
    # pivot longer to format for bar plots
    pivot_longer(cols = ends_with("avg"), names_to = "group", values_to = "avg")
  
  return( results )
}

# summarize drug class results - for bar plot
subgroup1_tnfi <- drugClassSummaries(prefer_tnfi) 
subgroup2_both <- drugClassSummaries(prefer_tnfi_il12) 
subgroup3_il12 <- drugClassSummaries(prefer_il12)

################################################################################
# Plots
################################################################################

# plot titles
titles <- c(
  sprintf("Prefer Anti-TNF Only (N = %d)",             nrow(prefer_tnfi)),
  sprintf("Prefer Anti-TNF or Anti-IL-12/23 (N = %d)", nrow(prefer_tnfi_il12)),
  sprintf("Prefer Anti-IL-12/23 Only (N = %d)",        nrow(prefer_il12))
)

# subgroup labels
labels <- c('Prefer Anti-TNF', 
            'Prefer Anti-TNF or Anti-IL-12/23', 
            'Prefer Anti-IL-12/23')

# subgroup colors
colors <- brewer.pal(3, "Set1")

#----------------------------- 
# Col 1 - Subgroup Bar Plots
#-----------------------------

# bar plot function
subgroupBarPlot <- function(subgroup_df, title = '', ylab = '', xlab = '', drug_color = 'black') {
  
  # add placebo.avg + drug.avg - print total at top of bar plot
  totals <- subgroup_df %>% 
    group_by(subgroup) %>% 
    summarize(total = round(sum(avg),1))
  
  # bar plot
  p <- subgroup_df %>% 
    group_by(subgroup) %>% 
    
    # add vars for stacked bar label (P, D) and height (middle of bar) 
    mutate(bar_label = ifelse(group == 'placebo.avg', 'P' ,'D')) %>%
    mutate(label_y = cumsum(avg) - 0.5 * avg) %>% 
    
    ggplot(aes(x = reorder(subgroup, -avg), y = avg, fill = group)) + 
    # bar plot
    geom_bar(stat = 'identity', width = 0.5, color = 'black') + 
    
    # add total average cdai reduction to top of each bar
    geom_text(data = totals, aes(x = subgroup, y = total, label = total, fill = NULL), nudge_y = 10) + 
    
    # add bar label (P, D) in middle of stacked bar
    geom_text(aes(y = label_y, label = bar_label), size = 4, color = 'grey') +
    
    # add title, labels
    ggtitle(title) + xlab(xlab) + ylab(ylab) +
    
    # manually fill bar colors 
    scale_fill_manual(
      name = NULL,
      labels = c("Drug Attributable", "Placebo Attributable"),
      values = c(drug_color,"white")) +
    
    # scale y axis - constant throughout all bar plots
    scale_y_continuous(limits = c(0,150)) + 
    
    # theme specifications
    theme_bw() + 
    theme(
      plot.title = element_text(size = 12),
      panel.grid.major.x = element_blank()) + 
    theme(legend.position = "none")
  
  p
}

# row 1 (r1) bar plots
r1.tnfi <- subgroupBarPlot(subgroup1_tnfi, 
                           title      = titles[1],
                           drug_color = colors[1]) + 
  geom_vline(xintercept = 1.5, linewidth = 1.5) # separate preferred vs not-preferred

r1.both <- subgroupBarPlot(subgroup2_both, 
                           title      = titles[2],
                           drug_color = colors[2]) + 
  geom_vline(xintercept = 2.5, linewidth = 1.5) # separate preferred vs not-preferred

r1.il12 <- subgroupBarPlot(subgroup3_il12, 
                           title      = titles[3],
                           drug_color = colors[3]) + 
  geom_vline(xintercept = 1.5, linewidth = 1.5) # separate preferred vs not-preferred

#----------------------------- 
# Col 2 - Binary Covariates (Bar)
#-----------------------------

# combine subgroups, add label column
combined_df <- rbind( 
    prefer_tnfi %>% mutate(label=labels[1]),
    prefer_tnfi_il12 %>% mutate(label=labels[2]),
    prefer_il12 %>% mutate(label=labels[3])
  ) %>% 
  # select covariates (continuous, binary) + label
  dplyr::select(CDAI_baseline:CRP, HxOfTNFi:Ileal, label) %>% 
  mutate(label = factor(label, levels = labels))

r2.binary <- combined_df %>% 
  
  # select binary covariates, label
  dplyr::select(HxOfTNFi:Ileal, label) %>% 
  
  # format (melt) data for bar plots
  group_by(label) %>%
  summarize(Female                    = 100 * ( 1 - sum(Sex_Male)/n() ),
            "Prior Anti-TNF Use"      = 100 * ( sum(HxOfTNFi)/n() ),
            "Immunomodulator Use"     = 100 * ( sum(ImmUse)/n() ),
            "Steroid Use"             = 100 * ( sum(SteroidUse)/n() ),
            "Disease Location: Ileum" = 100 * ( sum(Ileal)/n() ) ) %>% 
  reshape2::melt() %>% 
  
  # ggplot
  ggplot(aes(x = as.factor(variable), y = value,  fill = as.factor(label) )) + 
  # bar plot
  geom_bar(stat = 'identity', position = "dodge2", width = 0.4) + 
  scale_fill_manual(values = colors ) +
  
  # axis labels
  xlab("") + ylab("Percentage (%)") +
  
  # theme specifications  
  theme_bw() + 
  theme(
    legend.position = "none", 
    panel.grid.major.x = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  
  # dashed line at 50%
  geom_hline(yintercept=50, linetype = "dashed") + 
  facet_wrap(vars(label))

r2.female <- combined_df %>% 
  group_by(label) %>%
  summarize(Female                    = 100 * ( 1 - sum(Sex_Male)/n() )) %>% 
  reshape2::melt() %>% 
  ggplot(aes(x = as.factor(variable), y = value,  fill = as.factor(label) )) + 
  # bar plot
  geom_bar(stat = 'identity', position = "dodge2", width = 0.4) + 
  scale_fill_manual(values = colors ) +
  ggtitle("Female") + xlab("") + ylab("") +
  ylim(0,100) + 
  theme_bw() + 
  theme(
    legend.position = "none", 
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_blank()) + 
  geom_hline(yintercept=50, linetype = "dashed")

r2.tnfi <- combined_df %>% 
  group_by(label) %>%
  summarize("Prior Anti-TNF Use"      = 100 * ( sum(HxOfTNFi)/n() )) %>% 
  reshape2::melt() %>% 
  ggplot(aes(x = as.factor(variable), y = value,  fill = as.factor(label) )) + 
  # bar plot
  geom_bar(stat = 'identity', position = "dodge2", width = 0.4) + 
  scale_fill_manual(values = colors ) +
  ggtitle("Prior Anti-TNF Use") + xlab("") + ylab("") +
  ylim(0,100) + 
  theme_bw() + 
  theme(
    legend.position = "none", 
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_blank()) + 
  geom_hline(yintercept=50, linetype = "dashed")

r2.immune <- combined_df %>% 
  group_by(label) %>%
  summarize("Immunomodulator Use"     = 100 * ( sum(ImmUse)/n() )) %>% 
  reshape2::melt() %>% 
  ggplot(aes(x = as.factor(variable), y = value,  fill = as.factor(label) )) + 
  # bar plot
  geom_bar(stat = 'identity', position = "dodge2", width = 0.4) + 
  scale_fill_manual(values = colors ) +
  ggtitle("Current Immunomodulator Use") + xlab("") + ylab("") +
  ylim(0,100) + 
  theme_bw() + 
  theme(
    legend.position = "none", 
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_blank()) + 
  geom_hline(yintercept=50, linetype = "dashed")

r2.steroid <- combined_df %>% 
  group_by(label) %>%
  summarize("Steroid Use"             = 100 * ( sum(SteroidUse)/n() )) %>% 
  reshape2::melt() %>% 
  ggplot(aes(x = as.factor(variable), y = value,  fill = as.factor(label) )) + 
  # bar plot
  geom_bar(stat = 'identity', position = "dodge2", width = 0.4) + 
  scale_fill_manual(values = colors ) +
  ggtitle("Current Steroid Use") + xlab("") + ylab("") +
  ylim(0,100) + 
  theme_bw() + 
  theme(
    legend.position = "none", 
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_blank()) + 
  geom_hline(yintercept=50, linetype = "dashed")

r2.ileal <- combined_df %>% 
  group_by(label) %>%
  summarize("Disease Location: Ileum" = 100 * ( sum(Ileal)/n() )) %>% 
  reshape2::melt() %>% 
  ggplot(aes(x = as.factor(variable), y = value,  fill = as.factor(label) )) + 
  # bar plot
  geom_bar(stat = 'identity', position = "dodge2", width = 0.4) + 
  scale_fill_manual(values = colors ) +
  ggtitle("Disease Locatio: Ileum") + xlab("") + ylab("") +
  ylim(0,100) + 
  theme_bw() + 
  theme(
    legend.position = "none", 
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_blank()) + 
  geom_hline(yintercept=50, linetype = "dashed")


#----------------------------- 
# Col 3 - Continuous Covariates (Violin)
#-----------------------------

createViolinPlot <- function(combined_df, var = '', 
                             title='', xlab='', ylab=''){
  
  violin <- combined_df %>% 
    # ggplot
    ggplot(aes(x = label, y = .data[[var]], fill = label)) + 
    
    # violin plot
    geom_violin(alpha = 0.8) +
    
    # box-quantile plot (within violin plot)
    geom_boxplot(width = 0.1, color = "black", alpha = 0.2) +
    
    # title, labels
    ggtitle(title) + xlab(xlab) + ylab(ylab) + 
    
    # colors
    scale_fill_manual(values = colors) +
    
    # theme specifications
    theme_bw() + 
    theme(
      legend.position    = "none", 
      plot.title         = element_text(size = 12),
      axis.text.x        = element_blank(), 
      panel.grid.major.x = element_blank()) + 
    
    # rotate x axis labels 90 degrees
    scale_x_discrete(guide = guide_axis(angle = 90))
  
  violin
  
}


r3.cdai <- createViolinPlot(combined_df, 
                            var   = 'CDAI_baseline', 
                            title = 'CDAI Baseline')

r3.age <- createViolinPlot(combined_df, 
                           var   = 'Age', 
                           title = 'Age (yrs)')

r3.bmi <- createViolinPlot(combined_df, 
                           var   = 'BMI', 
                           title = 'BMI (kg/m2)')

r3.crp <- createViolinPlot(combined_df, 
                           var   = 'CRP', 
                           title = 'CRP (mg/L)')

################################################################################
# Combine Plots - cowplot
################################################################################

#----------------------------- 
# Legend
#-----------------------------

# dummy plot - just need legend details
p.legend <- data.frame(a = seq(1, length(labels)),
                       b = factor(labels, levels = labels) ) %>%
  ggplot(aes(y = a, fill = b)) +
  geom_bar() +
  scale_fill_manual(
    name = NULL,
    values = c(colors) ) + 
  theme(legend.position = "bottom")

# extract legend
overall.legend <- cowplot::get_legend(p.legend)

#----------------------------- 
# Plots by columns
#-----------------------------

# row 1
overall.c1 <- grid.arrange(
  # concat plots into row
  cowplot::plot_grid(r1.tnfi, r1.both, r1.il12, ncol = 1),
  left = textGrob("Average CDAI Reduction", gp=gpar(fontsize=11), rot=90, vjust = 1.5), 
  bottom = textGrob("Treatment Received", gp=gpar(fontsize=11), vjust = -1) )

# row 2
overall.c2 <- grid.arrange(
  # concat plots into row
  cowplot::plot_grid(r2.female, r2.tnfi, r2.immune, r2.steroid, ncol=1), 
  left = textGrob("Percentage (%)", gp=gpar(fontsize=11), rot=90, vjust = 1.5))

# row 3
overall.c3 <- grid.arrange(
  # concat plots into row
  cowplot::plot_grid(r3.age, r3.cdai, r3.bmi, r3.crp, ncol=1))

#-----------------------------
# Final Plot
#-----------------------------

# combined plot
p.final <- cowplot::plot_grid(
  cowplot::plot_grid(overall.c1, overall.c2, overall.c3, 
                     ncol = 3, nrow = 1, rel_widths = c(2, 1, 1)),
  overall.legend, 
  ncol = 1, nrow = 2, rel_heights = c(1, 0.1)
)

# save 
ggsave(filename = '01-Figures/F1-Sungroup-Comparisons.jpeg', 
       plot = p.final, 
       width = 15, height = 9, units='in')


