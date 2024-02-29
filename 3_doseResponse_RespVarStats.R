## USAGE: This script plots response variables for exploratory analysis 
## from the dose response bleaching database

# written by Harmony Martell 2023

# Loads in dataset, subsets variables
# Plots the density distribution and boxplots of each response variable by heated and control groupings, with points colored by study id
# Performs the statistical analyses of response variables

###############
rm(list = ls())
#load libraries
library(ggplot2)
library(tidyverse)
library(stats)
library(data.table)
library(lme4)
library(hrbrthemes)
library(gridExtra)
library(car)
library(rstatix)
library(dplyr)
library(psych)

# Set working directory
setwd("~/Desktop/UBC/_Projects/doseResponse_analysis/writing/FinalDocs_ERL/polishedAnalyses/")

# Load in datasets
  S = read.csv("symbiodiniaceae.csv", header=TRUE, stringsAsFactors=TRUE)
  C = read.csv("chlorophylla.ug.cm2.csv", header=TRUE, stringsAsFactors=TRUE)
Cpc = read.csv("chlorophylla.pg.cell1.csv", header=TRUE, stringsAsFactors=TRUE)
 Cc = read.csv("chlorophyllc2.ug.cm2.csv", header=TRUE, stringsAsFactors=TRUE)
 tC = read.csv("totalchlorophyll.csv", header=TRUE, stringsAsFactors=TRUE)
 Fv = read.csv("fvfm.csv", header=TRUE, stringsAsFactors=TRUE)

#################################################
### Response Variables Exploratory Statistics ###
#################################################
# 1. Summary stats for all & by heated and control treatment

# 2. Assumption checking (Shapiro Wilk and Variance Test)

# 3. Kolmogorov–Smirnov test tests for difference in shape (CDF)
 # Null: CDF of x does not lie below that of y (CDFs are not diff, shapes are the same)
 # Alt: CDF of x (control) lies below that of y (heated) (CDFs are different, shapes are not the same)

# 4. Wilcoxon Rank Sum test tests for difference in median distribution
 # Null: distribution of x is not shifted to the right of y (median distribution is not less)
 # Alt: distrubution of x (control) is shifted to the right of y (heated) (the median dist is less)

# 5. Kruskal Wallis test tests for explicit difference in the medians
 # Null: No difference in the medians of x and y
 # Alt: There is a difference in the medians of x (control) and y (heated)

library(plyr)
#SYMBIONT DENSITY
# 1. Summary stats & Parsing Treatments
summary(S)
sumStats_Sym <- describeBy(S$resp.var, S$trt); sumStats_Sym
fwrite(sumStats_Sym, "control.heated_symStats.csv")
heat <- S %>% filter(trt == "HEAT")
summary(heat)
control <- S %>% filter(trt == "CONTROL")
summary(control)

# 2. Assumption Checking
shapiro.test(heat$resp.var) # samples are non-normal, use non parametric, median comparisons
shapiro.test(control$resp.var) # samples are non-normal, use non parametric, median comparisons

par(mfrow=c(1,2))
qqPlot(heat$resp.var) # samples are non-normal
qqPlot(control$resp.var) # samples are non-normal

var.test(heat$resp.var, control$resp.var) # samples are not of equal variances

# 3. Kolmogorov–Smirnov test
S_ksresult <- ks.test(control$resp.var, heat$resp.var, alternative = "less", exact = TRUE); S_ksresult # the CDF distributions are the same
kstest <- c(S_ksresult$statistic, S_ksresult$p.value); kstest

# 4. Wilcoxon Rank Sum test 
S_wilcoxresult <- wilcox.test(control$resp.var, heat$resp.var, alternative = "greater", conf.int=TRUE, conf.level = .95, na.action = na.exclude ); S_wilcoxresult
wilcoxtest <- c(S_wilcoxresult$statistic, S_wilcoxresult$p.value); wilcoxtest
median(control$resp.var) - median(heat$resp.var)

# 5. Kruskal Wallis test
S_kwresult <- kruskal.test(resp.var ~ trt, data = S); S_kwresult
kwtest <- c(S_kwresult$method, S_kwresult$statistic, S_kwresult$parameter, S_kwresult$p.value); kwtest

write.csv(c(kwtest, kstest, wilcoxtest), file= "stats_Sym_control-heated.csv", na = " ")

### CHLOROPHYLL A (ug.cm-2)
# 1. Summary stats & Parsing Treatments
summary(C)
sumStats_Chla <- describeBy(C$resp.var, C$trt); sumStats_Chla
fwrite(sumStats_Chla, "control.heated_chlaStats.csv")
heat <- C %>% filter(trt == "HEAT")
summary(heat)
control <- C %>% filter(trt == "CONTROL")
summary(control)

# 2. Assumption Checking
shapiro.test(heat$resp.var) # samples are normal
shapiro.test(control$resp.var) # samples are non-normal, use non parametric, median comparisons

par(mfrow=c(1,2))
qqPlot(heat$resp.var) # samples are normal
qqPlot(control$resp.var) # samples are non-normal

var.test(heat$resp.var, control$resp.var) # samples are not of equal variances

# 3. Kolmogorov–Smirnov test
C_ksresult <- ks.test(control$resp.var, heat$resp.var, alternative = "less", exact = TRUE); C_ksresult # the CDF distributions are the same
kstest <- c(C_ksresult$statistic, C_ksresult$p.value); kstest

# 4. Wilcoxon Rank Sum test 
C_wilcoxresult <- wilcox.test(control$resp.var, heat$resp.var, alternative = "greater", conf.int=TRUE, conf.level = .95, na.action = na.exclude ); C_wilcoxresult
wilcoxtest <- c(C_wilcoxresult$statistic, C_wilcoxresult$p.value); wilcoxtest
median(control$resp.var) - median(heat$resp.var)

# 5. Kruskal Wallis test
C_kwresult <- kruskal.test(resp.var ~ trt, data = C); C_kwresult
kwtest <- c(C_kwresult$method, C_kwresult$statistic, C_kwresult$parameter, C_kwresult$p.value); kwtest

write.csv(c(kwtest, kstest, wilcoxtest), file= "stats_Chla.ugcm2_control-heated.csv", na = " ")


### CHLOROPHYLL A (pg.cell-1)
# 1. Summary stats & Parsing Treatments
summary(Cpc)
sumStats_ChlpCell <- describeBy(Cpc$resp.var, Cpc$trt); sumStats_ChlpCell
fwrite(sumStats_ChlpCell, "control.heated_chlapercellStats.csv")
heat <- Cpc %>% filter(trt == "HEAT")
summary(heat)
control <- Cpc %>% filter(trt == "CONTROL")
summary(control)

# 2. Assumption Checking
shapiro.test(heat$resp.var) # samples are normal
shapiro.test(control$resp.var) # samples are non-normal, use non parametric, median comparisons

par(mfrow=c(1,2))
qqPlot(heat$resp.var) # samples are non-normal
qqPlot(control$resp.var) # samples are non-normal

var.test(heat$resp.var, control$resp.var) # samples are not of equal variances

# 3. Kolmogorov–Smirnov test
Cpc_ksresult <- ks.test(control$resp.var, heat$resp.var, alternative = "less", exact = TRUE); Cpc_ksresult # the CDF distributions are the same
kstest <- c(Cpc_ksresult$statistic, Cpc_ksresult$p.value); kstest

# 4. Wilcoxon Rank Sum test 
Cpc_wilcoxresult <- wilcox.test(control$resp.var, heat$resp.var, alternative = "greater", conf.int=TRUE, conf.level = .95, na.action = na.exclude ); Cpc_wilcoxresult
wilcoxtest <- c(Cpc_wilcoxresult$statistic, Cpc_wilcoxresult$p.value); wilcoxtest
median(control$resp.var) - median(heat$resp.var)

1# 5. Kruskal Wallis test
Cpc_kwresult <- kruskal.test(resp.var ~ trt, data = Cpc); Cpc_kwresult
kwtest <- c(Cpc_kwresult$method, Cpc_kwresult$statistic, Cpc_kwresult$parameter, Cpc_kwresult$p.value); kwtest

write.csv(c(kwtest, kstest, wilcoxtest), file= "stats_ChlperCell_control-heated.csv", na = " ")


### CHLOROPHYLL C2 (ug.cm-2)
# 1. Summary stats & Parsing Treatments
summary(Cc)
sumStats_Chlc <- describeBy(Cc$resp.var, Cc$trt); sumStats_Chlc
fwrite(sumStats_Chlc, "control.heated_chlcStats.csv")
heat <- Cc %>% filter(trt == "HEAT")
summary(heat)
control <- Cc %>% filter(trt == "CONTROL")
summary(control)

# 2. Assumption Checking
shapiro.test(heat$resp.var) # samples are normal
shapiro.test(control$resp.var) # samples are non-normal, use non parametric, median comparisons

par(mfrow=c(1,2))
qqPlot(heat$resp.var) # samples are non-normal
qqPlot(control$resp.var) # samples are non-normal

var.test(heat$resp.var, control$resp.var) # samples are not of equal variances

# 3. Kolmogorov–Smirnov test
Cc_ksresult <- ks.test(control$resp.var, heat$resp.var, alternative = "less", exact = TRUE); Cc_ksresult # the CDF distributions are the same
kstest <- c(Cc_ksresult$statistic, Cc_ksresult$p.value); kstest

# 4. Wilcoxon Rank Sum test 
Cc_wilcoxresult <- wilcox.test(control$resp.var, heat$resp.var, alternative = "greater", conf.int=TRUE, conf.level = .95, na.action = na.exclude ); Cc_wilcoxresult
wilcoxtest <- c(Cc_wilcoxresult$statistic, Cc_wilcoxresult$p.value); wilcoxtest
median(control$resp.var) - median(heat$resp.var)

# 5. Kruskal Wallis test
Cc_kwresult <- kruskal.test(resp.var ~ trt, data = Cpc); Cc_kwresult
kwtest <- c(Cc_kwresult$method, Cc_kwresult$statistic, Cc_kwresult$parameter, Cc_kwresult$p.value); kwtest

write.csv(c(kwtest, kstest, wilcoxtest), file= "stats_Chlc_control-heated.csv", na = " ")


### TOTAL CHLOROPHYLL (ug.cm-2)
# 1. Summary stats & Parsing Treatments
summary(tC)
sumStats_tC <- describeBy(tC$resp.var, tC$trt); sumStats_tC
fwrite(sumStats_tC, "control.heated_totalChlStats.csv")
heat <- tC %>% filter(trt == "HEAT")
summary(heat)
control <- tC %>% filter(trt == "CONTROL")
summary(control)

# 2. Assumption Checking
shapiro.test(heat$resp.var) # samples are normal
shapiro.test(control$resp.var) # samples are normal

par(mfrow=c(1,2))
qqPlot(heat$resp.var) # samples are normal
qqPlot(control$resp.var) # samples are normal

var.test(heat$resp.var, control$resp.var) # samples are not of equal variances

# 3. Kolmogorov–Smirnov test
tC_ksresult <- ks.test(control$resp.var, heat$resp.var, alternative = "less", exact = TRUE); tC_ksresult # the CDF distributions are the same
kstest <- c(tC_ksresult$statistic, tC_ksresult$p.value); kstest

# 4. Wilcoxon Rank Sum test 
tC_wilcoxresult <- wilcox.test(control$resp.var, heat$resp.var, alternative = "greater", conf.int=TRUE, conf.level = .95, na.action = na.exclude ); tC_wilcoxresult
wilcoxtest <- c(tC_wilcoxresult$statistic, tC_wilcoxresult$p.value); wilcoxtest
median(control$resp.var) - median(heat$resp.var)

# 5. Kruskal Wallis test
tC_kwresult <- kruskal.test(resp.var ~ trt, data = tC); tC_kwresult
kwtest <- c(tC_kwresult$method, tC_kwresult$statistic, tC_kwresult$parameter, tC_kwresult$p.value); kwtest

write.csv(c(kwtest, kstest, wilcoxtest), file= "stats_totalChl_control-heated.csv", na = " ")


#FV/FM
# 1. Summary stats & Parsing Treatments
summary(Fv)
sumStats_Fv <- describeBy(Fv$resp.var, Fv$trt); sumStats_Fv
fwrite(sumStats_Fv, "control.heated_fvfmStats.csv")
heat <- Fv %>% filter(trt == "HEAT")
summary(heat)
control <- Fv %>% filter(trt == "CONTROL")
summary(control)

# 2. Assumption Checking
shapiro.test(heat$resp.var) # samples are non-normal, use non parametric, median comparisons
shapiro.test(control$resp.var) # samples are non-normal, use non parametric, median comparisons

par(mfrow=c(1,2))
qqPlot(heat$resp.var) # samples are non-normal
qqPlot(control$resp.var) # samples are non-normal

var.test(heat$resp.var, control$resp.var) # samples are not of equal variances

# 3. Kolmogorov–Smirnov test
Fv_ksresult <- ks.test(control$resp.var, heat$resp.var, alternative = "less", exact = TRUE); Fv_ksresult # the CDF distributions are the same
kstest <- c(Fv_ksresult$statistic, Fv_ksresult$p.value); kstest

# 4. Wilcoxon Rank Sum test 
Fv_wilcoxresult <- wilcox.test(control$resp.var, heat$resp.var, alternative = "greater", conf.int=TRUE, conf.level = .95, na.action = na.exclude ); Fv_wilcoxresult
wilcoxtest <- c(Fv_wilcoxresult$statistic, Fv_wilcoxresult$p.value); wilcoxtest
median(control$resp.var) - median(heat$resp.var)

# 5. Kruskal Wallis test
Fv_kwresult <- kruskal.test(resp.var ~ trt, data = Fv); Fv_kwresult
kwtest <- c(Fv_kwresult$method, Fv_kwresult$statistic, Fv_kwresult$parameter, Fv_kwresult$p.value); kwtest

write.csv(c(kwtest, kstest, wilcoxtest), file= "stats_FvFm_control-heated.csv", na = " ")


## FV/FM SUBSET
library(dplyr)
Match <- S %>% 
  filter(S$study.id %in% Fv$study.id)
ids <- unique(Match$study.id)
Fv_sym <- subset(Fv, sapply(study.id,  \(Fv) any(ids %in% Fv))) # keeping the studies that also reported symbiont densities

# 1. Summary stats & Parsing Treatments
summary(Fv_sym)
sumStats_Fv_sym <- describeBy(Fv_sym$resp.var, Fv_sym$trt); sumStats_Fv_sym
fwrite(sumStats_Fv_sym, "control.heated_fv.fm_subsetStats.csv")
heat <- Fv_sym %>% filter(trt == "HEAT")
summary(heat)
control <- Fv_sym %>% filter(trt == "CONTROL")
summary(control)

# 2. Assumption Checking
shapiro.test(heat$resp.var) # samples are non-normal, use non parametric, median comparisons
shapiro.test(control$resp.var) # samples are non-normal, use non parametric, median comparisons

par(mfrow=c(1,2))
qqPlot(heat$resp.var) # samples are non-normal
qqPlot(control$resp.var) # samples are non-normal

var.test(heat$resp.var, control$resp.var) # samples are not of equal variances

# 3. Kolmogorov–Smirnov test
Fv_ksresult <- ks.test(control$resp.var, heat$resp.var, alternative = "less", exact = TRUE); Fv_ksresult # the CDF distributions are the same
kstest <- c(Fv_ksresult$statistic, Fv_ksresult$p.value); kstest

# 4. Wilcoxon Rank Sum test 
Fv_wilcoxresult <- wilcox.test(control$resp.var, heat$resp.var, alternative = "greater", conf.int=TRUE, conf.level = .95, na.action = na.exclude ); Fv_wilcoxresult
wilcoxtest <- c(Fv_wilcoxresult$statistic, Fv_wilcoxresult$p.value); wilcoxtest
median(control$resp.var) - median(heat$resp.var)

# 5. Kruskal Wallis test
Fv_kwresult <- kruskal.test(resp.var ~ trt, data = Fv_sym); Fv_kwresulkwtest <- c(Fv_kwresult$method, Fv_kwresult$statistic, Fv_kwresult$parameter, Fv_kwresult$p.value); kwtest

write.csv(c(kwtest, kstest, wilcoxtest), file= "stats_FvFmsubset_control-heated.csv", na = " ")


#############
## FIGURES ##
#############

library(ggpubr)

# Set plotting parameters
black.bold.text<-element_text(face="bold",color="black", size=12)
small.text<-element_text(face="bold",color="black", size=10)
heatPalette <- c("#1A85FF","#D41159")

# 6. FIGURES
# SYMBIODINIACEAE
sym_all <- ggplot(data = S, aes(x = trt, y = resp.var), fill=study.id) +
  geom_boxplot() +
  geom_jitter(alpha = 0.7, aes(color = study.id), width = 0.32) +
  theme_bw() +
  labs(colour = "Study ID", x="", y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*') x 10'^6*"")))) +
  theme(legend.title = black.bold.text, axis.title=black.bold.text, axis.text = small.text)
sym_dist <- ggdensity(S, x = "resp.var",
                      add = "median", rug = TRUE,
                      color = "trt", fill = "trt",
                      palette = heatPalette) +
  theme_bw() +
  theme(legend.position=c(.85,.88), legend.title=element_blank(), strip.text.x = element_text(size=12, colour="black", face="bold"), strip.background = element_rect(colour="black", fill="white")) +
  labs(x=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*') x 10'^6*""))), y = "Density") +
  theme(axis.title=black.bold.text, axis.text = small.text)
S.plot <-grid.arrange(sym_dist, sym_all, nrow = 1)

#CHLOROPHYLL A
chla_all <- ggplot(data = C, aes(x = trt, y = resp.var), fill=study.id) +
  geom_boxplot() +
  geom_jitter(alpha = 0.7, aes(color = study.id), width = 0.32) +
  theme_bw() +
  labs(colour = "Study ID", x="", y=expression(bold(paste('Chlorophyll a (' ~mu *'g Chl cm'^-2*") ")))) +
  theme(axis.title=black.bold.text, axis.text = small.text, legend.title=black.bold.text)
chla_dist <- ggdensity(C, x = "resp.var",
                       add = "median", rug = TRUE,
                       color = "trt", fill = "trt",
                       palette = heatPalette) +
  theme_bw() +
  theme(legend.position="none", legend.title=element_blank(), strip.text.x = element_text(size=12, colour="black", face="bold"), strip.background = element_rect(colour="black", fill="white")) +
  labs(x=expression(bold(paste('Chlorophyll a (' ~mu *'g Chl cm'^-2*") "))), y = "Density") +
  theme(axis.title=black.bold.text, axis.text = small.text)
C.plot<-grid.arrange(chla_dist, chla_all, nrow = 1)

#CHLOROPHYLL A PER CELL
Cpc.plot.data <- Cpc[Cpc$resp.var < 100,] # one extreme outlier removed (S0067 Obs so you could see the values
chla_all <- ggplot(data = Cpc.plot.data, aes(x = trt, y = resp.var), fill=study.id) +
  geom_boxplot() +
  geom_jitter(alpha = 0.7, aes(color = study.id), width = 0.32) +
  theme_bw() +
  labs(colour = "Study ID", x="", y=expression(bold(paste('Chlorophyll a (pg Chl cell'^-1*") ")))) +
  theme(axis.title=black.bold.text, axis.text = small.text, legend.title=black.bold.text)
chla_dist <- ggdensity(Cpc.plot.data, x = "resp.var",
                       add = "median", rug = TRUE,
                       color = "trt", fill = "trt",
                       palette = heatPalette) +
  theme_bw() +
  theme(legend.position="none", legend.title=element_blank(), strip.text.x = element_text(size=12, colour="black", face="bold"), strip.background = element_rect(colour="black", fill="white")) +
  labs(x=expression(bold(paste('Chlorophyll a (pg Chl cell'^-1*") "))), y = "Density") +
  theme(axis.title=black.bold.text, axis.text = small.text)
Cpc.plot<-grid.arrange(chla_dist, chla_all, nrow = 1)

#CHLOROPHYLL C
chlc_all <- ggplot(data = Cc, aes(x = trt, y = resp.var), fill=study.id) +
  geom_boxplot() +
  geom_jitter(alpha = 0.7, aes(color = study.id), width = 0.32) +
  theme_bw() +
  labs(colour = "Study ID", x="", y=expression(bold(paste("Chlorophyll c"[2], " ("~mu *"g Chl cm"^-2* ")" )))) +
  theme(axis.title=black.bold.text, axis.text = small.text, legend.title=black.bold.text)
chlc_dist <- ggdensity(Cc, x = "resp.var",
                       add = "median", rug = TRUE,
                       color = "trt", fill = "trt",
                       palette = heatPalette) +
  theme_bw() +
  theme(legend.position="none", legend.title=element_blank(), strip.text.x = element_text(size=12, colour="black", face="bold"), strip.background = element_rect(colour="black", fill="white")) +
  labs(x=expression(bold(paste("Chlorophyll c"[2], " ("~mu *"g Chl cm"^-2* ")"))), y = "Density") +
  theme(axis.title=black.bold.text, axis.text = small.text)
Cc.plot<-grid.arrange(chlc_dist, chlc_all, nrow = 1)


#TOTAL CHLOROPHYLL
chl_all <- ggplot(data = tC, aes(x = trt, y = resp.var), fill=study.id) +
  geom_boxplot() +
  geom_jitter(alpha = 0.7, aes(color = study.id), width = 0.32) +
  theme_bw() +
  labs(colour = "Study ID", x="", y=expression(bold(paste("Total Chlorophyll ("~mu *"g Chl cm"^-2* ")" )))) +
  theme(axis.title=black.bold.text, axis.text = small.text, legend.title=black.bold.text)
chl_dist <- ggdensity(tC, x = "resp.var",
                       add = "median", rug = TRUE,
                       color = "trt", fill = "trt",
                       palette = heatPalette) +
  theme_bw() +
  theme(legend.position="none", legend.title=element_blank(), strip.text.x = element_text(size=12, colour="black", face="bold"), strip.background = element_rect(colour="black", fill="white")) +
  labs(x=expression(bold(paste("Total Chlorophyll ("~mu *"g Chl cm"^-2* ")" ))), y = "Density") +
  theme(axis.title=black.bold.text, axis.text = small.text)
tC.plot<-grid.arrange(chl_dist, chl_all, nrow = 1)


### FVFM
fvfm_all <- ggplot(data = Fv, aes(x = trt, y = resp.var), fill=study.id) +
  geom_boxplot() +
  geom_jitter(alpha = 0.5, aes(color = study.id), width = 0.32) +
  theme_bw() +
  labs(colour = "Study ID", x="", y=expression(bolditalic(paste('F'[V]*'/F'[M]*"" )))) +
  theme(axis.title=black.bold.text, axis.text = small.text, legend.title = black.bold.text)
fvfm_dist <- ggdensity(Fv, x = "resp.var",
                       add = "median", rug = TRUE,
                       color = "trt", fill = "trt",
                       palette = heatPalette) +
  theme_bw() +
  theme(legend.position="none", legend.title=element_blank(), strip.text.x = element_text(size=12, colour="black", face="bold"), strip.background = element_rect(colour="black", fill="white")) +
  labs(x=expression(bolditalic(paste('F'[V]*'/F'[M]*"" ))), y = "Density") +
  theme(axis.title=black.bold.text, axis.text = small.text)
Fv.plot<-grid.arrange(fvfm_dist, fvfm_all, nrow = 1)

fvfm_all <- ggplot(data = Fv_sym, aes(x = trt, y = resp.var), fill=study.id) +
  geom_boxplot() +
  geom_jitter(alpha = 0.7, aes(color = study.id), width = 0.32) +
  theme_bw() +
  labs(colour = "Study ID", x="", y=expression(bolditalic(paste('F'[V]*'/F'[M]*"" )))) +
  theme(axis.title=black.bold.text, axis.text = small.text, legend.title = black.bold.text)
fvfm_dist <- ggdensity(Fv_sym, x = "resp.var",
                       add = "median", rug = TRUE,
                       color = "trt", fill = "trt",
                       palette = heatPalette) +
  theme_bw() +
  theme(legend.position="none", legend.title=element_blank(), strip.text.x = element_text(size=12, colour="black", face="bold"), strip.background = element_rect(colour="black", fill="white")) +
  labs(x=expression(bolditalic(paste('F'[V]*'/F'[M]*"" ))), y = "Density") +
  theme(axis.title=black.bold.text, axis.text = small.text)
Fv.subsetplot<-grid.arrange(fvfm_dist, fvfm_all, nrow = 1)

## Big Grid Plots
grid.arrange(S.plot)
grid.arrange(C.plot)
grid.arrange(Cpc.plot)
grid.arrange(Fv.plot)
grid.arrange(Fv.subsetplot)
