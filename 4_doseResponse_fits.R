##USAGE This script determines the relation between response variables and temperature variables
# from the dose response bleaching database using a linear and linear mixed effects models

# written by Harmony Martell 2023

# Loads in dataset
# Plots relation of response variables with heating temperature, duration, dhd & rate
    # symbiodiniaceae density (x 10^6 cm^-2)
    # chla (ug.cm^2 and pg.cell^-1)
    # fvfm (full dataset and subset)
    # NOTE: Total Chl and Chl c were excluded for lack of sufficient number of observations & studies
# Fits a curve to each relation
# performs an lme with study_id as the random effect

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
library(ggpubr)
source("~/Documents/Documents - Harmony’s MacBook Pro - 1/ODU/BARSHIS_LAB_ALL/Chapter1_Recent Thermal History/Data/RsquaredGLMM.R")

# Load in response vars
setwd("~/Desktop/UBC/_Projects/doseResponse_analysis/writing/FinalDocs_ERL/polishedAnalyses/")
  sym = read.csv("symbiodiniaceae_heated.csv", header=TRUE, stringsAsFactors=TRUE)
  chl.cm = read.csv("chlorophylla.ug.cm2_heated.csv", header=TRUE, stringsAsFactors=TRUE)
  chl.cell = read.csv("chlorophylla.pg.cell1_heated.csv", header=TRUE, stringsAsFactors=TRUE)
  fvfm.all = read.csv("fvfm_heated.csv", header=TRUE, stringsAsFactors=TRUE)
  fvfm.sub = read.csv("fvfm.sub_heated.csv", header=TRUE, stringsAsFactors=TRUE)

  summary(fvfm.all$time.of.day)
  summary(fvfm.all$dark.adapt.time)
  boxplot(fvfm.all$dark.adapt.time)
  unique(fvfm.all$name)
  
## Tidy each response variable df
  #keep only "study.id","resp.var","temp","duration","dhd","rate"
  sym = sym[,c(1,4:7,9)] 
  head(sym)
  summary(sym)
  
  chl.cm = chl.cm[,c(1,4:7,9)]
  head(chl.cm)
  summary(chl.cm)
  
  chl.cell = chl.cell[,c(1,4:7,9)]
  head(chl.cell)
  summary(chl.cell)
  
  fvfm.all = fvfm.all[,c(1,4:7,9)]
  head(fvfm.all)
  summary(fvfm.all)
  
  fvfm.sub = fvfm.sub[,c(1,4:7,9)]
  head(fvfm.sub)
  summary(fvfm.sub)
  
## Subset the data by MHW values 

 ## SYM
  # by the maximum temperature of any MHW
  library(dplyr)
  sym.mhw.temp <- sym %>%
    filter(max.temp.C < 35.84) # n = 66/67

  # by the minimum duration of any MHW
  sym.mhw.time <- sym %>% 
    filter(duration.days > 2) # n = 49/67
  
  # by temp and time of any MHW
  sym.mhw.temptime <- sym %>%
    filter(max.temp.C < 35.84 & duration.days > 2) # n = 48/67
  
  # by the maximum heating rate of any MHW
  sym.mhw.rate <- sym %>%
    filter(rate.C.day.1  < 1.52) # n = 40/67
  
  # by the maximum temp, time & heating rate of any MHW
  sym.mhw.all <- sym %>%
    filter(max.temp.C < 35.84) %>%
    filter(duration.days > 2) %>%
    filter(rate.C.day.1  < 1.52) # n = 34/67
    
 ## CHLA ug.cm^-2
  # by the maximum temperature of any MHW
  chlcm.mhw.temp <- chl.cm %>%
    filter(max.temp.C < 35.84) # n = 15/15
  
  # by the minimum duration of any MHW
  chlcm.mhw.time <- chl.cm %>% 
    filter(duration.days > 2) # n = 15/15
  
  # by temp and time of any MHW
  chlcm.mhw.temptime <- chl.cm %>%
    filter(max.temp.C < 35.84 & duration.days > 2) # n = 15/15
  
  # by the maximum heating rate of any MHW
  chlcm.mhw.rate <- chl.cm %>%
    filter(rate.C.day.1  < 1.52) # n = 13/15
  
  # by the maximum temp, time & heating rate of any MHW
  chlcm.mhw.all <- chl.cm %>%
    filter(max.temp.C < 35.84) %>%
    filter(duration.days > 2) %>%
    filter(rate.C.day.1  < 1.52) # n = 13/15
  
 ## CHLA pg.cell^-1
  # by the maximum temperature of any MHW
  chlcell.mhw.temp <- chl.cell %>%
    filter(max.temp.C < 35.84) # n = 19/19
  
  # by the minimum duration of any MHW
  chlcell.mhw.time <- chl.cell %>% 
    filter(duration.days > 2) # n = 19/19
  
  # by temp and time of any MHW
  chlcell.mhw.temptime <- chl.cell %>%
    filter(max.temp.C < 35.84 & duration.days > 2) # n = 19/19
  
  # by the maximum heating rate of any MHW
  chlcell.mhw.rate <- chl.cell %>%
    filter(rate.C.day.1  < 1.52) # n = 13/19
  
  # by the maximum temp, time & heating rate of any MHW
  chlcell.mhw.all <- chl.cell %>%
    filter(max.temp.C < 35.84) %>%
    filter(duration.days > 2) %>%
    filter(rate.C.day.1  < 1.52) # n = 13/19
  
 ## FV/FM full dataset
  # by the maximum temperature of any MHW
  fvfm.mhw.temp <- fvfm.all %>%
    filter(max.temp.C < 35.84) # n = 287/316
  
  # by the minimum duration of any MHW
  fvfm.mhw.time <- fvfm.all %>% 
    filter(duration.days > 2) # n = 291/316
  
  # by temp and time of any MHW
  fvfm.mhw.temptime <- fvfm.all %>%
    filter(max.temp.C < 35.84 & duration.days > 2) # n = 262/316
  
  # by the maximum heating rate of any MHW
  fvfm.mhw.rate <- fvfm.all %>%
    filter(rate.C.day.1 < 1.52) # n = 272/316
  
  # by the maximum temp, time & heating rate of any MHW
  fvfm.mhw.all <- fvfm.all %>%
    filter(max.temp.C < 35.84) %>%
    filter(duration.days > 2) %>%
    filter(rate.C.day.1 < 1.52) # n = 230/316
 
## FV/FM subset
  # by the maximum temperature of any MHW
  fvfmsub.mhw.temp <- fvfm.sub %>%
    filter(max.temp.C < 35.84) # n = 73/73
  
  # by the minimum duration of any MHW
  fvfmsub.mhw.time <- fvfm.sub %>% 
    filter(duration.days > 2) # n = 60/73
  
  # by temp and time of any MHW
  fvfmsub.mhw.temptime <- fvfm.sub %>%
    filter(max.temp.C < 35.84 & duration.days > 2) # n = 60/73 
  
  # by the maximum heating rate of any MHW
  fvfmsub.mhw.rate <- fvfm.sub %>%
    filter(rate.C.day.1 < 1.52) # n = 34/73

  fvfmsub.mhw.all <- fvfm.sub %>%
    filter(max.temp.C < 35.84) %>%
    filter(duration.days > 2) %>%
    filter(rate.C.day.1 < 1.52) # n = 33/73
  
  ########################
  #POLYNOMIAL REGRESSIONS#
  ########################
  # General equation: y = ax^2 - bx + c
  ### y is response variable (i.e., sym, chl.cm, chl.cell, fvfm)
  ### x is the predictor variable (i.e., max.temp, duration, dhd, rate)
  ### a-d are constants, the coefficients
  
  # First, set x and y
  
  ########################## SYMBIONTS #########################   
  ### y = Sym
  names(sym)
  y <- sym$resp.var # chl.cm$resp.var, chl.cell$resp.var, fvfm.all$resp.var
  g <- sym$study.id # set the random effect variable, study.id
  
  ### x = Temp
  x <- sym$max.temp.C # set the x (temp, dhd, duration, rate)
  # Sym =  -0.01603 x² + 0.89818 x - 11.37860
  #    R² = 0.06 marginal, 0.67 conditional
  #    p = (x²) 0.19, (x) 0.22 (NS)
  
  par(mfrow=c(1,1))
  plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g))
  # Examine regression summary details
  model
  print(summary(model))
  plot(model)
  # Examine the R squared values
  ##  - marginal (entirety of the variance explained by the fixed effects)
  ##  - conditional (portion of variance explained by the random effect of Study ID)
  Rsq <-rsquared.glmm(model); Rsq
  # Examine the names of the summary output table
  names(summary(model))
  # Examine the summary output table
  print(summary(model))
  # Perform the anova
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  # Write out the Anova Table to a file
  fwrite(Rsq,file = "sym_v_temp_Rsq.csv")
  write.csv(AovTbl, file= "sym_v_temp_Aov.csv", na = " ")
  
  ### x = Duration
  x <- sym$duration.days # set the x (temp, dhd, duration, rate)
  
  par(mfrow=c(1,1))
  plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<- rsquared.glmm(model); Rsq
  print(summary(model))
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "sym_v_duration_Rsq.csv")
  write.csv(AovTbl, file= "sym_v_duration_Aov.csv", na = " ")
  
  ### x = DHD
  x <- sym$dhd.C.days
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq <-rsquared.glmm(model); Rsq
  print(summary(model))
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "sym_v_dhd_Rsq.csv")
  write.csv(AovTbl, file= "sym_v_dhd_Aov.csv", na = " ")
  
  ### x = Heating Rate
  x <- sym$rate.C.day.1
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<- rsquared.glmm(model); Rsq
  print(summary(model))
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "sym_v_rate_Rsq.csv")
  write.csv(AovTbl, file= "sym_v_rate_Aov.csv", na = " ")
  
  # set plot params
  black.bold.text<-element_text(face="bold",color="black", size=12)
  small.text<-element_text(face="bold",color="black", size=10)
  legend.text<-element_text(face="bold",color="black", size=10)

  # 4 panel figure of Sym vs. Heating Metrics
  p1<- ggplot(sym, aes(x=max.temp.C, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    #theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), legend.text = legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=legend.text, legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
   # theme(axis.title=black.bold.text, axis.text=small.text, legend.position="right", legend.title=legend.text, legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(colour = "Study ID",title = "", x=expression(bold(paste("Maximum Temperature (",degree,"C)"))), y= "") # expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  p2<- ggplot(sym, aes(x=duration.days, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), legend.text = legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("Duration (days)"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
 
  p3<- ggplot(sym, aes(x=dhd.C.days, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("Degree Heating Days (",degree,"C days)"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  p4<- ggplot(sym, aes(x=rate.C.day.1, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title = "", x=expression(bold(paste("Heating Rate (",degree,"C day"^-1*")"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  S.plot.unf<- grid.arrange(p1,p2,p3,p4,nrow = 2)

  ########################## CHL A PER CM² #########################   
  ### y = Chlorophyll a (ug cm²)
  y <- chl.cm$resp.var # chl.cm$resp.var, chl.cell$resp.var, fvfm.all$resp.var
  g <- chl.cm$study.id # set the random effect variable, study.id
  
  ### x = Temp
  x <- chl.cm$max.temp.C # set the x
  
  par(mfrow=c(1,1))
  plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq <- rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cm_v_temp_Rsq.csv")
  write.csv(AovTbl, file= "chl.cm_v_temp_Aov.csv", na = " ")
  
  ### x = Duration
  x <- chl.cm$duration.days
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cm_v_duration_Rsq.csv")
  write.csv(AovTbl, file= "chl.cm_v_duration_Aov.csv", na = " ")
  
  ### x = DHD
  x <-chl.cm$dhd.C.days
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cm_v_dhd_Rsq.csv")
  write.csv(AovTbl, file= "chl.cm_v_dhd_Aov.csv", na = " ")
  
  ### x = Heating Rate
  x <- chl.cm$rate.C.day.1
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cm_v_rate_Rsq.csv")
  write.csv(AovTbl, file= "chl.cm_v_rate_Aov.csv", na = " ")

  # 4 panel figure of Chla ug cm^-2 vs. Heating Metrics
  p1<- ggplot(chl.cm, aes(x=max.temp.C, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    #theme(axis.title=black.bold.text, axis.text=small.text, legend.position="right", legend.title=legend.text, legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(colour = "Study ID", title="", x=expression(bold(paste("Maximum Temperature (",degree,"C)"))), y="" ) #expression(bold(paste('Chlorophyll a ('~mu *'g Chl cm'^-2*") "))) ) 
  
  p2<- ggplot(chl.cm, aes(x=duration.days, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("Duration (days)"))), y="") #y=expression(bold(paste('Chlorophyll a ('~mu *'g Chl cm'^-2*") ")))
  
  p3<- ggplot(chl.cm, aes(x=dhd.C.days, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("Degree Heating Days (",degree,"C days)"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  p4<- ggplot(chl.cm, aes(x=rate.C.day.1, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title = "", x=expression(bold(paste("Heating Rate (",degree,"C day"^-1*")"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  Chlcm.plot.unf<- grid.arrange(p1,p2,p3,p4,nrow = 2)
  
  
  ########################## CHL A PER CELL #########################  
  ### y = Chlorophyll a (pg cell-1)
  y <- chl.cell$resp.var
  g <- chl.cell$study.id # set the random effect variable, study.id
  
  ### x = Temp
  x <- chl.cell$max.temp.C # singularity
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cell_v_temp_Rsq.csv")
  write.csv(AovTbl, file= "chl.cell_v_temp_Aov.csv", na = " ")
  
  ### x = Duration
  x <- chl.cell$duration.days # singularity

  par(mfrow=c(1,1))
  plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cell_v_duration_Rsq.csv")
  write.csv(AovTbl, file= "chl.cell_v_duration_Aov.csv", na = " ")

  ### x = DHD
  x <- chl.cell$dhd.C.days # singularity
  
  par(mfrow=c(1,1))
  plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cell_v_dhd_Rsq.csv")
  write.csv(AovTbl, file= "chl.cell_v_dhd_Aov.csv", na = " ")
  # *** Is Singular fit, likely because one study has several entries very close to each other. Could simply do a polyfit to see? ***
  
  ### x = Heating Rate
  x <- chl.cell$rate.C.day.1 # singularity
  
  par(mfrow=c(1,1))
  plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cell_v_rate_Rsq.csv")
  write.csv(AovTbl, file= "chl.cell_v_rate_Aov.csv", na = " ")
 
  # 4 panel figure of Chla pg. cell-1 vs. Heating Metrics
  p1<- ggplot(chl.cell, aes(x=max.temp.C, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    #theme(axis.title=black.bold.text, axis.text=small.text, legend.position="right", legend.title=legend.text, legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(colour = "Study ID", title="", x=expression(bold(paste("Maximum Temperature (",degree,"C)"))), y= "" ) #expression(bold(paste('Chlorophyll a (pg Chl cm'^-2*") "))) ) 
  
  p2<- ggplot(chl.cell, aes(x=duration.days, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("Duration (days)"))), y="") #y=expression(bold(paste('Chlorophyll a ('~mu *'g Chl cm'^-2*") ")))
  
  p3<- ggplot(chl.cell, aes(x=dhd.C.days, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("Degree Heating Days (",degree,"C days)"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  p4<- ggplot(chl.cell, aes(x=rate.C.day.1, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title = "", x=expression(bold(paste("Heating Rate (",degree,"C day"^-1*")"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  Chlcell.plot.unf<- grid.arrange(p1,p2,p3,p4,nrow = 2)

   
########################## FV/FM #########################  
  ### y = Fv/Fm (all samples included)
  y <- fvfm.all$resp.var
  g <- fvfm.all$study.id
  
  ### x = Temp
  x <- fvfm.all$max.temp.C 
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "fvfm.all_v_temp_Rsq.csv")
  write.csv(AovTbl, file= "fvfm.all_v_temp_Aov.csv", na = " ")
  
  ### x = Duration
  x <- fvfm.all$duration.days # set the x
  
  par(mfrow=c(1,1));plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "fvfm.all_v_duration_Rsq.csv")
  write.csv(AovTbl, file= "fvfm.all_v_duration_Aov.csv", na = " ")
  
  ### x = DHD
  x <- fvfm.all$dhd.C.days
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "fvfm.all_v_dhd_Rsq.csv")
  write.csv(AovTbl, file= "fvfm.all_v_dhd_Aov.csv", na = " ")
  
  ### x = Heating Rate
  x <- fvfm.all$rate.C.day.1
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "fvfm.all_v_rate_Rsq.csv")
  write.csv(AovTbl, file= "fvfm.all_v_rate_Aov.csv", na = " ")

  # 4 panel figure of Chla ug cm^-2 vs. Heating Metrics
  p1<- ggplot(fvfm.all, aes(x=max.temp.C, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    #theme(axis.title=black.bold.text, axis.text=small.text, legend.position="right", legend.title=legend.text, legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(colour = "Study ID", title="", x=expression(bold(paste("Maximum Temperature (",degree,"C)"))), y= " ") #expression(bolditalic(paste('F'[V]*'/F'[M]*"" ))) )
  
  p2<- ggplot(fvfm.all, aes(x=duration.days, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("Duration (days)"))), y="") #y=expression(bold(paste('Chlorophyll a ('~mu *'g Chl cm'^-2*") ")))
  
  p3<- ggplot(fvfm.all, aes(x=dhd.C.days, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("Degree Heating Days (",degree,"C days)"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  p4<- ggplot(fvfm.all, aes(x=rate.C.day.1, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title = "", x=expression(bold(paste("Heating Rate (",degree,"C day"^-1*")"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  FvFm.plot.unf<- grid.arrange(p1,p2,p3,p4,nrow = 2)
  
  ## PLOT EVERYTHING
  grid.arrange(S.plot.unf, Chlcm.plot.unf, Chlcell.plot.unf, FvFm.plot.unf, nrow = 2)
   
 #########################################################
 ###    RERUN ANALYSES WITH MHW FILTERED DATASETS      ### 
 ######################################################### 

  ################### MHW-FILTERED SYMBIONTS #########################   
  ### y = Sym
  names(sym.mhw.all)
  y <- sym.mhw.all$resp.var # chl.cm$resp.var, chl.cell$resp.var, fvfm.all$resp.var
  g <- sym.mhw.all$study.id # set the random effect variable, study.id
  
  ### x = Temp
  x <- sym.mhw.all$max.temp.C # set the x (temp, dhd, duration, rate)
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g))
  # Examine regression summary details
  model
  print(summary(model))
  plot(model)
  # Examine the R squared values
  ##  - marginal (entirety of the variance explained by the fixed effects)
  ##  - conditional (portion of variance explained by the random effect of Study ID)
  Rsq <- rsquared.glmm(model); Rsq
  # Examine the names of the summary output table
  names(summary(model))
  # Examine the summary output table
  print(summary(model))
  # Perform the anova
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  # Write out the stuff to a file
fwrite(Rsq,file = "sym.mhw.all_v_temp_Rsq.csv")
write.csv(AovTbl, file= "sym.mhw.all_v_temp_AovTbl.csv", na = " ")
  
  ### x = Duration
  x <- sym.mhw.all$duration.days # set the x (temp, dhd, duration, rate)
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  print(summary(model))
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "sym.mhw.all_v_duration_Rsq.csv")
  write.csv(AovTbl, file= "sym.mhw.all_v_duration_Aov.csv", na = " ")
  
  ### x = DHD
  x <- sym.mhw.all$dhd.C.days

  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  print(summary(model))
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "sym.mhw.all_v_dhd_Rsq.csv")
  write.csv(AovTbl, file= "sym.mhw.all_v_dhd_Aov.csv", na = " ")
  
  ### x = Heating Rate
  x <- sym.mhw.all$rate.C.day.1
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  print(summary(model))
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "sym.mhw.all_v_rate_Rsq.csv")
  write.csv(AovTbl, file= "sym.mhw.all_v_rate_Aov.csv", na = " ")
  
# set plot params
  black.bold.text<-element_text(face="bold",color="black", size=18)
  small.text<-element_text(face="bold",color="black", size=14)
  legend.text<-element_text(face="bold",color="black", size=14)
  
  # 4 panel figure of mhw-filtered Sym vs. Heating Metrics
  p1<- ggplot(sym.mhw.all, aes(x=max.temp.C, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    #theme(axis.title=black.bold.text, axis.text=small.text, legend.position="right", legend.title=legend.text, legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(colour= "Study ID", title="", x=expression(bold(paste("Maximum Temperature (",degree,"C)"))), y = "" ) #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  p2<- ggplot(sym.mhw.all, aes(x=duration.days, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("Duration (days)"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  p3<- ggplot(sym.mhw.all, aes(x=dhd.C.days, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("Degree Heating Days (",degree,"C days)"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  p4<- ggplot(sym.mhw.all, aes(x=rate.C.day.1, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title = "", x=expression(bold(paste("Heating Rate (",degree,"C day"^-1*")"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  grid.arrange(p1,p2,p3,p4,nrow = 2)
  
  ############## MHW-FILTERED CHL A PER CM² #########################   
  ### y = Chlorophyll a (ug cm²)
  y <- chlcm.mhw.all$resp.var # chl.cm$resp.var, chl.cell$resp.var, fvfm.all$resp.var
  g <- chlcm.mhw.all$study.id # set the random effect variable, study.id
  
  ### x = Temp
  x <- chlcm.mhw.all$max.temp.C # set the x

  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cm.mhw.all_v_temp_Rsq.csv")
  write.csv(AovTbl, file= "chl.cm.mhw.all_v_temp_Aov.csv", na = " ")
  
  ### x = Duration
  x <- chlcm.mhw.all$duration.days
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cm.mhw.all_v_duration_Rsq.csv")
  write.csv(AovTbl, file= "chl.cm.mhw.all_v_duration_Aov.csv", na = " ")
  
  ### x = DHD
  x <-chlcm.mhw.all$dhd.C.days
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cm.mhw.all_v_dhd_Rsq.csv")
  write.csv(AovTbl, file= "chl.cm.mhw.all_v_dhd_Aov.csv", na = " ")
  
  ### x = Heating Rate
  x <- chlcm.mhw.all$rate.C.day.1 
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cm.mhw.all_v_rate_Rsq.csv")
  write.csv(AovTbl, file= "chl.cm.mhw.all_v_rate_Aov.csv", na = " ")
  
  # 4 panel figure of mhw-filtered Chla ug cm^-2 vs. Heating Metrics
  p1<- ggplot(chlcm.mhw.all, aes(x=max.temp.C, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    #theme(axis.title=black.bold.text, axis.text=small.text, legend.position="right", legend.title=element_blank(), legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("Maximum Temperature (",degree,"C)"))), y = "") #y=expression(bold(paste('Chlorophyll a ('~mu *'g Chl cm'^-2*") "))) ) 
  
  p2<- ggplot(chlcm.mhw.all, aes(x=duration.days, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("Duration (days)"))), y="") #y=expression(bold(paste('Chlorophyll a ('~mu *'g Chl cm'^-2*") ")))
  
  p3<- ggplot(chlcm.mhw.all, aes(x=dhd.C.days, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("Degree Heating Days (",degree,"C days)"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  p4<- ggplot(chlcm.mhw.all, aes(x=rate.C.day.1, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title = "", x=expression(bold(paste("Heating Rate (",degree,"C day"^-1*")"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  grid.arrange(p1,p2,p3,p4,nrow = 2)
  
  
  ############### MHW-FILTERED CHL A PER CELL #########################  
  ### y = Chlorophyll a (pg cell-1)
  y <- chlcell.mhw.all$resp.var
  g <- chlcell.mhw.all$study.id # set the random effect variable, study.id
  
  ### x = Temp
  x <- chlcell.mhw.all$max.temp.C
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cell.mhw.all_v_temp_Rsq.csv")
  write.csv(AovTbl, file= "chl.cell.mhw.all_v_temp_Aov.csv", na = " ")

  ### x = Duration
  x <- chlcell.mhw.all$duration.days
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cell.mhw.all_v_duration_Rsq.csv")
  write.csv(AovTbl, file= "chl.cell.mhw.all_v_duration_Aov.csv", na = " ")
  
  ### x = DHD
  x <- chlcell.mhw.all$dhd.C.days 
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cell.mhw.all_v_dhd_Rsq.csv")
  write.csv(AovTbl, file= "chl.cell.mhw.all_v_dhd_Aov.csv", na = " ")
  # *** Is Singular fit, likely because one study has several entries very close to each other. Could simply do a polyfit to see? ***
  ## DHD of high values not expected to be where we see a hormetic bump, actually...
  
  ### x = Heating Rate
  x <- chlcell.mhw.all$rate.C.day.1
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cell.mhw.all_v_rate_Rsq.csv")
  write.csv(AovTbl, file= "chl.cell.mhw.all_v_rate_Aov.csv", na = " ")
  
  # 4 panel figure of  mhw-filtered Chla pg. cell-1 vs. Heating Metrics
  p1<- ggplot(chlcell.mhw.all, aes(x=max.temp.C, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    #theme(axis.title=black.bold.text, axis.text=small.text, legend.position="right", legend.title=element_blank(), legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("Maximum Temperature (",degree,"C)"))), y="") #expression(bold(paste('Chlorophyll a (pg Chl cm'^-2*") "))) ) 
  
  p2<- ggplot(chlcell.mhw.all, aes(x=duration.days, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("Duration (days)"))), y="") #y=expression(bold(paste('Chlorophyll a ('~mu *'g Chl cm'^-2*") ")))
  
  p3<- ggplot(chlcell.mhw.all, aes(x=dhd.C.days, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("Degree Heating Days (",degree,"C days)"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  p4<- ggplot(chlcell.mhw.all, aes(x=rate.C.day.1, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title = "", x=expression(bold(paste("Heating Rate (",degree,"C day"^-1*")"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  grid.arrange(p1,p2,p3,p4,nrow = 2)
  
  
  ################# MHW-FILTERED FV/FM #########################  
  ### y = Fv/Fm (all samples included)
  y <- fvfm.mhw.all$resp.var
  g <- fvfm.mhw.all$study.id
  
  ### x = Temp
  x <- fvfm.mhw.all$max.temp.C
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "fvfm.mhw.all_v_temp_Rsq.csv")
  write.csv(AovTbl, file= "fvfm.mhw.all_v_temp_Aov.csv", na = " ")
  
  ### x = Duration
  x <- fvfm.mhw.all$duration.days
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "fvfm.mhw.all_v_duration_Rsq.csv")
  write.csv(AovTbl, file= "fvfm.mhw.all_v_duration_Aov.csv", na = " ")
  
  ### x = DHD
  x <- fvfm.mhw.all$dhd.C.days
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "fvfm.mhw.all_v_dhd_Rsq.csv")
  write.csv(AovTbl, file= "fvfm.mhw.all_v_dhd_Aov.csv", na = " ")
  
  ### x = Heating Rate
  x <- fvfm.mhw.all$rate.C.day.1

  par(mfrow=c(1,1));plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "fvfm.mhw.all_v_rate_Rsq.csv")
  write.csv(AovTbl, file= "fvfm.mhw.all_v_rate_Aov.csv", na = " ")
  
  # 4 panel figure of mhw-filtered FvFm vs. Heating Metrics
  p1<- ggplot(fvfm.mhw.all, aes(x=max.temp.C, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    #theme(axis.title=black.bold.text, axis.text=small.text, legend.position="right", legend.title=element_blank(), legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("Maximum Temperature (",degree,"C)"))), y="") #expression(bolditalic(paste('F'[V]*'/F'[M]*"" ))) )
  
  p2<- ggplot(fvfm.mhw.all, aes(x=duration.days, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("Duration (days)"))), y="") #y=expression(bold(paste('Chlorophyll a ('~mu *'g Chl cm'^-2*") ")))
  
  p3<- ggplot(fvfm.mhw.all, aes(x=dhd.C.days, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("Degree Heating Days (",degree,"C days)"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  p4<- ggplot(fvfm.mhw.all, aes(x=rate.C.day.1, y=resp.var, colour = study.id)) +
    stat_smooth(method="lm", se=TRUE, fill=NA, formula=y~poly(x,2,raw=TRUE), colour="black" ) +
    geom_point(size=5) +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title = "", x=expression(bold(paste("Heating Rate (",degree,"C day"^-1*")"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  grid.arrange(p1,p2,p3,p4,nrow = 2)
  
  
  