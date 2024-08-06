## Martell & Donner 2024: Headwinds to Understanding Stress Response Physiology: A Systematic Review Reveals Mismatch between Real and Simulated Marine Heatwaves

## USAGE: Database Filtering Script: This script takes the complete coral bleaching database

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
#library(remotes)
#install_version("Matrix", "1.6")
#oo <- options(repos = "https://cran.r-project.org/")
#utils::install.packages("Matrix")
#utils::install.packages("lme4")
#options(oo)
#install.packages("Matrix")
#install.packages("lme4")
#utils::install.packages("lme4", type = "source")

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
setwd("~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/_6_Fitting/")
  sym = read.csv("symbiodiniaceae_heated.csv", header=TRUE, stringsAsFactors=TRUE)
  chl.cm = read.csv("chlorophylla.ug.cm2_heated.csv", header=TRUE, stringsAsFactors=TRUE)
  chl.cell = read.csv("chlorophylla.pg.cell1_heated.csv", header=TRUE, stringsAsFactors=TRUE)
  fvfm.all = read.csv("fvfm_heated.csv", header=TRUE, stringsAsFactors=TRUE)
  fvfm.sub = read.csv("fvfm.sub_heated.csv", header=TRUE, stringsAsFactors=TRUE)

  summary(fvfm.all$time.of.day)
  summary(fvfm.all$dark.adapt.time)
  par(mfrow=c(1,1))
  boxplot(fvfm.all$dark.adapt.time)
  unique(fvfm.all$name)
  
## Tidy each response variable df
  #keep only "study.id","resp.var","temp","duration","dhd","rate"
  sym = sym[,c(1:2,4:7,9)] 
  head(sym)
  summary(sym)
  
  chl.cm = chl.cm[,c(1:2,4:7,9)]
  head(chl.cm)
  summary(chl.cm)
  
  chl.cell = chl.cell[,c(1:2,4:7,9)]
  head(chl.cell)
  summary(chl.cell)
  
  fvfm.all = fvfm.all[,c(1:2,4:7,9)]
  head(fvfm.all)
  summary(fvfm.all)
  
  fvfm.sub = fvfm.sub[,c(1:2,4:7,9)]
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
  
  # create a dataset for plotting to separate the relevant from the not
  sym.new <- sym
  sym.new$relev <- "NA"
  relev.list <- data.frame(unique(sym.mhw.all$obs.id)); relev.list # get the list of Obs that are ecologically relevant
  sym.new$relev <- as.factor(ifelse(sym.new$obs %in% sym.mhw.all$obs.id, "realistic", "unrealistic"))
  summary(sym.new)
  
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
  
  # create a dataset for plotting to separate the relevant from the not
  chlcm.new <- chl.cm
  chlcm.new$relev <- "NA"
  relev.list <- data.frame(unique(chlcm.mhw.all$obs.id)); relev.list # get the list of Obs that are ecologically relevant
  chlcm.new$relev <- as.factor(ifelse(chlcm.new$obs %in% chlcm.mhw.all$obs.id, "realistic", "unrealistic"))
  summary(chlcm.new)
  
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
 
  # create a dataset for plotting to separate the relevant from the not
  chlcell.new <- chl.cell
  chlcell.new$relev <- "NA"
  relev.list <- data.frame(unique(chlcell.mhw.all$obs.id)); relev.list # get the list of Obs that are ecologically relevant
  chlcell.new$relev <- as.factor(ifelse(chlcell.new$obs %in% chlcell.mhw.all$obs.id, "realistic", "unrealistic"))
  summary(chlcell.new)
  
  
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
 
  # create a dataset for plotting to separate the relevant from the not
  fvfm.new <- fvfm.all
  fvfm.new$relev <- "NA"
  relev.list <- data.frame(unique(fvfm.mhw.all$obs.id)); relev.list # get the list of Obs that are ecologically relevant
  fvfm.new$relev <- as.factor(ifelse(fvfm.new$obs %in% fvfm.mhw.all$obs.id, "realistic", "unrealistic"))
  summary(fvfm.new)
  
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
  
  # create a dataset for plotting to separate the relevant from the not
  fvfmsub.new <- fvfm.sub
  fvfmsub.new$relev <- "NA"
  relev.list <- data.frame(unique(fvfmsub.mhw.all$obs.id)); relev.list # get the list of Obs that are ecologically relevant
  fvfmsub.new$relev <- as.factor(ifelse(fvfmsub.new$obs %in% fvfmsub.mhw.all$obs.id, "realistic", "unrealistic"))
  summary(fvfmsub.new)
  
  
  
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
  sym.temp.model <-model
  
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
  sym.duration.model <-model
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "sym_v_duration_Rsq.csv")
  write.csv(AovTbl, file= "sym_v_duration_Aov.csv", na = " ")
  
  ### x = DHD
  x <- log(sym$dhd.C.days)
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  Rsq <-rsquared.glmm(model); Rsq
  print(summary(model))
  sym.dhd.model <-model
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
  sym.rate.model <-model
  plot(model)
  Rsq<- rsquared.glmm(model); Rsq
  print(summary(model))
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "sym_v_rate_Rsq.csv")
  write.csv(AovTbl, file= "sym_v_rate_Aov.csv", na = " ")
  
# set plot params
  black.bold.text<-element_text(face="bold",color="black", size=18)
  small.text<-element_text(face="bold",color="black", size=14)
  legend.text<-element_text(face="bold",color="black", size=14)
  
  # 4 panel figure of Sym vs. Heating Metrics
  p1<- ggplot() +
    geom_point(data= sym.new, size=5, aes(x=max.temp.C, y=resp.var, color = study.id, shape = relev)) +
    scale_shape_manual(values = c(19,21)) + 
    stat_smooth(data= sym.new, aes(x=max.temp.C, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "dotdash") +
    stat_smooth(data= sym.mhw.all, aes(x=max.temp.C, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "solid") +
    theme_bw() +
    xlim(22,37) +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=legend.text, legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(colour = "Study ID",title = "", x=expression(bold(paste("Maximum Temperature (",degree,"C)"))), y= "") # expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) +
    
  p2<- ggplot() +
    geom_point(data= sym.new, size=5, aes(x=duration.days, y=resp.var, color = study.id, shape = relev)) +
    scale_shape_manual(values = c(19,21)) + 
    stat_smooth(data= sym.new, aes(x=duration.days, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "dotdash") +
    stat_smooth(data= sym.mhw.all, aes(x=duration.days, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "solid") +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=legend.text, legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(colour = "Study ID",title = "", shape = "Relevance", x=expression(bold(paste("Duration (days)"))), y= "") # expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) +
    #theme(axis.title=black.bold.text, axis.text=small.text, legend.position="right", legend.title= legend.text, legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    #labs(title = "", colour = "Study ID", shape = "Relevance", x=expression(bold(paste("Duration (days)"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
    #labs(colour = "Study ID",title = "", x=expression(bold(paste("Maximum Temperature (",degree,"C)"))), y= expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")"))))  # expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  p3<- ggplot() +
    geom_point(data= sym.new, size=5, aes(x=log(dhd.C.days), y=resp.var, color = study.id, shape = relev)) +
    scale_shape_manual(values = c(19,21)) + 
    stat_smooth(data= sym.new, aes(x=log(dhd.C.days), y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "dotdash") +
    stat_smooth(data= sym.mhw.all, aes(x=log(dhd.C.days), y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "solid") +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("log Degree Heating Days (",degree,"C days)"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 

  p4<- ggplot() +
    geom_point(data= sym.new, size=5, aes(x=rate.C.day.1, y=resp.var, color = study.id, shape = relev)) +
    scale_shape_manual(values = c(19,21)) + 
    stat_smooth(data= sym.new, aes(x=rate.C.day.1, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "dotdash") +
    stat_smooth(data= sym.mhw.all, aes(x=rate.C.day.1, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "solid") +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title = "", x=expression(bold(paste("Heating Rate (",degree,"C day"^-1*")"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  S.plot <- grid.arrange(p1,p2,p3,p4,nrow = 2)

  
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
  chlcm.temp.model <-model
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
  chlcm.time.model <- model
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cm_v_duration_Rsq.csv")
  write.csv(AovTbl, file= "chl.cm_v_duration_Aov.csv", na = " ")

  ### x = DHD
  x <-log(chl.cm$dhd.C.days)

  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3)

  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  chlcm.dhd.model <- model
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
  chlcm.rate.model <- model
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cm_v_rate_Rsq.csv")
  write.csv(AovTbl, file= "chl.cm_v_rate_Aov.csv", na = " ")

# 4 panel figure of Chl.cm vs. Heating Metrics
  p1<- ggplot() +
    geom_point(data= chlcm.new, size=5, aes(x=max.temp.C, y=resp.var, color = study.id, shape = relev)) +
    scale_shape_manual(values = c(19,21)) + 
    stat_smooth(data= chlcm.new, aes(x=max.temp.C, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "dotdash") +
    stat_smooth(data= chlcm.mhw.all, aes(x=max.temp.C, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "solid") +
    theme_bw() +
    #xlim(22,37) +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=legend.text, legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(colour = "Study ID",title = "", x=expression(bold(paste("Maximum Temperature (",degree,"C)"))), y= "") 
  
  p2<- ggplot() +
    geom_point(data= chlcm.new, size=5, aes(x=duration.days, y=resp.var, color = study.id, shape = relev)) +
    scale_shape_manual(values = c(19,21)) + 
    stat_smooth(data= chlcm.new, aes(x=duration.days, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "dotdash") +
    stat_smooth(data= chlcm.mhw.all, aes(x=duration.days, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "solid") +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=legend.text, legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(colour = "Study ID",title = "", shape = "Relevance", x=expression(bold(paste("Duration (days)"))), y= "") # expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) +
  #theme(axis.title=black.bold.text, axis.text=small.text, legend.position="right", legend.title= legend.text, legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
  #labs(title = "", colour = "Study ID", shape = "Relevance", x=expression(bold(paste("Duration (days)"))), y="") +
  #labs(colour = "Study ID",title = "", x=expression(bold(paste("Maximum Temperature (",degree,"C)"))), y= expression(bold(paste('Chlorophyll a ('~mu *'g Chl cm'^-2*") "))) ) 
  
  p3<- ggplot() +
    geom_point(data= chlcm.new, size=5, aes(x=log(dhd.C.days), y=resp.var, color = study.id, shape = relev)) +
    scale_shape_manual(values = c(19,21)) + 
    stat_smooth(data= chlcm.new, aes(x=log(dhd.C.days), y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "dotdash") +
    stat_smooth(data= chlcm.mhw.all, aes(x=log(dhd.C.days), y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "solid") +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("log Degree Heating Days (",degree,"C days)"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  p4<- ggplot() +
    geom_point(data= chlcm.new, size=5, aes(x=rate.C.day.1, y=resp.var, color = study.id, shape = relev)) +
    scale_shape_manual(values = c(19,21)) + 
    stat_smooth(data= chlcm.new, aes(x=rate.C.day.1, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "dotdash") +
    stat_smooth(data= chlcm.mhw.all, aes(x=rate.C.day.1, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "solid") +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title = "", x=expression(bold(paste("Heating Rate (",degree,"C day"^-1*")"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  Chl.cm.plot <- grid.arrange(p1,p2,p3,p4,nrow = 2)

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
  chlcell.temp.model <- model
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
  chlcell.time.model <- model
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cell_v_duration_Rsq.csv")
  write.csv(AovTbl, file= "chl.cell_v_duration_Aov.csv", na = " ")

  ### x = DHD
  x <- log(chl.cell$dhd.C.days) # singularity

  par(mfrow=c(1,1))
  plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3)

  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  chlcell.dhd.model <- model
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
  chlcell.rate.model <- model
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "chl.cell_v_rate_Rsq.csv")
  write.csv(AovTbl, file= "chl.cell_v_rate_Aov.csv", na = " ")
  
  # 4 panel figure of Chl p cell vs. Heating Metrics
  p1<- ggplot() +
    geom_point(data= chlcell.new, size=5, aes(x=max.temp.C, y=resp.var, color = study.id, shape = relev)) +
    scale_shape_manual(values = c(19,21)) + 
    stat_smooth(data= chlcell.new, aes(x=max.temp.C, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "dotdash") +
    stat_smooth(data= chlcell.mhw.all, aes(x=max.temp.C, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "solid") +
    theme_bw() +
    #xlim(22,37) +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=legend.text, legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(colour = "Study ID",title = "", x=expression(bold(paste("Maximum Temperature (",degree,"C)"))), y= "") 
  
  p2<- ggplot() +
    geom_point(data= chlcell.new, size=5, aes(x=duration.days, y=resp.var, color = study.id, shape = relev)) +
    scale_shape_manual(values = c(19,21)) + 
    stat_smooth(data= chlcell.new, aes(x=duration.days, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "dotdash") +
    stat_smooth(data= chlcell.mhw.all, aes(x=duration.days, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "solid") +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=legend.text, legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(colour = "Study ID",title = "", shape = "Relevance", x=expression(bold(paste("Duration (days)"))), y= "") # expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) +
  #theme(axis.title=black.bold.text, axis.text=small.text, legend.position="right", legend.title= legend.text, legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
  #labs(title = "", colour = "Study ID", shape = "Relevance", x=expression(bold(paste("Duration (days)"))), y="") +
  #labs(colour = "Study ID",title = "", x=expression(bold(paste("Maximum Temperature (",degree,"C)"))), y= expression(bold(paste('Chlorophyll a ('~mu *'g Chl cm'^-2*") "))) ) 
  
  p3<- ggplot() +
    geom_point(data= chlcell.new, size=5, aes(x=log(dhd.C.days), y=resp.var, color = study.id, shape = relev)) +
    scale_shape_manual(values = c(19,21)) + 
    stat_smooth(data= chlcell.new, aes(x=log(dhd.C.days), y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "dotdash") +
    stat_smooth(data= chlcell.mhw.all, aes(x=log(dhd.C.days), y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "solid") +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("log Degree Heating Days (",degree,"C days)"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  p4<- ggplot() +
    geom_point(data= chlcell.new, size=5, aes(x=rate.C.day.1, y=resp.var, color = study.id, shape = relev)) +
    scale_shape_manual(values = c(19,21)) + 
    stat_smooth(data= chlcell.new, aes(x=rate.C.day.1, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "dotdash") +
    stat_smooth(data= chlcell.mhw.all, aes(x=rate.C.day.1, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "solid") +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title = "", x=expression(bold(paste("Heating Rate (",degree,"C day"^-1*")"))), y="") #y=expression(bold(paste('Symbiodiniaceae Density (cells cm'^-2*' x 10'^6*")")))) 
  
  Chl.cell.plot <- grid.arrange(p1,p2,p3,p4,nrow = 2)

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
  fvfm.temp.model <- model
  plot(model)
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "fvfm.all_v_temp_Rsq.csv")
  write.csv(AovTbl, file= "fvfm.all_v_temp_Aov.csv", na = " ")
  
  ### x = Duration
  x <- log(fvfm.all$duration.days)
  
  par(mfrow=c(1,1));plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  fvfm.duration.model <- model
  plot(model)
  fvfm.duration.model <- model
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "fvfm.all_v_duration_Rsq.csv")
  write.csv(AovTbl, file= "fvfm.all_v_duration_Aov.csv", na = " ")
  
  ### x = DHD
  x <- log(fvfm.all$dhd.C.days)
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  fvfm.dhd.model <- model
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
  fvfm.rate.model <- model
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "fvfm.all_v_rate_Rsq.csv")
  write.csv(AovTbl, file= "fvfm.all_v_rate_Aov.csv", na = " ")

  # 4 panel figure of Fv/Fm vs. Heating Metrics
  p1<- ggplot() +
    geom_point(data= fvfm.new, size=5, aes(x=max.temp.C, y=resp.var, color = study.id, shape = relev)) +
    scale_shape_manual(values = c(19,21)) + 
    stat_smooth(data= fvfm.new, aes(x=max.temp.C, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "dotdash") +
    stat_smooth(data= fvfm.mhw.all, aes(x=max.temp.C, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "solid") +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=legend.text, legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(colour = "Study ID", shape = "Relevance", title = "", x=expression(bold(paste("Maximum Temperature (",degree,"C)"))), y= "") # expression(bolditalic(paste('F'[V]*'/F'[M]*"" ))) ) +
  
  p2<- ggplot() +
    geom_point(data= fvfm.new, size=5, aes(x=log(duration.days), y=resp.var, color = study.id, shape = relev)) +
    scale_shape_manual(values = c(19,21)) + 
    stat_smooth(data= fvfm.new, aes(x=log(duration.days), y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "dotdash") +
    stat_smooth(data= fvfm.mhw.all, aes(x=log(duration.days), y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "solid") +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=legend.text, legend.text=legend.text, strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(colour = "Study ID",title = "", shape = "Relevance", x=expression(bold(paste("log Duration (days)"))), y= "")
  
  p3<- ggplot() +
    geom_point(data= fvfm.new, size=5, aes(x=log(dhd.C.days), y=resp.var, color = study.id, shape = relev)) +
    scale_shape_manual(values = c(19,21)) + 
    stat_smooth(data= fvfm.new, aes(x=log(dhd.C.days), y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "dotdash") +
    stat_smooth(data= fvfm.mhw.all, aes(x=log(dhd.C.days), y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "solid") +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title="", x=expression(bold(paste("log Degree Heating Days (",degree,"C days)"))), y="")
  
  p4<- ggplot() +
    geom_point(data= fvfm.new, size=5, aes(x=rate.C.day.1, y=resp.var, color = study.id, shape = relev)) +
    scale_shape_manual(values = c(19,21)) + 
    stat_smooth(data= fvfm.new, aes(x=rate.C.day.1, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "dotdash") +
    stat_smooth(data= fvfm.mhw.all, aes(x=rate.C.day.1, y=resp.var), method="lm", se=TRUE, fill = NA, formula=y~poly(x,2,raw=TRUE), colour="black", linetype = "solid") +
    theme_bw() +
    theme(axis.title=black.bold.text, axis.text=small.text, legend.position="none", legend.title=element_blank(), strip.text.x=element_text(size=12, colour="black", face="bold"), strip.background=element_rect(colour="black", fill="white")) +
    labs(title = "", x=expression(bold(paste("Heating Rate (",degree,"C day"^-1*")"))), y="") 
  
  F.plot <- grid.arrange(p1,p2,p3,p4,nrow = 2)
  
 
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
  sym.temp.mhw.model <- model
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
  sym.duration.mhw.model <- model
  Rsq<-rsquared.glmm(model); Rsq
  print(summary(model))
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "sym.mhw.all_v_duration_Rsq.csv")
  write.csv(AovTbl, file= "sym.mhw.all_v_duration_Aov.csv", na = " ")
  
  ### x = DHD
  x <- log(sym.mhw.all$dhd.C.days)

  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  # Perform the regression
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  sym.dhd.mhw.model <- model
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
  sym.rate.mhw.model <- model
  Rsq<-rsquared.glmm(model); Rsq
  print(summary(model))
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "sym.mhw.all_v_rate_Rsq.csv")
  write.csv(AovTbl, file= "sym.mhw.all_v_rate_Aov.csv", na = " ")
 
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
  x <-log(chlcm.mhw.all$dhd.C.days)

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
  x <- log(chlcell.mhw.all$dhd.C.days)

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
  fvfm.temp.mhw.model <- model
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "fvfm.mhw.all_v_temp_Rsq.csv")
  write.csv(AovTbl, file= "fvfm.mhw.all_v_temp_Aov.csv", na = " ")
  
  ### x = Duration
  x <- log(fvfm.mhw.all$duration.days)
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  fvfm.duration.mhw.model <- model
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "fvfm.mhw.all_v_duration_Rsq.csv")
  write.csv(AovTbl, file= "fvfm.mhw.all_v_duration_Aov.csv", na = " ")
  
  ### x = DHD
  x <- log(fvfm.mhw.all$dhd.C.days)
  
  par(mfrow=c(1,1)); plot(x,y,col=rgb(0.4,0.4,0.8,0.6),pch=16 ,cex=1.3) 
  
  model <- lmer(y ~ I(x^2) + I(x) + (1|g)); model
  print(summary(model))
  plot(model)
  fvfm.dhd.mhw.model <- model
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
  fvfm.rate.mhw.model <- model
  Rsq<-rsquared.glmm(model); Rsq
  j <-anova(model); j
  AovTbl<-Anova(model, test.statistic = "F"); AovTbl
  fwrite(Rsq,file = "fvfm.mhw.all_v_rate_Rsq.csv")
  write.csv(AovTbl, file= "fvfm.mhw.all_v_rate_Aov.csv", na = " ")

########################################################################
##### Write out the model summaries for +- CI comparisons of Coeffs ####
########################################################################
df1<-summary(sym.temp.model)
df2<-summary(sym.temp.mhw.model)
sym.temp<-cbind(df1$coefficients,df2$coefficients)  
sym.temp <- sym.temp[,c(1,2,4,5)]

df1<-summary(sym.duration.model)
df2<-summary(sym.duration.mhw.model)
sym.time<-cbind(df1$coefficients,df2$coefficients)  
sym.time <- sym.time[,c(1,2,4,5)]

df1<-summary(sym.dhd.model)
df2<-summary(sym.dhd.mhw.model)
sym.dhd<-cbind(df1$coefficients,df2$coefficients)  
sym.dhd <- sym.dhd[,c(1,2,4,5)]

df1<-summary(sym.rate.model)
df2<-summary(sym.rate.mhw.model)
sym.rate<-cbind(df1$coefficients,df2$coefficients)  
sym.rate <- sym.rate[,c(1,2,4,5)]

sym.all.model.coefs <- rbind(sym.temp, sym.time, sym.dhd, sym.rate)
rbind(sym.temp, sym.time, sym.dhd, sym.rate)
rownames(sym.all.model.coefs) <- c("Intercept Temp", "Ix^2 Temp", "Ix Temp","Intercept Time", "Ix^2 Time", "Ix Time","Intercept DHD", "Ix^2 DHD", "Ix DHD","Intercept Rate", "Ix^2 Rate", "Ix Rate")
fwrite(sym.all.model.coefs, "sym.coefs.csv")  

df1<-summary(fvfm.temp.model)
df2<-summary(fvfm.temp.mhw.model)
fvfm.temp<-cbind(df1$coefficients,df2$coefficients)  
fvfm.temp <- fvfm.temp[,c(1,2,4,5)]

df1<-summary(fvfm.duration.model)
df2<-summary(fvfm.duration.mhw.model)
fvfm.time<-cbind(df1$coefficients,df2$coefficients)  
fvfm.time <- fvfm.time[,c(1,2,4,5)]

df1<-summary(fvfm.dhd.model)
df2<-summary(fvfm.dhd.mhw.model)
fvfm.dhd<-cbind(df1$coefficients,df2$coefficients)  
fvfm.dhd <- fvfm.dhd[,c(1,2,4,5)]

df1<-summary(fvfm.rate.model)
df2<-summary(fvfm.rate.mhw.model)
fvfm.rate<-cbind(df1$coefficients,df2$coefficients)  
fvfm.rate <- fvfm.rate[,c(1,2,4,5)]

fvfm.all.model.coefs <- rbind(fvfm.temp, fvfm.time, fvfm.dhd, fvfm.rate)
rbind(fvfm.temp, fvfm.time, fvfm.dhd, fvfm.rate)
rownames(fvfm.all.model.coefs) <- c("Intercept Temp", "Ix^2 Temp", "Ix Temp","Intercept Time", "Ix^2 Time", "Ix Time","Intercept DHD", "Ix^2 DHD", "Ix DHD","Intercept Rate", "Ix^2 Rate", "Ix Rate")
fwrite(fvfm.all.model.coefs, "fvfm.coefs.csv")  
