## USAGE: This script plots temperature variables from the dose response bleaching 
## database and marine heatwaves, then performs statistical comparison of temperature
## metrics between the database and MHW

# written by Harmony Martell 2023

# Plots the flat violin plots (i.e., distributions) and boxplots of each heat variable
# Plots the boxplots of heating metrics between database and MHW values

###############
rm(list = ls())

# load libraries
library(ggplot2)
library(tidyverse)
library(stats)
library(data.table)
library(hrbrthemes)
library(car)
library(rstatix)
library(dplyr)
library(psych)
library(ggpubr)
library(grid)
library(gridExtra)

# Load in dataset
setwd("~/Desktop/UBC/_Projects/doseResponse_analysis/writing/FinalDocs_ERL/polishedAnalyses/")
  db=read.csv("heatedFilteredDatabase.csv", header=TRUE, stringsAsFactors=TRUE)

# Set some plot text stuff
black.bold.text<-element_text(face="bold",color="black", size=12)
small.text<-element_text(face="bold",color="black", size=10)
heatPalette <- c("#D41159","#1A85FF")

####################################################
###  Comparison of DB heating variables to MHWs  ###
####################################################
## Load in & tidy MHW & databse datasets (remove zeros and NAs, assign source)
## bind them together
## subset 95th %ile
## plot 95th %ile
## perform hypothesis tests

## Maximum Temperature
# loads in hot spot peak data, adds MMM to get temperature, then tidies the dataset (removes zero hot spot values and NAs)
HSpeak <- read.csv("~/Desktop/UBC/_Projects/doseResponse_analysis/Analysis/MHWtemps_fromXinru/HSpeak_HSY20092018.csv", header = TRUE)
HSpeak[HSpeak==0] <- NA # replace 0 values with NA so they won't be added
# load in the MMM data that corresponds to each row (grid cell)
mmm <- read.csv("~/Desktop/UBC/_Projects/doseResponse_analysis/Analysis/MHWtemps_fromXinru/MMM_HSY20092018.csv")
mhwMaxTemp <- HSpeak + mmm$MMM_sel #sum the MMM to the HSpeak to get max temperature during the marine heatwave

mhwMaxTemp <- data.frame(value=unlist(mhwMaxTemp)) # put all the columns in a single row
maxtemp.mhw <- data.frame(group = "marine heatwave", value = mhwMaxTemp) # add a column for the data source: marine heatwave
maxtemp.mhw <- na.omit(maxtemp.mhw) # removes NAs
head(maxtemp.mhw)
summary(maxtemp.mhw) # take a look

# sets the database max temp 
dbMaxTemp <- data.frame(value=db$max.temp.C)
rownames(dbMaxTemp) <- db$observation_id
maxtemp.lab <- data.frame(group = "laboratory", value = dbMaxTemp)
maxtemp.lab <- na.omit(maxtemp.lab)
head(maxtemp.lab)
summary(maxtemp.lab)

# Combine the max temp data from both marine heatwaves and database
maxtemp <- rbind(maxtemp.mhw, maxtemp.lab)
library(psych)
headTail(maxtemp)
summary(maxtemp)
sumStats_maxtemp <- describeBy(maxtemp$value, maxtemp$group)
library(data.table)
fwrite(sumStats_maxtemp, "mhw_v_lab_summaryStats_maxTemp.csv")

# Find range and 95th %ile for database and MHW observations
quantile(maxtemp.mhw$value, c(0.05,0.95)) # find quantiles for MHW
maxtemp.mhw.95 <- maxtemp.mhw[maxtemp.mhw$value > 28.42 & maxtemp.mhw$value < 31.57, ]
summary(maxtemp.mhw.95)
summary(maxtemp.mhw)

quantile(maxtemp.lab$value, c(0.05,0.95)) # find quantiles for lab observations
maxtemp.lab.95 <- maxtemp.lab[maxtemp.lab$value > 27.50 & maxtemp.lab$value < 39.85, ]
summary(maxtemp.lab.95)
summary(maxtemp.lab)

maxtemp.95 <- rbind(maxtemp.mhw.95, maxtemp.lab.95)
headTail(maxtemp.95)
summary(maxtemp.95)

# plot 95th %ile data
plotmaxtemp95 <- ggplot(data = maxtemp.95, aes(x = group, y = value), fill=group) +
  geom_boxplot(fill=heatPalette) +
  theme_bw() +
  labs(x="", y= expression(bold(paste("Maximum Temperature (",degree,"C)"))) ) +
  theme(axis.title=black.bold.text, axis.text = small.text, legend.title = black.bold.text)
#plotmaxtemp95

# perform the stats
shap <- shapiro_test(maxtemp.lab,value); shap
qqnorm(maxtemp.lab$value)
# SW statistic = 0.952
# p = 4.28e-10 (significant result, not drawn from normal distribution)

# kruskal wallis test
maxtemp_kwresult <- kruskal.test(value ~ group, data = maxtemp); maxtemp_kwresult
kwtest <- c(maxtemp_kwresult$method, maxtemp_kwresult$statistic, maxtemp_kwresult$parameter, maxtemp_kwresult$p.value); kwtest
write.csv(kwtest, file= "kw.test_max.temp_MHWcomparison.csv", na = " ")
maxtemp%>%
  group_by(group)%>% 
  summarise(GroupMed=median(value))

### Duration
duration <- read.csv("~/Desktop/UBC/_Projects/doseResponse_analysis/Analysis/MHWtemps_fromXinru/Dc_HSY20092018.csv", header = TRUE)
duration[duration==0] <- NA # replace 0 values with NA so they won't be added

mhwDuration <- data.frame(value=unlist(duration)) # put all the columns in a single row
duration.mhw <- data.frame(group = "marine heatwave", value = mhwDuration) # add a column for the data source: marine heatwave
duration.mhw <- na.omit(duration.mhw) # removes NAs
head(duration.mhw)
summary(duration.mhw) # take a look

# sets the database max temp 
dbDuration <- data.frame(value=db$duration.days)
dbDuration[dbDuration==0] <- NA
rownames(dbDuration) <- db$observation_id
duration.lab <- data.frame(group = "laboratory", value = dbDuration)
duration.lab <- na.omit(duration.lab)
head(duration.lab)
summary(duration.lab)

# Combine the max temp data from both marine heatwaves and database
duration <- rbind(duration.mhw, duration.lab)
headTail(duration)
summary(duration)
sumStats_duration <- describeBy(duration$value, duration$group)
fwrite(sumStats_duration, "mhw_v_lab_summaryStats_duration.csv")

# Find 95th %ile for database and MHW observations
quantile(duration.mhw$value, c(0.05,0.95)) # find quantiles for MHW
duration.mhw.95 <- duration.mhw[duration.mhw$value > 6 & duration.mhw$value < 101, ]
summary(duration.mhw.95)
summary(duration.mhw)

quantile(duration.lab$value, c(0.05,0.95)) # find quantiles for lab observations
duration.lab.95 <- duration.lab[duration.lab$value > 0.48 & duration.lab$value < 298.39, ]
summary(duration.lab.95)
summary(duration.lab)

duration.95 <- rbind(duration.mhw.95, duration.lab.95)
headTail(duration.95)
summary(duration.95)

# plot 95th %ile data
plotduration95 <- ggplot(data = duration.95, aes(x = group, y = value), fill=group) +
  geom_boxplot(fill=heatPalette) +
  theme_bw() +
  labs(x="", y= expression(bold(paste("Duration (days)"))) ) +
  theme(axis.title=black.bold.text, axis.text = small.text, legend.title = black.bold.text)

# perform the stats
shap <- shapiro_test(duration.lab,value); shap
qqnorm(duration.lab$value)
# SW statistic = 0.697
# p = 3.05e-26 (Significant, not from a normal distribution)

# kruskal wallis test
duration_kwresult <- kruskal.test(value ~ group, data = duration); duration_kwresult
kwtest <- c(duration_kwresult$method, duration_kwresult$statistic, duration_kwresult$parameter, duration_kwresult$p.value); kwtest
write.csv(kwtest, file= "kw.test_duration_MHWcomparison.csv", na = " ")
duration%>%
  group_by(group)%>% 
  summarise(GroupMed=median(value))

### Degree Heating Days
# Consider filtering by DHD < 28 (Bleaching Warning) and < 56 (CRW Bleaching Alert)
dhds <- read.csv("~/Desktop/UBC/_Projects/doseResponse_analysis/Analysis/MHWtemps_fromXinru/DHD_hsy20092018.csv", header = TRUE)
dhds[dhds==0] <- NA # replace 0 values with NA so they won't be added

mhwDhds <- data.frame(value=unlist(dhds)) # put all the columns in a single row
dhds.mhw <- data.frame(group = "marine heatwave", value = mhwDhds) # add a column for the data source: marine heatwave
dhds.mhw <- na.omit(dhds.mhw) # removes NAs
head(dhds.mhw)
summary(dhds.mhw) # take a look

# sets the database
dbDhds <- data.frame(value=db$dhd.C.days)
dbDhds[dbDhds==0] <- NA
rownames(dbDhds) <- db$observation_id
dhds.lab <- data.frame(group = "laboratory", value = dbDhds)
dhds.lab <- na.omit(dhds.lab)
head(dhds.lab)
summary(dhds.lab)

# Combine the max temp data from both marine heatwaves and database
dhd <- rbind(dhds.mhw, dhds.lab)
headTail(dhd)
summary(dhd)
sumStats_dhd <-describeBy(dhd$value, dhd$group)
fwrite(sumStats_dhd, "mhw_v_lab_summaryStats_dhd.csv")

# Find 95th %ile for database and MHW observations
quantile(dhds.mhw$value, c(0.05,0.95)) # find quantiles for MHW
dhds.mhw.95 <- dhds.mhw[dhds.mhw$value > 2.10 & dhds.mhw$value < 65.55, ]
summary(dhds.mhw.95)
summary(dhds.mhw)

quantile(dhds.lab$value, c(0.05,0.95)) # find quantiles for lab observations
dhds.lab.95 <- dhds.lab[dhds.lab$value > 5.42 & dhds.lab$value < 1026.52, ]
summary(dhds.lab.95)
summary(dhds.lab)

dhd.95 <- rbind(dhds.mhw.95, dhds.lab.95)
headTail(dhd.95)
summary(dhd.95)

# plot 95th %ile data
plotdhd95 <- ggplot(data = dhd.95, aes(x = group, y = value), fill=group) +
  geom_boxplot(fill=heatPalette) +
  theme_bw() +
  labs(x="", y= expression(bold(paste("Degree Heating Days (",degree,"C days)"))) ) +
  theme(axis.title=black.bold.text, axis.text = small.text, legend.title = black.bold.text)
#plotdhd95

# perform the stats
shap <- shapiro_test(dhds.lab,value); shap
qqnorm(dhds.lab$value)
# SW statistic = 0.671
# p = 3.59e-27 (significant, not normal)

# kruskal wallis test
dhd_kwresult <- kruskal.test(value ~ group, data = dhd); dhd_kwresult
kwtest <- c(dhd_kwresult$method, dhd_kwresult$statistic, dhd_kwresult$parameter, dhd_kwresult$p.value); kwtest
write.csv(kwtest, file= "kw.test_dhd_MHWcomparison.csv", na = " ")

dhd%>%
  group_by(group)%>% 
  summarise(GroupMed=median(value))

### Heating Rate
rate <- read.csv("~/Desktop/UBC/_Projects/doseResponse_analysis/Analysis/MHWtemps_fromXinru/HRc_HSY20092018.csv", header = TRUE)
rate[rate==0] <- NA # replace 0 values with NA so they won't be added

mhwRate <- data.frame(value=unlist(rate)) # put all the columns in a single row
head(mhwRate)
rate.mhw <- data.frame(group = "marine heatwave", value = mhwRate) # add a column for the data source: marine heatwave
rate.mhw <- na.omit(rate.mhw) # removes NAs
head(rate.mhw)
summary(rate.mhw) # take a look

# sets the database max temp 
dbRate <- data.frame(value=db$rate.C.day.1)
dbRate[dbRate==0] <- NA
rownames(dbRate) <- db$observation_id
rate.lab <- data.frame(group = "laboratory", value = dbRate)
rate.lab <- na.omit(rate.lab)
head(rate.lab)
summary(rate.lab)

# Combine the max temp data from both marine heatwaves and database
rate <- rbind(rate.mhw, rate.lab)
headTail(rate)
summary(rate)
sumStats_rate <- describeBy(rate$value, rate$group)
fwrite(sumStats_rate, "mhw_v_lab_summaryStats_rate.csv")

# Find 95th %ile for database and MHW observations
quantile(rate.mhw$value, c(0.05,0.95)) # find quantiles for MHW
rate.mhw.95 <- rate.mhw[rate.mhw$value > 0.021 & rate.mhw$value < 0.238, ]
summary(rate.mhw.95)
summary(rate.mhw)

quantile(rate.lab$value, c(0.05,0.95)) # find quantiles for lab observations
rate.lab.95 <- rate.lab[rate.lab$value > 0.04 & rate.lab$value < 16.86, ]
summary(rate.lab.95)
summary(rate.lab)

rate.95 <- rbind(rate.mhw.95, rate.lab.95)
headTail(rate.95)
summary(rate.95)

# plot 95th %ile data
plotrate95 <- ggplot(data = rate.95, aes(x = group, y = value), fill=group) +
  geom_boxplot(fill=heatPalette) +
  theme_bw() +
  labs(x="", y= expression(bold(paste("Heating Rate (",degree,"C day"^-1*")" ) ) ) ) +
  theme(axis.title=black.bold.text, axis.text = small.text, legend.title = black.bold.text)

# perform the stats
shap <- shapiro_test(rate.lab,value); shap
qqnorm(rate.lab$value)
# SW statistic = 0.358
# p = 3.63e-35 (Significant, not normal)

  # kw test
rate_kwresult <- kruskal.test(value ~ group, data = rate); rate_kwresult
kwtest <- c(rate_kwresult$method, rate_kwresult$statistic, rate_kwresult$parameter, rate_kwresult$p.value); kwtest
write.csv(kwtest, file= "kw.test_rate_MHWcomparison.csv", na = " ")

rate%>%
  group_by(group)%>% 
  summarise(GroupMed=median(value))

## Plot all boxplots on a single plot (95th %iles only)
grid.arrange(plotmaxtemp95,plotduration95,plotdhd95,plotrate95,ncol = 2, nrow = 2)

