## Martell & Donner 2024: Headwinds to Understanding Stress Response Physiology: A Systematic Review Reveals Mismatch between Real and Simulated Marine Heatwaves

## USAGE: Effect Size Calculations: This script takes symbiont density and fv/fm
# observations filtered database with outliers removed, and calculates the effect size 
# of all paired observations with a heat and corresponding control for each species in
# each study


#### Prepare filtered database for Effect Size Calculations ####
rm(list = ls())
setwd("~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/_5_EScalcs/")
db_noOut <- read.csv("~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/_4_filteredDatabase_outliersRemoved.csv", header = TRUE, stringsAsFactors = TRUE)

library(plyr)
library(dplyr)
library(psych)
library(data.table)
library(tidyr)

## gather and calculate the relevant stuff for resp var effect size calcs
    # symbiont density
df.sym <- ddply(db_noOut, c("study_id","species_name","heat.control"), summarise,  ### Symbiodiniaceae density (cells cm^-2)
                n = length(sym),
                mean = mean(sym),
                sd = sd(sym),
                maxtemp = mean(max.temp.C),
                duration = mean(duration.days),
                dhd = mean(dhd.C.days),
                rate = mean(rate.C.day.1)
)
df.sym
df.sym <- df.sym %>%
  filter(!is.na(mean)) # filter out NAs with no mean sym density
headTail(df.sym) # check it out
dim(df.sym)
fwrite(df.sym, "sym_means_for_ES_calcs.csv")
df.sym.wide <- df.sym %>% # put in wide format for ES calculation
  pivot_wider(id_cols = c(study_id, species_name), names_from = heat.control, values_from = c(mean, sd, n, maxtemp, duration, dhd, rate))
colnames(df.sym.wide) <- c("study","species", "mean.e", "mean.c","sd.e", "sd.c","n.e","n.c","maxtemp.e","maxtemp.c","duration.e","duration.c","dhd.e","dhd.c","rate.e","rate.c") # rename cols for consistency
df.sym.wide <- df.sym.wide[,c(1:9,11,13,15)]; #View(df.sym.wide) # remove control heating var cols, as these were not calculated, see n3_doseResponse_mergeCount.R usage statement
dim(df.sym.wide) # gives total n of effect sizes for resp vars
sym.wide = drop_na(df.sym.wide); sym.wide # drops NAs from any row in the DB for fitting (including observations that did not have either temp vars or resp vars)
dim(sym.wide) # gives total n of effect sizes for fits
fwrite(sym.wide,"sym_wide_forMeta.csv")

      # chl a (ug cm^-2)
df.chla.ug <- ddply(db_noOut, c("study_id","species_name","heat.control"), summarise,  ### Chlorophyll a (ug cm^-2)
                    n = length(chla),
                    mean = mean(chla),
                    sd = sd(chla),
                    maxtemp = mean(max.temp.C),
                    duration = mean(duration.days),
                    dhd = mean(dhd.C.days),
                    rate = mean(rate.C.day.1)
)
df.chla.ug <- df.chla.ug%>%
  filter(!is.na(mean))
df.chla.ug
fwrite(df.chla.ug, "chla.ug.cm_means_for_ES_calcs.csv")
df.chla.ug.wide  <- df.chla.ug %>% 
  pivot_wider(id_cols = c(study_id, species_name), names_from = heat.control, values_from = c(mean, sd, n, maxtemp, duration, dhd, rate))
names(df.chla.ug.wide)
colnames(df.chla.ug.wide) <- c("study","species", "mean.c", "mean.e","sd.c", "sd.e","n.c","n.e","maxtemp.c","maxtemp.e","duration.c","duration.e","dhd.c","dhd.e","rate.c","rate.e")
df.chla.ug.wide <- df.chla.ug.wide[,c(1:8,10,12,14,16)]; #View(df.chla.ug.wide) # use for ES calcs
chla.ug.wide = drop_na(df.chla.ug.wide) # leaves nothing for fits
#fwrite(chla.ug.wide,"chla.ug_wide_forMeta.csv")
### cannot fit Chl a (ug cm^-2) as bleaching response variable, too few paired observations (heat -- control)

      # chl a (pg cell^-1)
df.chla.pg <- ddply(db_noOut, c("study_id","species_name","heat.control"), summarise,  ### Chlorophyll a (pg cell^-1)
                    n = length(chlapcell),
                    mean = mean(chlapcell),
                    sd = sd(chlapcell),
                    maxtemp = mean(max.temp.C),
                    duration = mean(duration.days),
                    dhd = mean(dhd.C.days),
                    rate = mean(rate.C.day.1)
)
df.chla.pg <- df.chla.pg%>%
  filter(!is.na(mean))
df.chla.pg
fwrite(df.chla.pg, "chla.pg_means_for_ES_calcs.csv")
df.chla.pg.wide  <- df.chla.pg %>% 
  pivot_wider(id_cols = c(study_id, species_name), names_from = heat.control, values_from = c(mean, sd, n, maxtemp, duration, dhd, rate))
colnames(df.chla.pg.wide) <- c("study","species", "mean.e", "mean.c", "sd.e", "sd.c","n.e","n.c","maxtemp","maxtemp.c","duration","duration.c","dhd","dhd.c","rate","rate.c")
df.chla.pg.wide <- df.chla.pg.wide[,c(1:8,10,12,14,16)]; #View(df.chla.pg.wide) # use for ES calcs
chla.pg.wide = drop_na(df.chla.pg.wide) # ONE paired observation, cannot perform a fit
fwrite(chla.pg.wide,"chla.pg_wide_forMeta.csv")

      # Fv/Fm
df.fvfm <- ddply(db_noOut, c("study_id","species_name","heat.control"), summarise,  ### Fv/Fm
                 n = length(fvfm),
                 mean = mean(fvfm),
                 sd = sd(fvfm),
                 maxtemp = mean(max.temp.C),
                 duration = mean(duration.days),
                 dhd = mean(dhd.C.days),
                 rate = mean(rate.C.day.1)
)
df.fvfm <- df.fvfm%>%
  filter(!is.na(mean))
headTail(df.fvfm)
fwrite(df.fvfm, "fvfm_for_ES_calcs.csv")
df.fvfm.wide  <- df.fvfm %>% 
  pivot_wider(id_cols = c(study_id, species_name), names_from = heat.control, values_from = c(mean, sd, n, maxtemp, duration, dhd, rate))
colnames(df.fvfm.wide) <- c("study","species", "mean.e", "mean.c","sd.e", "sd.c","n.e","n.c","maxtemp.e","maxtemp.c","duration.e","duration.c","dhd.e","dhd.c","rate.e","rate.c") # rename cols for consistency
headTail(df.fvfm.wide)
df.fvfm.wide <- df.fvfm.wide[,c(1:9,11,13,15)]; headTail(df.fvfm.wide) # remove control heating var cols, as these were not calculated, see n3_doseResponse_mergeCount.R usage statement
dim(df.sym.wide) # gives total n of effect sizes for resp vars
fvfm.wide = drop_na(df.fvfm.wide); fvfm.wide # drops NAs from any row in the DB for fitting (including observations that did not have either temp vars or resp vars)
dim(fvfm.wide) # gives total n of effect sizes for fits
fwrite(fvfm.wide, "fvfm_forMeta.csv")

      # chl c
df.chlc.ug <- ddply(db_noOut, c("study_id","species_name","heat.control"), summarise,  ## Chlorophyll c2 (ug cm^-2)
                    n = length(chlc),
                    mean = mean(chlc),
                    sd = sd(chlc),
                    maxtemp = mean(max.temp.C),
                    duration = mean(duration.days),
                    dhd = mean(dhd.C.days),
                    rate = mean(rate.C.day.1)
)
df.chlc.ug <- df.chlc.ug%>%
  filter(!is.na(mean))
headTail(df.chlc.ug) ## 2 studies, 3 obs, only one paired study for ES calc without replication
fwrite(df.chlc.ug, "chlc.ug_for_ES_calcs.csv")
df.chlc.ug.wide  <- df.chlc.ug %>%  
  pivot_wider(id_cols = c(study_id, species_name), names_from = heat.control, values_from = c(mean, sd, n, maxtemp, duration, dhd, rate))
colnames(df.chlc.ug.wide) <- c("study","species", "mean.c", "mean.e","sd.c", "sd.e","n.c","n.e","maxtemp.c","maxtemp.e","duration.c","duration.e","dhd.c","dhd.e","rate.c","rate.e")
head(df.chlc.ug.wide)
df.chlc.ug.wide <- df.chlc.ug.wide[,c(1:8,10,12,14,16)]; #View(df.chla.pg.wide) # use for ES calcs
chlc.ug.wide = drop_na(df.chlc.ug.wide) ### No Data Remaining
#fwrite(chlc.ug.wide, "chlc.ug_forMeta.csv")

      # Total Chl
df.totalchl <- ddply(db_noOut, c("study_id","species_name","heat.control"), summarise,  ### Total Chlorophyll (ug cm^-2)
                     n = length(totalchl),
                     mean = mean(totalchl),
                     sd = sd(totalchl),
                     maxtemp = mean(max.temp.C),
                     duration = mean(duration.days),
                     dhd = mean(dhd.C.days),
                     rate = mean(rate.C.day.1)
)
df.totalchl <- df.totalchl%>%
  filter(!is.na(mean))
headTail(df.totalchl) ## No Data Remaining
#View(df.totalchl) 


      # Total Chl per Cell
df.totalchlpcell <- ddply(db_noOut, c("study_id","species_name","heat.control"), summarise,  ### Total Chlorophyll (pg cell^-1)
                          n = length(totalchlpcell),
                          mean = mean(totalchlpcell),
                          sd = sd(totalchlpcell),
                          maxtemp = mean(max.temp.C),
                          duration = mean(duration.days),
                          dhd = mean(dhd.C.days),
                          rate = mean(rate.C.day.1)
)
df.totalchlpcell <- df.totalchlpcell%>%
  filter(!is.na(mean))
headTail(df.totalchlpcell) ## No Data Remaining
#View(df.totalchlpcell) 

#### Calculate effect sizes and variance for each response variable with temperature data, excluding chlc, totalchl and totalchlpcell because insufficient samples
# Harrer, M., Cuijpers, P., Furukawa, T.A., & Ebert, D.D. (2021). Doing Meta-Analysis with R: A Hands-On Guide. Boca Raton, FL and London: Chapmann & Hall/CRC Press. ISBN 978-0-367-61007-4.

library(esc)
library(tidyverse)

## Calculate effect size (Hedge's 'g') for each paired heat-control response variable in the database

### Symbiont density
ESsym <- pmap_dfr(sym.wide, 
                  function(mean.e, sd.e, n.e, mean.c,
                           sd.c, n.c, study_id,species,...){
                    esc_mean_sd(grp1m = mean.e,
                                grp1sd = sd.e,
                                grp1n = n.e,
                                grp2m = mean.c,
                                grp2sd = sd.c,
                                grp2n = n.c,
                                study = study_id,
                                es.type = "g") %>% 
                      as.data.frame()}) 
ESsym
ESsym$species <- c(sym.wide$species)
glimpse(ESsym)
fwrite(sym.wide, "summaryStats_sym.csv")
fwrite(ESsym, "effectTable_sym.csv")

# create a df that includes the observed outcome, in this case SMDH (yi), the sample variance (vi), and the temp variable (xi = c(max temp, duration, or DHD) )
sym.fig.data <- ESsym[, c(1,10,2,6)] # get the study, the outcome (yi), and the sample variance (vi)
sym.fig.dat <- merge(sym.fig.data, sym.wide)
names(sym.fig.dat)
colnames(sym.fig.dat) <- c("study", "species", "yi", "vi", "heat.mean", "control.mean", "heated.sd", "control.sd", "heated.n", "control.n", "max.temp.C", "duration.d", "dhd.C.d", "rate.C.d.1")
#View(sym.fig.dat)
names(sym.fig.dat)
sym.fig.dat$es.id = as.factor(c(1:12))
str(sym.fig.dat)
sym.fig.dat2 <- sym.fig.dat[sym.fig.dat$rate.C.d.1<=2, ] # dropping observations that used absurd heating rates :)

par(mfrow=c(4,2))
plot(sym.fig.dat$max.temp.C, sym.fig.dat$yi, pch=21, col="black", bg=sym.fig.dat$study, cex=1.5/sqrt(sym.fig.dat$vi),
     las=1, bty="l", xlab=expression(paste("Maximum Temperature (",degree,"C)")), xlim = c(28, 35), ylim = c(-2,2), ylab="Mean Effect Size (SMD)", main="Symbiodiniaceae Density (all)")

plot(sym.fig.dat2$max.temp.C, sym.fig.dat2$yi, pch=21, col="black", bg=sym.fig.dat2$study, cex=1.5/sqrt(sym.fig.dat2$vi),
     las=1, bty="l", xlab=expression(paste("Maximum Temperature (",degree,"C)")), xlim = c(28, 35), ylim = c(-2,2),ylab="", main="Symbiodiniaceae Density (realistic)")

plot(sym.fig.dat$duration.d, sym.fig.dat$yi, pch=21, col="black", bg=sym.fig.dat$study, cex=1.5/sqrt(sym.fig.dat$vi),
     las=1, bty="l", xlab="Heating Duration (d)", xlim = c(0, 12), ylab="Mean Effect Size (SMD)", ylim = c(-2,2))

plot(sym.fig.dat2$duration.d, sym.fig.dat2$yi, pch=21, col="black", bg=sym.fig.dat2$study, cex=1.5/sqrt(sym.fig.dat2$vi),
     las=1, bty="l", xlab="Heating Duration (d)", xlim = c(0, 12), ylab= "", ylim = c(-2,2))

plot(sym.fig.dat$dhd.C.d, sym.fig.dat$yi, pch=21, col="black", bg=sym.fig.dat$study, cex=1.5/sqrt(sym.fig.dat$vi),
     las=1, bty="l", xlab=expression(paste("Degree Heating Days (",degree,"C days)")), xlim = c(0, 250), ylim = c(-2,2), ylab="Mean Effect Size (SMD)")

plot(sym.fig.dat2$dhd.C.d, sym.fig.dat2$yi, pch=21, col="black", bg=sym.fig.dat2$study, cex=1.5/sqrt(sym.fig.dat2$vi),
     las=1, bty="l", xlab=expression(paste("Degree Heating Days (",degree,"C days)")), xlim = c(0, 250), ylim = c(-2,2), ylab="")

plot(sym.fig.dat$rate.C.d.1, sym.fig.dat$yi, pch=21, col="black", bg=sym.fig.dat$study, cex=1.5/sqrt(sym.fig.dat$vi),
     las=1, bty="l", xlab=expression(paste("Heating Rate (",degree,"C day"^-1*")" )), xlim = c(0, 15), ylim = c(-2,2), ylab="Mean Effect Size (SMD)")

plot(sym.fig.dat2$rate.C.d.1, sym.fig.dat2$yi, pch=21, col="black", bg=sym.fig.dat2$study, cex=1.5/sqrt(sym.fig.dat2$vi),
     las=1, bty="l", xlab=expression(paste("Heating Rate (",degree,"C day"^-1*")" )), xlim = c(0, 15), ylab="", ylim = c(-2,2))

### Chl a per cell ###
ESchlapcell <- pmap_dfr(chla.pg.wide, 
                   function(mean.e, sd.e, n.e, mean.c,
                            sd.c, n.c, study_id,...){
                     esc_mean_sd(grp1m = mean.e,
                                 grp1sd = sd.e,
                                 grp1n = n.e,
                                 grp2m = mean.c,
                                 grp2sd = sd.c,
                                 grp2n = n.c,
                                 study = study_id,
                                 es.type = "g") %>% 
                       as.data.frame()}) 
glimpse(ESchlapcell)
fwrite(chla.pg.wide, "summaryStats_chla.pgcell.csv")
fwrite(ESchlapcell, "effectTable_chla.pgcell.csv")

#### FV/FM ####
ESfvfm <- pmap_dfr(fvfm.wide, 
                   function(mean.e, sd.e, n.e, mean.c,
                            sd.c, n.c, study_id,...){
                     esc_mean_sd(grp1m = mean.e,
                                 grp1sd = sd.e,
                                 grp1n = n.e,
                                 grp2m = mean.c,
                                 grp2sd = sd.c,
                                 grp2n = n.c,
                                 study = study_id,
                                 es.type = "g") %>% 
                       as.data.frame()}) 
glimpse(ESfvfm)
fwrite(fvfm.wide, "summaryStats_fvfm.csv")
fwrite(ESfvfm, "effectTable_fvfm.csv")

# create a df that include the observed outcome, in this case SMDH (yi), the sample variance (vi), and the temp variable (xi = c(max temp, duration, or DHD) )
fvfm.fig.data <- ESfvfm[, c(1, 2, 6)] # get the study, the outcome (yi), and the sample variance (vi)
fvfm.fig.dat <- merge(fvfm.fig.data, fvfm.wide, by = c("study"), all.x = TRUE, all.y = FALSE)
names(fvfm.fig.dat)
colnames(fvfm.fig.dat) <- c("study", "yi", "vi", "species", "heat.mean", "control.mean", "heated.sd", "control.sd", "heated.n", "control.n", "max.temp.C", "duration.d", "dhd.C.d","rate.C.d.1" )
#View(fvfm.fig.dat)

### Now recalculate sym and fvfm fits without absurd heating rates > 2 C d-1
names(fvfm.fig.dat)
fvfm.fig.dat2 <- fvfm.fig.dat[fvfm.fig.dat$rate.C.d.1<=2, ]

plot(fvfm.fig.dat$max.temp.C, fvfm.fig.dat$yi, pch=21, col="black", bg=fvfm.fig.dat$study, cex=1.5/sqrt(fvfm.fig.dat$vi),
     las=1, bty="l", xlab = expression(paste("Maximum Temperature (",degree,"C)")), xlim = c(28, 35), ylim = c(-2, 2), ylab="Mean Effect Size (SMD)", main=expression(bold(paste("F"[V]*"/F"[M]*" (all)"))) )

plot(fvfm.fig.dat2$max.temp.C, fvfm.fig.dat2$yi, pch=21, col="black", bg=fvfm.fig.dat2$study, cex=1.5/sqrt(fvfm.fig.dat2$vi),
     las=1, bty="l", xlab = expression(paste("Maximum Temperature (",degree,"C)")), xlim = c(28, 35), ylim = c(-2, 2),  ylab="", main=expression(bold(paste("F"[V]*"/F"[M]*" (realistic)"))) )

plot(fvfm.fig.dat$duration.d, fvfm.fig.dat$yi, pch=21, col="black", bg=fvfm.fig.dat$study, cex=1.5/sqrt(fvfm.fig.dat$vi),
     las=1, bty="l", xlab="Heating Duration (d)", xlim = c(0, 50),  ylim = c(-2, 2), ylab="Mean Effect Size (SMD)")

plot(fvfm.fig.dat2$duration.d, fvfm.fig.dat2$yi, pch=21, col="black", bg=fvfm.fig.dat2$study, cex=1.5/sqrt(fvfm.fig.dat2$vi),
     las=1, bty="l", xlab="Heating Duration (d)", xlim = c(0, 50),  ylim = c(-2, 2), ylab="")

plot(fvfm.fig.dat$dhd.C.d, fvfm.fig.dat$yi, pch=21, col="black", bg=fvfm.fig.dat$study, cex=1.5/sqrt(fvfm.fig.dat$vi),
     las=1, bty="l", xlab=expression(paste("Degree Heating Days (",degree,"C days)")), xlim = c(0, 150),  ylim = c(-2, 2), ylab="Mean Effect Size (SMD)")

plot(fvfm.fig.dat2$dhd.C.d, fvfm.fig.dat2$yi, pch=21, col="black", bg=fvfm.fig.dat2$study, cex=1.5/sqrt(fvfm.fig.dat2$vi),
     las=1, bty="l", xlab=expression(paste("Degree Heating Days (",degree,"C days)")), xlim = c(0, 150),  ylim = c(-2, 2), ylab="")

plot(fvfm.fig.dat$rate.C.d.1, fvfm.fig.dat$yi, pch=21, col="black", bg=fvfm.fig.dat$study, cex=1.5/sqrt(sym.fig.dat$vi),
     las=1, bty="l", xlab=expression(paste("Heating Rate (",degree,"C day"^-1*")" )), xlim = c(0, 15),  ylim = c(-2, 2), ylab="Mean Effect Size (SMD)")

plot(fvfm.fig.dat2$rate.C.d.1, fvfm.fig.dat2$yi, pch=21, col="black", bg=fvfm.fig.dat2$study, cex=1.5/sqrt(sym.fig.dat2$vi),
     las=1, bty="l", xlab=expression(paste("Heating Rate (",degree,"C day"^-1*")" )), xlim = c(0, 15),  ylim = c(-2, 2), ylab="")
