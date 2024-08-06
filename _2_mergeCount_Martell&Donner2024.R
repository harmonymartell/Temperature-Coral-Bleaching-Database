### Martell & Donner 2024: Headwinds to Understanding Stress Response Physiology: A Systemic Review Reveals Mismatch between Real and Simulated Marine Heatwaves

## USAGE: This script loads in the filtered database & heating parameters,
# counts all observations for each response variable,
# merges heated response variables with heating parameters on observation id,
# filters heated response variable observations & recounts,
# then writes the data to a file for MHW comparison

# removes one outlying study, and writes the data out to files for downstream analyses

# written by Harmony Martell 2024

###############
rm(list = ls())

#load libraries
#library(ggplot2)
#library(tidyverse)
#library(stats)
#library(lme4)
#library(hrbrthemes)
#library(gridExtra)
#library(car)
#library(rstatix)
library(dplyr)
library(psych)

# Load in full database without heating parameters
setwd("~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/")
fulldb <- read.csv("_2_filteredDatabase.csv", header = TRUE, stringsAsFactors = TRUE) #~/Desktop/UBC/_Projects/doseResponse_analysis/writing/FinalDocs_ERL/polishedAnalyses/
## the filteredDatabase.csv observation_id column was manually modified to remove the _ to enable merging
str(fulldb)
summary(fulldb)
unique(fulldb$sym_method)
unique(fulldb$chl_method) # needs to be tidied to fluorometer, spectrophotometric equations and chromatic equations
unique(fulldb$fvfm_name)
# Load in heating parameters
heatparams=read.csv("_2_heatingParams_HeatedObs_toMerge.csv", header = TRUE, stringsAsFactors = TRUE)
### the heating parameters were only calculated for HEATED treatment observations, to save computing time
### complete database is possible with calculation of CONTROL treatment observations, using heatingParameterCalculation.R
names(heatparams)
str(heatparams)
head(names(fulldb),5)
names(fulldb)
head(fulldb$observation_id)
head(heatparams$observation_id)

# Merge the db data by the heatparams dataframe
db_merged <- merge(fulldb, heatparams, by = c("study_id","observation_id"), all.x = TRUE, all.y = FALSE)
#db_merged <- merge(fulldb, heatparams, by = "observation_id")
head(db_merged)
names(db_merged)
dim(db_merged)
summary(db_merged$max.temp.C) # summary shows 546 observations do not have heating parameters (includes both trts)
summary(db_merged)
# Remove the control observations for downstream analyses & write to file
db_heat <- db_merged %>% 
  filter(heat.control == "HEAT") %>%
  filter(!max.temp.C == "NA")
summary(db_heat$max.temp.C)
library(data.table)
fwrite(db_heat, "~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/_3_MHWcomparison/_3_heatedFilteredDatabase.csv")

# Counting all observations of each response variable (includes control and heated observations) before any outlier removal
# sym
dim( db_merged %>% filter(!is.na(sym)) ) #245
# chla
dim( db_merged %>% filter(!is.na(chla)) ) #67
# chlapcell
dim( db_merged %>% filter(!is.na(chlapcell))) # 122
# totalchl
dim( db_merged %>% filter(!is.na(totalchl))) # 14
# chlc
dim( db_merged %>% filter(!is.na(chlc)) ) # 13
# fvfm
dim( db_merged %>% filter(!is.na(fvfm)) ) # 778
fvcount <- db_merged %>% 
  filter(!is.na(fvfm)) %>%
  filter(dark_adapt == "YES") %>% 
  filter(!dark_adapt_time == "NR"); dim(fvcount) # 675

## the database was filtered by _1_dbFiltering_Martell&Donner2024.R and written out as filteredDatabase.csv
## the fvfm_names in unfilteredDatabase.csv were manually tidied to fix typos, resulting in ~16 different name descriptions
## Fv/Fm observations were counted after excluding those that were NaNs, were not dark adapted and reported the dark adapting time

## Remove outlying study for recount 
db_noOut <- db_merged %>%
  filter(!study_id == "S0093") #Remove S0093 from database, because entries are those that were derived by means other than regular air brushing or water piking (using decalcified tissue to cut a 1 cm^2 tissue to homogenize and then derive chl from a similar process using methanol)
summary(db_noOut) # examine database values

# Recount each response variable with S0093 removed
## Change the response variable to get a count! ##
respvar <- db_noOut %>% filter(!is.na(chlapcell)) # change the response variable from sym to what else you want to calculate
#respvar <- fvcount %>% filter(!is.na(fvfm)) # uncomment to use the filtered FvFm database
heat <- respvar %>% filter(heat.control == "HEAT")
cont <- respvar %>% filter(heat.control == "CONTROL")
count(cont)
count(heat)

db_merged_NA <- db_merged %>% filter(is.na(max.temp.C));  # count the NAs in the merged db (obs. without both resp. var + temp)
dim(db_merged_NA) # 546 observations did not have temp data
length(db_merged$max.temp.C) - length(db_merged_NA$max.temp.C) # count the total number of observations 

library(data.table)
fwrite(db_noOut, "~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/_5_EScalcs/_4_filteredDatabase_outliersRemoved.csv")

### Parse each response variable into smaller files for downstream Response Variable & Fitting analysis ###

# Symbiondiniaceae Density
S <- db_noOut[ ,c("study_id", "observation_id", "heat.control", "sym", "max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "sym_method")]; dim(S)
colnames(S) <- c("study.id", "obs.id", "trt", "resp.var","max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "method")
S <- S %>% filter(!is.na(resp.var)) # remove NAs
summary(S)
fwrite(S, "~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/_4_RespVarStats/symbiodiniaceae.csv") # write out data for response var stats

heat <- S %>% filter(trt == "HEAT") %>% filter(!is.na(resp.var)) %>% filter(!is.na(max.temp.C))
summary(heat)
  fwrite(heat, "~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/_6_Fitting/symbiodiniaceae_heated.csv") # write out data for dose fitting

# Chlorophyll a (ug cm-2)
C <- db_noOut[ ,c("study_id", "observation_id", "heat.control", "chla", "max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "chl_method")]; dim(C)
colnames(C) <- c("study.id", "obs.id", "trt", "resp.var","max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "method")
C <- C %>% filter(!is.na(resp.var)) # remove NAs
summary(C)
fwrite(C, "~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/_4_RespVarStats/chlorophylla.ug.cm2.csv") # write out data for response var stats

heat <- C %>% filter(trt == "HEAT") %>% filter(!is.na(resp.var)) %>% filter(!is.na(max.temp.C))
summary(heat)
fwrite(heat, "~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/_6_Fitting/chlorophylla.ug.cm2_heated.csv") # write out data for dose fitting

# Chlorophyll a (pg cell-1)
Cpc <- db_noOut[ ,c("study_id", "observation_id", "heat.control", "chlapcell", "max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "chl_method")]; dim(Cpc)
colnames(Cpc) <- c("study.id", "obs.id", "trt", "resp.var","max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "method")
Cpc <- Cpc %>% filter(!is.na(resp.var)) # remove NAs
summary(Cpc)
fwrite(Cpc, "~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/_4_RespVarStats/chlorophylla.pg.cell1.csv") # write out data for response var stats

heat <- Cpc %>% filter(trt == "HEAT") %>% filter(!is.na(resp.var)) %>% filter(!is.na(max.temp.C))
summary(heat)
fwrite(heat, "~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/_6_Fitting/chlorophylla.pg.cell1_heated.csv") # write out data for dose fitting

# Chlorophyll c2 (ug cm-2)
Cc <- db_noOut[ ,c("study_id", "observation_id", "heat.control", "chlc", "max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "chl_method")]; dim(Cc)
colnames(Cc) <- c("study.id", "obs.id", "trt", "resp.var","max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "method")
Cc <- Cc %>% filter(!is.na(resp.var)) # remove NAs
summary(Cc)
fwrite(Cc, "~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/_4_RespVarStats/chlorophyllc2.ug.cm2.csv") # write out data for response var stats

heat <- Cc %>% filter(trt == "HEAT") %>% filter(!is.na(resp.var)) %>% filter(!is.na(max.temp.C))
summary(heat) ### n = 3 studies, n = 6 obs -- remove for singularity
fwrite(heat, "~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/_6_Fitting/chlorophyllc2.ug.cm2_heated.csv") # write out data for dose fitting

# Total Chlorophyll (ug cm-2)
tC <- db_noOut[ ,c("study_id", "observation_id", "heat.control", "totalchl", "max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "chl_method")]; dim(tC)
colnames(tC) <- c("study.id", "obs.id", "trt", "resp.var","max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "method")
tC <- tC %>% filter(!is.na(resp.var)) # remove NAs
summary(tC)
fwrite(tC, "~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/_4_RespVarStats/totalchlorophyll.csv") # write out data for response var stats

heat <- tC %>% filter(trt == "HEAT") %>% filter(!is.na(resp.var)) %>% filter(!is.na(max.temp.C))
summary(heat) ### n = 1 study, n = 4 obs -- remove for singularity
fwrite(heat, "~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/_6_Fitting/chlorophyllc2.ug.cm2_heated.csv") # write out data for dose fitting

# Fv/Fm
Fv <- db_noOut[ ,c("study_id", "observation_id", "heat.control", "fvfm", "max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "fvfm_name", "fvfm_time","dark_adapt", "dark_adapt_time")]; dim(Fv)
colnames(Fv) <- c("study.id", "obs.id", "trt", "resp.var","max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "name", "time.of.day","dark.adapt","dark.adapt.time")
Fv <- Fv %>% filter(!is.na(resp.var)) %>% # 778
  filter(dark.adapt == "YES") %>% # n = 680
  filter(!dark.adapt.time == "NR") # n = 675
summary(Fv) #
fwrite(Fv, "~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/_4_RespVarStats/fvfm.csv") # write out data for response var stats

heat <- Fv %>% filter(trt == "HEAT") %>% filter(!is.na(resp.var)) %>% filter(!is.na(max.temp.C))
summary(heat)
fwrite(heat, "~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/_6_Fitting/fvfm_heated.csv") # write out data for dose fitting

Match <- S %>% 
  filter(S$study.id %in% Fv$study.id)
ids <- unique(Match$study.id)
Fv_sym <- subset(Fv, sapply(study.id,  \(Fv) any(ids %in% Fv))) # keeping the studies that also reported symbiont densities
heat <- Fv_sym %>% filter(trt == "HEAT") %>% filter(!is.na(resp.var)) %>% filter(!is.na(max.temp.C))
fwrite(heat, "~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/_6_Fitting/fvfm.sub_heated.csv") # write out data for dose fitting
