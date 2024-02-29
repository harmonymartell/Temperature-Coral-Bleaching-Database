## USAGE: This script merges heated response variables with heating parameters
## from the dose response bleaching database then counts and parses each response variable

# written by Harmony Martell 2023

# Loads in filtered database & heating parameters
# Merges on observation id
# Counts heated and control trts of each response variable before removing outliers (resp var n)
# Removes outlying studies, recounts (resp var with no outliers n)
# Parses only heated treatment observations
# Writes each heated response variable, recounts (fitting n)
# Writes heated treatment heating parameters for MHW comparison
#------------------------------------------------------------------#

rm(list = ls())

# Loads in filtered database (no heat parameters yet) & heating parameters
setwd("~/Desktop/UBC/_Projects/doseResponse_analysis/writing/FinalDocs_ERL/polishedAnalyses/")
fulldb <- read.csv("filteredDatabase.csv", header = TRUE, stringsAsFactors = TRUE) # the filteredDatabase.csv observation_id column was manually modified to remove the _ to enable merging
str(fulldb)

heatparams=read.csv("heatingParams_HeatedObs_toMerge.csv", header = TRUE, stringsAsFactors = TRUE) # the heating parameters were only calculated for HEATED treatment observations, to save computing time
    ### complete database is possible with calculation of CONTROL treatment observations, using heatingParameterCalculation.R
names(heatparams)
dim(heatparams)
head(names(fulldb),3)
names(fulldb)
head(fulldb$observation_id)
head(heatparams$observation_id)

# Merge the db by the heatparams on observation id
db_merged <- merge(fulldb, heatparams, by = c("study_id","observation_id"), all.x = TRUE, all.y = FALSE)
head(db_merged)
names(db_merged)
dim(db_merged)
summary(db_merged$max.temp.C) # summary shows 546 observations do not have heating parameters (includes both trts)

# Count heated and control trts of each response variable before removing outliers (resp var n)
# sym
db_merged %>% count(is.na(sym)) # 245
# chla
db_merged %>% count(is.na(chla)) # 67
# chlapcell
db_merged %>% count(is.na(chlapcell)) # 122
# totalchl
db_merged %>% count(is.na(totalchl)) # 14
# chlc
db_merged %>% count(is.na(chlc)) # 13
# fvfm
db_merged %>% count(!is.na(fvfm)) # 778
fvcount <- db_merged %>% 
  filter(!is.na(fvfm)) %>%
  filter(dark_adapt == "YES") %>% 
  filter(!dark_adapt_time == "NR"); count(fvcount) # 675

## the database was filtered by doseResponse_dbFiltering.R and written out as filteredDatabase.csv
## the fvfm_names in filteredDatabase.csv were manually tidied to fix typos, resulting in 11 names
## Fv/Fm observations were counted after excluding those that were NaNs, were not dark adapted and reported the dark adapting time

# Remove outlying study S0093, recount (resp var with no outliers n)
db_noOut <- db_merged %>%
  filter(!study_id == "S0093") #Remove S0093 from database, because entries are those that were derived by means other than regular air brushing or water piking (using decalcified tissue to cut a 1 cm^2 tissue to homogenize and then derive chl from a similar process using methanol)
summary(db_noOut) # examine database values

# Count heated and control trts of each response variable before removing outliers (resp var n)
# sym
db_noOut %>% count(is.na(sym)) # 209
# chla
db_noOut %>% count(is.na(chla)) # 67
# chlapcell
db_noOut %>% count(is.na(chlapcell)) # 86
# totalchl
db_noOut %>% count(is.na(totalchl)) # 14
# chlc
db_noOut %>% count(is.na(chlc)) # 13
# fvfm
db_noOut %>% count(!is.na(fvfm)) # 778
fvcount <- db_noOut %>% 
  filter(!is.na(fvfm)) %>%
  filter(dark_adapt == "YES") %>% 
  filter(!dark_adapt_time == "NR"); count(fvcount) # 675

### Count each response variable by treatment (heat & control) after removing outliers (fitting n)
### Write out all data for Response Variable Statistics & Heated Data for Fitting

  # Symbiondiniaceae Density
S <- db_noOut[ ,c("study_id", "observation_id", "heat.control", "sym", "max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "sym_method")]; dim(S)
colnames(S) <- c("study.id", "obs.id", "trt", "resp.var","max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "method")
S <- S %>% filter(!is.na(resp.var)) # remove NAs
summary(S)
fwrite(S, "symbiodiniaceae.csv") # write out data for response var stats

heat <- S %>% filter(trt == "HEAT") %>% filter(!is.na(resp.var)) %>% filter(!is.na(max.temp.C))
summary(heat)
fwrite(heat, "symbiodiniaceae_heated.csv") # write out data for dose fitting

  # Chlorophyll a (ug cm-2)
C <- db_noOut[ ,c("study_id", "observation_id", "heat.control", "chla", "max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "chl_method")]; dim(C)
colnames(C) <- c("study.id", "obs.id", "trt", "resp.var","max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "method")
C <- C %>% filter(!is.na(resp.var)) # remove NAs
summary(C)
fwrite(C, "chlorophylla.ug.cm2.csv") # write out data for response var stats

heat <- C %>% filter(trt == "HEAT") %>% filter(!is.na(resp.var)) %>% filter(!is.na(max.temp.C))
summary(heat)
fwrite(heat, "chlorophylla.ug.cm2_heated.csv") # write out data for dose fitting

# Chlorophyll a (pg cell-1)
Cpc <- db_noOut[ ,c("study_id", "observation_id", "heat.control", "chlapcell", "max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "chl_method")]; dim(Cpc)
colnames(Cpc) <- c("study.id", "obs.id", "trt", "resp.var","max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "method")
Cpc <- Cpc %>% filter(!is.na(resp.var)) # remove NAs
summary(Cpc)
fwrite(Cpc, "chlorophylla.pg.cell1.csv") # write out data for response var stats

heat <- Cpc %>% filter(trt == "HEAT") %>% filter(!is.na(resp.var)) %>% filter(!is.na(max.temp.C))
summary(heat)
fwrite(heat, "chlorophylla.pg.cell1_heated.csv") # write out data for dose fitting

# Chlorophyll c2 (ug cm-2)
Cc <- db_noOut[ ,c("study_id", "observation_id", "heat.control", "chlc", "max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "chl_method")]; dim(Cc)
colnames(Cc) <- c("study.id", "obs.id", "trt", "resp.var","max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "method")
Cc <- Cc %>% filter(!is.na(resp.var)) # remove NAs
summary(Cc)
fwrite(Cc, "chlorophyllc2.ug.cm2.csv") # write out data for response var stats

heat <- Cc %>% filter(trt == "HEAT") %>% filter(!is.na(resp.var)) %>% filter(!is.na(max.temp.C))
summary(heat) ### n = 3 studies, n = 6 obs -- remove for singularity
fwrite(heat, "chlorophyllc2.ug.cm2_heated.csv") # write out data for dose fitting

# Total Chlorophyll (ug cm-2)
tC <- db_noOut[ ,c("study_id", "observation_id", "heat.control", "totalchl", "max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "chl_method")]; dim(tC)
colnames(tC) <- c("study.id", "obs.id", "trt", "resp.var","max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "method")
tC <- tC %>% filter(!is.na(resp.var)) # remove NAs
summary(tC)
fwrite(tC, "totalchlorophyll.csv") # write out data for response var stats

heat <- tC %>% filter(trt == "HEAT") %>% filter(!is.na(resp.var)) %>% filter(!is.na(max.temp.C))
summary(heat) ### n = 1 study, n = 4 obs -- remove for singularity
fwrite(heat, "chlorophyllc2.ug.cm2_heated.csv") # write out data for dose fitting

# Fv/Fm
Fv <- db_noOut[ ,c("study_id", "observation_id", "heat.control", "fvfm", "max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "fvfm_name", "fvfm_time","dark_adapt", "dark_adapt_time")]; dim(Fv)
colnames(Fv) <- c("study.id", "obs.id", "trt", "resp.var","max.temp.C", "duration.days", "dhd.C.days", "range.C", "rate.C.day.1", "name", "time.of.day","dark.adapt","dark.adapt.time")
Fv <- Fv %>% filter(!is.na(resp.var)) %>% # 778
  filter(dark.adapt == "YES") %>% # n = 680
  filter(!dark.adapt.time == "NR") # n = 675
summary(Fv) #
fwrite(Fv, "fvfm.csv") # write out data for response var stats

heat <- Fv %>% filter(trt == "HEAT") %>% filter(!is.na(resp.var)) %>% filter(!is.na(max.temp.C))
summary(heat)
fwrite(heat, "fvfm_heated.csv") # write out data for dose fitting

Match <- S %>% 
  filter(S$study.id %in% Fv$study.id)
ids <- unique(Match$study.id)
Fv_sym <- subset(Fv, sapply(study.id,  \(Fv) any(ids %in% Fv))) # keeping the studies that also reported symbiont densities
heat <- Fv_sym %>% filter(trt == "HEAT") %>% filter(!is.na(resp.var)) %>% filter(!is.na(max.temp.C))
fwrite(heat, "fvfm.sub_heated.csv") # write out data for dose fitting

# Remove the control observations for downstream analyses (MHW Comparison) & write to file
db_heat <- db_merged %>% 
  filter(heat.control == "HEAT") %>%
  filter(!max.temp.C == "NA")
summary(db_heat$max.temp.C)
fwrite(db_heat, "heatedFilteredDatabase.csv")
fwrite(db_merged, "filteredDatabase_outliersRemoved.csv") # write out filtered db with outliers removed for 

