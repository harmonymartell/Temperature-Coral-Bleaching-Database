## Martell & Donner 2024: Headwinds to Understanding Stress Response Physiology: A Systematic Review Reveals Mismatch between Real and Simulated Marine Heatwaves

## USAGE: Database Filtering Script: This script takes the complete coral bleaching database
# and filters observations that were from 
# 1. studies of non-coral organisms 
# 2. studies with non adult coral lifestages
# 3. studies of multiple stressors
# 4. studies of field-based measurements of bleaching
# 5. studies with cold water bleaching observations
#    (with temperatures below 15C, the lower thermal threshold of Scleractinian corals (Kemp et al. 2012) &
#    with decreasing or decreasing variable thermal exposures)
# 6. studies on facultatively symbiotic corals (Oculina spp.)
# 7. studies that used invalid Fv/Fm (i.e., Yi instead of Fv/Fm)
# 8. studies of non-tropical Scleractinians (i.e., Mediterranean corals)

# written by Harmony Martell February 2024

rm(list = ls())

setwd("~/Desktop/UBC/_Projects/doseResponse_analysis/writing/_PLoSOne_Submission/") # set the working directory

## load the data
database=read.csv("_1_unfilteredDatabase.csv", header=TRUE, stringsAsFactors=TRUE)
database$sym = database$sym/10^6
database$genus <- substr(database$species_name, 1, 5) # grabs first four characters of species names

### Examine the complete dataset
names(database) # what variables are there
str(database) # what is the structure
dim(database) # n = 2998
length(unique(database$study_id))

## Count the number of observations from each filtered category
detach("package:dplyr", unload = TRUE)
library(plyr)

count(database$organism) # n = 10 orgs, 700 non coral
count(database$organism == "coral") # n = 2298
count(database$multiple_stressor) # n = 319
count(database$source) # n = 923 field, n = 24 both
count(database$life_stage) # n = 94
count(database$temp_exposure == "DECREASED") # n = 37
count(database$temp_exposure == "DECREASED VARIABLE") # n = 54
count(database$observation_temp < 15) # n = 10
  # excluded by cold temperature exposures (n = 37 + 54 + 10 = 101)
count(database$study_id == "S0202") # n = 54 # no valid response variables (used Yi instead of Fv/Fm)
count(database$study_id == "S0094") # n = 8 # facultatively symbiotic orgs
count(database$basin == "Mediterranean") # n = 10 # non-tropical Med coral

library(dplyr)
## 1. Filter out 700 non coral observations, n = 700/2998 observations
unique(database$organism)
library(dplyr)
db1 <- database %>% 
  filter(organism == "coral")
dim(db1) # n = 2,298, 700 non coral obs removed
detach("package:dplyr", unload = TRUE)
library(plyr)
count(db1$genus)
count(db1$species_name)

## 2. Filter out 87 non-adult life stage observations, n = 87/2998 observations
library(dplyr)
db2 <- db1 %>% filter(life_stage == "adult") # removes non-adult organism observations
dim(db2) # n = 2204

## 3. Filter out 228 multiple stressor observations, n = 319/2998 observations
db3 <-db2 %>% filter(multiple_stressor == "NO")
unique(db3$multiple_stressor)
dim(db3) # n = 1976

## 4. Filter out 888 field based observations, n = 923/2998 observations
db4 <- db3 %>% filter(!source == "F") # removes field only observations
dim(db4) # n = 1088

## 5. Filter out 81 cold water bleaching observations
  # includes 6 with temps < 15C bc below Scleractinian coral survival, n = 10/2998 observations
  # includes 75 decreasing and decreasing variable temp exposures, n = 75/2998 observations
db5 <- db4 %>%
  filter(!observation_temp < 15) %>%
  filter(!temp_exposure == "DECREASED") %>%
  filter(!temp_exposure == "DECREASED VARIABLE")
dim(db5) # n = 1007

## 5. Filter out 2 facultatively-symbiotic organisms (e.g., Oculina spp.), n = 8/2998 observations
db6 <- db5 %>% filter(!study_id == "S0094") # removes this single study with n = 2 obs in the remaining database
dim(db6) # n = 1005

## 6. Filter out 54 observations with clearly invalid Fv/Fm measurements
db7 <- db6 %>%
  filter(!study_id == "S0202") # removes 54 obs from Study 202 which measured only Yi with a pseudo dark period at the start of a rapid light curve, not Fv/Fm
dim(db7) # n = 951

## 7. Filter out 4 observations from single study from the Mediterranean which does not represent the majority of tropical Scleractinian corals
db8 <- db7 %>%
  filter(!study_id == "S0017")
dim(db8) # n = 947

library(data.table)
fwrite(db8,"_2_filteredDatabase.csv")
