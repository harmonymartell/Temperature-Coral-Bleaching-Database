## Martell & Donner 2024: Headwinds to Understanding Stress Response Physiology: A Systematic Review Reveals Mismatch between Real and Simulated Marine Heatwaves

## USAGE: Heat parameter calculating script

# This script takes datetime and temperature data from a processed & formatted .csv file,
# reads it into R, converts datetime from ISO POSIXct standardized format to unix timestamp, 
# identifies the maximum temperature in the series (max_temp);
# calculates the duration of heating and converts from seconds to days
# calculates the accummulated heating (DHD) sum of magnitudes divided by duration
# calculates the temp range between start and max, and then the rate (using the time of the last min temp before the max, and the timing of the first max temp)

## Load in packages
library(readr)
library(psych)
#library(ggplot2)
library(dplyr)

rm(list = ls())

####################################
#set working directory
setwd("~/Desktop/UBC/_Projects/doseResponse_analysis/Analysis/heatrateRecalc/")

# data <- read_csv("filesToProcess/S0282_OBO4804_temp_processed.csv", col_names = TRUE)
#  data <- read_csv("filesToProcess/S0033_OBO5101_temp_processed.csv", col_names = TRUE)
files <- list.files("S0128_S0148files/", pattern="*.csv",full.names=TRUE) # specifies the directory with the files in it
files

# initalise  
heatparams <- data.frame()

for (file in files)
  {
  data <- read_csv(file, col_names = TRUE) # read in the file
  data$unixdate <- as.numeric(data$datetime) # create a col of unix dates from datetime
  
  ## Filename
  filename = file
  
    pdf(file = paste0(filename, ".pdf", sep="") ) # create a pdf file with the filename containing the study ID and obs ID
    plot(data$datetime,data$temperature) # plot the data
    #ggplot(data = data, aes(x = datetime, y = temperature) ) + geom_point() + ggtitle(filename)
    dev.off() # close the pdf file
  
  ## max_temp (DEGREES C)
  max_temp = max(data$temperature); max_temp # maximum temperature of entire temperature series
  
  ## duration (DAYS)
  duration_data = data[data$temperature > data$temperature[1], ] # capture only durations greater than starting temp (ambient)
  duration = round( (duration_data$unixdate[length(duration_data$unixdate)] - duration_data$unixdate[1]) /86400 , 2); duration # gets the difference in seconds, divided by the numbers of seconds in a day, spits out the duration in days

  ## dhd (DEGREES C DAYS)
  tempdiffs = (data$temperature - data$temperature[1]); tempdiffs # calculate the difference from the starting temp, which is the assumed ambient temp
  greaterthan0 = sum(round(tempdiffs[tempdiffs>0],2)); greaterthan0 # returns a vector of all values greater than 0, which are accummulated heating amounts from each tempdiff greater than 
  dhd = round(greaterthan0/duration,1); dhd
  
  ## heating rate (DEGREES C/DAY)
  
  # get time when max temp first reached
  maxsubsetdata <- data[data$temperature == max_temp, ]; #maxsubsetdata # subsets only data that have max temps
  time_firstMaxTemp = maxsubsetdata$unixdate[1]; time_firstMaxTemp # get the unix timestamp of first max temperature
  
  # get min temp and time when last start temp reached
  min_to_max <- data[data$unixdate <= time_firstMaxTemp, ]; headTail(min_to_max) # get only the data before the max temp is reached
  plot(min_to_max$datetime, min_to_max$temperature)
  
  min_temp = min(min_to_max$temperature); min_temp # gets the minimum temperature value from all data before the max is reached
  
  mintemps <- min_to_max %>% # keeps only minimum temps, and only the ones before first max temp
    filter(temperature == min_temp) %>% # keeps only minimum temperature values from subsetted series before max temp is reached
    filter(unixdate <= time_firstMaxTemp) # keeps only min temp values on the date before the first max temp is reached, should be redundant
  mintemps
  
  time_lastMinTemp = tail(mintemps$unixdate, n = 1); time_lastMinTemp  # get the timestamp of the last Start Temp
  
  range = max_temp - min_temp; range # calculates the difference in temperature between max and start temps
  
  ratetime = (time_firstMaxTemp - time_lastMinTemp)/86400; ratetime # gets time (in seconds) between last start temp and first max temp, then converts to DAYS
  
  rate = range/ratetime; rate # this is the heating rate in DEGREES C PER DAY
  
  output <- c(filename, max_temp, duration, dhd, range, rate)
  
  ## Write Filename, max_temp, duration, dhd, rate to file
  heatparams <- rbind(heatparams, output)
  }

names(heatparams)
colnames(heatparams) <- c("filename", "max temp", "duration", "dhd", "range", "rate")
headTail(heatparams)

write_csv(heatparams, "heatingParams_HeatedObs_S0128S0148.csv") 
