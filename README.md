# Temperature-Coral-Bleaching-Database

The .xlsx file named databaseREADME includes a Column Name, Description and Instructions on how each column of the database was populated for all entries.

The .R files are the scripts used to complete the analyses.

The _1_unfilteredDatabase.csv file includes all entries before _a priori_ filtering for use in _1_dbFiltering_Martell&Donner2024.R

_1_dbFiltering_Martell&Donner2024.R takes the unfiltered database (__1_unfilteredDatabase.csv_) and systematically filters it, generating __2_filteredDatabase.csv_ (provided here, written out by this script).

_2_mergeCount_Martell&Donner2024.R script merges the filtered database (__2_filteredDatabase.csv_) with the heating metrics of the heated observations (__2_heatingParams_HeatedObs_toMerge.csv_); the heating metrics are the output from the heatingParameterCalculation.R; the script also removes outlying observations, counts the remaining observations and parses into smaller response variable datasets for downstream analyses. Note that the control heated metrics can be calculated by modification of this script. Please email drharmonymartell@gmail.com for the raw heated or control temperature series, if desired.          

_3_MHWcomparison_Martell&Donner2024.R performs a comparison of laboratory vs. marine heatwave observations for figures and stats of control vs. heated observation of response variables
   *** NOTE: MHW datafiles are large, please email drharmonymartell@gmail.com for file transfer ***
   
_4_RespVarStats_Martell&Donner2024.R includes response variable statistics and figures for the control vs heated comparisons

_5_EScalcs_Martell&Donner2024.R performs effect size calculations (SMD, Hedge's g) and bubble plots found in the Suppplement (Figs S9, S10)

_6_fitting_Martell&Donner2024.R includes statistics and figures of hormetic fits of each response variable against each heating metric, with all observations and excluding those that were not ecologically relevant, that is, outside of the parameters observated during recent MHWs (2010-2018) on coral reefs.
