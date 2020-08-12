# ---------------------------------------------------------------------------------------------------------------------- #

###### ------------------------- Code to calculate Metadata Values for Corridor Mapping -------------------------- ###### 
# -------------------------- Emily Gelzer & Jerod Merkle & Holly Copeland & Matthew Cuzzocreo -------------------------- #
# ---------------- Contact: egelzer@uwyo.edu | jmerkle@uwyo.edu | hcopelan@uwyo.edu | mcuzzocr@uwyo.edu ---------------- #

# ---------------------------------------------------------------------------------------------------------------------- #

####  Load Required Packages ####
library(lubridate)
library(stringr)
library(foreign)
library(tidyselect)
# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
#### CHANGE THESE VALUES TO MATCH YOUR DATA ####
#(no need to change anything below this section, UNLESS you changed output folder/files names!!!)

#Set Location of Migration Mapper Output File ----CHANGE THESE FILES AND FILEPATH!!!!!!!!!!!!
# ALSO CHECK YOUR DATE FORMAT AND MODIFY LINES 75 and 76, 149 & 150 TO MATCH DATE FORMAT IF NEEDED!!!!!!
setwd("G:/My Drive/Analysis/Deer/SFPeaks/SFPeaksMM060920") 
pointsout = "pointsOut.dbf"
metadata_migration = "metadata_migration.csv"
winter_metadata = "metadata_winter.csv"
# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
#### ANIMAL CAPTURE AND DATA COLLECTION ####

# Code to Calculate how many fall and spring squences are used in analyses -----------------------------------------------
# Get a character string of the files in the specified folder and have it print the number of squences within that folder
directory  <-  paste0(getwd(),"/UDs")
UDs <- vars_select(dir(path = directory), ends_with(".asc")) #only .asc files

directory  <-  paste0(getwd(),"/UDsWinter")
UDsW <- vars_select(dir(path = directory), ends_with(".asc")) #only .asc files

directory  <-  paste0(getwd(),"/UDs_pop")
UDs_pop <- dir(path = directory,pattern = "ASCII")

directory  <-  paste0(getwd(),"/UDs_popWinter")
UDs_popW <- dir(path = directory, pattern = "ASCII")

directory <-  paste0(getwd(),"/UDs")
fls_UDs_fall <- vars_select(UDs, contains("fa")) #select from the UDs list created in line 32 rather than the directory to avoid duplicate fall UDs
fls_UDs_spring <- vars_select(UDs, contains("sp")) #select from the UDs list created in line 32 rather than the directory to avoid duplicate spring UDs

print(paste0("You have ", length(UDs), " migration sequences from " ,length(UDs_pop), " individuals"))
print(paste0("You have ", length(fls_UDs_fall) , " fall sequences."))
print(paste0("You have ", length(fls_UDs_spring), " spring sequences."))
print(paste0("You have ", length(UDsW), " winter sequences from " ,length(UDs_popW), " individuals"))


# Sample Size ------------------------------------------------------------------------------------------------------------
UDs_pop <- as.list (UDs_pop)
UDs_popw <- as.list (UDs_popW)
allUDspop <- c(UDs_pop, UDs_popW)
unique1 <- unique(allUDspop)
print(paste0("Sample size = ", length(unique1), " individuals "))

rm(UDs, fls_UDs_fall, fls_UDs_spring, UDsW, UDs_popW)


# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
#### CORRIDOR AND STOPOVER SUMMARY ####

# Code to Calculate Median Migration Start and End Dates (Spring and Fall) -----------------------------------------------

# Bring in metadata.csv
mig.metadata <- read.csv(paste0(getwd(),"/metadata/", metadata_migration), stringsAsFactors=TRUE)
str(mig.metadata)

mig.metadata$Season <- ifelse(grepl(pattern = "sp", mig.metadata$input.file) == TRUE, "sp", "fa")

# Convert start_season and end_season columns to "date" format
mig.metadata$Start.Date_1 <- as.Date (mig.metadata$Start.Date, format = "%Y-%m-%d")
mig.metadata$End.Date_1 <- as.Date (mig.metadata$End.Date, format = "%Y-%m-%d")

# Make a month column because the order matters here
mig.metadata$start.month <- month(mig.metadata$Start.Date_1)
mig.metadata$end.month <- month(mig.metadata$End.Date_1)


# Need to organize months by 11, 12, 01, 02 which requires changing the years
year(mig.metadata$Start.Date_1) <- ifelse(mig.metadata$start.month > 8, 2016, 2017)   # MAY NEED TO CHANGE THE 8 IF SEASONS ARE WEIRD!!!!!!!!!!!!!!!!
year(mig.metadata$End.Date_1) <- ifelse(mig.metadata$end.month > 8, 2016, 2017)    # MAY NEED TO CHANGE THE 8 IF SEASONS ARE WEIRD!!!!!!!!!!!!!!!!

# Get median values!
sp.start <- sort(mig.metadata$Start.Date_1[mig.metadata$Season == "sp"])
print(paste0("The Median Spring Migration Start Time is ", sp.start[round(length(sp.start)/2,0)], ". Ignore the Year"))

sp.end <- sort(mig.metadata$End.Date_1[mig.metadata$Season == "sp"])
print(paste0("The Median Spring Migration End Time is ", sp.end[round(length(sp.end)/2,0)], ". Ignore the Year"))


fa.start <- sort(mig.metadata$Start.Date_1[mig.metadata$Season == "fa"])
print(paste0("The Median Fall Migration Start Time is ", fa.start[round(length(fa.start)/2,0)], ". Ignore the Year"))

fa.end <- sort(mig.metadata$End.Date_1[mig.metadata$Season == "fa"])
print(paste0("The Median Fall Migration End Time is ", fa.end[round(length(fa.end)/2,0)], ". Ignore the Year"))


rm(fa.end, fa.start, sp.end, sp.start, mig.metadata)

# Calculate average migration length by season ---------------------------------------------------------------------------
# Bring in metadata.csv
metadata <- read.csv(paste0(getwd(),"/metadata/", metadata_migration), stringsAsFactors=TRUE)
metadata$input.file<- as.character(metadata$input.file)

# split season string and calculate length
sp.length<-round(mean(metadata$num.days[!is.na(str_extract(metadata$input.file, pattern = "sp"))]))
fa.length<-round(mean(metadata$num.days[!is.na(str_extract(metadata$input.file, pattern = "fa"))]))

print(paste0("The mean spring migration length is ", sp.length, " days"))
print(paste0("The mean fall migration length is ", fa.length, " days"))

rm(sp.length, fa.length)

# To calculate Migration corridor length (mean and max) please refer to Holly Copeland's python script on GitHub! --------

# Calculate the High Use Migration Corridor Areas and Stopover Areas (in acres) ------------------------------------------
# Calculate the Migration Corridor Areas and Stopover Areas (in acres) ------------------------------------------

migarea <- read.csv("finalProducts/output_areas.csv")
migarea$migacreslow <- round(migarea$GRIDCODE_1_area * 247.10)
migarea$migacres2 <- round(migarea$GRIDCODE_2_area * 247.10)
migarea$migacresmed <- round(migarea$GRIDCODE_10_area * 247.10)
migarea$migacreshigh <- round(migarea$GRIDCODE_20_area * 247.10)

migarea$stopacres <- round(migarea$stopover_area * 247.10)

print (paste0("Migration Corridor Area (low use) = ", migarea$migacreslow, " acres"))
print (paste0("Migration Corridor Area (2 or more) = ", migarea$migacres2, " acres"))
print (paste0("Migration Corridor Area (med use) = ", migarea$migacresmed, " acres"))
print (paste0("Migration Corridor Area (high use) = ", migarea$migacreshigh, " acres"))

print (paste0("Stopover Area = ", migarea$stopacres, " acres"))

# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
#### WINTER RANGE SUMMARY ####

# Winter Median Start and End Date ---------------------------------------------------------------------------------------

# Bring in metadata.csv
winter.metadata <- read.csv(paste0(getwd(),"/metadata/",winter_metadata), stringsAsFactors=TRUE)
str(winter.metadata)

# Convert start_season and end_season columns to "date" format
winter.metadata$Start.Date_1 <- as.Date (winter.metadata$Start.Date, format = "%Y-%m-%d")
winter.metadata$End.Date_1 <- as.Date (winter.metadata$End.Date, format = "%Y-%m-%d")

# Make a month column because the order matters here
winter.metadata$start.month <- month(winter.metadata$Start.Date_1)
winter.metadata$end.month <- month(winter.metadata$End.Date_1)

# Need to organize months by 11, 12, 01, 02 which requires changing the years
year(winter.metadata$Start.Date_1) <- ifelse(winter.metadata$start.month > 8, 2016, 2017)   # MAY NEED TO CHANGE THE 8 IF SEASONS ARE WEIRD!!!!!!!!!!!!!!!!
year(winter.metadata$End.Date_1) <- ifelse(winter.metadata$end.month > 8, 2016, 2017)    # MAY NEED TO CHANGE THE 8 IF SEASONS ARE WEIRD!!!!!!!!!!!!!!!!

# Get median values!
winter.start <- sort(winter.metadata$Start.Date_1)
print(paste0("The Median Winter Start Time is ", winter.start[round(length(winter.start)/2,0)], ". Ignore the Year"))

winter.end <- sort(winter.metadata$End.Date_1)
print(paste0("The Median Winter End Time is ", winter.end[round(length(winter.end)/2,0)], ". Ignore the Year"))


rm(winter.end, winter.start, winter.metadata)

# Calculate Winter Length (mean) -----------------------------------------------------------------------------------------
# Bring in metadata.csv
metadata_winter <- read.csv(paste0(getwd(),"/metadata/",winter_metadata), stringsAsFactors=TRUE)
metadata_winter$input.file<- as.character(metadata_winter$input.file)

# split season string and calculate length
winter.length<-round(mean(metadata_winter$num.days))

print(paste0("The mean winter length is ", winter.length, " days"))

# Calculate Winter Range (50% countour) Area -----------------------------------------------------------------------------
winterarea <- read.csv("finalProductsWinter/output_areas_winter.csv")
winterarea$winter05acres <- round(winterarea$GRIDCODE_0.5_area * 247.10)
print (paste0("The Winter Range Area= ", winterarea$winter05acres, " acres"))






