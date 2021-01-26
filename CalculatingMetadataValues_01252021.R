# ---------------------------------------------------------------------------------------------------------------------- #

###### ------------------------- Code to calculate Metadata Values for Corridor Mapping -------------------------- ###### 
###### ----------------------------------------- Updated January 2021 -------------------------------------------- ###### 
# --------------------- Emily Gelzer & Jerod Merkle & Holly Copeland & Matthew Cuzzocreo & Lucas Olson -------------------------- #
# ---------------- Contact: egelzer@uwyo.edu | jmerkle@uwyo.edu | hcopelan@uwyo.edu | lolson@azgfd.gov ---------------- #

# ---------------------------------------------------------------------------------------------------------------------- #

####  Load Required Packages ####
library(lubridate)
library(stringr)
library(foreign)
library(tidyselect)
library(dplyr)
# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
#### CHANGE THESE VALUES TO MATCH YOUR DATA ####
#(no need to change anything below this section, UNLESS you changed output folder/files names!!!)

#Set Location of Migration Mapper Output File ----CHANGE THESE FILES AND FILEPATH!!!!!!!!!!!!
setwd("G:/My Drive/Analysis/Deer/Analysis/Paunsaugunt_Kaibab_2020") 
metadata_migration = "metadata_migration.csv"
# ------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------
#### ANIMAL CAPTURE AND DATA COLLECTION ####

# Code to Calculate how many fall and spring squences are used in analyses -----------------------------------------------
# Get a character string of the files in the specified folder and have it print the number of squences within that folder
directory  <-  paste0(getwd(),"/UDs")
UDs <- vars_select(dir(path = directory), ends_with(".asc")) #only .asc files

directory  <-  paste0(getwd(),"/UDsW")#May need to change this directory!
UDsW <- vars_select(dir(path = directory), ends_with(".asc")) #only .asc files

directory  <-  paste0(getwd(),"/UDs_pop")
UDs_pop <- dir(path = directory,pattern = "ASCII")

directory  <-  paste0(getwd(),"/UDs_popW")#May need to change this directory!
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
mig.metadata <- read.csv(paste0(getwd(),'/',metadata_migration), stringsAsFactors=TRUE)

mig.metadata$Season <- ifelse(grepl(pattern = "sp", mig.metadata$input.file) == TRUE, "sp", "fa")
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
metadata <- read.csv(paste0(getwd(),"/", metadata_migration), stringsAsFactors=TRUE)
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

migarea <- read.csv("final_products/output_areas.csv")#May need to change this directory!
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
# Calculate Winter Length(mean), Start Date(Median), End Date(Median) -----------------------------------------------------------------------------------------
# Bring in migration metadata.csv
df <- metadata
df$seas <- substr(str_split_fixed(df$input.file, "_",2)[,2],1,2)
df$Start.Date <- as.POSIXct(strptime(df$Start.Date,format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
df$End.Date <- as.POSIXct(strptime(df$End.Date,format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
df$id <- str_split_fixed(df$input.file, "_",2)[,1]

df<-df %>%
  group_by(id) %>%
  arrange(Start.Date, .by_group = TRUE)

for(i in 1:nrow(df)){
  if(df$seas[i] == 'sp'){
    df$End.Date[i] <- NA
  }
}

for(i in 1:nrow(df)){
  if(df$seas[i] == 'fa'){
    df$Start.Date[i] <- NA
  }
}

#### Calculate winter start/end dates
df2<-df
df2$W.Start.Date <- as.Date (df2$End.Date, format = "%Y-%m-%d")
df2$W.End.Date <- as.Date (df2$Start.Date, format = "%Y-%m-%d")

# Need to organize months by 11, 12, 01, 02 which requires changing the years
year(df2$W.Start.Date) <- ifelse(month(df2$W.Start.Date) > 8, 2016, 2017)   # MAY NEED TO CHANGE THE 8 IF SEASONS ARE WEIRD!!!!!!!!!!!!!!!!
year(df2$W.End.Date) <- ifelse(month(df2$W.End.Date) > 8, 2016, 2017)    # MAY NEED TO CHANGE THE 8 IF SEASONS ARE WEIRD!!!!!!!!!!!!!!!!

# Get median winter start/end values
winter.end <- median(df2$W.End.Date, na.rm = T)
winter.start <- median(df2$W.Start.Date, na.rm = T)

## Calculate mean winter length
df2<-df %>%
  mutate(Start.Date = replace(Start.Date, row_number() == 1, NA))
df2<-ungroup(df2)
df2<-mutate(df2, days = lead(df2$Start.Date)-df2$End.Date)

winter.length<-round(mean(df2$days, na.rm = TRUE))

print(paste0("The Median Winter Start Time is ", winter.start, ". Ignore the Year."))
print(paste0("The Median Winter End Time is ", winter.end, ". Ignore the Year."))
print(paste0("The mean winter length is ", winter.length, " days."))

# Calculate Winter Range (50% countour) Area -----------------------------------------------------------------------------
winterarea <- read.csv("final_productsW/output_areas_winter.csv") #May need to change this directory!
winterarea$winter05acres <- round(winterarea$GRIDCODE_0.5_area * 247.10)
print (paste0("The Winter Range Area= ", winterarea$winter05acres, " acres"))





