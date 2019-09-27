###### Code to calculate Metadata Values for Corridor Mapping
###### Emily Gelzer & Jerod Merkle & Holly Copeland
###### Contact: egelzer@uwyo.edu


### Load Packages
library(lubridate)
library(stringr)
library(foreign)

###CHANGE THESE VALUES TO MATCH YOUR DATA (no need to change anything below this section)
#Set Location of Migration Mapper Output File ----CHANGE THESE FILES AND FILEPATH!!!!!!!!!!!!
setwd("C:/Holly/Work/Projects/EasternYellowstoneMuleDeer2016/Data/CMT_Dubois_8262019") 
migtime = "migtime_20190701131444.csv"
pointsout = "pointsOut_20190701131444.dbf"

###Get Sample Size

# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------

### Code to Calculate Median Migration Start and End Dates (Spring and Fall)

# Bring in migtime_number.csv
migtime <- read.csv(migtime, stringsAsFactors=TRUE)
str(migtime)
# Convert start_season and end_season columns to "date" format

migtime$startSpring <- ymd(migtime$startSpring)
migtime$endSpring <- ymd(migtime$endSpring)

migtime$startFall <- ymd(migtime$startFall)
migtime$endFall <- ymd(migtime$endFall)

# Create a column for Julien Day for each start_season and end_Season

migtime$startSpring_J <- as.numeric(format(migtime$startSpring, "%j"))
migtime$endSpring_J <- as.numeric(format(migtime$endSpring, "%j"))
migtime$startFall_J <- as.numeric(format(migtime$startFall, "%j"))
migtime$endFall_J <- as.numeric(format(migtime$endFall, "%j"))

# Create data frame with just Julien Days
JD_startSpring <- data.frame(migtime$startSpring_J)
JD_endSpring <- data.frame(migtime$endSpring_J)
JD_startFall <- data.frame(migtime$startFall_J)
JD_endFall <- data.frame(migtime$endFall_J)

# Rename for convience 
names(JD_startSpring) <- c("startSpring_J")
names(JD_endSpring)<-c ("endSpring_J")
names(JD_startFall) <- c("startFall_J")
names(JD_endFall) <- c("endFall_J")

# Sort the Julien Day columns to find the median Julien Day for each start_season and end_season

median_startSpring_J <- median(JD_startSpring$startSpring_J, na.rm = TRUE)
median_startSpring <-as.Date(median_startSpring_J, origin=as.Date("1899-12-31"))
medianstartSpring <- paste(month(median_startSpring), day(median_startSpring), sep ="-")

median_endSpring_J <- median(JD_endSpring$endSpring_J, na.rm = TRUE)
median_endSpring <-as.Date(median_endSpring_J, origin=as.Date("1899-12-31"))
medianendSpring <- paste(month(median_endSpring), day(median_endSpring), sep ="-")

JD_startFall$startFall_J<- ifelse(JD_startFall$startFall_J<100,JD_startFall$startFall_J+365,JD_startFall$startFall_J)
median_startFall_J <- median(JD_startFall$startFall_J, na.rm = TRUE)
median_startFall_J <- ifelse(median_startFall_J>365, median_startFall_J-365,median_startFall_J)
median_startFall <- as.Date(median_startFall_J, origin=as.Date("1899-12-31"))
medianstartFall <- paste(month(median_startFall), day(median_startFall), sep ="-")

median_endFall_J <- median(JD_endFall$endFall_J, na.rm = TRUE)
median_endFall <- as.Date(median_endFall_J, origin=as.Date("1899-12-31"))
medianendFall <- paste(month(median_endFall), day(median_endFall), sep ="-")

print(paste0("The Median Spring Migration Start Time is ", medianstartSpring))
print(paste0("The Median Spring Migration End Time is ", medianendSpring))
print(paste0("The Median Fall Migration Start Time is ", medianstartFall))
print(paste0("The Median Fall Migration End Time is ", medianendFall))

rm(median_startSpring, median_startSpring_J, medianstartSpring,median_endSpring, median_endSpring_J, medianendSpring,
   median_startFall, median_startFall_J, medianstartFall, median_endFall, median_endFall_J, medianendFall)


# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------

###### Code to Find Begining and End Migration Dates

summary(migtime)

SpringMig_MIN <-min(migtime$startSpring, na.rm = TRUE)


FallMig_MIN <-min(migtime$startFall, na.rm = TRUE)


BeginingMig <- ifelse(FallMig_MIN < SpringMig_MIN, FallMig_MIN, SpringMig_MIN)

BeginingMig <- as.Date(BeginingMig,origin=as.Date("1899-12-31"))

FallMig_MAX <-max(migtime$endFall, na.rm = TRUE)

SpringMig_MAX <-max(migtime$endSpring, na.rm = TRUE)

EndMig <- ifelse(SpringMig_MAX > FallMig_MAX, SpringMig_MAX, FallMig_MAX)

EndMig <- as.Date(EndMig,origin=as.Date("1899-12-31"))


print(paste0("The Begining Migration Date is ", BeginingMig))
print(paste0("The End Migration Date is ", EndMig))

# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------

###### Code to Calculate Average fix rate (in hour) for corridor and winter range analyses 

# Bring in the pointsout DBF file ---- CHANGE THIS FILEPATH
pointsout <- read.dbf(pointsout, as.is = FALSE)

head(pointsout)

table(pointsout$seas)

# -----------------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------------

###### Calculate average migration length by season
###### Matthew Cuzzocreo  
###### Contact: mcuzzocr@uwyo.edu

### Bring in metadata.csv
metadata <- read.csv("metadata.csv")
metadata$input.file<- as.character(metadata$input.file)

#split season string and calculate length
sp.length<-round(mean(metadata$num.days[!is.na(str_extract(metadata$input.file, pattern = "sp"))]))
fa.length<-round(mean(metadata$num.days[!is.na(str_extract(metadata$input.file, pattern = "fa"))]))

print(paste0("The mean spring migration length is ", sp.length, " days"))
print(paste0("The mean fall migration length is ", fa.length, " days"))


# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------

###### Code to Calculate how many fall and spring squences are used in analyses 

### Get a character string of the files in the specified folder and have it print the number of squences within that folder
directory  <-  paste0(getwd(),"/UDs")
UDs <- dir(path = directory, pattern = "*.asc")

directory  <-  paste0(getwd(),"/UDsW")
UDsW <- dir(path = directory,pattern = "*.asc")

directory  <-  paste0(getwd(),"/UDs_pop")
UDs_pop <- dir(path = directory,pattern = "*.asc")

directory  <-  paste0(getwd(),"/UDs_popW")
UDs_popW <- dir(path = directory, pattern = "*.asc")

directory <-  paste0(getwd(),"/UDs")
fls_UDs_fall <- list.files(directory, pattern="[0-9]{1,3}[_]{1}[fa]{1,2}[0-9]{1,2}[_]{1}[ASCII]{1,5}\\.asc")
fls_UDs_spring <- list.files(directory, pattern="[0-9]{1,3}[_]{1}[sp]{1,2}[0-9]{1,2}[_]{1}[ASCII]{1,5}\\.asc")

print(paste0("You have ", length(UDs), " migration sequences from " ,length(UDs_pop), " individuals"))
print(paste0("You have ", length(fls_UDs_fall) , " fall sequences."))
print(paste0("You have ", length(fls_UDs_spring), " spring sequences."))
print(paste0("You have ", length(UDsW), " winter sequences from " ,length(UDs_popW), " individuals"))

##
#Sample Size
UDs_pop <- as.list (UDs_pop)
UDs_pop <- as.list (UDs_popW)
allUDspop <- c(UDs_pop, UDs_popW)
unique1 <- unique(allUDspop)
print(paste0("Sample size = ", length(unique1), " individuals "))

rm(UDs, fls_UDs_fall, fls_UDs_spring)
# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------

###### Calculate average winter length
###### Adopted from:  Calculate average migration length by season, by Matthew Cuzzocreo, mcuzzocr@uwyo.edu 

### Bring in metadata.csv
metadata_winter <- read.csv("metadata_winter.csv")
metadata_winter$input.file<- as.character(metadata_winter$input.file)

#split season string and calculate length
winter.length<-mean(metadata_winter$num.days)

print(paste0("The mean winter length is ", winter.length, " days"))

rm(winter.length)

# --------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------

### Code to Calculate Winter Median Start and End Dates 

# Bring in migtime_number.csv ---- CHANGE THIS FILEPATH!!!!!!!!!!!
metadata_winter <- read.csv("metadata_winter.csv")

# Convert start_season and end_season columns to "date" format

metadata_winter$startWinter <- ymd_hms(metadata_winter$Start.Date)
metadata_winter$endWinter <- ymd_hms(metadata_winter$End.Date)

# Isolate just the date

metadata_winter$startWinter <- date(metadata_winter$Start.Date)
metadata_winter$endWinter <- date(metadata_winter$End.Date)

# Create a column for Julien Day for each start_season and end_Season

metadata_winter$startWinter_J <- as.numeric(format(metadata_winter$startWinter, "%j"))
metadata_winter$endWinter_J <- as.numeric(format(metadata_winter$endWinter, "%j"))


# Create data frame with just Julien Days
JD_startWinter <- data.frame(metadata_winter$startWinter_J)
JD_endWinter <- data.frame(metadata_winter$endWinter_J)

# Rename for convience 
names(JD_startWinter) <- c("startWinter_J")
names(JD_endWinter)<-c ("endWinter_J")

# Sort the Julien Day columns to find the median Julien Day for each start_season and end_season

JD_startWinter$startWinter_J<- ifelse(JD_startWinter$startWinter_J<100,JD_startWinter$startWinter_J+365,JD_startWinter$startWinter_J)
median_startWinter_J <- median(JD_startWinter$startWinter_J, na.rm = TRUE)
median_startWinter_J <- ifelse(median_startWinter_J>365, median_startWinter_J-365,median_startWinter_J)
median_startWinter <- as.Date(median_startWinter_J, origin=as.Date("1899-12-31"))
medianstartWinter <- paste(month(median_startWinter), day(median_startWinter), sep ="-")

median_endWinter_J <- median(JD_endWinter$endWinter_J, na.rm = TRUE)
median_endWinter <- as.Date(median_endWinter_J, origin=as.Date("1899-12-31"))
medianendWinter <- paste(month(median_endWinter), day(median_endWinter), sep ="-")

print(paste0("The Median Winter Start Time is ", medianstartWinter))
print(paste0("The Median Winter End Time is ", medianendWinter))

####Calculate areas
migarea <- read.csv("final_products/output_areas.csv")
migarea$migacres <- round(migarea$GRIDCODE_1_area * 247.10)
migarea$stopacres <- round(migarea$stopover_area * 247.10)
print (paste0("Migration Corridor Area =", migarea$migacres, " acres"))
print (paste0("Stopover Area =", migarea$stopacres, " acres"))

winterarea <- read.csv("final_productsW/output_areas_winter.csv")
winterarea$winter05acres <- round(winterarea$GRIDCODE_0.5_area * 247.10)
print (paste0("The Winter Range Area=", winterarea$winter05acres, " acres"))

rm(metadata_winter, JD_startWinter, JD_endWinter, median_endWinter, median_endWinter_J, medianendWinter,
   median_startWinter, median_startWinter_J)
