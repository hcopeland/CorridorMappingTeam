# This function takes the outputs from Migration Mapper tab 6 and creates
# dbf files of each individual animal's migration sequences
# written by Jerod Merkle, 7 Feb 2019
# do not add '.shp' to the end of the shpfl_name

create.seqs.W <- function(shpfl_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output", 
                          shpfl_name= "pointsOut_20190207093941" ,
                          idname="newUid",  #name of the column representing animal ID
                          datename="nwMstrD",   #name of the column representing date in POSIX format
                          mig.metadata.file="C:/Users/jmerkle/Desktop/Mapp2/tab6output/metadata.csv",  # metadata file from migration part of analysis
                          qtl.end.fall.mig=0.95, 
                          qtl.start.spring.mig=0.05,
                          out_fldr="C:/Users/jmerkle/Desktop/Mapp2/tab6output/sequencesW",
                          out_proj="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"){
  #manage packages
  if(all(c("rgdal","foreign","stringr","sf") %in% installed.packages()[,1])==FALSE)
    stop("You must install the following packages: rgdal, foreign, sf, and stringr")
  require(rgdal)
  require(foreign)
  require(stringr)
  require(sf)
  
  #check teh new directory
  if(dir.exists(out_fldr)==FALSE){
    dir.create(out_fldr)
  }
  if(length(dir(out_fldr))> 0)
    stop("Your out_fldr Has something in it. It should be empty!")
  
  fl <- read.csv(mig.metadata.file)
  fl$seas <- substr(str_split_fixed(fl$input.file, "_",2)[,2],1,2)
  fl$yr <- substr(str_split_fixed(fl$input.file, "_",2)[,2],3,4)
  fl$yr <- as.numeric(ifelse(as.numeric(fl$yr) > 80, paste0("19",fl$yr), paste0("20",fl$yr)))
  fl$Start.Date <- as.POSIXct(strptime(fl$Start.Date,format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
  fl$End.Date <- as.POSIXct(strptime(fl$End.Date,format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
  # fl$Start.Date <- as.POSIXct(strptime(paste0("2018", substr(fl$Start.Date,5,19)),format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
  # fl$End.Date <- as.POSIXct(strptime(paste0("2018", substr(fl$End.Date,5,19)),format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
  
  sp <- fl[fl$seas == "sp",]
  fa <- fl[fl$seas == "fa",]
  
  yr_rng <- unique(fl$yr)
  
  dts <- do.call(rbind, lapply(0:length(yr_rng), function(i){
    if(i == 0){
      wintr <- paste(yr_rng[i+1]-1, yr_rng[i+1], sep="_")
      wint.start <- round(quantile(fa$End.Date[fa$yr == yr_rng[i+1]-1], probs = qtl.end.fall.mig),0)
      wint.end <- round(quantile(sp$Start.Date[sp$yr == yr_rng[i+1]], probs = qtl.start.spring.mig),0)
    }else{
      wintr <- paste(yr_rng[i], yr_rng[i]+1, sep="_")
      wint.start <- round(quantile(fa$End.Date[fa$yr == yr_rng[i]], probs = qtl.end.fall.mig),0)
      wint.end <- round(quantile(sp$Start.Date[sp$yr == yr_rng[i]+1], probs = qtl.start.spring.mig),0)
    }
    return(data.frame(winter=wintr, wint.start=substr(wint.start,1,10), 
                      wint.end=substr(wint.end,1,10)))
  }))
  dts$wint.start <- as.character(dts$wint.start)
  dts$wint.end <- as.character(dts$wint.end)
  mnths <- as.numeric(substr(na.omit(dts$wint.start), 6,7))
  mn_strt <- mean(as.Date(strptime(paste0(ifelse(mnths > 8, "2017","2018"),"-",substr(na.omit(as.character(dts$wint.start)),6,10)),format = "%Y-%m-%d")))
  if(as.numeric(substr(mn_strt,6,7)) > 8){
    yrs <- substr(dts$winter[is.na(dts$wint.start)==TRUE],1,4)
  }else{
    yrs <- substr(dts$winter[is.na(dts$wint.start)==TRUE],6,9)
  }
  dts$wint.start[is.na(dts$wint.start)==TRUE] <- paste0(yrs, "-", substr(mn_strt,6,10))
  
  mnths <- as.numeric(substr(na.omit(dts$wint.end), 6,7))
  mn_end <- mean(as.Date(strptime(paste0(ifelse(mnths > 8, "2017","2018"),"-",substr(na.omit(as.character(dts$wint.end)),6,10)),format = "%Y-%m-%d")))
  if(as.numeric(substr(mn_end,6,7)) > 8){
    yrs <- substr(dts$winter[is.na(dts$wint.end)==TRUE],1,4)
  }else{
    yrs <- substr(dts$winter[is.na(dts$wint.end)==TRUE],6,9)
  }
  dts$wint.end[is.na(dts$wint.end)==TRUE] <- paste0(yrs, "-", substr(mn_end,6,10))

  print("Here are the winter range start and end dates by year. Average winter start and end dates were added where corresponding migration dates were missing.")
  print(dts)
  
  loca <- str_locate_all(mig.metadata.file,"/")
  loca <- loca[[1]][nrow(loca[[1]]),1]
  write.csv(dts, file=paste0(substr(mig.metadata.file,1,loca),"winter_dates.csv"), row.names = FALSE)
  
  #change columns back to posix
  dts2 <- dts
  dts$wint.start <- as.POSIXct(strptime(paste(dts$wint.start," 00:00:01",sep=""), format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
  dts$wint.end <- as.POSIXct(strptime(paste(dts$wint.end," 23:59:59",sep=""), format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
  
  # read in the shapefile
  print("Loading the shapefile (this can take a few minutes)...")
  # d <- readOGR(shpfl_fldr, shpfl_name)
  d <- st_read(shpfl_fldr, shpfl_name)
  d <- as(d, "Spatial")
  
  print(paste0("Your shapefiles has ", nrow(d), " rows."))
  if(all(c(idname,datename) %in% names(d)) == FALSE) 
    stop("There is an issue with the columns in your sequences. See Error 1.")
  
  #reproject to new projection, if necessary
  proj <- proj4string(d)
  if(proj != out_proj){    
    d <- spTransform(d, CRS(out_proj))
  }
  
  #reduce the dataset to columns of interest
  d <- d[,c(idname,datename)]
  names(d) <- c("id","date")
  
  #fix the dates  (it is OK to specify GMT, since all dates will be in GMT!)
  d$date <- as.POSIXct(strptime(d$date,format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
  
  if(any(is.na(d$date)))
    stop("You have a problem with your date column! It does not seem to be a character that was previously in POSIX format.")
  
  #keep only the data within the date ranges
  d$year <- as.numeric(strftime(d$date, format = "%Y", tz = "GMT"))
  d$month <- as.numeric(strftime(d$date, format = "%m", tz = "GMT"))

  d <- lapply(1:nrow(dts), function(i){
    tmp <- d[d$date >= dts$wint.start[i] & d$date <= dts$wint.end[i],]
    if(nrow(tmp)==0){
      return(NULL)
    }
    tmp$winter <- dts$winter[i]
    #create a unique winter period ID (so the year of teh winter means the year at the end of the winter period)
    if(as.numeric(substr(dts$wint.start[i],6,7)) <6){ #if you specified start of winter after 1 January, instead of in the previous year
      tmp$id_yr <- paste0(tmp$id,"_wi", substr(tmp$year,3,4))
    }else{
      tmp$id_yr <- paste0(tmp$id,"_wi", substr(ifelse(tmp$month>6,tmp$year+1,tmp$year),3,4))
    }

    return(tmp)
  })
  d[sapply(d, is.null)] <- NULL
  d <- do.call(rbind, d)
  
  if(any(is.na(d$id_yr)))
    stop("You have an issue. See error 2.")
  
  #identify the unique sequences
  seqs <- unique(d$id_yr)

  print("Identifying sequences...")
  # loop through the sequences
  
  for(i in 1:length(seqs)){
    tmp <- d[d$id_yr == seqs[i],]
    tmp <- cbind(as.data.frame(tmp)[,c("id","date")],coordinates(tmp))    # get it out of sp object
    tmp$date <- as.character(tmp$date)    #need to switch this back to character for dbf files
    names(tmp) <- c("id","date","x","y")   # rename columns
    write.dbf(tmp, file = paste0(out_fldr,"/",seqs[i],".dbf"))    #write as dbf
  }
  print(paste0("Finished creating ",length(dir(out_fldr))," sequences. Check your sequences folder."))
  return(dts2)
}
