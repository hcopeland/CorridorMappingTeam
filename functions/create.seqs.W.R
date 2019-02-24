# This function takes the outputs from Migration Mapper tab 6 and creates
# dbf files of each individual animal's migration sequences
# written by Jerod Merkle, 7 Feb 2019
# do not add '.shp' to the end of the shpfl_name

create.seqs.W <- function(shpfl_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output", 
                          shpfl_name= "pointsOut_20190207093941" ,
                          idname="newUid",  #name of the column representing animal ID
                          datename="nwMstrD",   #name of the column representing date in POSIX format
                          winterStart="12-01",  #start of winter %m-%d  
                          winterEnd="02-28",  #end of winter  %m-%d
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
  
  #remove data from outside of the date range
  d$year <- as.numeric(strftime(d$date, format = "%Y", tz = "GMT"))
  d$month <- as.numeric(strftime(d$date, format = "%m", tz = "GMT"))
  d$jul <- as.numeric(strftime(d$date, format = "%j", tz = "GMT"))
  
  if(as.numeric(substr(winterStart,1,2)) <6){ #if you specified start of winter after 1 January, instead of in the previous year
    d$WinterS <- as.POSIXct(strptime(paste(d$year,"-",winterStart," 00:00:01",sep=""), format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
    d$WinterE <- as.POSIXct(strptime(paste(d$year,"-",winterEnd," 23:59:59",sep=""), format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
  }else{
    d$WinterS <- as.POSIXct(strptime(paste(ifelse(d$month > 6, d$year, d$year-1),"-",winterStart," 00:00:01",sep=""), format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
    d$WinterE <- as.POSIXct(strptime(paste(ifelse(d$month > 6, d$year+1, d$year),"-",winterEnd," 23:59:59",sep=""), format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
  }
  
  #remove all data except for 
  d <- d[d$date >= d$WinterS & d$date <= d$WinterE,]

  #create a unique winter period ID (so the year of teh winter means the year at the end of the winter period)
  if(as.numeric(substr(winterStart,1,2)) <6){ #if you specified start of winter after 1 January, instead of in the previous year
    d$id_yr <- paste0(d$id,"_wi", substr(d$year,3,4))
  }else{
    d$id_yr <- paste0(d$id,"_wi", substr(ifelse(d$month>6,d$year+1,d$year),3,4))
  }
  
  if(any(is.na(d$id_yr)))
    stop("You have an issue. See error 2.")
  
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
  return(paste0("Finished creating ",length(dir(out_fldr))," sequences. Check your sequences folder."))
}
