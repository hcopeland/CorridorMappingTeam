# This function takes the outputs from Migration Mapper tab 6 and creates
# dbf files of each individual animal's migration sequences
# written by Jerod Merkle, 7 Feb 2019
# do not add '.shp' to the end of the shpfl_name

create.seqs <- function(shpfl_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output", 
                        shpfl_name= "pointsOut_20190207093941" ,
                        migtbl_name="migtime_20190207093941.csv",
                        out_fldr="C:/Users/jmerkle/Desktop/Mapp2/tab6output/sequences",
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
  
  mt <- read.csv(paste(shpfl_fldr, migtbl_name,sep="/"))
  
  # read in the shapefile
  print("Loading the shapefile (this can take a few minutes)...")
  # d <- readOGR(shpfl_fldr, shpfl_name)
  d <- st_read(shpfl_fldr, shpfl_name)
  d <- as(d, "Spatial")
    
  
  if(all(c("newUid","nwMstrD") %in% names(d)) == FALSE) 
    stop("There is an issue with the columns in your sequences. See Error 1.")
  
  #reproject to new projection, if necessary
  proj <- proj4string(d)
  if(proj != out_proj){    
    d <- spTransform(d, CRS(out_proj))
  }
  
  #reduce the dataset to columns of interest
  d <- d[,c("newUid","nwMstrD")]
  names(d) <- c("id","date")
  
  #fix the dates  (it is OK to specify GMT, since all dates will be in GMT!)
  d$date <- as.POSIXct(strptime(d$date,format = "%Y-%m-%d %H:%M:%S"), tz="GMT")  
  mt$startSpring <- as.POSIXct(strptime(paste(mt$startSpring, "00:00:01"),format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
  mt$endSpring <- as.POSIXct(strptime(paste(mt$endSpring, "23:59:59"),format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
  mt$startFall <- as.POSIXct(strptime(paste(mt$startFall, "00:00:01"),format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
  mt$endFall <- as.POSIXct(strptime(paste(mt$endFall, "23:59:59"),format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
  
  print("Identifying sequences...")
  # loop through spring migrations
  for(i in 1:nrow(mt)){
    if(is.na(mt$startSpring[i])==FALSE & mt$shouldSpringRun[i] == TRUE & is.na(mt$shouldSpringRun[i])==FALSE){
      tmp <- d[d$id == mt$newUid[i] & d$date >= mt$startSpring[i] & d$date <= mt$endSpring[i],]
      tmp <- as.data.frame(tmp)    # get it out of sp object
      tmp$date <- as.character(tmp$date)    #need to switch this back to character for dbf files
      names(tmp) <- c("id","date","x","y")   # rename columns
      write.dbf(tmp, file = paste(out_fldr,"/",mt$newUid[i], "_sp",substr(mt$year_bio[i],3,4),".dbf",sep=""))    #write as dbf
    }else{
      next
    }
  }
  
  # loop through spring migrations
  for(i in 1:nrow(mt)){
    if(is.na(mt$startFall[i])==FALSE & mt$shouldFallRun[i] == TRUE & is.na(mt$shouldFallRun[i])==FALSE){
      tmp <- d[d$id == mt$newUid[i] & d$date >= mt$startFall[i] & d$date <= mt$endFall[i],]
      tmp <- as.data.frame(tmp)    # get it out of sp object
      tmp$date <- as.character(tmp$date)    #need to switch this back to character for dbf files
      names(tmp) <- c("id","date","x","y")   # rename columns
      write.dbf(tmp, file = paste(out_fldr,"/",mt$newUid[i], "_fa",substr(mt$year_bio[i],3,4),".dbf",sep=""))    #write as dbf
    }else{
      next
    }
  }
  print("Finished creating the sequences.")
  return("Success!")
}
