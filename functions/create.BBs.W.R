# This function takes a folder of migration sequences (in .dbf format)
# and creates regular brownian bridges from them
# location.error is in meters, cell size is in meters. max.lag is in hours
# written by Jerod Merkle, 24 Feb 2019 (but based on Hall Sawyer's code)

create.BBs.W <- function(seqs_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/sequencesW",
                       BBs_out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/UDsW",
                       metadata_fldr="C:/Users/jmerkle/Desktop/Mapp2/tab6output",
                       mindays=30,   #if an individual animal has less than this many days in a given year's winter data, it will be removed
                       cores=11, location.error=20, cell.size=50, max.lag=8, time.step=5,mult4buff=0.2,
                       proj_of_dbfs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"){
  
  #manage packages
  if(all(c("rgdal","foreign","R.utils","stringr","BBMM","snowfall","raster") %in% installed.packages()[,1])==FALSE)
    stop("You must install the following packages: rgdal, foreign, BBMM, snowfall, raster, and stringr")
  require(rgdal)
  require(foreign)
  require(stringr)
  require(BBMM)
  require(raster)
  require(snowfall)
  require(R.utils)
  
  #check the new directories
  if(dir.exists(BBs_out_fldr)==FALSE){
    dir.create(BBs_out_fldr)
  }
  if(length(dir(BBs_out_fldr))> 0)
    stop("Your BBs_out_fldr Has something in it. It should be empty!")
  if(dir.exists(metadata_fldr)==FALSE){
    dir.create(metadata_fldr)
  }
  
  print(paste0("Start time: ", Sys.time()))
  
  #load up teh data into a single database
  fls <- dir(seqs_fldr)
  print(paste0("You have ", length(fls), " sequences."))
  d <- do.call(rbind, lapply(1:length(fls), function(i){
    db <- read.dbf(paste(seqs_fldr, fls[i],sep="/"), as.is=TRUE)
    db$wint <- sub(".dbf","",fls[i])
    return(db)
  }))
  #check and make sure the columns are correct. 
  if(all(c("date","x","y") %in% names(d)) == FALSE) 
    stop("There is an issue with the columns in your sequences. See Error 2.")
  
  d$date <- as.POSIXct(strptime(d$date,format = "%Y-%m-%d %H:%M:%S"), tz="GMT") 
  if(any(is.na(d$date))) 
    stop("There is an issue with your date columns. See Error 2a.")
  
  # remove any sequences that do not have > mindays
  jul <- as.numeric(strftime(d$date, format = "%j", tz = "GMT"))
  u <- unique(d$wint)
  ids2keep <- do.call(c, lapply(1:length(u), function(i){
    if(length(unique(jul[d$wint == u[i]])) >= mindays){
      return(u[i])
    }else{
      return(NULL)
    }
  }))
  print(paste0("You are removing ", length(u)-length(ids2keep), " sequences because there are less than ", mindays, " days worth of data!"))
  
  d <- d[d$wint %in% ids2keep,]
  
  #Create a grid to calculate BBs over
  dsp <- d
  coordinates(dsp) <- c("x","y")
  proj4string(dsp) <- proj_of_dbfs
  ext <- raster::extent(dsp)
  multiplyers <- c((ext[2]-ext[1])*mult4buff, (ext[4]-ext[3])*mult4buff)   # add about mult4buff around the edges of your extent (you can adjust this if necessary)
  ext <- raster::extend(ext, multiplyers)
  grd <- raster(ext)
  res(grd) <- cell.size     
  projection(grd) <- proj4string(dsp)
  
  print(paste0("The number of cells in your population grid is ", ncell(grd), "."))
  
  u <- unique(d$wint)
  
  sfInit(parallel = T, cpus = cores)   #must change the cpus
  sfExport("d", "grd", "u", "max.lag", "location.error",
           "BBs_out_fldr","time.step","mult4buff")
  sfLibrary(BBMM)
  sfLibrary(raster)
  sfLibrary(R.utils)
  regBB <- do.call(rbind, sfClusterApplyLB(1:length(u), function(i){
    start.time <- Sys.time()
    
    temp <- d[d$wint==u[i],]
    temp <- temp[order(temp$date),]
    
    #prepare only the cells to run BB over
    ext2 <- raster::extent(temp)
    multiplyers <- c((ext2[2]-ext2[1])*mult4buff, (ext2[4]-ext2[3])*mult4buff)   # add about mult4buff around the edges of your extent (you can adjust this if necessary)
    ext2 <- raster::extend(ext2, multiplyers)
    cels <- raster::cellsFromExtent(grd, ext2)

    # this is the function to calculate the regular BB
    bb <- R.utils::withTimeout({
      try(BBMM::brownian.bridge(x=temp$x,
                                y=temp$y,
                                time.lag=diff(as.numeric(temp$date)/60),
                                area.grid=coordinates(grd)[cels,],
                                max.lag=max.lag*60,
                                time.step=time.step,
                                location.error=location.error), #this is the location error of your collars
          silent=TRUE)
    }, envir=environment(), timeout = 10800, onTimeout = "warning")   #it'll stop if it runs for more than 10800 seconds (3 hours)
    
    
    #write out results to file too, so if there is an error you don't loose all your work!
    if(class(bb)=="try-error" | class(bb)== "character"){
      return(data.frame(input.file=u[i],
                        brownian.motion.variance=NA,
                        grid.size=NA,
                        grid.cell.size=NA,
                        date.created=Sys.time(),
                        execution_time=paste(round(difftime(Sys.time(), start.time, units="min"),2)," minutes",sep=""),
                        num.locs=nrow(temp),
                        Start.Date=min(temp$date),
                        End.Date=max(temp$date),
                        num.days=length(unique(strftime(temp$date, format = "%j", tz = "GMT"))),
                        errors=ifelse(class(bb)=="try-error", attr(bb, "condition")$message, "Ran too long")))    
    }else{
      
      
      #set to 0 any values that are outside of the < 0.9999 contour
      cutoff <- sort(bb$probability, decreasing=TRUE)
      vlscsum <- cumsum(cutoff)
      cutoff <- cutoff[vlscsum > .9999][1]
      bb$probability[bb$probability < cutoff] <- 0
      
      #rescale probabilities so they equal 1
      bb$probability <- bb$probability/sum(bb$probability)
      
      #output ASCII file
      m <- data.frame(x=coordinates(grd)[,1], y=coordinates(grd)[,2],z=0)
      m$z[cels] <- bb$probability
      m <- SpatialPixelsDataFrame(points = m[c("x", "y")], data=m)
      m <- as(m, "SpatialGridDataFrame")
      write.asciigrid(m, paste(BBs_out_fldr,"/",u[i],"_ASCII.asc",sep=""), attr=3)
      
      #gather summary info
      return(data.frame(input.file=u[i],
                        brownian.motion.variance=round(bb[[1]],2),
                        grid.size=length(bb$x),
                        grid.cell.size=abs(bb$x[1]-bb$x[2]),
                        date.created=Sys.time(),
                        execution_time=paste(round(difftime(Sys.time(), start.time, units="min"),2)," minutes",sep=""),
                        num.locs=nrow(temp),
                        Start.Date=min(temp$date),
                        End.Date=max(temp$date),
                        num.days=length(unique(strftime(temp$date, format = "%j", tz = "GMT"))),
                        errors="None"))
    }
  }))
  sfStop()
  if(length(u) != nrow(regBB))
    print("There is a problem (see error 1).")
  if(all(is.na(regBB$brownian.motion.variance))){
    print("All of your BBs failed!")
  }else{
    if(any(is.na(regBB$brownian.motion.variance))){
      print(paste0("WARNING: You have ",nrow(regBB)-length(na.omit(regBB$brownian.motion.variance))," BBs that failed. Take a look at the errors column in your metadata file."))
    }else{
      print(paste0("All ", nrow(regBB), " sequences produced a successful BB."))
    }
    if(any(na.omit(regBB$brownian.motion.variance) > 8000))
      print("WARNING: You have BB motion variances that are > 8,000!!! Take a look at your metadata file.")
  }
  write.csv(regBB, file=paste(metadata_fldr,"/metadata_winter.csv",sep=""), row.names=FALSE)
  print(paste0("End time: ", Sys.time()))
  return("Done. Check your folders.")
}#end function

