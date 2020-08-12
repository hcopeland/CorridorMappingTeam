# This function takes a folder of migration sequences (in .dbf format)
# and creates regular brownian bridges from them
# location.error is in meters, cell size is in meters. max.lag is in hours
# written by Jerod Merkle, 24 Feb 2019 (but based on Hall Sawyer's code)

create.BBs.W <- function(seqs_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/sequencesW",
                         BBs_out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/UDsW",
                         metadata_fldr="C:/Users/jmerkle/Desktop/Mapp2/tab6output",
                         mindays=30,   #if an individual animal has less than this many days in a given year's winter data, it will be removed
                         BMvar=NULL, # if Null, will run regular BB and calculate motion variance. If a number is specified, it will invoke the Forced motion variance method
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
  
  #make the custom brownian bridge function... FOR FMV
  BrownianBridgeCustom<-function (x, y, time.lag, location.error, area.grid = NULL, cell.size = NULL, 
                                  time.step = 10, max.lag = NULL, BMvar=BMvar) 
  {
    if (is.null(x) | is.null(y) | (length(x) != length(y))) {
      stop("data is missing or unequal number of x and y coordinates")
    }
    if (is.null(location.error)) 
      stop("must specify 'location.error'")
    if (is.null(area.grid) & is.null(cell.size)) {
      stop("'area.grid' or 'cell.size' must be specified")
    }
    if (!is.null(area.grid) & is.null(cell.size)) {
      cell.size <- abs(area.grid[1, 1] - area.grid[2, 1])
    }
    if (is.null(area.grid) & !is.null(cell.size)) {
      range.x <- range(x)
      range.y <- range(y)
      min.grid.x <- round(range.x[1] - 1 * sd(x))
      max.grid.x <- round(range.x[2] + 1 * sd(x))
      min.grid.y <- round(range.y[1] - 1 * sd(y))
      max.grid.y <- round(range.y[2] + 1 * sd(y))
      x. <- seq(min.grid.x, max.grid.x, cell.size)
      y. <- seq(min.grid.y, max.grid.y, cell.size)
      area.grid <- merge(x., y.)
    }
    if (is.null(max.lag)) {
      max.lag = max(time.lag) + 1
    }
    if (length(location.error) == 1) {
      location.error <- rep(location.error, length(x))
    }
    n.locs <- length(x)
    # BMvar <- brownian.motion.variance(n.locs, time.lag, location.error, 
    #                                   x, y, max.lag)
    BMvar <- rep(BMvar, times = length(x))
    if (is.null(time.step)) 
      time.step <- 10
    grid.size <- nrow(area.grid)
    probability <- rep(0, grid.size)
    T.Total <- sum(time.lag)
    bbmm <- vector("list", 4)
    names(bbmm) <- c("Brownian motion variance", "x", "y", "probability")
    class(bbmm) <- "bbmm"
    probability <- NULL
    int <- 0
    for (i in 1:(n.locs - 1)) {
      if (time.lag[i] <= max.lag) {
        theta <- NULL
        tm <- 0
        while (tm <= time.lag[i]) {
          alpha <- tm/time.lag[i]
          mu.x <- x[i] + alpha * (x[i + 1] - x[i])
          mu.y <- y[i] + alpha * (y[i + 1] - y[i])
          sigma.2 <- time.lag[i] * alpha * (1 - alpha) * 
            BMvar[i] + ((1 - alpha)^2) * (location.error[i]^2) + 
            (alpha^2) * (location.error[i + 1]^2)
          ZTZ <- (area.grid[, 1] - mu.x)^2 + (area.grid[, 
                                                        2] - mu.y)^2
          theta <- (1/(2 * pi * sigma.2)) * exp(-ZTZ/(2 * 
                                                        sigma.2))
          int <- int + theta
          tm <- tm + time.step
        }
      }
    }
    probability <- int/T.Total
    probability <- probability/sum(probability)
    bbmm[[4]] <- probability
    bbmm[[1]] <- BMvar[1]
    bbmm[[2]] <- area.grid[, 1]
    bbmm[[3]] <- area.grid[, 2]
    return(bbmm)
  }
  #----------- END OF CUSTOM BB FUNCTION (for FMV)
  
  
  #load up teh data into a single database
  fls <- list.files(seqs_fldr, ".dbf$")
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
  
  sfInit(parallel = T, cpus = ifelse(length(u)<cores, length(u), cores))   #must change the cpus
  sfExport("d", "grd", "u", "max.lag", "location.error","BMvar","BrownianBridgeCustom",
           "BBs_out_fldr","time.step","mult4buff")
  sfLibrary(BBMM)
  sfLibrary(raster)
  sfLibrary(R.utils)
  regBB <- do.call(rbind, sfClusterApplyLB(1:length(u), function(i){
    start.time <- Sys.time()
    
    temp <- d[d$wint==u[i],]
    temp <- temp[order(temp$date),]
    
    jul <- as.numeric(strftime(temp$date, format = "%j", tz = "GMT"))
    if(length(unique(jul)) < mindays){
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
                        errors=paste0("Less than ",mindays, " days of data.")))
    }else{
      if(nrow(temp) <= length(unique(jul))){
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
                          errors="Only 1 point per day, on average."))
      }else{
        if(nrow(temp) < 4){
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
                            errors="Less than 4 points."))
        }else{
          
          #prepare only the cells to run BB over
          ext2 <- raster::extent(temp)
          multiplyers <- c((ext2[2]-ext2[1])*mult4buff, (ext2[4]-ext2[3])*mult4buff)   # add about mult4buff around the edges of your extent (you can adjust this if necessary)
          ext2 <- raster::extend(ext2, multiplyers)
          cels <- raster::cellsFromExtent(grd, ext2)
          
          # this is the function to calculate the regular BB
          if(class(BMvar)=="NULL"){
            bb <- R.utils::withTimeout({
              try(BBMM::brownian.bridge(x=temp$x,
                                        y=temp$y,
                                        time.lag=diff(as.numeric(temp$date)/60),
                                        area.grid=coordinates(grd)[cels,],
                                        max.lag=max.lag*60,
                                        time.step=time.step,
                                        location.error=location.error), #this is the location error of your collars
                  silent=TRUE)
            }, envir=environment(), timeout = 43200, onTimeout = "warning")
          }else{ #THIS IS FMV code
            bb <- R.utils::withTimeout({
              try(BrownianBridgeCustom(x=temp$x,
                                       y=temp$y,
                                       time.lag=diff(as.numeric(temp$date)/60),
                                       area.grid=coordinates(grd)[cels,],
                                       max.lag=max.lag*60,
                                       time.step=time.step,
                                       BMvar=BMvar,
                                       location.error=location.error), #this is the location error of your collars
                  silent=TRUE)
            }, envir=environment(), timeout = 43200, onTimeout = "warning")
          }
          
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
                              errors=ifelse(class(bb)=="try-error", attr(bb, "condition")$message, "Ran too long or other issue")))    
          }else{
            if(length(bb$probability) == 1 | all(is.na(bb$probability))==TRUE){
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
                                errors="BB Failed for unknown reason"))    
              
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
          }
        }
      }
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

