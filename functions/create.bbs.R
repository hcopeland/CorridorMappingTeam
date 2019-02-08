# This function takes a folder of migration sequences (in .dbf format)
# and creates regular brownian bridges from them
# location.error is in meters, cell size is in meters. max.lag is in hours
# contour is the % contour around UD for the footprint. Usually 99
# written by Jerod Merkle, 7 Feb 2019 (but based on Hall Sawyer's code)

create.bbs <- function(seqs_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/sequences",
                       BBs_out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/UD",
                       footprint_out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/Footprint",
                       metadata_fldr="C:/Users/jmerkle/Desktop/Mapp2/tab6output/",
                       cores=11, location.error=20, cell.size=50, max.lag=8, contour=99,
                       proj_of_dbfs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"){
  
  #manage packages
  if(all(c("rgdal","foreign","stringr","BBMM","snowfall","raster") %in% installed.packages()[,1])==FALSE)
    stop("You must install the following packages: rgdal, foreign, BBMM, snowfall, raster, and stringr")
  require(rgdal)
  require(foreign)
  require(stringr)
  require(BBMM)
  require(raster)
  require(snowfall)
  
  #check the new directories
  if(dir.exists(BBs_out_fldr)==FALSE){
    dir.create(BBs_out_fldr)
  }
  if(length(dir(BBs_out_fldr))> 0)
    stop("Your BBs_out_fldr Has something in it. It should be empty!")
  if(dir.exists(footprint_out_fldr)==FALSE){
    dir.create(footprint_out_fldr)
  }
  if(length(dir(footprint_out_fldr))> 0)
    stop("Your footprint_out_fldr Has something in it. It should be empty!")
  if(dir.exists(metadata_fldr)==FALSE){
    dir.create(metadata_fldr)
  }
  
  #load up teh data into a single database
  fls <- dir(seqs_fldr)
  d <- do.call(rbind, lapply(1:length(fls), function(i){
    db <- read.dbf(paste(seqs_fldr, fls[i],sep="/"), as.is=TRUE)
    db$mig <- sub(".dbf","",fls[i])
    return(db)
  }))
  #check and make sure the columns are correct. 
  if(all(c("date","x","y") %in% names(d)) == FALSE) 
    stop("There is an issue with the columns in your sequences. See Error 2.")
  
  d$date <- as.POSIXct(strptime(d$date,format = "%Y-%m-%d %H:%M:%S"), tz="GMT")  
  
  #Create a grid to calculate BBs over
  dsp <- d
  coordinates(dsp) <- c("x","y")
  proj4string(dsp) <- proj_of_dbfs
  ext <- extent(dsp)
  multiplyers <- c((ext[2]-ext[1])*0.2, (ext[4]-ext[3])*0.2)   # add about 20% around the edges of your extent (you can adjust this if necessary)
  ext <- extend(ext, multiplyers)
  grd <- raster(ext)
  res(grd) <- cell.size     
  projection(grd) <- proj4string(dsp)
  
  
  u <- unique(d$mig)
  
  sfInit(parallel = T, cpus = cores)   #must change the cpus
  sfExport("d", "grd", "u", "max.lag", "location.error","contour",
           "BBs_out_fldr","footprint_out_fldr")
  sfLibrary(BBMM)
  sfLibrary(raster)
  regBB <- do.call(rbind, sfClusterApplyLB(1:length(u), function(i){
    start.time <- Sys.time()
    
    temp <- d[d$mig==u[i],]
    temp <- temp[order(temp$date),]
    # this is the function to calculate the regular BB
    bb <- try(brownian.bridge(x=temp$x,
                              y=temp$y,
                              time.lag=diff(as.numeric(temp$date)/60),
                              area.grid=coordinates(grd),
                              max.lag=max.lag*60,
                              location.error=location.error), #this is the location error of your collars
              silent=TRUE)    
    #write out results to file too, so if there is an error you don't loose all your work!
    if(class(bb)=="try-error"){
      return(data.frame(input.file=u[i],
                        brownian.motion.variance=NA,
                        grid.size=NA,
                        grid.cell.size=NA,
                        date.created=Sys.time(),
                        execution_time=paste(round(difftime(Sys.time(), start.time, units="min"),2)," minutes",sep=""),
                        num.locs=nrow(temp),
                        Start.Date=min(temp$date),
                        End.Date=max(temp$date),
                        num.days=length(unique(strftime(temp$date, format = "%j", tz = "GMT")))))    
    }else{
      
      #set to 0 any values that have prob<0.00000001
      bb$probability[bb$probability < 0.00000001] <- 0
      
      # bb$x<-bb$x[bb$probability >= 0.00000001]
      # bb$y<-bb$y[bb$probability >= 0.00000001]
      # bb$probability<-bb$probability[bb$probability >= 0.00000001]
      
      #rescale probabilities
      bb$probability <- bb$probability/sum(bb$probability)
      
      #output ASCII file
      m <- data.frame(x=bb$x, y=bb$y,z=bb$probability)
      m <- SpatialPixelsDataFrame(points = m[c("x", "y")], data=m)
      m <- as(m, "SpatialGridDataFrame")
      write.asciigrid(m, paste(BBs_out_fldr,"/",u[i],"_ASCII.asc",sep=""), attr=3)
      
      #output 99% contour
      contours <- bbmm.contour(bb, levels=contour, plot=F)
      
      # Create data.frame indicating cells within the each contour and export as Ascii Grid
      contour.99 <- data.frame(x = bb$x, y = bb$y, probability = bb$probability)
      # contour.99 <- contour.99[contour.99$probability >= contours$Z[1],]
      contour.99$in.out <- ifelse(contour.99$probability >= contours$Z[1], 1, 0)
      
      m <- SpatialPixelsDataFrame(points = contour.99[c("x", "y")], data=contour.99)
      m <- as(m, "SpatialGridDataFrame")
      write.asciigrid(m, paste(footprint_out_fldr,"/",u[i],"_99pct_contour.asc",sep=""), attr=ncol(m))
      
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
                        num.days=length(unique(strftime(temp$date, format = "%j", tz = "GMT")))))
    }
  }))
  sfStop()
  if(length(u) != nrow(regBB))
    print("There is a problem (see error 1).")
  write.csv(regBB, file=paste(metadata_fldr,"/metadata.csv",sep=""), row.names=FALSE)
  print("Success! Check your folders.")
  return("Done.")
}#end function

