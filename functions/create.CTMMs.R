# This function takes a folder of migration sequences (in .dbf format)
# and creates regular brownian bridges from them
# location.error is in meters, cell size is in meters. max.lag is in hours
# contour is the % contour around UD for the footprint. Usually 99
# written by Jerod Merkle, 7 Feb 2019 (but based on Hall Sawyer's code)

create.CTMMs <- function(seqs_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/sequences",
                       CTMMs_out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/UD",
                       metadata_fldr="C:/Users/jmerkle/Desktop/Mapp2/tab6output",
                       cores=11, cell.size=50, mult4buff=0.2,
                       Information_Criteria="LOOCV", #Information criteria used for selection. Can be "AICc", "AIC", "BIC", "LOOCV" or none (NA). Use LOOCV for 12 hour data. Use something else for more frequent fixes
                       proj_of_dbfs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"){
  
  #manage packages
  if(all(c("rgdal","fields","foreign","snowfall","raster","data.table","ctmm") %in% installed.packages()[,1])==FALSE)
    stop("You must install the following packages: rgdal, fields, foreign, snowfall, raster, data.table, and ctmm")
  require(snowfall)
  require(raster)
  require(data.table)
  require(rgdal)
  require(ctmm)
  require(fields)
  require(foreign)
  
  #check the new directories
  if(dir.exists(CTMMs_out_fldr)==FALSE){
    dir.create(CTMMs_out_fldr)
  }
  if(length(dir(CTMMs_out_fldr))> 0)
    stop("Your CTMMs_out_fldr Has something in it. It should be empty!")
  if(dir.exists(metadata_fldr)==FALSE){
    dir.create(metadata_fldr)
  }
  
  print(paste0("Start time: ", Sys.time()))
  
  #load up teh data into a single database
  fls <- dir(seqs_fldr)
  print(paste0("You have ", length(fls), " sequences."))
  d <- do.call(rbind, lapply(1:length(fls), function(i){
    db <- read.dbf(paste(seqs_fldr, fls[i],sep="/"), as.is=TRUE)
    db$mig <- sub(".dbf","",fls[i])
    return(db)
  }))
  #check and make sure the columns are correct. 
  if(all(c("date","x","y") %in% names(d)) == FALSE) 
    stop("There is an issue with the columns in your sequences. See Error 2.")
  
  # create migration distances
  u <- unique(d$mig)
  
  mig_dists <- do.call(rbind, lapply(1:length(u), function(i){
    return(data.frame(mig=u[i], max_dist=max(rdist(d[d$mig==u[i],c("x","y")]))/1000))
  }))
  
  head(mig_dists)
  hist(mig_dists$max_dist)
  mig_dists$mean_max_dist <- mean(mig_dists$max_dist)
  mig_dists$sd_max_dist <- sd(mig_dists$max_dist)
  mig_dists$min_max_dist <- min(mig_dists$max_dist)
  mig_dists$max_max_dist <- max(mig_dists$max_dist)
  write.csv(mig_dists, file=paste(metadata_fldr,"/migration_distance_info.csv",sep=""), row.names=FALSE)
  
  #prepare the data for ctmm analysis
  d$date <- as.POSIXct(strptime(d$date,format = "%Y-%m-%d %H:%M:%S"), tz="GMT")  
  coordinates(d) <- c("x","y")
  proj4string(d) <- proj_of_dbfs
  spr_tr <- spTransform(d, CRS("+proj=longlat"))
  spr_tr$longitude <-  coordinates(spr_tr)[,1]
  spr_tr$latitude <-  coordinates(spr_tr)[,2]
  proj_new <- CRS(proj4string(spr_tr))
  d <- as.data.frame(spr_tr)
  
  # set extent
  dsp <- d # used data from all individs. here to prevent extent problems in pop-level averaging
  coordinates(dsp) <- c("longitude", "latitude")
  proj4string(dsp) <- CRS("+proj=longlat")
  dsp_proj <- spTransform(dsp, crs(proj_of_dbfs))
  dsp_proj$x <-  coordinates(dsp_proj)[,1]
  dsp_proj$y <-  coordinates(dsp_proj)[,2]
  ext <- extent(dsp_proj)
  multiplyers <- c((ext[2]-ext[1])*mult4buff, (ext[4]-ext[3])*mult4buff)
  ext <- extend(ext, multiplyers)   #this will be used inside the loop to extend each UD so it can be summed up
  
  u <- unique(d$mig)
  
  sfInit(parallel = T, cpus = cores)   #must change the cpus
  sfExport("d", "ext", "u", "CTMMs_out_fldr","Information_Criteria")
  sfLibrary(ctmm)
  sfLibrary(raster)
  sfLibrary(R.utils)
  sfLibrary(data.table)
  
  results <- do.call(rbind, sfClusterApplyLB(1:length(u), function(i){
    start.time <- Sys.time()

    df <- d[d$mig == u[i],]
    
    if(nrow(df) <4){   #throw out an NA line if there are less than 4 points.
      return(data.frame(Animal_ID=u[i], model=NA, deltaLOOCV=NA, deltaMSPE=NA,
                        DOF=NA, grid.size=NA, grid.cell.size=NA, date.created=Sys.time(),
                        execution_time=NA, Start.Date=min(df$date), End.Date=max(df$date),
                        num.days=length(unique(strftime(df$date, format = "%j", tz = "GMT"))),
                        errors="Less than 4 points"))
      
    }else{
      df <- df[order(df$date), ]
      
      telem <- as.telemetry(df, projection = CRS(proj_of_dbfs))
      #plot(telem)
      
      # store guesstimated parameters of models (not including BM)
      GUESS <- ctmm.guess(telem, interactive = FALSE)
      
      # run model selection
      FITS <- ctmm.select(telem, CTMM = GUESS, IC = Information_Criteria, 
                          verbose = TRUE, method = "pHREML")
      
      # summary - sorted by LOOCV
      #summary(FITS, IC = "LOOCV", MSPE = "velocity")
      
      # extract top 3 model fits
      mod1 <- FITS[[1]]
      mod2 <- FITS[[2]]
      mod3 <- FITS[[3]]
      LOOCV = c(mod1$LOOCV, mod2$LOOCV, mod3$LOOCV)
      #AICc = c(mod1$AICc, mod2$AICc, mod3$AICc)
      temp1 <- setDT(data.frame(summary(FITS[1:3], IC = "LOOCV", MSPE = "velocity")), 
                     keep.rownames = "model")
      colnames(temp1) <- c("model", "deltaLOOCV", "deltaMSPE", "DOF")
      temp2 <- data.frame(Animal_ID = rep(u[i], nrow(temp1)))
      temp3 <- cbind(temp2, temp1)
      
      # run BM model
      #inits.bm <- ctmm(tau = Inf)
      ##bm <- ctmm.fit(telem, CTMM = inits.bm)
      #temp4 <- data.frame(Animal_ID = id, AICc = bm$AICc, model = "BM", deltaAICc = NA)
      
      # Krige the data given the BM model
      #     out.krige <- occurrence(telem, CTMM = bm, grid = list(dr = c(cell.size,cell.size), align.to.origin = TRUE))
      
      # Krige the data given one of the other models
      out.krige <- occurrence(telem, CTMM = mod1, 
                              grid = list(dr = c(cell.size,cell.size), 
                                          align.to.origin = TRUE))
      
      
      #rm(dsp, dsp_proj)
      
      # write out raster file
      rkrig <- raster(out.krige, DF = "PDF", values = TRUE)#, template = ext) #, ext = ext)
      #rm(out.krige)
      rkrig[is.na(rkrig[])] <- 0   # make sure there are no NAs
      
      # set to 0 any values that are outside of the < 0.9999 contour
      cutoff <- sort(rkrig@data@values, decreasing=TRUE)
      vlscsum <- cumsum(cutoff)
      cutoff <- cutoff[vlscsum > .9999][1]
      rkrig@data@values[rkrig@data@values < cutoff] <- 0
      rm(vlscsum)
      
      # rescale probabilities so they equal 1
      prt <- rkrig/cellStats(rkrig, sum)
      rkrig2 <- extend(prt, ext, value = 0)
      
      # write out ascii file
      m <- as(rkrig2, "SpatialGridDataFrame")
      write.asciigrid(m, paste(CTMMs_out_fldr,"/",u[i],"_ASCII.asc",sep=""))
      
      temp3$grid.size <- ncell(prt)
      temp3$grid.cell.size <- res(rkrig2)[1]
      temp3$date.created <- Sys.time()
      temp3$execution_time <- paste(round(difftime(Sys.time(), start.time, units="min"),2)," minutes",sep="")
      temp3$Start.Date <- min(df$date)
      temp3$End.Date <- max(df$date)
      temp3$num.days <- length(unique(strftime(df$date, format = "%j", tz = "GMT")))
      temp3$errors <- "None"
      # some other stuff to output????
      return(temp3)   #return the model results 
    }
  }))
  sfStop()
  
  write.csv(results, file=paste(metadata_fldr,"/metadata.csv",sep=""), row.names=FALSE)
  print(paste0("End time: ", Sys.time()))
  return("Done. Check your folders.")
}#end function

