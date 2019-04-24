#code to create a spatial lines data frame from a sequences folder
# by Jerod Merkle, 25 Feb 2019

seqs_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/sequences"
out_fldr="C:/Users/jmerkle/Desktop/Mapp2/tab6output/hack_corridors"
proj_of_dbfs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
buff=200
cell.size=50

hack.corridors <- function(seqs_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/sequences",
                            out_fldr="C:/Users/jmerkle/Desktop/Mapp2/tab6output/hack_corridors",
                            proj_of_dbfs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0",
                            buff=200, cell.size=50){   
  if(all(c("rgdal","foreign","stringr","rgeos","raster") %in% installed.packages()[,1])==FALSE)
    stop("You must install the following packages: rgdal, foreign,rgeos, raster, and stringr")
  require(rgdal)
  require(foreign)
  require(stringr)
  require(rgeos)
  require(raster)
  
  #check the new directories
  if(dir.exists(out_fldr)==FALSE){
    dir.create(out_fldr)
  }
  if(length(dir(out_fldr))> 0)
    stop("Your out_fldr Has something in it. It should be empty!")
  
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
  
  d$date <- as.POSIXct(strptime(d$date,format = "%Y-%m-%d %H:%M:%S"), tz="GMT")  
  d <- d[order(d$mig, d$date),]
  coordinates(d) <- c("x","y")
  proj4string(d) <- proj_of_dbfs
  
  u <- unique(d$mig)
  lns <- lapply(1:length(u), function(e){
    return(Lines(list(Line(coordinates(d[d$mig==u[e],]))),ID=u[e]))
  })
  lns <- SpatialLines(lns)
  df <- data.frame(mig=u, 
                   firstdate=do.call(c, lapply(u, function(e){min(d$date[d$mig==e], na.rm=TRUE)})),
                   lastdate=do.call(c, lapply(u, function(e){max(d$date[d$mig==e], na.rm=TRUE)})))
  rownames(df) <- df$mig
  df$id <- str_split_fixed(df$mig, "_", 2)[,1]
  df$season <- substr(str_split_fixed(df$mig, "_", 2)[,2],1,2)
  df$year <- as.numeric(paste0("20",substr(str_split_fixed(df$mig, "_", 2)[,2],3,4)))
  lns <- SpatialLinesDataFrame(lns, df)
  proj4string(lns) <- proj_of_dbfs
  
  head(lns@data)
  
  ids <- lns$id
  ids_unique <- unique(ids)
  
  id_buffs <- do.call(rbind, lapply(1:length(ids_unique), function(i){
    tmp <- lns[lns$id %in% ids_unique[i],]
    tmpb <- gBuffer(tmp, byid=TRUE, width=buff)
    tmpb <- gUnaryUnion(tmpb, id=NULL)
    tmpb <- SpatialPolygonsDataFrame(tmpb, data=data.frame(id=ids_unique[i]))
  }))
  nrow(id_buffs)

  # library(sf)
  # id_buffs2 <- as(id_buffs, "sf")
  # 
  # u <- st_combine(id_buffs2)

  ext <- extent(id_buffs)
  grd <- raster(ext)
  res(grd) <- cell.size     
  projection(grd) <- proj4string(id_buffs)
  
  grd <- rasterize(id_buffs[1,], grd, background=0)
  for(i in 2:nrow(id_buffs)){
    grd <- grd+rasterize(id_buffs[i,], grd, background=0)
  }
  plot(grd)
  writeRaster(grd, filename=paste0(out_fldr, "/footprints_hack.img"), format="HFA")
  # writeOGR(lns, out_fldr, "migration_lines", driver="ESRI Shapefile")
  return("Finished. Check your out folder.")
}
