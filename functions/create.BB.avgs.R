# This function takes in individual BBs and creates averaged individual BBs and footprints
# and population UDs and population footprints
# written by Jerod Merkle, 8 Jan 2019
# based on code from Hall Sawyer

create.BB.avgs <- function(BBs_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/UD",
                        pop_BBs_out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/UD_pop",
                        pop_footprint_out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/Footprint_pop",
                        contour=99,
                        proj_of_ascs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"){
  
  #manage packages
  if(all(c("stringr", "raster", "BBMM") %in% installed.packages()[,1])==FALSE)
    stop("You must install the following packages: raster, BBMM, and stringr")
  require(stringr)
  require(raster)
  require(BBMM)
  
  #check the new directories
  if(dir.exists(pop_BBs_out_fldr)==FALSE){
    dir.create(pop_BBs_out_fldr)
  }
  if(length(dir(pop_BBs_out_fldr))> 0)
    stop("Your pop_BBs_out_fldr Has something in it. It should be empty!")
  if(dir.exists(pop_footprint_out_fldr)==FALSE){
    dir.create(pop_footprint_out_fldr)
  }
  if(length(dir(pop_footprint_out_fldr))> 0)
    stop("Your pop_footprint_out_fldr Has something in it. It should be empty!")
  
  print(paste0("Start time: ", Sys.time()))
  
  fls <- dir(BBs_fldr)
  print(paste("You have ", length(fls), " sequences with successful BBs.", sep=""))
  ids <- str_split_fixed(fls, "_",3)[,1]
  ids_unique <- unique(ids)
  print(paste("You have ", length(ids_unique), " unique individuals with successful BBs.", sep=""))
  
  for(i in 1:length(ids_unique)){
    if(i == round(length(ids_unique)/2,0))
      print("Your are half done.")
    if(length(ids[ids == ids_unique[i]])==1){
      bb <- raster(paste(BBs_fldr, "/", fls[ids == ids_unique[i]], sep=""))   # when there is just 1 individual
      bbtemp <- list("Brownian motion variance" = 0, "x" = coordinates(bb)[,1], "y" = coordinates(bb)[,2], "probability" = values(bb))
      qtl <- bbmm.contour(bbtemp, levels = contour, plot = FALSE)
      bb <- reclassify(bb, rcl=matrix(c(-1,qtl$Z,0),2,3, byrow=T))
      bb <- bb/sum(values(bb))   #verify that they add up to 1
    }else{
      fls2 <- fls[ids == ids_unique[i]]
      bb <- raster(paste(BBs_fldr, "/", fls2[1], sep=""))
      bbtemp <- list("Brownian motion variance" = 0, "x" = coordinates(bb)[,1], "y" = coordinates(bb)[,2], "probability" = values(bb))
      qtl <- bbmm.contour(bbtemp, levels = contour, plot = FALSE)
      bb <- reclassify(bb, rcl=matrix(c(-1,qtl$Z,0),2,3, byrow=T))
      bb <- bb/sum(values(bb))   #verify that they add up to 1
      for(e in 2:length(fls2)){
        bb2 <- raster(paste(BBs_fldr, "/", fls2[e], sep=""))
        bbtemp <- list("Brownian motion variance" = 0, "x" = coordinates(bb2)[,1], "y" = coordinates(bb2)[,2], "probability" = values(bb2))
        qtl <- bbmm.contour(bbtemp, levels = contour, plot = FALSE)
        bb2 <- reclassify(bb2, rcl=matrix(c(-1,qtl$Z,0),2,3, byrow=T))
        bb2 <- bb2/sum(values(bb2))   #verify that they add up to 1
        bb <- addLayer(bb, bb2)
      }
      if(nlayers(bb) != length(fls2))
        stop("You have a problem. See error 1.")
      bb <- mean(bb)
      bb <- bb/sum(values(bb))   #verify that they add up to 1
    }
    #output averaged individual ASCII file
    m <- as(bb, "SpatialGridDataFrame")
    
    # make sure the cell size is teh same in x and y direction. 
    cs <- slot(slot(m, "grid"), "cellsize") 
    slot(slot(m, "grid"), "cellsize") <- rep(mean(cs), 2) 
    write.asciigrid(m, paste(pop_BBs_out_fldr,"/",ids_unique[i],"_ASCII.asc",sep=""), attr=1)
    
    #99% contours
    bb <- list('Brownian motion variance'=0,x=coordinates(bb)[,1],y=coordinates(bb)[,2],probability=values(bb))
    contours <- bbmm.contour(bb, levels=contour, plot=F)
    # Create data.frame indicating cells within the each contour and export as Ascii Grid
    contour.99 <- data.frame(x = bb$x, y = bb$y, probability = bb$probability)
    # contour.99 <- contour.99[contour.99$probability >= contours$Z[1],]
    contour.99$in.out <- ifelse(contour.99$probability >= contours$Z[1], 1, 0)
    #write out footprint for individual
    m <- SpatialPixelsDataFrame(points = contour.99[c("x", "y")], data=contour.99)
    m <- as(m, "SpatialGridDataFrame")
    cs <- slot(slot(m, "grid"), "cellsize") 
    slot(slot(m, "grid"), "cellsize") <- rep(mean(cs), 2) 
    write.asciigrid(m, paste(pop_footprint_out_fldr,"/",ids_unique[i],"_99pct_contour.asc",sep=""), attr=ncol(m))
  }
  
  #create overall averages
  fls <- dir(pop_BBs_out_fldr)
  bb <- raster(paste(pop_BBs_out_fldr, "/", fls[1], sep=""))
  for(i in 2:length(fls)){
    bb <- addLayer(bb, raster(paste(pop_BBs_out_fldr, "/", fls[i], sep="")))
  }
  if(nlayers(bb) != length(fls))
    stop("You have a problem. See error 3.")
  bb <- mean(bb)
  bb <- bb/sum(values(bb))   #verify that they add up to 1
  #output averaged individual ASCII file
  m <- as(bb, "SpatialGridDataFrame")
  cs <- slot(slot(m, "grid"), "cellsize") 
  slot(slot(m, "grid"), "cellsize") <- rep(mean(cs), 2) 
  write.asciigrid(m, paste(pop_BBs_out_fldr,"/","averageUD.asc",sep=""), attr=1)
  projection(bb) <- proj_of_ascs
  writeRaster(bb, filename = paste(pop_BBs_out_fldr,"/","averageUD.img",sep=""), format="HFA")
  
  #create overall population footprint
  fls <- dir(pop_footprint_out_fldr)
  bb <- raster(paste(pop_footprint_out_fldr, "/", fls[1], sep=""))
  for(i in 2:length(fls)){
    bb <- addLayer(bb, raster(paste(pop_footprint_out_fldr, "/", fls[i], sep="")))
  }
  if(nlayers(bb) != length(fls))
    stop("You have a problem. See error 4.")
  bb <- sum(bb)
  #output averaged individual ASCII file
  m <- as(bb, "SpatialGridDataFrame")
  cs <- slot(slot(m, "grid"), "cellsize") 
  slot(slot(m, "grid"), "cellsize") <- rep(mean(cs), 2) 
  write.asciigrid(m, paste(pop_footprint_out_fldr,"/","popFootprint.asc",sep=""), attr=1)
  projection(bb) <- proj_of_ascs
  writeRaster(bb, filename = paste(pop_footprint_out_fldr,"/","popFootprint.img",sep=""), format="HFA")
  print(paste0("End time: ", Sys.time()))
  return("Done. Check your folders.")
} #end of function

