# This function takes in individual BBs and creates averaged individual BBs and footprints
# and population UDs and population footprints
# written by Jerod Merkle, 8 Jan 2019
# based on code from Hall Sawyer

create.BB.avgs.W <- function(BBs_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/UDsW",
                             pop_BBs_out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/UD_popW",
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
  
  print(paste0("Start time: ", Sys.time()))
  
  fls <- list.files(BBs_fldr, ".asc$")
  if(length(fls)==0)
    stop("You have no sequences with successful BBs!")
  print(paste("You have ", length(fls), " sequences with successful BBs.", sep=""))
  ids <- str_split_fixed(fls, "_",3)[,1]
  ids_unique <- unique(ids)
  print(paste("You have ", length(ids_unique), " unique individuals with successful BBs.", sep=""))
  
  for(i in 1:length(ids_unique)){
    if(i == round(length(ids_unique)/2,0))
      print("Your are half done.")
    if(length(ids[ids == ids_unique[i]])==1){
      bb <- raster(paste(BBs_fldr, "/", fls[ids == ids_unique[i]], sep=""))   # when there is just 1 individual
    }else{
      fls2 <- fls[ids == ids_unique[i]]
      bb <- raster(paste(BBs_fldr, "/", fls2[1], sep=""))
      for(e in 2:length(fls2)){
        bb <- addLayer(bb, raster(paste(BBs_fldr, "/", fls2[e], sep="")))
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
    
  }
  
  #create overall averages
  fls <- list.files(pop_BBs_out_fldr, ".asc$")

  if(length(fls)==1){
    bb <- raster(paste(pop_BBs_out_fldr, "/", fls[1], sep=""))
    m <- as(bb, "SpatialGridDataFrame")
    cs <- slot(slot(m, "grid"), "cellsize") 
    slot(slot(m, "grid"), "cellsize") <- rep(mean(cs), 2) 
    write.asciigrid(m, paste(pop_BBs_out_fldr,"/","averageUD_winter.asc",sep=""), attr=1)
    projection(bb) <- proj_of_ascs
    writeRaster(bb, filename = paste(pop_BBs_out_fldr,"/","averageUD_winter.img",sep=""), format="HFA")
    print(paste0("End time: ", Sys.time()))
    return("Done. Check your folders. Note that you only have 1 individual, so no averaging occured!")
  }else{
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
    write.asciigrid(m, paste(pop_BBs_out_fldr,"/","averageUD_winter.asc",sep=""), attr=1)
    projection(bb) <- proj_of_ascs
    writeRaster(bb, filename = paste(pop_BBs_out_fldr,"/","averageUD_winter.img",sep=""), format="HFA")
    print(paste0("End time: ", Sys.time()))
    return("Done. Check your folders.")
  }
  
  
} #end of function

