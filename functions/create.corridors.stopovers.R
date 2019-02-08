# This function creates final stopover and corridor files
# written by Jerod Merkle, 8 Jan 2019
# based on code from Hall Sawyer

create.corridors.stopovers <- function(PopUD_asc = "C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6/UDs_pop/averageUD.asc",
                                       PopFootprint_asc = "C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6/Footprints_pop/popFootprint.asc",
                                       pop_BBs_fldr = "C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6/Footprints_pop",
                                       out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6/final_products",
                                       stopover_percent=10, corridor_percents=c(10, 20),min_area = 20000,
                                       proj_of_ascs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"){
  
  #manage packages
  if(all(c("raster","move","igraph","rgeos") %in% installed.packages()[,1])==FALSE)
    stop("You must install the following packages: raster,igraph, rgeos, move")
  require(raster)
  require(move)
  require(igraph)
  require(rgeos)

  #check the new directories
  if(dir.exists(out_fldr)==FALSE){
    dir.create(out_fldr)
  }
  if(length(dir(out_fldr))> 0)
    stop("Your out_fldr Has something in it. It should be empty!")

  popUD <- raster(PopUD_asc)
  projection(popUD) <- proj_of_ascs
  
  popUD <- getVolumeUD(as(popUD, Class="DBBMM"))
  qtl <- quantile(values(popUD)[values(popUD)!=1], probs=stopover_percent/100)
  stopovers <- reclassify(popUD, rcl=matrix(c(0,qtl,1,qtl,Inf,NA),2,3, byrow=T))
  plot(stopovers)
  
  #remove patches that are smaller than min_area
  clmps <- clump(stopovers)
  clmpvals <- na.omit(values(clmps))
  clmpvals <- data.frame(table(clmpvals))
  clmpvals$Freq <- clmpvals$Freq*res(stopovers)[1]*res(stopovers)[2]   #this turns Freq into actuall m^2
  clmpvals <- as.numeric(as.character(clmpvals$clmpvals)[clmpvals$Freq < min_area])
  clmps[is.na(values(clmps))==TRUE] <- -5
  stopovers[values(clmps) %in% clmpvals == TRUE] <- NA
  
  #convert final raster to polygons
  stopovers <- rasterToPolygons(stopovers, dissolve=TRUE)
  writeOGR(stopovers, out_fldr, "stopovers", driver="ESRI Shapefile")
  
  # now for corridors
  popFootprint <- raster(PopFootprint_asc)
  projection(popFootprint) <- proj_of_ascs
  
  # for the low use corridors
  lowusecorridor <- reclassify(popFootprint, rcl=matrix(c(0,0.000001,0,0.000001,Inf,1),2,3, byrow=T))
  lowusecorridor[values(lowusecorridor)==0] <- NA
  lowusecorridor <- rasterToPolygons(lowusecorridor, dissolve=TRUE)
  writeOGR(lowusecorridor, out_fldr, "corridor_lowuse1ormore", driver="ESRI Shapefile")
  
  lowuse2ormore <- reclassify(popFootprint, rcl=matrix(c(0,1.9,0,1.9,Inf,1),2,3, byrow=T))
  lowuse2ormore [values(lowuse2ormore)==0] <- NA
  lowuse2ormore <- rasterToPolygons(lowuse2ormore, dissolve=TRUE)
  writeOGR(lowuse2ormore, out_fldr, "corridor_lowuse2ormore", driver="ESRI Shapefile")
  
  # for the other corridors based on % of individuals collared
  numb_ids <- length(list.files(pop_BBs_fldr, pattern=".asc"))-1
  print(paste("The number of ids in your populatoin is ", numb_ids, ". If that doesn't look right, then do some investigation.", sep=""))
  popFootprint <- popFootprint/numb_ids
  
  for(i in 1:length(corridor_percents)){
    medusecorridor <- reclassify(popFootprint, rcl=matrix(c(0,corridor_percents[i]/100,0,corridor_percents[i]/100,1,1),2,3, byrow=T))
    medusecorridor[values(medusecorridor)==0] <- NA
    medusecorridor <- rasterToPolygons(medusecorridor, dissolve=TRUE)
    writeOGR(medusecorridor, out_fldr, paste("corridor_percent",corridor_percents[i],sep=""), driver="ESRI Shapefile")
  }
  print("Success. Check your folders")
  return("Done.")
}