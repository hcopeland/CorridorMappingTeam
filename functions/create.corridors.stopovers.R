# This function creates final stopover and corridor files
# written by Jerod Merkle, 8 Jan 2019
# based on code from Hall Sawyer

create.corridors.stopovers <- function(PopUD_asc = "C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6/UDs_pop/averageUD.asc",
                                       PopFootprint_asc = "C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6/Footprints_pop/popFootprint.asc",
                                       pop_BBs_fldr = "C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6/Footprints_pop",
                                       out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6/final_products",
                                       stopover_percent=10, corridor_percents=c(10, 20),min_area = 20000,
                                       simplify = TRUE, tolerance = 25, # unites are meters
                                       proj_of_ascs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"){
  
  #manage packages
  if(all(c("raster","BBMM","igraph","rgdal","rgeos","maptools") %in% installed.packages()[,1])==FALSE)
    stop("You must install the following packages: raster,igraph, rgeos, BBMM, maptools, and rgdal.")
  require(raster)
  require(BBMM)
  require(rgdal)
  require(igraph)
  require(rgeos)
  require(maptools)
  
  #check the new directories
  if(dir.exists(out_fldr)==FALSE){
    dir.create(out_fldr)
  }
  if(length(dir(out_fldr))> 0)
    stop("Your out_fldr Has something in it. It should be empty!")
  
  stopovers <- raster(PopUD_asc)
  projection(stopovers) <- proj_of_ascs
  
  bb <- list("Brownian motion variance" = 0, "x" = coordinates(stopovers)[,1], "y" = coordinates(stopovers)[,2], "probability" = values(stopovers))
  qtl <- bbmm.contour(bb, levels = stopover_percent, plot = FALSE)
  stopovers <- reclassify(stopovers, rcl=matrix(c(-1,qtl$Z,NA,qtl$Z,Inf,1),2,3, byrow=T))

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
  stopoverarea <- gArea(stopovers)/1000000
  if(simplify==TRUE){
    stopovers<- thinnedSpatialPoly(stopovers, tolerance = tolerance, topologyPreserve=TRUE)
  }
  writeOGR(stopovers, out_fldr, "stopovers", driver="ESRI Shapefile")
  
  #-------------------#
  # now for corridors #
  #-------------------#
  popFootprint <- raster(PopFootprint_asc)
  projection(popFootprint) <- proj_of_ascs
  
  # for the other corridors based on % of individuals collared
  numb_ids <- length(list.files(pop_BBs_fldr, pattern=".asc"))-1
  print(paste0("The number of IDs in your populatoin is ", numb_ids, "."))
  maxval <- max(values(popFootprint))
  print(paste0("The max value of your population footprint grid is ",maxval,"."))
  
  popFootprint <- popFootprint/numb_ids
  
  corridor_percents <- corridor_percents[(corridor_percents/100) > (2/numb_ids)]
  
  thresholdQuantiles = c(0, c(2/numb_ids)-0.0001, (corridor_percents/100)-0.0001, .99)  #this starts with low 1 or more, then low 2 or more, then the other percents
  thresholdQuantiles_names <- c(1,2,corridor_percents)
  
  print(paste0("GRIDCODEs of ", paste(thresholdQuantiles_names,collapse = ", "), 
               " represent the following ranges of individuals: 0 - ", 
               paste(floor(thresholdQuantiles*numb_ids)+1,collapse=" - "),"."))
  
  # compute the contours
  classifiedRaster = cut(popFootprint,breaks=thresholdQuantiles)
  
  # extract the contours as polygons
  polygonList = rasterToPolygons(classifiedRaster,dissolve=T)
  
  # Methods - Data Munging -------------------------------------------
  # reads and updates the uniqueID global variable to assign 
  # identifiers to each polygon (unique across shapefiles).
  makeIndividualPolygons<-function(multiPoly,attr,minArea=0){
    #polygon.list = polygonList.acc[[1]]@polygons[[1]]@Polygons
    polygon.list = multiPoly@Polygons
    area = sapply(polygon.list,function(x) x@area)
    isHole = sapply(polygon.list,function(x) x@hole)
    area[which(isHole)]=-1*area[which(isHole)]
    
    holeIndices = which(!isHole)
    
    diffArea=numeric()
    startID=uniqueID+1
    if(length(holeIndices)>1){
      polygon.list.wrap = lapply(1:(length(holeIndices)-1),function(x) {
        uniqueID<<-uniqueID+1
        temp = holeIndices[x]:(holeIndices[x+1]-1)
        diffArea[x]<<-sum(area[temp])
        Polygons(polygon.list[temp],uniqueID)
      })
    }else{
      polygon.list.wrap = list()
    }
    
    uniqueID<<-uniqueID+1
    temp = holeIndices[length(holeIndices)]:length(polygon.list)
    polygon.list.wrap[[length(polygon.list.wrap)+1]] = Polygons(polygon.list[temp],uniqueID)
    diffArea[length(holeIndices)]=sum(area[temp])
    endID=uniqueID
    id=startID:endID
    
    temp = which(diffArea>minArea)
    diffArea=diffArea[temp]
    polygon.list.wrap=polygon.list.wrap[temp]  
    id=id[temp]
    
    #sum(sapply(1:length(polygon.list.wrap), function(x) sapply(polygon.list.wrap[[x]]@Polygons, function(y) y@hole)))
    #sapply(polygon.list.wrap[[2]]@Polygons,function(x) x@hole)
    #print(sapply(1:length(index), function(x) sapply(polygon.list.wrap[[x]]@Polygons, function(y) y@hole=isHole[x])))
    
    attrTable = data.frame(ID=id, GRIDCODE=attr)
    row.names(attrTable) = id
    
    x=SpatialPolygonsDataFrame(SpatialPolygons(polygon.list.wrap), data=attrTable)
    
    #print(length(sapply(sapply(x@polygons,function(y) sapply(y@Polygons, function(z) z@hole)),sum)))
    #print(sum(sapply(sapply(x@polygons,function(y) sapply(y@Polygons, function(z) z@hole)),sum)))
    x
  }
  
  
  
  # NOTE: The simplification is set to preserve topology; although not 
  # strictly necessary for visualizaiton, this avoids the scenario where
  # large polygons were mysteriously getting excluded from the output.
  # each polygon contains only the area within a single cut
  # join the polygons with any in higher cuts (cumulative area for each level)
  n=length(polygonList@polygons)
  
  polygonList.acc <- lapply(1:n, function(x){
    idList = c(rep(NA,x-1),rep(1,n-x+1))
    #print(idList)
    unionSpatialPolygons(SpatialPolygons(polygonList@polygons),idList)
  })
  
  # line simplification
  if(simplify==TRUE){
    polygonList.acc <- lapply(polygonList.acc,thinnedSpatialPoly,
                              tolerance,topologyPreserve=T)
  }
  
  # split up multi-polygons into separate polygons
  uniqueID<-0 # reset before calling makeIndividual polygons
  spDataFrame.list <- sapply(1:length(polygonList.acc),function(x) makeIndividualPolygons(
    polygonList.acc[[x]]@polygons[[1]],thresholdQuantiles_names[x],min_area))
  #sapply((sapply(spDataFrame.list,function(x) sapply(x@polygons,function(y) sapply(y@Polygons, function(z) z@hole)))),sum)
  
  for(x in 1:length(spDataFrame.list)){
    proj4string(spDataFrame.list[[x]])<-proj_of_ascs
  }
  
  bigData.df = data.frame(ID=numeric(), GRIDCODE=numeric())
  bigData.sp = list()
  # accumulate the data
  for(i in 1:length(spDataFrame.list)){
    bigData.df = rbind(bigData.df,spDataFrame.list[[i]]@data)
    bigData.sp = append(bigData.sp,spDataFrame.list[[i]]@polygons)
  }
  # create the spatial data frame with a projection (same as above)
  bigData = SpatialPolygonsDataFrame(SpatialPolygons(bigData.sp),bigData.df)
  proj4string(bigData) <- proj_of_ascs
  
  bigData <- aggregate(bigData, "GRIDCODE")   #dissolve the polygons based on the gridcode
  bigData <- bigData[order(bigData$GRIDCODE),]
  
  # write to a shapefile
  writeOGR(bigData,out_fldr,"corridors", driver="ESRI Shapefile")
  print(paste0("The area of your stopovers is ", round(stopoverarea,1), " square KM."))
  
  toreport <- data.frame(numb_ids=numb_ids, max_value_footprint_grid=maxval, stopover_area=round(stopoverarea,1))
  
  grdcds <- unique(bigData$GRIDCODE)
  for(i in 1:length(grdcds)){
    toreport$temp <- round(gArea(bigData[bigData$GRIDCODE >= grdcds[i],], byid = FALSE)/1000000,1)
    names(toreport)[names(toreport)=="temp"] <- paste0("GRIDCODE_",grdcds[i],"_area")
    print(paste0("The area of GRIDCODE ", grdcds[i], " is ", round(gArea(bigData[bigData$GRIDCODE >= grdcds[i],], byid = FALSE)/1000000,1), " square KM."))
  }
  write.csv(toreport, file = paste0(out_fldr,"/output_areas.csv"), row.names = FALSE)
   
  print("Success. Check your folders")
  return("Done.")
}