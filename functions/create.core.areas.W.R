# This function creates final stopover and corridor files
# written by Jerod Merkle, 8 Jan 2019
# based on code from Hall Sawyer

create.core.areas.W <- function(PopUD_asc = "C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6/UDs_popW/averageUD_winter.asc",
                                       out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6/final_productsW",
                                       core_contours=c(10,20,30,40,50,60,70,80,90), min_area = 20000,
                                       simplify = TRUE, tolerance = 25, # unites are meters
                                       proj_of_ascs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"){
  
  #manage packages
  if(all(c("raster","move","BBMM","igraph","rgdal","rgeos","maptools") %in% installed.packages()[,1])==FALSE)
    stop("You must install the following packages: raster,igraph, move, rgeos, BBMM, maptools, and rgdal.")
  require(raster)
  require(BBMM)
  require(rgdal)
  require(igraph)
  require(rgeos)
  require(maptools)
  require(move)
  
  #check the new directories
  if(dir.exists(out_fldr)==FALSE){
    dir.create(out_fldr)
  }
  if(length(dir(out_fldr))> 0)
    stop("Your out_fldr Has something in it. It should be empty!")
  
  #these are the functions needed:
  
  # NOTE: The simplification is set to preserve topology; although not 
  # strictly necessary for visualizaiton, this avoids the scenario where
  # large polygons were mysteriously getting excluded from the output.
  processPolygons<-function(){
    # each polygon contains only the area within a single cut
    # join the polygons with any in higher cuts (cumulative area for each level)
    n=length(polygonList@polygons)
    
    
    polygonList.acc <<- lapply(1:n, 
                               function(x){
                                 ##################################
                                 # --- Updated January 2016
                                 # accommodate the layers being out of order
                                 idList = rep(NA_integer_,n)
                                 idList[which(polygonList@data$layer>=x)]=1
                                 #print(idList)
                                 #print(polygonList@data$layer[which(idList==1)])
                                 
                                 unionSpatialPolygons(
                                   SpatialPolygons(polygonList@polygons),idList)
                               })
    
    # line simplification
    if(simplify){
      polygonList.acc <<- lapply(polygonList.acc,thinnedSpatialPoly,
                                 tolerance,topologyPreserve=T)
    }
    
    # split up multi-polygons into separate polygons
    uniqueID<<-0 # reset before calling makeIndividual polygons
    
    spDataFrame.list <<- sapply(1:length(polygonList.acc),
                                function(x) makeIndividualPolygons(
                                  polygonList.acc[[x]]@polygons[[1]],core_contours[x],min_area))
    #sapply((sapply(spDataFrame.list,function(x) sapply(x@polygons,function(y) sapply(y@Polygons, function(z) z@hole)))),sum)
  }
  
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
    
    attrTable = data.frame(GRIDCODE=attr,area=diffArea)
    row.names(attrTable) = id
    
    x=SpatialPolygonsDataFrame(SpatialPolygons(polygon.list.wrap), 
                               data=attrTable)
    
    #print(length(sapply(sapply(x@polygons,function(y) sapply(y@Polygons, function(z) z@hole)),sum)))
    #print(sum(sapply(sapply(x@polygons,function(y) sapply(y@Polygons, function(z) z@hole)),sum)))
    x
  }
  
  
  # OK, now bring in the popUD_asc for teh winter
  wintUDs <- raster(PopUD_asc)
  projection(wintUDs) <- proj_of_ascs
  # wintUDs[wintUDs == 0] <- NA
  
  #identify the contour
  breakValues <- raster::quantile(wintUDs[wintUDs != 0], probs = core_contours)
  #*** FLIP 'thresholdQuantiles' SO THE POLYGONS PERTAIN TO THE AREA HS EXPECTS - STARTING AT 2ND ELEMENT AS CODE APPEARS TO DISMISS LAST ELEMENT
  core_contours <- core_contours[(length(core_contours) - 1):1]
  
  # compute the contours
  classifiedRaster = cut(wintUDs,breaks=breakValues)
  
  # extract the contours as polygons
  polygonList = rasterToPolygons(classifiedRaster,dissolve=T)
  
  # transform the multipolygons into simplified individual polygons
  # accounting for the holes and filtering as needed
  processPolygons()
  
  # write all of the output to a single shapefile
  # variables to hold the accumulated data
  bigData.df = data.frame(contour=numeric(),area=numeric())
  bigData.sp = list()
  # accumulate the data
  for(i in 1:length(spDataFrame.list)){
    bigData.df = rbind(bigData.df,spDataFrame.list[[i]]@data)
    bigData.sp = append(bigData.sp,spDataFrame.list[[i]]@polygons)
  }
  # create the spatial data frame with a projection (same as above)
  bigData = SpatialPolygonsDataFrame(SpatialPolygons(bigData.sp),bigData.df)
  proj4string(bigData) <- proj_of_ascs
  
  bigData$area <- NULL   #remove the area column

  # bigData$GRIDCODE <- as.numeric(as.character(bigData$GRIDCODE))
  # bigData <- aggregate(bigData, "GRIDCODE")   #dissolve the polygons based on the gridcode
  # bigData <- bigData[order(bigData$GRIDCODE, decreasing = FALSE),]
  
  # write to a shapefile
  writeOGR(bigData,out_fldr,"winter_ranges", driver="ESRI Shapefile")

  toreport <- data.frame(results="yes")
  grdcds <- unique(bigData$GRIDCODE)
  for(i in 1:length(grdcds)){
    toreport$temp <- round(gArea(bigData[bigData$GRIDCODE <= grdcds[i],], byid = FALSE)/1000000,1)
    names(toreport)[names(toreport)=="temp"] <- paste0("GRIDCODE_",grdcds[i],"_area")
    print(paste0("The area of GRIDCODE ", grdcds[i], " is ", round(gArea(bigData[bigData$GRIDCODE <= grdcds[i],], byid = FALSE)/1000000,1), " square KM."))
  }
  
  write.csv(toreport, file = paste0(out_fldr,"/output_areas_winter.csv"), row.names = FALSE)
  
  return("Success. Check your folders.")
}