# This function creates final stopover and corridor files
# written by Jerod Merkle, 8 Jan 2019
# based on code from Hall Sawyer

create.core.areas.W <- function(PopUD_asc = "C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6/UDs_popW/averageUD_winter.asc",
                                       out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6/final_productsW",
                                       core_contours=c(10,20,30,40,50,60,70,80,90,99), min_area = 20000,
                                       simplify = TRUE, tolerance = 100, # unites are meters
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
  
  wintUDs <- raster(PopUD_asc)
  projection(wintUDs) <- proj_of_ascs
  
  bb <- list("Brownian motion variance" = 0, "x" = coordinates(wintUDs)[,1], "y" = coordinates(wintUDs)[,2], "probability" = values(wintUDs))
  qtl <- bbmm.contour(bb, levels = core_contours, plot = FALSE)
  
  thresholdQuantiles = c(0, qtl$Z)  
  thresholdQuantiles_names <- sub("%","",qtl$Contour)
  
  if(length(thresholdQuantiles)!= length(unique(thresholdQuantiles))){
    thresholdQuantiles_names <- thresholdQuantiles_names[duplicated(thresholdQuantiles)[2:length(thresholdQuantiles)]==FALSE]
    thresholdQuantiles <- thresholdQuantiles[duplicated(thresholdQuantiles)==FALSE]
  }
  
  # compute the contours
  classifiedRaster = cut(wintUDs,breaks=thresholdQuantiles)
  
  
  #remove patches that are smaller than min_area
  clmps <- clump(classifiedRaster)
  clmpvals <- na.omit(values(clmps))
  clmpvals <- data.frame(table(clmpvals))
  clmpvals$Freq <- clmpvals$Freq*res(classifiedRaster)[1]*res(classifiedRaster)[2]   #this turns Freq into actuall m^2
  clmpvals <- as.numeric(as.character(clmpvals$clmpvals)[clmpvals$Freq < min_area])
  clmps[is.na(values(clmps))==TRUE] <- -5
  classifiedRaster[values(classifiedRaster) %in% clmpvals == TRUE] <- NA
  
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
  
  bigData$GRIDCODE <- as.numeric(as.character(bigData$GRIDCODE))
  
  # write to a shapefile
  writeOGR(bigData,out_fldr,"winter_ranges", driver="ESRI Shapefile")

  grdcds <- unique(bigData$GRIDCODE)
  for(i in 1:length(grdcds)){
    print(paste0("The area of GRIDCODE ", grdcds[i], " is ", round(gArea(bigData[bigData$GRIDCODE >= grdcds[i],], byid = FALSE)/1000000,1), " square KM."))
  }
  return("Success. Check your folders.")
}