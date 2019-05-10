#code to create a spatial lines data frame from a sequences folder
# by Jerod Merkle, 25 Feb 2019

hack.corridors <- function(lines_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/migration_lines",
                           lines_name = "migration_lines",
                            out_fldr="C:/Users/jmerkle/Desktop/Mapp2/tab6output/hack_corridors",
                           corridor_percents=c(10, 20),min_area = 20000,
                           simplify = TRUE, tolerance = 25, # unites are meters
                            buff=200, cell.size=50){   
  if(all(c("rgdal","maptools","stringr","rgeos","raster") %in% installed.packages()[,1])==FALSE)
    stop("You must install the following packages: rgdal, maptools,rgeos, raster, and stringr")
  require(rgdal)
  require(maptools)
  require(stringr)
  require(rgeos)
  require(raster)
  
  #check the new directories
  if(dir.exists(out_fldr)==FALSE){
    dir.create(out_fldr)
  }
  if(length(dir(out_fldr))> 0)
    stop("Your out_fldr Has something in it. It should be empty!")
  
  print("Loading lines file")
  
  lns <- readOGR(lines_fldr, lines_name)
  
  if("AID" %in% names(lns)){
    lns$id <- str_split_fixed(lns$AID, "_", 2)[,1]
  }
  
  ids <- lns$id
  ids_unique <- unique(ids)
  
  print("Buffering lines")
  
  id_buffs <- do.call(rbind, lapply(1:length(ids_unique), function(i){
    tmp <- lns[lns$id %in% ids_unique[i],]
    tmpb <- gBuffer(tmp, byid=TRUE, width=buff)
    tmpb <- gUnaryUnion(tmpb, id=NULL)
    tmpb <- SpatialPolygonsDataFrame(tmpb, data=data.frame(id=ids_unique[i]))
  }))
  nrow(id_buffs)

  ext <- extent(id_buffs)
  popFootprint <- raster(ext)
  res(popFootprint) <- cell.size     
  projection(popFootprint) <- proj4string(id_buffs)
  
  print("Rasterizing buffered lines (this can take a while)")
  
  popFootprint <- rasterize(id_buffs[1,], popFootprint, background=0)
  for(i in 2:nrow(id_buffs)){
    popFootprint <- popFootprint+rasterize(id_buffs[i,], popFootprint, background=0)
  }
  # plot(popFootprint)
  writeRaster(popFootprint, filename=paste0(out_fldr, "/popFootprint.img"), format="HFA")
  numb_ids <- length(ids_unique)
  popFootprint <- popFootprint/numb_ids
  # plot(popFootprint)
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
    proj4string(spDataFrame.list[[x]])<-proj4string(id_buffs)
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
  proj4string(bigData) <- proj4string(id_buffs)
  
  bigData <- aggregate(bigData, "GRIDCODE")   #dissolve the polygons based on the gridcode
  bigData <- bigData[order(bigData$GRIDCODE),]
  
  # write to a shapefile
  writeOGR(bigData,out_fldr,"corridors", driver="ESRI Shapefile")
  return("Finished. Check your out folder.")
}
