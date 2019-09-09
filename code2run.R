#-----------------------------------#
# Corridor and stopover analysis ####
#-----------------------------------#

# The following code finishes teh BB analysis of migration corridors after one has 
# completed tab 5 in Migration Mapper and clicked the export button on tab 6.

#source the functions you will need
source("C:/Users/jmerkle/Documents/GitHub/CorridorMappingTeam/functions/create.seqs.R")
source("C:/Users/jmerkle/Documents/GitHub/CorridorMappingTeam/functions/create.BBs.R")
source("C:/Users/jmerkle/Documents/GitHub/CorridorMappingTeam/functions/create.CTMMs.R")
source("C:/Users/jmerkle/Documents/GitHub/CorridorMappingTeam/functions/create.BB.avgs.R")
source("C:/Users/jmerkle/Documents/GitHub/CorridorMappingTeam/functions/create.corridors.stopovers.R")
source("C:/Users/jmerkle/Documents/GitHub/CorridorMappingTeam/functions/create.lns.file.R")

#step 1. Create sequences from the output of Migration Mapper (tab 6)
create.seqs(shpfl_fldr = "C:/Users/jmerkle/Desktop/tab6output",  #this is the folder where files from tab 6 are located
             shpfl_name= "pointsOut_20190207093941" ,  #name of the actual shapefile without the .shp 
             migtbl_name="migtime_20190207093941.csv",   #name of the actual csv file
             out_fldr="C:/Users/jmerkle/Desktop/tab6output/sequences",   #it will make this folder for you
             out_proj="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")   #name a projection you want the output to be in. Then carry this proj through the rest of the steps

# step 1b. Create a lines file from the sequences.
create.lns.file(seqs_fldr = "C:/Users/jmerkle/Desktop/tab6output/sequences",
                out_fldr="C:/Users/jmerkle/Desktop/tab6output/migration_lines",
                proj_of_dbfs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")  # this is the proj4string of your data. (should be carried through from previous functions)

#step 2. Conduct BB analyses. You will want use parrallel processing for this one. 
# This also spits out a metadata file of the results of the BB analysis
create.BBs(seqs_fldr = "C:/Users/jmerkle/Desktop/tab6output/sequences",    #this is the folder where all the sequences are saved
           BBs_out_fldr = "C:/Users/jmerkle/Desktop/tab6output/UDs",  #it will make this folder for you
           footprint_out_fldr = "C:/Users/jmerkle/Desktop/tab6output/Footprints",  #it will make this folder for you
           metadata_fldr="C:/Users/jmerkle/Desktop/tab6output",
           cores=11,    ## number of cores/threads to use for parrallel processing
           location.error=20,   ## location error of your GPS data (in meters)
           cell.size=50,    ## cell size of the raster you'd like to fit the BBs over 
           max.lag=8,     ## maximum amoung of time (in hours) to allow any two sequential points to be connected when fitting BBs
           contour=99,    ## Contour level used to create footprints for corridor analysis
           time.step=5,   ## how often (in minutes) that BB integrates between sequential points
           BMvar=NULL, ## if NULL, a motion variance will be calculated for each BB. If a number (in meters) is specified, that motion variance will be invoked for each BB.
           mult4buff=0.2, ## proportion of space around your gps data (i.e., a buffer) that is used to create the grid for conducting BBs
           proj_of_dbfs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")    # this is the proj4string of your data. (should be carried through from previous functions)


# step 2 ALTERNATIVE. (This code does the ctmm alternative method). You'll want to parralellize this one, fo sho.
# This also spits out a metadata file of the results of the CTMM analysis
# devtools::install_github("ctmm-initiative/ctmm")   # need to download specific packages for what we're doing
create.CTMMs(seqs_fldr = "C:/Users/jmerkle/Desktop/tab6output/sequences",  #this is the folder where all the sequences are saved
             CTMMs_out_fldr = "C:/Users/jmerkle/Desktop/tab6output/UDs2", #it will make this folder for you
             metadata_fldr="C:/Users/jmerkle/Desktop/tab6output",
             cores=11,   ## number of cores/threads to use for parrallel processing
             cell.size=50,  ## cell size of the raster you'd like to fit the BBs over 
             mult4buff=0.2, ## proportion of space around your gps data (i.e., a buffer) that is used to create the grid for conducting BBs
             Information_Criteria="LOOCV", ## information criteria used for model selection. Can be "AICc", "AIC", "BIC", "LOOCV" or none (NA). Use LOOCV for 12 hour data. Use something else like AIC for more frequent fix intervals.
             proj_of_dbfs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")


#Step 3. Calculate average BBs for each individual and and a population UD and population footprint.
create.BB.avgs(BBs_fldr = "C:/Users/jmerkle/Desktop/tab6output/UDs",    #this is the folder where all the UDs are saved.
               pop_BBs_out_fldr = "C:/Users/jmerkle/Desktop/tab6output/UDs_pop",  #it will make this folder for you
               pop_footprint_out_fldr = "C:/Users/jmerkle/Desktop/tab6output/Footprints_pop",  #it will make this folder for you
               contour=99,  ## Contour level used to create footprints for corridor analysis
               proj_of_ascs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")  # this is the proj4string of your data. (should be carried through from previous functions)

# Step 4. Calculate the stopover files and the low, medium, high use corridors as shapefiles
create.corridors.stopovers(PopUD_asc = "C:/Users/jmerkle/Desktop/tab6output/UDs_pop/averageUD.asc",   #this is the file path for the POPud ascii file
                           PopFootprint_asc = "C:/Users/jmerkle/Desktop/tab6output/Footprints_pop/popFootprint.asc",  #this is the file path for the POfootprint ascii file
                           pop_BBs_fldr = "C:/Users/jmerkle/Desktop/tab6output/Footprints_pop",   #this is the pop BBs folder
                           out_fldr = "C:/Users/jmerkle/Desktop/tab6output/final_products",     #this is an empty folder where you want the results to be saved
                           stopover_percent=10, ## the quantile of the UD values for which to call stopovers (larger values results in larger stopover areas)
                           corridor_percents=c(10, 20),  ## corridor percents that should be calculated (these are in addition to 1 or more, and 2 or more individuals)
                           min_area = 20000, ## if there are unique polygons smaller than this value (in squared meters), they will be removed
                           simplify = TRUE, ## if TRUE, polygons will be simplified (i.e., the square edges of the grid will be smoothed)
                           tolerance = 25, ## tolerance value when simplifying polygons (unites are meters)
                           proj_of_ascs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")   # this is the proj4string of your data. (should be carried through from previous functions)

#--------------------------#
# Winter range analysis ####
#--------------------------#

# The following code finishes teh BB analysis of winter ranges after one has 
# completed tab 5 in Migration Mapper and clicked the export button on tab 6.

#source the functions you will need
source("C:/Users/jmerkle/Documents/GitHub/CorridorMappingTeam/functions/create.seqs.W.R")
source("C:/Users/jmerkle/Documents/GitHub/CorridorMappingTeam/functions/create.BBs.W.R")
source("C:/Users/jmerkle/Documents/GitHub/CorridorMappingTeam/functions/create.BB.avgs.W.R")
source("C:/Users/jmerkle/Documents/GitHub/CorridorMappingTeam/functions/create.core.areas.W.R")

#step 1. Create winter sequences from the output of Migration Mapper (tab 6)
create.seqs.W(shpfl_fldr = "C:/Users/jmerkle/Desktop/tab6output", 
              shpfl_name= "pointsOut_20190207093941" ,
              idname="newUid",  #name of the column representing animal ID
              datename="nwMstrD",   #name of the column representing date in POSIX format
              mig.metadata.file="C:/Users/jmerkle/Desktop/tab6output/metadata_migration.csv",  # metadata file from migration part of analysis
              use_same_annual_dates=TRUE, ## if TRUE, then will identify the same winter period for each year based on the qtl.end.fall.mig and qtl.start.spr.mig. If FALSE, winter is defined for each id-yr as the end of fall migration to start of spring migration
              qtl.end.fall.mig=0.95, ## quantile of end of fall migration dates, which serves to start the winter period for each year
              qtl.start.spring.mig=0.05,   ## quantile of start of spring migration dates, which serves to end the winter period for each year
              out_fldr="C:/Users/jmerkle/Desktop/tab6output/sequencesW",   #this is where you want to sequences saved to
              out_proj="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")   #name a projection you want the output to be in. Then carry this proj through the rest of the steps


#step 2. Conduct BB analyses. You will want use parrallel processing for this one. 
# This also spits out a metadata file of the results of the BB analysis
create.BBs.W(seqs_fldr = "C:/Users/jmerkle/Desktop/tab6output/sequencesW",    #this is the folder where all the sequences are saved
             BBs_out_fldr = "C:/Users/jmerkle/Desktop/tab6output/UDsW",  #it will make this folder for you
             metadata_fldr="C:/Users/jmerkle/Desktop/tab6output",
             mindays=30,   ## if an individual animal has less than this many days in a given sequence of winter data, it will be removed
             cores=11,    ## number of cores/threads to use for parrallel processing
             location.error=20,   ## location error of your GPS data (in meters)
             cell.size=50,    ## cell size of the raster you'd like to fit the BBs over 
             max.lag=8,     ## maximum amoung of time (in hours) to allow any two sequential points to be connected when fitting BBs
             time.step=5,   ## how often (in minutes) that BB integrates between sequential points
             BMvar=NULL, ## if NULL, a motion variance will be calculated for each BB. If a number (in meters) is specified, that motion variance will be invoked for each BB.
             mult4buff=0.2, ## proportion of space around your gps data (i.e., a buffer) that is used to create the grid for conducting BBs
             proj_of_dbfs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")    # this is the proj4string of your data. (should be carried through from previous functions)


#Step 3. Calculate average BBs for each individual and and a population UD.
create.BB.avgs.W(BBs_fldr = "C:/Users/jmerkle/Desktop/tab6output/UDsW",    #this is the folder where all the UDs are saved.
                 pop_BBs_out_fldr = "C:/Users/jmerkle/Desktop/tab6output/UDs_popW",  #it will make this folder for you
                 proj_of_ascs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")    # this is the proj4string of your data. (should be carried through from previous functions)


# Step 4. Calculate the winter contours
create.core.areas.W(PopUD_asc = "C:/Users/jmerkle/Desktop/tab6output/UDs_popW/averageUD_winter.asc",   #this is the file path of the averageUD_winter.asc
                    out_fldr = "C:/Users/jmerkle/Desktop/tab6output/final_productsW",   #this is an empty folder where you want the final products to go
                    core_contours=c(.01,0.10,0.20,0.30,0.40,0.50,0.60,0.70,0.80,0.90,0.99,1),   # You must include a 1. Function will spit out polygons for each of these contours (0.01 = the top 1% of the values/area; 0.99 = the top 99% of the values/area)
                    min_area = 20000, ## if there are unique polygons smaller than this value (in squared meters), they will be removed
                    simplify = TRUE, ## if TRUE, polygons will be simplified (i.e., the square edges of the grid will be smoothed)
                    tolerance = 25, ## tolerance value when simplifying polygons (unites are meters)
                    proj_of_ascs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")    # this is the proj4string of your data. (should be carried through from previous functions)






# These are Extra bits of code

# Method for teh hack corridors (based solely on buffers around the lines)
source("C:/Users/jmerkle/Documents/GitHub/CorridorMappingTeam/functions/hack.corridors.R")
hack.corridors(lines_fldr = "C:/Users/jmerkle/Desktop/tab6output/tab6output/migration_lines",  #this is the folder where your migration lines shapefile is stored
               lines_name = "migration_lines",   # This is the name of the migration lines shapefile (do not include .shp)
               out_fldr="C:/Users/jmerkle/Desktop/tab6output/tab6output/hack_corridors",   #an empty folder where you want the outputs to be stored
               corridor_percents=c(10, 20),  #the corridor percents that are provided (these are in addition to 1 or more, and 2 or more corridors)
               min_area = 20000, #if there are polygons smaller than this (in squared meters), they will be removed
               simplify = TRUE, #should polygons be simplified?
               tolerance = 25, # how to polygons are simplified (unites are meters)
               buff=200, # how many meters do you want to buffer the lines?
               cell.size=50) #this is the cell size of the output (should be 50m)
