#-----------------------------------#
# Corridor and stopover analysis ####
#-----------------------------------#

# The following code finishes teh BB analysis of migration corridors after one has 
# completed tab 5 in Migration Mapper and clicked the export button on tab 6.

#source the functions you will need
source("C:/Users/jmerkle/Documents/GitHub/CorridorMappingTeam/functions/create.seqs.R")
source("C:/Users/jmerkle/Documents/GitHub/CorridorMappingTeam/functions/create.BBs.R")
source("C:/Users/jmerkle/Documents/GitHub/CorridorMappingTeam/functions/create.BB.avgs.R")
source("C:/Users/jmerkle/Documents/GitHub/CorridorMappingTeam/functions/create.corridors.stopovers.R")
source("C:/Users/jmerkle/Documents/GitHub/CorridorMappingTeam/functions/create.lns.file.R")

#step 1. Create sequences from the output of Migration Mapper (tab 6)
create.seqs(shpfl_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output",  #this is the folder where files from tab 6 are located
             shpfl_name= "pointsOut_20190207093941" ,  #name of the actual shapefile without the .shp 
             migtbl_name="migtime_20190207093941.csv",   #name of the actual csv file
             out_fldr="C:/Users/jmerkle/Desktop/Mapp2/tab6output/sequences",   #it will make this folder for you
             out_proj="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")   #name a projection you want the output to be in. Then carry this proj through the rest of the steps

# step 1b. Create a lines files from the sequences.
create.lns.file(seqs_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/sequences",
                out_fldr="C:/Users/jmerkle/Desktop/Mapp2/tab6output/migration_lines",
                proj_of_dbfs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

#step 2. Conduct BB analyses. You will want use parrallel processing for this one. 
# This also spits out a metadata file of the results of the BB analysis
create.BBs(seqs_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/sequences",    #this is the folder where all the sequences are saved
                  BBs_out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/UDs",  #it will make this folder for you
                  footprint_out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/Footprints",  #it will make this folder for you
                  metadata_fldr="C:/Users/jmerkle/Desktop/Mapp2/tab6output",
                  cores=11, location.error=20, cell.size=50, max.lag=8, contour=99, time.step=5,
                  proj_of_dbfs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

#Step 3. Calculate average BBs for each individual and and a population UD and population footprint.

create.BB.avgs(BBs_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/UDs",    #this is the folder where all the UDs are saved.
               pop_BBs_out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/UDs_pop",  #it will make this folder for you
               pop_footprint_out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/Footprints_pop",  #it will make this folder for you
               contour=99,
               proj_of_ascs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

# Step 4. Calculate the stopover files and the low, medium, high use corridors as shapefiles
create.corridors.stopovers(PopUD_asc = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/UDs_pop/averageUD.asc",
                           PopFootprint_asc = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/Footprints_pop/popFootprint.asc",
                           pop_BBs_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/Footprints_pop",
                           out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/final_products",
                           stopover_percent=10, corridor_percents=c(10, 20), min_area = 20000, #this is in squared meters
                           simplify = TRUE, tolerance = 25, # how to polygons are simplified (unites are meters)
                           proj_of_ascs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

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
create.seqs.W(shpfl_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output", 
              shpfl_name= "pointsOut_20190207093941" ,
              idname="newUid",  #name of the column representing animal ID
              datename="nwMstrD",   #name of the column representing date in POSIX format
              mig.metadata.file="C:/Users/jmerkle/Desktop/Mapp2/tab6output/metadata.csv",  # metadata file from migration part of analysis
              qtl.end.fall.mig=0.95, 
              qtl.start.spring.mig=0.05,
              out_fldr="C:/Users/jmerkle/Desktop/Mapp2/tab6output/sequencesW",
              out_proj="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")   #name a projection you want the output to be in. Then carry this proj through the rest of the steps

#step 2. Conduct BB analyses. You will want use parrallel processing for this one. 
# This also spits out a metadata file of the results of the BB analysis
create.BBs.W(seqs_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/sequencesW",    #this is the folder where all the sequences are saved
             BBs_out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/UDsW",  #it will make this folder for you
             metadata_fldr="C:/Users/jmerkle/Desktop/Mapp2/tab6output",
             mindays=30,   #if an individual animal has less than this many days in a given sequence of winter data, it will be removed
             cores=11, location.error=20, cell.size=50, max.lag=8, time.step=5,
             proj_of_dbfs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

#Step 3. Calculate average BBs for each individual and and a population UD.
create.BB.avgs.W(BBs_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/UDsW",    #this is the folder where all the UDs are saved.
                 pop_BBs_out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/UDs_popW",  #it will make this folder for you
                 proj_of_ascs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")


# Step 4. Calculate the winter contours
create.core.areas.W(PopUD_asc = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/UDs_popW/averageUD_winter.asc",
                    out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/tab6output/final_productsW",
                    core_contours=c(10,20,30,40,50,60,70,80,90,99), min_area = 20000,  #this is in squared meters
                    simplify = TRUE, tolerance = 25, # unites are meters
                    proj_of_ascs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")



