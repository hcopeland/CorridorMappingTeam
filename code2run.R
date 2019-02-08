
# The following code finishes teh BB analysis of migration corridors after one has 
# completed tab 5 in Migration Mapper and clicked the export button on tab 6.

#source the functions you will need
source("C:/Users/jmerkle/Documents/GitHub/CorridorMappingTeam/functions/create.seqs.R")
source("C:/Users/jmerkle/Documents/GitHub/CorridorMappingTeam/functions/create.BBs.R")
source("C:/Users/jmerkle/Documents/GitHub/CorridorMappingTeam/functions/create.BB.avgs.R")

#step 1. Create sequences from the output of Migration Mapper (tab 6)
create.seqs(shpfl_fldr = "C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6",  #this is the folder where files from tab 6 are located
             shpfl_name= "pointsOut_20190208080234" ,   
             migtbl_name="migtime_20190208080234.csv",
             out_fldr="C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6/sequences",   #it will make this folder for you
             out_proj="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")   #name a projection you want the output to be in. Then carry this proj through the rest of the steps

#step 2. Conduct BB analyses. You will want use parrallel processing for this one. 
# This also spits out a metadata file of the results of the BB analysis
create.BBs(seqs_fldr = "C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6/sequences",
                  BBs_out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6/UDs",  #it will make this folder for you
                  footprint_out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6/Footprints",  #it will make this folder for you
                  metadata_fldr="C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6",
                  cores=11, location.error=20, cell.size=500, max.lag=8, contour=99,
                  proj_of_dbfs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

#Step 3. Calculate average BBs for each individual and and a population UD and population footprint.
create.BB.avgs(BBs_fldr = "C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6/UDs",
               pop_BBs_out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6/UDs_pop",  #it will make this folder for you
               pop_footprint_out_fldr = "C:/Users/jmerkle/Desktop/Mapp2/Elk_BenchCorral_Tab6/Footprints_pop",  #it will make this folder for you
               contour=99,
               proj_of_ascs="+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

# Step 4. Calculate the stopover files and the low, medium, high use corridors as shapefiles


