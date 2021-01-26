# Name: FixFilesforWGFD.py
# Bring in corridor file
# Separate by gridcode - select and export
# Smooth all separate files
# erase each grid code - next highest 2 from 1, 10 from 2, 20 from 10. that leaves 3 and then
# merge output from those three plus just 20. if loop - some don't have 10 - so if empty - then erase high from 2 or more.

# Import system modules
import arcpy, os
from arcpy import env
import arcpy.cartography as CA

arcpy.env.overwriteOutput = True
import pandas

# Set environment settings CHANGE THIS LOCATION
arcpy.env.workspace = "C:/Users/jmerkle_local/Desktop/corridors"


#Corridor, Winter Range and Stopover names CHANGE THESE #####################
corridors = "PH_WY_Sublette_Corridors_Ver1_2020.shp"   # file has GRIDCODE Column with 1, 10, 20
corridorOutput = "PH_WY_Sublette_CorridorsSmoothed_Edited_2020.shp"   # new name for smoothed output
wgfdcorridors ="PH_WY_Sublette_WGFDDesignatedCorridor_Ver1_2020.shp"   # output file for 

stopovers = "Stopovers_Edited.shp"   
stopoversmoothed = "StopoversSmoothed_Edited.shp"   # name for smoothed file

# winterrange = "MD_WY_UpperShoshone_WinterRange50Core_Ver1_2019.shp"
# winterrangesmoothed = "MD_WY_UpperShoshone_WinterRange50CoreSmoothed_Ver1_2019.shp"

####Change names above this line ##############
eraseOutput10 = "eraseOutput10.shp"
eraseOutput1 = "eraseOutput1.shp"
temp = "temp.shp"
corridors1 = "temp_1.shp"
corridors2 = "temp_2.shp"
corridors10 = "temp_10.shp"
corridors20 = "temp_20.shp"
corridorsunion = "temp_union.shp"

try:
    # Smooth corridor, stopover, winter range
    # query individual corridors and smooth
    expression1 = "\"GRIDCODE\" = 1"
    arcpy.Select_analysis(corridors, temp, expression1 )
    CA.SmoothPolygon(temp, corridors1, "PAEK", "200 Meters", "", "FLAG_ERRORS")

    expression1 = "\"GRIDCODE\" = 2"
    arcpy.Select_analysis(corridors, temp, expression1 )
    CA.SmoothPolygon(temp, wgfdcorridors, "PAEK", "200 Meters", "", "FLAG_ERRORS")

    expression1 = "\"GRIDCODE\" = 10"
    arcpy.Select_analysis(corridors, temp, expression1 )
    CA.SmoothPolygon(temp, corridors10, "PAEK","200 Meters", "", "FLAG_ERRORS")
    
    expression1 = "\"GRIDCODE\" = 20"
    arcpy.Select_analysis(corridors, temp, expression1 )
    CA.SmoothPolygon(temp, corridors20, "PAEK", "200 Meters", "", "FLAG_ERRORS")
    
    CA.SmoothPolygon(stopovers, stopoversmoothed, "PAEK","200 Meters", "", "FLAG_ERRORS")

    #expression1 = "\"GRIDCODE\" = 0.5"
    #arcpy.Select_analysis(winterrange, temp, expression1 )
    #CA.SmoothPolygon(temp, winterrangesmoothed, "PAEK","200 Meters", "", "FLAG_ERRORS")

    #erase features sequentially from the one above
    arcpy.Erase_analysis(corridors1, corridors10, eraseOutput1, "1 Meters")
    arcpy.Erase_analysis(corridors10, corridors20, eraseOutput10, "1 Meters")
    
    #merge/union back together
    arcpy.Union_analysis([eraseOutput1, eraseOutput10,corridors20], corridorsunion)
    arcpy.MakeFeatureLayer_management(corridorsunion, "corridor_lyr")
    #fix gridcode attribute
    arcpy.SelectLayerByAttribute_management("corridor_lyr", 'NEW_SELECTION', "\"GRIDCODE_1\" = 20")
    arcpy.CalculateField_management("corridor_lyr", "GRIDCODE", '!GRIDCODE_1!', "PYTHON_9.3")
    arcpy.SelectLayerByAttribute_management("corridor_lyr", 'NEW_SELECTION', "\"GRIDCODE_2\" = 10")
    arcpy.CalculateField_management("corridor_lyr", "GRIDCODE", '!GRIDCODE_2!', "PYTHON_9.3")
    arcpy.SelectLayerByAttribute_management("corridor_lyr", 'CLEAR_SELECTION')

    #delete all extra fields
    FCfields = [f.name for f in arcpy.ListFields("corridor_lyr")]
    DontDeleteFields = ['FID','Shape','GRIDCODE']
    fields2Delete = list(set(FCfields) - set(DontDeleteFields))
    arcpy.DeleteField_management("corridor_lyr", fields2Delete)
    arcpy.CopyFeatures_management("corridor_lyr", corridorOutput)
                                            
    arcpy.Delete_management("temp.shp")
    arcpy.Delete_management("temp_1.shp")
    arcpy.Delete_management("temp_2.shp")
    arcpy.Delete_management("temp_10.shp")
    arcpy.Delete_management("temp_20.shp")
    arcpy.Delete_management("temp_union.shp")
    arcpy.Delete_management("eraseOutput1.shp")
    arcpy.Delete_management("eraseOutput10.shp")
    print "Script finished"        
except arcpy.ExecuteError:
    print "Script error occurred"
    print arcpy.GetMessages(2)
