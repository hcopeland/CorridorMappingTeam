# Name: CalcMigrationDistance.py
# Description: Calculates migration distances for metadata
# Requirements: none

# Import system modules
import arcpy, os
from arcpy import env
arcpy.env.cellSize = 100
arcpy.env.overwriteOutput = True
import pandas

# Set environment settings
#Location of line file and local variables - CHANGE THIS DIRECTORY to where your data is!!!
arcpy.env.workspace = "C:/Holly/Work/Projects/EasternYellowstoneMuleDeer2016/Data/CMT_Dubois_8262019/migration_lines/"

#Execute GeneratePointsAlongLines by percentage and convert back to lines
in_features = 'migration_lines.shp'
arcpy.GeneratePointsAlongLines_management(in_features, 'migration_points_simplify.shp', 'PERCENTAGE',
                                          Percentage=5, Include_End_Points='END_POINTS')

arcpy.PointsToLine_management('migration_points_simplify.shp', 'migration_lines_simplify.shp', 'mig')

#add area field
arcpy.AddField_management('migration_lines_simplify.shp','MILES','DOUBLE')
arcpy.CalculateField_management("migration_lines_simplify.shp",'MILES','!shape.length@miles!','PYTHON_9.3')


#Statistics to get mean, max, min
out_table = r"in_memory\stats_table"
stat_fields = [['MILES', 'MEAN'], ['MILES', 'MAX'], ['MILES', 'MIN']]

stats = arcpy.Statistics_analysis('migration_lines_simplify.shp', out_table, stat_fields)

# Get a list of field names to display
field_names = [i.name for i in arcpy.ListFields(out_table) if i.type != 'OID']

# Open a cursor to extract results from stats table
cursor = arcpy.da.SearchCursor(out_table, field_names)

# Create a pandas dataframe to display results
df = pandas.DataFrame(data=[row for row in cursor],
                      columns=field_names)

print(df)
print("Check your output in ArcGIS to be sure you have captured the longest migration!")

#    print ("Max migration length = {} miles".format(max(cursor)))
