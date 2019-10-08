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
#Location of line file - CHANGE THIS DIRECTORY to where your data is!!!
arcpy.env.workspace = "C:/Holly/Work/Projects/EasternYellowstoneMuleDeer2016/Data/CMT_Dubois_8262019/migration_lines/"

# Replace a layer/table view name with a path to a dataset (which can be a layer file) or create the layer/table view within the script
# The following inputs are layers or table views: "migration_lines"
arcpy.SimplifyLine_cartography(in_features="migration_lines.shp", out_feature_class="migration_lines_simplify.shp",
                               algorithm="POINT_REMOVE", tolerance="500 Meters",
                               error_resolving_option="RESOLVE_ERRORS", collapsed_point_option="NO_KEEP",
                               error_checking_option="CHECK", in_barriers="")

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

#    print ("Max migration length = {} miles".format(max(cursor)))
