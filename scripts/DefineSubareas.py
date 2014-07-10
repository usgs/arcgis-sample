# Name:         Define subareas
# Purpose:      Set up raster study area for equal-area sample selection
# Author:       cprice
# Created:      3/25/2013 10:47:10 AM
# Environment:  ArcGIS 10.x, Python 2.6,2.7 arcpy
# -------------------------------------------------------------------

import os
import sys
import traceback
import arcpy
from arcpy import env
import samputil
from samputil import GPMsg, MsgError, ScratchFolder

def DefineSubareas(study_area, procRaster,
                   subarea_meth="NUMBER",subarea_val=1e5,
                   clip_extent=True):
    """Set up raster study area for sample selection.

    study_area   polygon study area
    procRaster   Subareas raster for use for sampling.
                 All non-target cells are set to NoData.
    subarea_meth subarea method: "NUMBER"|"DISTANCE" (opt)
    subarea_val  value for subarea_meth (opt)
    clip_extent  subareas defined based on selected features
    """
    try:
        # process arguments
        outWS = os.path.dirname(procRaster)
        if not arcpy.Exists(outWS):
            raise Exception, "Output workspace not found: " + outWS

        subarea_val = float(subarea_val)

        tv, Rows = None, None



        msg = "Defining subarea characteristics"
        GPMsg(msg)

        # define subarea creation extent
        D = arcpy.Describe(study_area)
        xyUnits = D.spatialReference.linearUnitName
        ext = D.extent

        fmtI = "  {0:<35s}{1:>8}"        # format for integer/string values
        fmtF = "  {0:<35s}{1:>8.1f} {2}" # format for float values w/ units

        # set up grid cell resolution
        subarea_val = float(subarea_val)
        if subarea_meth == "NUMBER":
            subarea_val = int(subarea_val)
            GPMsg(fmtI.format('Approximate number of subareas: ', subarea_val))
            subarea_dist = (ext.width * ext.height / float(subarea_val)) ** 0.5
        elif subarea_meth == "WIDTH":
            subarea_dist = float(subarea_val)
        else:
            raise Exception(
                  "Invalid value for subarea method: " + subarea_meth)

        GPMsg(fmtF.format('Subarea width (cell size):', subarea_dist,
                          xyUnits[0].lower()[0]))

        # set up raster processing environment
        env.workspace = env.workspace = ScratchFolder()
        env.cellSize = subarea_dist

        # expand extent by one cell (to avoid unintended clipping)
        ext = arcpy.Extent(ext.XMin - subarea_dist,
                           ext.YMin - subarea_dist,
                           ext.XMax + subarea_dist,
                           ext.YMax + subarea_dist)
        env.extent = ext

        # add a single-value field and convert
        # polygons to raster
        arcpy.AddField_management(study_area, "XXVAL", "LONG")
        arcpy.CalculateField_management(study_area, "XXVAL", "1", "PYTHON_9.3")
        arcpy.FeatureToRaster_conversion(study_area, "XXVAL",
                                         procRaster, subarea_dist)

        # verify that the raster has data
        try:
            arcpy.GetRasterProperties_management(procRaster, "MINIMUM")
        except:
            raise MsgError("Subarea creation failed")

        # get cell count
        tv = "tv"
        arcpy.MakeTableView_management(procRaster, tv)
        Rows = arcpy.SearchCursor(tv, "", "", "", "COUNT")
        s = Rows.next().getValue("COUNT")
        GPMsg(fmtI.format("Number of subareas created:", s))


    except MsgError, xmsg:
        GPMsg("e",str(xmsg))
    except arcpy.ExecuteError:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        GPMsg("e",tbinfo.strip())
        numMsg = arcpy.GetMessageCount()
        for i in range(0, numMsg):
            GPMsg("return",i)
    except Exception, xmsg:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        GPMsg("e",tbinfo + str(xmsg))
    finally:
        del Rows
        for f in [tv]:
            if f:
                try:
                    arcpy.Delete_management(f)
                except:
                    pass

if __name__ == "__main__":
    # ArcGIS Script tool interface
    # Get arguments, call tool
    argv = tuple(arcpy.GetParameterAsText(i) \
        for i in range(arcpy.GetArgumentCount()))
    DefineSubareas(*argv)

