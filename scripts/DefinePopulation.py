# Name:         Define Population
# Purpose:      Prepare sample sites for a study area
# Author:       cprice
# Created:      3/25/2013 10:47:10 AM
# Modified:     4/2/2013 3:56:35 PM
# Environment:  ArcGIS 10.x, Python 2.6,2.7 arcpy
# -------------------------------------------------------------------

import os
import sys
import traceback
import random
import arcpy
from arcpy import env
from samputil import GPMsg, MsgError, ScratchFolder

def DefinePopulation(study_area, out_points,
                     newsite_meth="NUMBER", newsite_val=10000, site_points=None):
    """Define population

    study_area    input study area polygons
    out_points    output points
    newsite_meth  new site generation parameter: "NUMBER","DISTANCE"
    newsite_val   value to use for above method
    site_points   existing site points
    """

    try:

        # initialize temp file variables
        lyrStudy, lyrSites, tmpFC, tmpRas, tmpPoints, numSites = [None] * 6
        rasWK = None

        fmtI = "  {0:<35s}{1:>8}"        # format to report integer/string values
        fmtF = "  {0:<35s}{1:>8.1f} {2}" # format to report float values w/ units

        lyrStudy = "lyr1"
        arcpy.MakeFeatureLayer_management(study_area, lyrStudy)

        # set processing environment
        arcpy.env.workspace = os.path.dirname(out_points)
        D = arcpy.Describe(study_area)
        env.extent = ext = D.extent
        env.outputCoordinateSystem = D.spatialReference
        xyUnits = D.spatialReference.linearUnitName
        arcpy.ClearEnvironment("snapRaster")
        rasWK = ScratchFolder()

        procLabel = "Defining population characteristics"
        GPMsg(procLabel)

        if not site_points:

            GPMsg("  Creating points...")

            # Prepare a population of inside study area

            if newsite_meth == "NUMBER":
                newsite_val = int(newsite_val)
                GPMsg(fmtI.format('Approximate number of sites:', newsite_val))
                samp_dist = ((ext.width * ext.height) / newsite_val) ** 0.5
            elif newsite_meth == "DISTANCE":
                samp_dist = float(newsite_val)
            else:
                raise Exception("Invalid new site method " + newsite_meth)

            GPMsg(fmtF.format(
                'Sample distance:', samp_dist, xyUnits.lower()[0]))

            # randomize the lattice origin
            xmin = ext.XMin - samp_dist * random.random()
            ymin = ext.YMin - samp_dist * random.random()
            env.extent = arcpy.Extent(xmin, ymin, ext.XMax, ext.YMax)

            # Report number sites
            n = int((env.extent.width * env.extent.height) /
                    (samp_dist ** 2))
            GPMsg(fmtI.format(
                "Building a population with", n) + " sites")

            # Create a raster covering the the study area
            tmpRas = arcpy.CreateScratchName("saras", "", "raster", rasWK)
            arcpy.FeatureToRaster_conversion(lyrStudy, D.OIDFieldName,
                                             tmpRas, samp_dist)

            # check raster - are there data cells?
            try:
                arcpy.GetRasterProperties_management(tmpRas, "MINIMUM")
            except:
                GPMsg()
                raise MsgError("No points created")

            # Generate a point lattice from raster cell centroids
            tmpPoints = arcpy.CreateScratchName("pt", "",
                                                "featureclass", rasWK)
            arcpy.RasterToPoint_conversion(tmpRas, tmpPoints, "VALUE")
            lyrSites = "lyrSites"
            arcpy.MakeFeatureLayer_management(tmpPoints, lyrSites)
            arcpy.DeleteField_management(lyrSites, "GRID_CODE;GRIDCODE")

            # count points
            numSites = int(arcpy.GetCount_management(lyrSites).getOutput(0))
            GPMsg(fmtI.format("Points inside study area:", numSites))

        else:

            # Select points from an existing point feature class

            lyrSites = "lyrSites"
            arcpy.MakeFeatureLayer_management(site_points, lyrSites)
            numSites = int(arcpy.GetCount_management(lyrSites).getOutput(0))
            # select points within study area
            arcpy.SelectLayerByLocation_management(lyrSites, "WITHIN", lyrStudy)

            # check number of sites selected
            numSelected = int(arcpy.GetCount_management(lyrSites).getOutput(0))
            if not numSelected:
                raise MsgError("No points selected")
            nsel = "{0}/{1}".format(numSelected, numSites)
            GPMsg(fmtI.format("Points inside study area:", nsel))
            numSites = numSelected

        # copy points to output
        arcpy.CopyFeatures_management(lyrSites, out_points)

    except MsgError, xmsg:
        GPMsg("e", str(xmsg))
    except arcpy.ExecuteError:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        GPMsg("e", tbinfo.strip())
        numMsg = arcpy.GetMessageCount()
        for i in range(0, numMsg):
            GPMsg("return", i)
    except Exception, xmsg:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        GPMsg("e", tbinfo + str(xmsg))
    finally:
        for f in [lyrStudy, lyrSites, tmpFC, tmpRas, tmpPoints, rasWK]:
            if f:
                try:
                    arcpy.Delete_management(f)
                except:
                    pass
        arcpy.ClearEnvironment("extent")

    return numSites

if __name__ == "__main__":
    # ArcGIS Script tool interface
    # Get arguments, call tool
    argv = tuple(arcpy.GetParameterAsText(i)
        for i in range(arcpy.GetArgumentCount()))
    DefinePopulation(*argv)

