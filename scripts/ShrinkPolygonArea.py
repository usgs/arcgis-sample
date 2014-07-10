"""
ArcGIS Script Tool - Shrink Polygon Area
"""
# Name:         ShrinkPolygonArea.py
# Purpose:      Modifies a study area to ensure selected sample points
#               are fully inside the areas of interest
# Author:       Curtis Price, U.S. Geological Survey, cprice@usgs.gov
# Modified:     03/22/2013 10:11:13 AM
#
# ---------------------------------------------------------------------------
# This is a ArcGIS/Python implementation of SHRINK.AML, described in:
#
# Squillace, P.J., and Price, C.V., 1996, Urban land-use study plan for NAWQA,
# U.S.Geological Survey Open-File Report 96-217, 19 p. (Appendix B, SHRINK.AML)
# http://pubs.er.usgs.gov/publication/ofr96217
#
# Environment:  ArcGIS 10.x, Python 2.6,2.7 arcpy
# Requires ArcGIS Desktop Spatial Analyst extension
# ---------------------------------------------------------------------------

import sys
import os
import traceback
import arcpy
from arcpy import env
from arcpy.sa import *
import samputil
from samputil import GPMsg, MsgError, ScratchFolder

def ShrinkPolygonArea(in_polys, where_expr, out_polys, buf="1000 Meters",
                      pct_inside=75, shrink_dist="200 Meters",
                      proc_cell=None,
                      ps_polys=None,buf100_polys=None,simplify=True):
    """Modifies study area polygons to ensure selected sample points
       are fully inside the area of interest

       arguments:

       in_polys
       where_expr
       out_polys
       buf
       pct_inside
       shrink_dist
       proc_cell
       ps_polys
       buf100_polys
       simplify

       This is a ArcGIS/Python implementation of SHRINK.AML, described in:

       Squillace, P.J., and Price, C.V., 1996, Urban land-use study plan
       for NAWQA, U.S.Geological Survey Open-File Report 96-217, 19 p.
       (Appendix B, SHRINK.AML)
       http://pubs.er.usgs.gov/publication/ofr96217
       """

    try:

        # temp dataset variables
        lyrPolys, tmpFC, tmpPoly, tmpClip, \
            tmpMask, tmpSum, tmpPct, tmpShrink,\
            = [None] * 8

        # check for Spatial Analyst license
        if arcpy.CheckOutExtension("spatial") == "NotLicensed":
            raise MsgError, "Spatial Analyst license is required"

        # prepare parameters

        buf_dist = float(buf.split()[0])

        # Percent inside to select (default = 75)
        try:
            pct_inside = float(pct_inside)
            if pct_inside <= 0 or pct_inside >= 100: raise Exception
        except:
            GPMsg("w","Percent_inside must be between 1 and 100")
            GPMsg("w","Using value: 75.0")
            pct_inside = 75.0

        # "Shrink" distance (optional "final shrink" from study area edge)
        SD = shrink_dist.split()
        try:
            shrink_dist = float(SD[0])
            try:
                shrink_units = SD[1]
            except:
                shrink_units = "Unknown"
        except:
            shrink_dist = None

        # Set processing cell size
        try:
            # use user-specified cell size
            proc_cell = float(proc_cell)
        except:
            try:
                # else use environment setting
                proc_cell = float(env.cellSize)
                GPMsg("w","Using environment cell size "\
                      "{0:.1f}".format(proc_cell))
            except:
                # else estimate a value from buffer
                fact = 25.0
                proc_cell = buf_dist / fact

                GPMsg("w","Using cell size {0:.1f} ".format(proc_cell) + \
                      "(Buffer radius / {0}".format(fact))

        # support boolean or string input (True, "SIMPLIFY", "true")
        if str(simplify).lower()[0] in "ts":
            simplify = "SIMPLIFY"
        else:
            simplify = "NO_SIMPLIFY"

        # Set up environment

        env.overwriteOutput = 1
        env.cellSize = proc_cell

        # set up output workspace
        outWS = os.path.realpath(os.path.dirname(out_polys))
        # use a folder for raster processing
        env.workspace = env.scratchWorkspace = ScratchFolder()

        # first create mask area (select/project)
        tmpFC = arcpy.CreateScratchName("xxshrink","","featureclass",outWS)
        lyrPolys = "lyrPolys"
        arcpy.MakeFeatureLayer_management(in_polys,lyrPolys,where_expr)
        arcpy.CopyFeatures_management(lyrPolys,tmpFC)

        # set processing extent to input polygons + 1/2 buffer radius
        ext = arcpy.Describe(tmpFC).extent
        ext.XMin -= buf_dist * 0.5
        ext.YMin -= buf_dist * 0.5
        ext.XMax += buf_dist * 0.5
        ext.YMax += buf_dist * 0.5
        env.extent = ext

        # Check extent for problems
        ExtDiag = ( ext.width ** 2 + ext.height ** 2 ) ** 0.5
        if buf_dist > ExtDiag / 2.0:
            raise MsgError, "Buffer size too large for extent of input features"

        # do raster processing in a folder workspace (grids are fastest)
        if arcpy.Describe(outWS).workspaceType == "LocalDatabase":
            tmpWK = os.path.dirname(outWS)
            env.workspace = env.scratchWorkspace = outWS

        # total up "inside cells"
        GPMsg("Counting neighborhood cells...")
        tmpMask = CreateConstantRaster(1)
        tmpClip = arcpy.CreateScratchName("","","raster")
        arcpy.Clip_management(tmpMask,"#",tmpClip,tmpFC,"#","ClippingGeometry")

        Buf = NbrCircle(buf_dist,"MAP")
        tmpSum = Int(FocalStatistics(tmpClip, Buf, "SUM", True))

        # estimate total number of cells in a circular buffer
        bufCells = math.pi * ( buf_dist ** 2 ) / ( proc_cell ** 2 )

        GPMsg("Calculating percent inside areas...")
        # int ( 100 * data cells in buffer / total cells in buffer)
        tmpPct = Int(100 * tmpSum / bufCells)
        tmpShrink = SetNull(tmpPct <= pct_inside,1)

        if buf100_polys:
            try:
                buf100 = SetNull(tmpPct < 99, 1)
                arcpy.RasterToPolygon_conversion(buf100, buf100_polys, simplify)
            except:
                GPMsg("w", arcpy.GetMessages(0))

        if shrink_dist:
            # additional shrink (distance from edge of area)

            # save pre-shrink raster if specified
            if ps_polys:
                arcpy.RasterToPolygon_conversion(tmpShrink, ps_polys, simplify)
            GPMsg("Shrinking edges by %.1f ..." % shrink_dist)
            tmpShrink = SetNull( \
                EucDistance(SetNull(IsNull(tmpShrink) == 0, 1)) \
                  <= shrink_dist, 1)
            env.cellSize = None
            arcpy.RasterToPolygon_conversion(tmpShrink, out_polys, simplify)
        else:
            env.cellSize = None
            arcpy.RasterToPolygon_conversion(tmpShrink, out_polys, simplify)

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
        ##GPMsg("cleaning up")
        for f in [lyrPolys, tmpFC, tmpClip, tmpPoly]:
            try:
                if f: arcpy.Delete_management(f)
            except:
                pass

if __name__ == "__main__":
    # ArcGIS Script tool interface
    # Get arguments, call tool
    argv = tuple(arcpy.GetParameterAsText(i)
        for i in range(arcpy.GetArgumentCount()))
    ShrinkPolygonArea(*argv)