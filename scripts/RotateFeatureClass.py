# Name:         RotateFC.py
# Purpose:      Rotates a feature class
# Author:       Curtis Price, cprice@usgs.gov
# Created:      09/17/2013 11:49:25 AM
# Environment:  ArcGIS 10.x
# -------------------------------------------------------------------

import os
import sys
import traceback
import arcpy
from arcpy import env
from samputil import GPMsg, MsgError
import samputil

def RotateXY(x, y, xc=0, yc=0, angle=0, units="DEGREES"):
    """Rotate an xy cooordinate about a specified origin

    x,y      xy coordinates
    xc,yc   center of rotation
    angle   angle
    units    "DEGREES" (default) or "RADIANS"
    """
    import math
    x = x - xc
    y = y - yc
    # make angle clockwise (like Rotate_management)
    angle = angle * -1
    if units == "DEGREES": angle = math.radians(angle)
    xr = (x * math.cos(angle)) - (y * math.sin(angle)) + xc
    yr = (x * math.sin(angle)) + (y * math.cos(angle)) + yc
    return xr, yr


def RotateFeatureClass(inputFC, outputFC,
                       angle=0, pivot_type="CENTER", pivot_point=None):
    """Rotate Feature Class

    inputFC     Input features
    outputFC    Output feature class
    angle       Angle to rotate, in degrees
    pivot_type  CENTER, LOWER_LEFT, LOWER_RIGHT, UPPER_LEFT, UPPER_RIGHT, XY
    pivot_point X,Y coordinates (as space-separated string)
                (only used if pivot_type = "XY")
    """

    # temp names for cleanup
    env_file, lyrTmp, tmpFC  = [None] * 3
    Row, Rows, oRow, oRows = [None] * 4 # cursors

    # warning thresholds for complex features
    warnPts = 10000
    warnPrt = 10

    try:
        # process parameters

        try:
            angle = float(angle)
        except:
            raise MsgError("Invalid value for rotation: {0}".format(angle))

        # determine pivot point

        ext = arcpy.Describe(inputFC).extent

        if pivot_type == "CENTER":
            pivot_point = (
                ext.XMin + ext.width * 0.5,
                ext.YMin + ext.height * 0.5)
        elif pivot_type == "LOWER_LEFT":
            pivot_point = ext.XMin, ext.YMin
        elif pivot_type == "UPPER_LEFT":
            pivot_point = ext.XMin, ext.YMax
        elif pivot_type == "UPPER_RIGHT":
            pivot_point = ext.XMax, ext.YMax
        elif pivot_type == "LOWER_RIGHT":
            pivot_point = ext.XMax, ext.YMin
        elif pivot_type == "XY":
            try:
                pivot_point = tuple([float(xy) for xy in pivot_point.split()])
            except:
                raise Exception(
                    "Invalid value for pivot point: %s" % pivot_point)
        xcen, ycen = pivot_point


        # msg = "Rotating {0} degrees around {1} ({2:.1f}, {3:.1f})"
        # GPMsg(msg.format(angle, pivot_type, xcen, ycen))

        # set up environment
        env_file = arcpy.CreateScratchName("xxenv",".xml","file",
                                           arcpy.GetSystemEnvironment("TEMP"))
        arcpy.SaveSettings(env_file)

        # Disable any GP environment clips or project on the fly
        arcpy.ClearEnvironment("extent")
        arcpy.ClearEnvironment("outputCoordinateSystem")

        WKS = env.workspace
        if not WKS:
            if os.path.dirname(outputFC):
                WKS = os.path.dirname(outputFC)
            else:
                WKS = os.path.dirname(
                    arcpy.Describe(inputFC).catalogPath)
        env.workspace = env.scratchWorkspace = WKS

        # get feature class properties
        dFC = arcpy.Describe(inputFC)
        shpField = dFC.shapeFieldName
        shpType = dFC.shapeType
        FID = dFC.OIDFieldName
        SR = dFC.spatialReference

        # create temp feature class
        tmpFC = arcpy.CreateScratchName("xxfc","","featureclass")
        arcpy.CreateFeatureclass_management(os.path.dirname(tmpFC),
                                            os.path.basename(tmpFC),
                                            shpType,
                                            spatial_reference=SR)
        # set up id field (used to join later)
        TFID = "ORIG_FID"
        arcpy.AddField_management(tmpFC, TFID, "LONG")
        arcpy.DeleteField_management(tmpFC, "ID")


        # rotate the feature class coordinates
        # only points, polylines, and polygons are supported

        ##GPMsg("t", "writing")
        # open read and write cursors

        Rows = arcpy.SearchCursor(inputFC, "", "",
                                  "%s;%s" % (shpField,FID))
        oRows = arcpy.InsertCursor(tmpFC)

        if shpType  == "Point":
            for Row in Rows:
                shp = Row.getValue(shpField)
                pnt = shp.getPart()
                pnt.X, pnt.Y = RotateXY(pnt.X, pnt.Y, xcen, ycen, angle)
                oRow = oRows.newRow()
                oRow.setValue(shpField, pnt)
                oRow.setValue(TFID, Row.getValue(FID))
                oRows.insertRow(oRow)
        elif shpType in ["Polyline","Polygon"]:

            # initialize total area / length
            totarea_in, totarea_out = 0.0, 0.0

            parts = arcpy.Array()
            rings = arcpy.Array()
            ring = arcpy.Array()
            for Row in Rows:
                shp = Row.getValue(shpField)
                if shpType == "Polygon":
                    totarea_in += shp.area
                p = 0
                for part in shp:
                    for pnt in part:
                        if pnt:
                            x, y = RotateXY(pnt.X, pnt.Y, xcen, ycen, angle)
                            ring.add(arcpy.Point(x, y, pnt.ID))
                        else:
                            # if we have a ring, save it
                            if len(ring) > 0:
                                rings.add(ring)
                                ring.removeAll()
                    # we have our last ring, add it
                    rings.add(ring)
                    ring.removeAll()
                    # if only one, remove nesting
                    if len(rings) == 1: rings = rings.getObject(0)
                    parts.add(rings)
                    rings.removeAll()
                    p += 1

                # if only one, remove nesting
                if len(parts) == 1: parts = parts.getObject(0)
                if dFC.shapeType == "Polyline":
                    shp = arcpy.Polyline(parts)
                else:
                    shp = arcpy.Polygon(parts)
                    totarea_out += shp.area
                if shp.pointCount > warnPts or shp.partCount > warnPrt:
                        GPMsg("w", ("Feature {0} contains {1} points, "
                                    "{2} parts").format(
                              Row.getValue(FID),
                              shp.pointCount,
                              shp.partCount))
                parts.removeAll()
                oRow = oRows.newRow()
                oRow.setValue(shpField, shp)
                oRow.setValue(TFID,Row.getValue(FID))
                oRows.insertRow(oRow)

        else:
            raise Exception, "Shape type {0} is not supported".format(shpType)

        del Row, Rows, oRow, oRows # close write cursor (ensure buffer written)
        Row, Rows, oRow, oRows = [None] * 4 # restore variables for cleanup

        if shpType == "Polygon":
            # check - did the area change more than 0.1 percent?
            diff = totarea_out - totarea_in
            diffpct = 100.0 * diff / totarea_in
            if abs(diffpct) > 1.0:
                GPMsg("w", "Please check output polygons")
                GPMsg("w", "Input area: {0:.1f}  Output area: {1:.1f}".format(
                    totarea_in, totarea_out))
                GPMsg("w", "({:.1f} percent change)".format(diffpct))




        # join attributes, and copy to output
        lyrTmp = "lyrTmp"
        arcpy.MakeFeatureLayer_management(tmpFC, lyrTmp)
        arcpy.AddJoin_management(lyrTmp, TFID, inputFC, FID)
        env.qualifiedFieldNames = False
        arcpy.CopyFeatures_management(lyrTmp, outputFC)
        env.qualifiedFieldNames = True

        # Fourth field [3] is a duplicate of TFID
        dropField = arcpy.ListFields(outputFC)[3].name
        arcpy.DeleteField_management(outputFC, dropField)

    except MsgError, xmsg:
        GPMsg("e",str(xmsg))
    except arcpy.ExecuteError:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        GPMsg("e",tbinfo.strip())
        arcpy.AddError(arcpy.GetMessages())
        numMsg = arcpy.GetMessageCount()
        for i in range(0, numMsg):
            GPMsg("return",i)
    except Exception, xmsg:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        GPMsg("e",tbinfo + str(xmsg))
    finally:
        # reset environment
        if env_file:
            arcpy.LoadSettings(env_file)
        # delete cursors
        try:
            for c in [Row, Rows, oRow, oRows]: del c
        except:
            pass
        # Clean up temp files
        for f in [lyrTmp, tmpFC, env_file]:
            try:
                if f: arcpy.Delete_management(f)
            except:
                GPMsg()
                pass


        # return pivot point
        try:
            pivot_point = "{0} {1}".format(*pivot_point)
        except:
            pivot_point = None

        return pivot_point


if __name__ == '__main__':
    """ArcGIS script tool interface

    Arguments:

    Input_features
    Output_features
    Rotation_angle
    Rotation_point ("CENTER","LOWER_LEFT",...,"XY")
    Pivot point (ignored unless Rotation_point is "XY")
    Pivot point value (derived)
    """
    # last argument is derived
    numArgs = arcpy.GetArgumentCount() - 1
    argv = tuple([arcpy.GetParameterAsText(i)
        for i in range(numArgs)])
    xy = RotateFeatureClass(*argv)
    # populate derived parameter (pivot point)
    arcpy.SetParameterAsText(5, xy)

