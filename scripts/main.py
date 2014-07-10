# -----------------------------------------------------------------------------
# Purpose:  Random site selection software
# Author:   Matthew W Collier, USGS, OKWSC; And Curtis V Price, USGS, SDWSC
# Date:     January to May in FY2013
# Modified:     7/23/2013
# Environment:  ArcGIS 10.0, Python 2.6.x; arcpy, numpy

# Science References:
# - Scott, J.S., 1991, USGS OFR 90-4101, available at
#   http://pubs.er.usgs.gov/publication/wri904101

# Programming References:
# - General ArcGIS
#   http://help.arcgis.com/en/arcgisdesktop/10.0/help/index.html
# - General Python
#   http://docs.python.org/2.6/
# - numpy module
#   http://docs.scipy.org/doc/numpy/user/

# Although this software program has been used by the U.S. Geological Survey
# (USGS), no warranty, expressed or implied, is made by the USGS or the
# U.S. Government as to the accuracy and functioning of the program and
# related program material nor shall the fact of distribution constitute
# any such warranty, and no responsibility is assumed by the USGS
# in connection therewith.

# system modules
import os
import random
import sys
import time
import traceback
from copy import deepcopy

# third-party modules
import numpy
import arcpy
from arcpy import env

# local modules
from samputil import GPMode, GPMsg, ScratchFolder
from DefinePopulation import DefinePopulation
from DefineSubareas import DefineSubareas
from RotateFeatureClass import RotateFeatureClass

# arcpy starting environment
env.overwriteOutput = True

# Euclidean Distance function used for exclusion
def EucDist(p1, p2):
    """2D Euclidean distance"""
    return ( (p2[0]-p1[0])*(p2[0]-p1[0]) + (p2[1]-p1[1])*(p2[1]-p1[1]) )**0.5

# function to convert x,y to numpy array coords
def XYtoJI(x, y, ext, cellsize):
    """return a tuple of array-space coordinates

    arguments

    x, y       x and y coordinates (numeric)
    ext       arcpy Extent object (array extent)
    cellsize  square array cell size (numeric)
    """
    jFloat = ((ext.YMin + ext.height - y) / cellsize)
    iFloat = ((x - ext.XMin) / cellsize)
    if jFloat < 0 or iFloat < 0: raise
    return jFloat, iFloat

def DefineArea(study_area, where_expr, out_poly):
    """Create polygon area from raster or polygon input

    study_area   input study area (polygon or raster)
    out_poly     output polygon feature class
    """

    # create study area polygon from input study area
    tmpLyr = "tmpLyr"
    tmpFC = arcpy.CreateScratchName("xxstudy", "", "featureclass", "in_memory")
    if "Raster" not in arcpy.Describe(study_area).dataType:
        # polygon input - copy selected features only
        arcpy.MakeFeatureLayer_management(study_area, tmpLyr, where_expr)
        arcpy.CopyFeatures_management(tmpLyr, tmpFC)
    else:
        # raster input, convert selected zones to polygons
        arcpy.MakeRasterLayer_management(study_area, tmpLyr, where_expr)
        arcpy.RasterToPolygon_conversion(study_area, tmpFC, "NO_SIMPLIFY")

    # Create study area polygons
    arcpy.AddField_management(tmpFC, "TVAL", "LONG")
    arcpy.CalculateField_management(tmpFC, "TVAL", "1", "PYTHON_9.3")
    arcpy.CopyFeatures_management(tmpFC, out_poly)

    # clean up
    arcpy.Delete_management(tmpLyr)
    arcpy.Delete_management(tmpFC)

    return


def RunSample():
    """Sample ArcGIS script tool interface"""
    # Get arguments
    args = [arcpy.GetParameterAsText(i) \
            for i in range(arcpy.GetArgumentCount())]

    #  0 method                 EQUAL_AREA | SIMPLE | ITERATING_GRIDS
    #  1 Study_area             study area (polygons or raster)
    #  2 {Select_expression}    select expression
    #  3 {Primary_sites}        number of primary sites to select
    #  4 {Alternate_sites}      number of alternate sites
    #  5 {Exclusion_distance}   exclusion distance  (linear unit)
    #  6 {Rotation_angle}       rotation angle
    #  7 {Output_points}        output site points
    #  8 {Output_subcells}      subcell polygons
    #  9 {NUMBER | DISTANCE}    new site creation method
    # 10 {New_site_value}       new site creation parameter
    # 11 {Existing_sites}       existing sites to select from (point dataset)
    # 12 {NUMBER | WIDTH}       subarea creation method
    # 13 {Subarea_value}        subarea creation parameter
    # 14 {Convergence_fraction} iterating grids parameter
    # 15 {Iterations}           iterating grids parameter
    # 16 {Damping_factor}       iterating grids parameter

    try:
        method = args[0].upper()
        study_area = args[1]
        where_expr = args[2]
        Np = args[3]
        nalt   = args[4]
        exdist = args[5]
        rotate_angle = args[6]
        out_points = args[7]
        out_polys = args[8]
        newsite_meth = args[9].upper()
        newsite_val = args[10]
        in_sites = args[11]
        subarea_meth = args[12].upper()
        subarea_val = args[13]
        iter_frac = args[14]
        iter_count = args[15]
        iter_damp = args[16]

    except:
        raise

    # temp datasets
    spoly, spolyTmp, ssite, tmpSites, \
        tmpOutSites, tmpOutPolys = [None] * 6

    # check parameters
    methods = ["SIMPLE", "EQUAL_AREA", "ITERATING_GRIDS"]
    titles = ["Simple Random Selection",
              "Random Selection With Equal-Area Distribution",
              "Random Selection With Iterating Grids"]
    methcode = "SEI".find(method[0].upper())
    if methcode == -1:
        raise Exception(
            "Invalid sampling method: \"{0}\"".format(method))
    else:
        method = methods[methcode]
        # display a banner
        GPMsg("\n{0}\n{1}".format(titles[methcode], "="*50,))

    # calculate rotate angle
    try:
        rotate_angle = float(rotate_angle)
        if rotate_angle > 0:
            GPMsg("Rotating data before sampling by {} degrees".format(
                rotate_angle))
    except:
        rotate_angle = 0

    # arcpy environment
    SR = arcpy.Describe(study_area).spatialReference
    arcpy.overwriteOutput = True
    arcpy.outputCoordinateSystem = SR
    outWS = os.path.dirname(out_points)
    env.workspace = env.scratchWorkspace = outWS

    try:
        spoly = arcpy.CreateScratchName("spoly","","featureclass")
        ssite = arcpy.CreateScratchName("ssite","","featureclass")
        if not rotate_angle:
            # define study area polygon
            DefineArea(study_area, where_expr, spoly)
            if in_sites:
                arcpy.CopyFeatures_management(in_sites, ssite)
            else:
                ssite = None
        else:
            # define study area polygon, AND rotate it (and sites -- if provided)
            spolyTmp = arcpy.CreateScratchName("spolyTmp", "",
                                             "featureclass", outWS)
            DefineArea(study_area, where_expr, spolyTmp)
            xy = RotateFeatureClass(spolyTmp, spoly, rotate_angle)
            if in_sites:
                ssite = arcpy.CreateScratchName("ssite", "",
                                                    "featureclass", outWS)
                RotateFeatureClass(in_sites, ssite, rotate_angle, "XY", xy)
            else:
                ssite = None

        # define sample population
        tmpSites = arcpy.CreateScratchName("xxsites", "", "featureclass", outWS)
        num_points = DefinePopulation(spoly, tmpSites,
                                   newsite_meth, newsite_val, ssite)

        if not arcpy.Exists(tmpSites):
            raise Exception(
                  "Temporary dataset {0} could not be created".format(
                      tmpSites))

        # set up temporary output FC paths
        tmpOutSites = arcpy.CreateScratchName("xxoutpt", "", "featureclass", outWS)
        tmpOutPolys = arcpy.CreateScratchName("xxoutpoly", "", "featureclass", outWS)
        # run sample method
        if method == "SIMPLE":
            SimpleRandomSelection(tmpSites, tmpOutSites,
                                  Np, nalt, exdist)
            out_polys = None
        elif method == "EQUAL_AREA":
            EqualAreaSelection(spoly, tmpSites, tmpOutSites, tmpOutPolys,
                               Np, nalt, exdist,
                               subarea_meth, subarea_val)

        elif method == "ITERATING_GRIDS":
            IteratingGridsSelection(spoly, tmpSites, tmpOutSites, tmpOutPolys,
                                    Np, nalt, exdist,
                                    subarea_meth, subarea_val,
                                    iter_frac, iter_count, iter_damp)

        # If something went wrong, the output will not be there
        # just skip, the error messages will tell us

        try:
            # if rotated the input, rotate it back into place

            if rotate_angle:
                rot = rotate_angle * -1
                RotateFeatureClass(tmpOutSites, out_points, rot, "XY", xy)
                arcpy.Delete_management(tmpOutSites)
                if out_polys:
                    RotateFeatureClass(tmpOutPolys, out_polys, rot, "XY", xy)
            else:
                try:
                    arcpy.CopyFeatures_management(tmpOutSites, out_points)
                    if out_polys:
                        arcpy.CopyFeatures_management(tmpOutPolys, out_polys)
                except:
                    pass
        except:
            GPMsg("e", "Could not create output")
            pass

    except:
        raise
    finally:
        for f in [spoly, spolyTmp, ssite, tmpSites,
        tmpOutSites, tmpOutPolys]:
            try:
                if f: arcpy.Delete_management(f)
            except:
                pass



##############################################################################
### SimpleRandomSelection ####################################################
##############################################################################

def SimpleRandomSelection(in_points, out_points,
                          Np="10", nalt="0", exdist="1000 Meters"):
    """Simple random selection of an ArcGIS feature layer or table view

    in_points    input sample points
    out_points   output sample points
    Np           number of primary sites to select
    nalt         number of alternate sites to select for each primary site
    exdist       exclusion distance   (linear unit)
    """
    try:

        lyr = None
        Row, Rows = None, None

        Np = int(Np)
        nalt = int(nalt)
        exdist = int(exdist.split()[0])

        import random as r

        outWS = os.path.dirname(out_points)
        lyr, tmpPt = None, None

        env.workspace = env.scratchWorkspace = outWS

        tmpPt = arcpy.CreateScratchName("xxspt", "", "featureclass")
        arcpy.CopyFeatures_management(in_points, tmpPt)
        lyr = "tmpLyr0"
        arcpy.MakeFeatureLayer_management(tmpPt, lyr)

        PopPrepLabel = "Preparing population"
        arcpy.SetProgressor("step", PopPrepLabel, 0, 5)
        arcpy.SetProgressorPosition(1)
        GPMsg(PopPrepLabel)

        # populate samplePop dictionary
        # elements are: j,i,chosen, SITENO
        samplePop = {}
        sRows = arcpy.SearchCursor(lyr)
        dLyr = arcpy.Describe(lyr)
        OIDField = dLyr.OIDFieldName
        ShapeField = dLyr.shapeFieldName
        CHOSEN = 0
        SITENO = 0
        sRows = arcpy.SearchCursor(lyr, "", "", OIDField + ";" + ShapeField)
        for sRow in sRows:
            pnt = sRow.shape.getPart()
            pid = sRow.getValue(OIDField)
            samplePop[pid] = [pnt.X, pnt.Y, CHOSEN, SITENO]

        # Count number of input sites
        popSize = len(samplePop.keys())

        # user messages
        fmt = "  {0:<26} {1:>6}"
        GPMsg(fmt.format("Number of input points:", popSize))
        GPMsg(fmt.format("Number of primary sites:", Np))
        GPMsg(fmt.format("Number of alternate sites:", nalt))

        SelSiteLabel = "Selecting sites"
        arcpy.SetProgressor("Step", SelSiteLabel, 0, 100)
        arcpy.SetProgressorPosition(20)
        GPMsg(SelSiteLabel)

        ##### Verify sufficient population before proceeding...
        # Nr is number of sites to select in each round
        Nr = nalt + 1
        nLeast = Np*Nr
        if popSize < nLeast:
            raise Exception(
                "Insufficient Population (%s) to select %s sites."
                  % (popSize, nLeast))

        # Select w/o replacement from the list of candidates
        candidates = samplePop.keys()
        count = 0
        ProgFlag = 5 # progressor update interval
        for samp_round in range(Np):
            GPMsg('  Selecting site #%d...' % (samp_round+1,))
            numSamples = len(candidates)
            if numSamples < Nr:
                sampleSize = numSamples
                GPMsg("w", ("Selection quantity requested (%s) "
                            "exceeds available population (%s)")
                             % (Nr, numSamples))
            else: sampleSize = Nr
            selection = random.sample(candidates, sampleSize)

            # assign CHOSEN in samplePop, and remove picked sites
            for n, pt in enumerate(selection):
                samplePop[pt][2] = n + 1            # set CHOSEN
                samplePop[pt][3] = samp_round + 1   # set SITENO (1-based)
                candidates.remove(pt)

            # remove distance excluded sites
            if exdist > 0.0:
                exlist = []
                for pt in selection:
                    p1 = [samplePop[pt][0], samplePop[pt][1]]
                    for other in candidates:
                        p2 = [samplePop[other][0], samplePop[other][1]]
                        if EucDist(p2, p1) < exdist:
                            exlist.append(other)
                # unique the list of sites exlist excluded
                exlist = list(set(exlist))
                # tag the sites excluded
                for expt in exlist:
                    # assign CHOSEN value to -1 (excluded)
                    samplePop[expt][2] = -1
                    # remove from candidate list
                    candidates.remove(expt)

            count += 1
            if (count % ProgFlag) == 0:
                arcpy.SetProgressorPosition((20 + int(70.0 * count / Np)))

        arcpy.SetProgressorLabel("Writing points...")

        arcpy.SetProgressorPosition(90)
        ProgFlag = max(popSize // 10, 1)

        # add output fields
        arcpy.AddField_management(lyr, "CHOSEN", "LONG")
        arcpy.AddField_management(lyr, "SITENO", "LONG")

        Rows = arcpy.UpdateCursor(lyr)
        k = 0
        for Row in Rows:
            # Set CHOSEN status for each point
            pid = Row.getValue(OIDField)
            Row.setValue("CHOSEN", samplePop[pid][2])
            Row.setValue("SITENO", samplePop[pid][3])
            Rows.updateRow(Row)
            if (k % ProgFlag) == 0:
                arcpy.SetProgressorPosition(90 + int(10.0 * k / popSize))
            k += 1
        # clear cursor
        del Row, Rows
        Row, Rows = None, None # for cleanup

        arcpy.SetProgressorPosition(100)

        # copy output
        arcpy.CopyFeatures_management(lyr, out_points)

    except arcpy.ExecuteError:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        arcpy.AddError(arcpy.GetMessages(0))
    except Exception, xmsg:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        arcpy.AddError(tbinfo + str(xmsg))
    finally:
        # Clean up here (delete cursors, temp files)
        del Row, Rows
        for f in [lyr, tmpPt]:
            if f:
                try:
                    arcpy.Delete_management(f)
                except:
                    pass


##############################################################################
### EqualAreaSelection #######################################################
##############################################################################

def EqualAreaSelection(studyarea, in_points, out_points, out_poly,
                       Np=10, nalt=0, exdist="1000 Meters",
                       subarea_meth="NUMBER", subarea_val="1e5"):
    """EqualAreaSelection

    studyarea    study area polygons
    in_points    input sites to select from
    out_points   output sample points
    out_poly     output sampling cells
    Np           number of primary sites to select
    nalt         number of alternate sites to select for each primary site
    exdist       exclusion distance (linear unit)
    subarea_meth method to create new sites NUMBER|WIDTH
    subarea_val  subarea creation parameter value
    """

    try:
        # Initialize variables
        Row, Rows, iRow, iRows = [None] * 4
        tmpFC,procRas,tmpRas = [None] * 3
        lyrCells = "lyrCells"
        lyrSites = "lyrSites"
        fmtI = "  {0:<35s}{1:>8}"           # format to report int/str values
        fmtF = "  {0:<35s}{1:>8.1f} {2}"    # format to report float values

        # convert arguments to numbers
        Np = int(Np)
        nalt = int(nalt)
        exdist = int(exdist.split()[0]) # "10 Meters" -> 10

        # arcpy environment
        dsc = arcpy.Describe(studyarea)
        arcpy.env.outputCoordinateSystem = dsc.spatialReference
        xyUnits = dsc.spatialReference.linearUnitName
        outWS = os.path.dirname(out_poly)
        arcpy.env.workspace = arcpy.env.scratchWorkspace = outWS

        lyrSites, lyrCells, tmpPt, tmpCells, procRas, \
            tmpPRas, tmpRas, rasWK = [None] * 8

        # temp processing raster
        procRas = arcpy.CreateScratchName("", "", "raster", ScratchFolder())

        # Define subareas

        # if we rotated, xcen,ycen is pivot point
        pivotPoint = DefineSubareas(studyarea, procRas,
                       subarea_meth, subarea_val)

        msg = "Defining Sampling Parameters"
        GPMsg(msg)
        arcpy.SetProgressor("step",msg,0,5)
        arcpy.SetProgressorPosition(1)

        # Report: primary and alternate sites, and exclusion distance
        GPMsg(fmtI.format("Number of primary sites:", Np))
        GPMsg(fmtI.format("Number of alternate sites:", nalt))
        GPMsg(fmtF.format("Exclusion distance", exdist,
                          xyUnits.lower()[0]))

        # shift format over to match table
        fmtI = "  {0:<35s}{1:>13}"          # format to report int/str values
        fmtF = "  {0:<35s}{1:>13.1f} {2}"   # format to report float values

        # prepare numpy array from sampling cells

        msg = "Creating {0} equal-area cells".format(Np)
        GPMsg(msg)
        arcpy.SetProgressor("step", msg, 0, 5)
        arcpy.SetProgressorPosition(2)

        numpyNoData = 0
        numpyData = 1
        R = arcpy.RasterToNumPyArray(procRas, nodata_to_value=numpyNoData)
        R = R.astype('int32')
        J = len(R)      # number of rows in studyarea raster
        I = len(R[0])   # number of columns in studyarea raster

        # get raster properties
        D = arcpy.Describe(procRas)
        proc_cell = D.meanCellHeight
        jiRadius = exdist / proc_cell       # assumes square cells!
        procExt = D.extent
        dx = procExt.width
        dy = procExt.height
        left = procExt.XMin
        bottom = procExt.YMin

        # read in the sampling points
        GPMsg("Reading sample points...")
        # temp feature classes
        tmpPt = arcpy.CreateScratchName("xxpt", "", "featureclass")
        arcpy.CopyFeatures_management(in_points, tmpPt)
        lyrSites = "lyrSites"
        arcpy.MakeFeatureLayer_management(tmpPt, lyrSites)

        # populate samplePop (a Python dictionary)
        # key is OIDField, values are a Python list of j,i,chosen, and cellno
        samplePop = {}
        candidates = []
        sRows = arcpy.SearchCursor(lyrSites)
        OIDField = arcpy.Describe(lyrSites).OIDFieldName
        CHOSEN = 0
        CELLNO = 0
        for sRow in sRows:
            pnt = sRow.shape.getPart()
            j,i = None,None
            try: jFloat,iFloat = XYtoJI(pnt.X, pnt.Y, procExt, proc_cell)
            except:  # x,y invalid for extent
                continue
            try:
                j, i = int(jFloat), int(iFloat)
                R[j,i] == numpyData   # edge points may otherwise be lost in discretization
                pid = sRow.getValue(OIDField)
                samplePop[pid] = [jFloat, iFloat, CHOSEN, CELLNO]
                candidates.append(pid)
            except: pass
        del sRow, sRows
        arcpy.env.scratchWorkspace = rasWK = ScratchFolder()
        tmpPRAS = arcpy.NumPyArrayToRaster(R, arcpy.Point(left, bottom),
                                           proc_cell, proc_cell, numpyNoData)

        # Get total number of subareas
        # ("data" cells inside study area)
        numcells = 0
        r = arcpy.BuildRasterAttributeTable_management(tmpPRAS)
        tvRas = r.getOutput(0) # tool creates a table view
        Rows = arcpy.SearchCursor(tvRas, "", "", "COUNT")
        for Row in Rows:
            numcells += Row.getValue("COUNT")
        del Row, Rows
        Row, Rows = None, None
        arcpy.Delete_management(tvRas)
        tmpPRAS = tmpPRAS.catalogPath  # for easier cleanup

        nstrips = int(Np ** 0.5)

        arcpy.SetProgressorPosition(3)

        # Create dictionary of category counts for each column
        colSums = {}
        for i in range(I):
            tmp = {}
            for j in range(J-1,-1,-1): tmp[R[j,i]] = tmp.get(R[j,i],0) + 1
            colSums[i] = tmp

        # Create dictionary of over-all counts
        allSums = {}
        colKeys = colSums.keys()
        for colKey in colKeys:
            catKeys = colSums[colKey].keys()
            for catKey in catKeys:
                try:
                    allSums[catKey] = colSums[colKey].get(catKey,0) + \
                                      allSums[catKey]
                except:
                    allSums[catKey] = colSums[colKey].get(catKey,0)

        # Generate vertical strip boundaries
        cat = 1
        numVerticalStrips = int(0.5 + Np**0.5)
        numCatCellsPerStrip = int(0.5 + allSums[cat]/float(numVerticalStrips))

        GPMsg(fmtI.format('Number of vertical strips: ', nstrips))
        GPMsg(fmtI.format('Total area:', int(numcells * proc_cell ** 2)))
        GPMsg(fmtI.format('Total number of subareas:', numcells))

        # track vertical boundaries to partitions
        verticalPartitionBounds = [0]
        cellCount = 0
        for i in range(I):
            cellCount = cellCount + colSums[i].get(cat,0)
            if cellCount > numCatCellsPerStrip:
                verticalPartitionBounds.append(i)
                cellCount = colSums[i][cat]
        verticalPartitionBounds[-1] = I

        # Generate areal unit numpy map
        areaMap = R[:,:]
#####
##        GPMsg('##### areaMapSum = {0}'.format(areaMap.sum()))
##        GPMsg('##### allSums[cat] = {0}'.format(allSums[cat]))
#####
        numCellsPerAreaFloat = allSums[cat]/float(Np) #####
        numCellsPerAreaInt = int(numCellsPerAreaFloat) #####
        numCellsPerAreaActual = numCellsPerAreaInt
        arealTruncationErrorPerArea = numCellsPerAreaFloat - numCellsPerAreaInt #####
        arealTruncationErrorSum = 0.0 #####
        cellCount = 0
        stripNum = 0
        areaIdNum = 1
        numVerticalStrips = len(verticalPartitionBounds) - 1
        for a in range(numVerticalStrips):
            up = (a+1)%2 # up, or down aggregation
            if up: (start, stop, increment) = (J-1,-1,-1)
            else: (start, stop, increment) = (0,J,1)
            for j in range(start,stop,increment):
                for i in range(verticalPartitionBounds[stripNum],
                               verticalPartitionBounds[stripNum+1]):
                    if R[j,i] != numpyNoData:
                        cellCount += 1
                    ## my alternative
                    areaIdNum = 1 + int(Np * (cellCount - 1) /  numcells)
##                    if cellCount > (numCellsPerAreaActual): # if True, new area construction done
##                        if areaIdNum < Np: areaIdNum += 1
##                        cellCount = 1
##                        arealTruncationErrorSum = arealTruncationErrorSum + arealTruncationErrorPerArea
##                        if arealTruncationErrorSum >= 1.0:
##                            numCellsPerAreaActual += 1
##                            arealTruncationErrorSum = 0.0
##                        else:
##                            numCellsPerAreaActual = numCellsPerAreaInt

                    areaMap[j,i] = areaIdNum * R[j,i]
            stripNum += 1

        # Create polygon feature class of sampling cells
        # do raster processing in a folder (make a grid)
        arcpy.env.scratchWorkspace = rasWK
        tmpRas = arcpy.NumPyArrayToRaster(areaMap, arcpy.Point(left, bottom),
                                           proc_cell, proc_cell, numpyNoData)

        # Done raster processing, set the workspaces to output location
        arcpy.env.workspace = arcpy.env.scratchWorkspace = outWS

#######/  dbg BEGIN
##        import pickle
##        label = 1
##        GPMsg('##### Vertical Partition Bound Length = {0}'.format(len(verticalPartitionBounds)))
##        GPMsg(verticalPartitionBounds)
##        GPMsg( '##### ending strimNum = {0}'.format(stripNum) )
##        pickle.dump(samplePop, open('d:/tmp/samplePop%s.p' % label,'wb'))
##        procRas.save('d:/tmp/procRas%s' % label)
##        GPMsg( '##### numCellsPerAreaFloat = {0}'.format(numCellsPerAreaFloat) )
##        GPMsg( '##### numCellsPerAreaInt = {0}'.format(numCellsPerAreaInt) )
##        GPMsg( '##### numCellsPerAreaActual = {0}'.format(numCellsPerAreaActual) )
##        GPMsg( '##### arealTruncationErrorPerArea = {0}'.format(arealTruncationErrorPerArea) )
##        GPMsg( '##### numCellsPerArea = {0}'.format(numCellsPerArea) )
##        GPMsg( '##### numVerticalStrips = {0}'.format(numVerticalStrips) )
##        GPMsg( '##### stripNum = {0}'.format(stripNum) )
##        GPMsg( '##### areaIdNum = {0}'.format(areaIdNum) )
#######\  dgb END

        # report subareas as we go
        arcpy.MakeTableView_management(tmpRas,"tv")
        Rows = arcpy.SearchCursor("tv")
        cellData = []
        s = 0
        for Row in Rows:
            cellData.append([Row.VALUE,Row.COUNT])
            s += Row.COUNT
        Row, Rows = None, None
        totArea = int(s * (proc_cell ** 2))
        targetCells = s // Np
        target = totArea // Np
        fmtt1 = "  {0:>6} {1:>10} {2:>14} {3:>8} {4:>6}"
        fmtt2 = "  {0:>6} {1:>10} {2:>14} {3:>8} {4:>6.1f}"
        GPMsg(fmtt1.format("","Subareas","Area" ,"Diff","Pct"))
        GPMsg(fmtt1.format("TARGET", targetCells, target, "----", "----"))

        for row in cellData:
            rowArea = int(row[1] * (proc_cell ** 2))
            rowPct = round((100.0 * row[1] / s), 1)
            GPMsg(fmtt2.format(row[0], row[1], rowArea,
                            (rowArea-target), rowPct))
        GPMsg(fmtt2.format("Total", s, totArea, "----", 100))

        # Convert raster to cell polygons

        GPMsg("  Creating equal-area cell polygons")
        tmpCells = arcpy.CreateScratchName("xxcell", "", "featureclass")
        arcpy.RasterToPolygon_conversion(tmpRas, tmpCells, "NO_SIMPLIFY", "VALUE")
        tmpRas = tmpRas.catalogPath # for easier deletion
        lyrCells = "lyrCells"
        arcpy.MakeFeatureLayer_management(tmpCells, lyrCells)
        if arcpy.ListFields(lyrCells, "GRIDCODE"): valField = "GRIDCODE"
        else: valField = "GRID_CODE"
        arcpy.AddField_management(lyrCells, "CELLNO", "LONG")
        expr = "!{0}!".format(valField)
        arcpy.CalculateField_management(lyrCells, "CELLNO", expr, "PYTHON_9.3")
        arcpy.DeleteField_management(lyrCells,["Id", valField])

        # Assign CELLNO from areaMap to the sampling points

        subCellLists = {}
        for CELLNO in range(Np): subCellLists[CELLNO] = []
        for pid in samplePop.keys():
            j, i = int(samplePop[pid][0]), int(samplePop[pid][1])
            # areaMap counts from one to handle NoData
            CELLNO = areaMap[j,i]
            samplePop[pid][3] = CELLNO - 1
            if CELLNO > 0:
                subCellLists[CELLNO - 1].append(pid)

        # Make Selections

        msg = ("  Selecting from {0} available points..."
               ).format(len(samplePop))
        GPMsg(msg)
        arcpy.SetProgressor("step", msg, 0, 5)
        arcpy.SetProgressorPosition(4)

        Nr = nalt + 1

        for cellno in range(Np):
            GPMsg('  Selecting sites from cell #{0}'.format(cellno + 1))
            numSamples = len(subCellLists[cellno])
            if numSamples < Nr:
                sampleSize = numSamples
                GPMsg("w", ("Selection quantity requested ({0}) "
                            "exceeds available population ({1})").format(
                            Nr, numSamples))
            else: sampleSize = Nr
            selection = random.sample(subCellLists[cellno],sampleSize)

            # assign CHOSEN value to samplePop, and remove any picked sites
            for n, pt in enumerate(selection):
                # assign CHOSEN value
                # 1: primary site;  2,3...: alternate sites
                samplePop[pt][2] = n + 1
                # remove from cell/site and global candidate lists
                try:
                    subCellLists[cellno].remove(pt)
                    candidates.remove(pt)
                except:
                    pass

            # remove distance excluded sites
            if exdist > 0.0:
                exlist = []
                for pt in selection:
                    p1 = [samplePop[pt][0], samplePop[pt][1]]
                    for other in candidates:
                        p2 = [samplePop[other][0], samplePop[other][1]]
                        if EucDist(p2, p1) < jiRadius:
                            exlist.append(other)
                # unique the list of sites exlist excluded
                exlist = list(set(exlist))
                # tag the sites excluded
                for expt in exlist:
                    # assign CHOSEN value to -1 (excluded)
                    samplePop[expt][2] = -1
                    # remove from cell list, and candidate list
                    try:
                        subCellLists[samplePop[expt][3]].remove(expt)
                        candidates.remove(expt)
                    except:
                        pass

            # remove all remaining sites in this cell
            # from general candidates list (for performance reasons)
            for pt in subCellLists[cellno]:
                try: candidates.remove(pt)
                except: pass

        # Record sample results
        GPMsg('  Recording sample results...')
        arcpy.AddField_management(lyrSites, "CHOSEN", "SHORT")
        arcpy.AddField_management(lyrSites, "CELLNO", "SHORT")

        iRows = arcpy.UpdateCursor(lyrSites,"","",
                                   OIDField + ";CHOSEN;CELLNO")
        count = 0
        for iRow in iRows:
            OID = iRow.getValue(OIDField)
            iRow.chosen = samplePop[OID][2]
            iRow.cellno = samplePop[OID][3] + 1 # zero to one-based
            iRows.updateRow(iRow)
            count += 1
##        GPMsg('count = %d' % count)
        del iRow, iRows
        iRow, iRows = None, None # for cleanup routine


        msg = "Writing output points and polygons"
        GPMsg(msg)
        arcpy.SetProgressor("step", msg, 0, 5)
        arcpy.SetProgressorPosition(5)

        # remove points marked <null> --
        # these are sites that did not discretize into
        # data cells in the areaMap
        arcpy.SelectLayerByAttribute_management(lyrSites, "",
                                                "CHOSEN IS NOT NULL")

        arcpy.CopyFeatures_management(lyrSites, out_points)
        arcpy.CopyFeatures_management(lyrCells, out_poly)

    except arcpy.ExecuteError:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        arcpy.AddError(tbinfo.strip())
        GPMsg(arcpy.GetMessages(0))
    except Exception, xmsg:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        arcpy.AddError(tbinfo + str(xmsg))
    finally: # Clean up here (delete cursors, temp files)
        for var in [Row,Rows,iRow,iRows]:
            try: del var
            except: pass
        for f in [lyrSites, lyrCells, tmpPt, tmpCells,
                  procRas, tmpPRAS, tmpRas, rasWK]:
            if f:
                try:
                    arcpy.Delete_management(f)
                except:
                    pass

##############################################################################
### IteratingGridsSelection ##################################################
##############################################################################

def IteratingGridsSelection(studyarea, in_points, out_points, out_poly,
                            Np=10, nalt=0, exdist="1000 Meters",
                            subarea_meth="NUMBER", subarea_val="1e5",
                            iter_frac=0.75, iter_count=5, iter_damp=0.9):
    """Iterating grids method

    studyarea       study area polygons
    in_points       input sites to select from (assumed to be inside studyarea)
    out_points      output sample points
    out_poly        output sampling cells
    Np              number of primary sites to select
    nalt            number of alternate sites to select for each primary site
    exdist          exclusion distance   (linear unit)
    subarea_meth    method to create new sites NUMBER|WIDTH
    subarea_val     subarea creation parameter value
    iter_frac       fraction of requested sample that must contain enough sites
    iter_count      maximum  number of iterations allowed
    iter_damp       damping factor to control rate of convergence
    """

    try:

        fmtI = "  {0:<35s}{1:>8}"           # format to report int/str values
        fmtF = "  {0:<35s}{1:>8.1f} {2}"    # format to report float values

        # convert arguments from arcpy.GetParameterAsText() to numbers
        Np = int(Np)
        nalt = int(nalt)
        exdist = float(exdist.split()[0]) # "10 Meters" -> 10.0
##        GPMsg('iter_frac %s' % iter_frac)
        iter_frac = float(iter_frac)
        iter_count = int(iter_count)
        iter_damp = float(iter_damp)

        # count points
        popSize = int(arcpy.GetCount_management(in_points).getOutput(0))

        # Verify sufficient population before proceeding...
        Nr = nalt + 1
        nLeast = Np*Nr
        if popSize < nLeast:
            raise Exception(
                ("Insufficient Population ({0}) to select {1} sites."
                 ).format(popSize, nLeast))

        # setup arcpy environment
        dsc = arcpy.Describe(studyarea)
        arcpy.env.outputCoordinateSystem = dsc.spatialReference
        xyUnits = dsc.spatialReference.linearUnitName
        outWS = os.path.dirname(out_poly)
        arcpy.env.workspace = outWS
        arcpy.env.scratchWorkspace = outWS

        # Report on user supplied parameters
        GPMsg('Sampling parameters:')
        GPMsg(fmtI.format("Number of primary sites:", Np))
        GPMsg(fmtI.format("Number of alternate sites:", nalt))
        GPMsg(fmtF.format("Exclusion distance:", exdist, xyUnits.lower()[0]))
        GPMsg(fmtF.format("Cell fraction:", iter_frac,""))
        GPMsg(fmtI.format("Maximum iterations:", iter_count))
        GPMsg(fmtF.format("Damping factor:", iter_damp, ""))

        # Initialize geoprocessing variables
        Row, Rows, iRow, iRows = [None] * 4
        lyrCells, lyrSites, \
            tmpPt, tmpCells, procRas, tmpRas,  rasWK = [None] * 7

        rasWK = ScratchFolder()
        procRas = arcpy.CreateScratchName("", "", "raster", rasWK)

        # Define subareas (this is xres and yres in ArcSpeak)
        # if we rotated, (xcen,ycen) is pivot point
        DefineSubareas(studyarea, procRas, subarea_meth, subarea_val)

        # Rasterize study area, find useful properties
        numpyData = 1
        numpyNoData = 0
        R = arcpy.RasterToNumPyArray(procRas,nodata_to_value=numpyNoData)
        J = len(R)      # number of rows in studyarea raster
        I = len(R[0])   # number of columns in studyarea raster
        D = arcpy.Describe(procRas)
        proc_cell = D.meanCellHeight
        jiRadius = exdist / proc_cell       # assumes square cells
        procExt = D.extent
        dx = procExt.width
        dy = procExt.height
        left = procExt.XMin
        bottom = procExt.YMin
        numcells = J*I
        GPMsg(fmtI.format('Total number of subareas: ', numcells))

        # Prepare the sampling points
        GPMsg("Preparing sample points for selection...")
        tmpPt = arcpy.CreateScratchName("xxpt", "", "featureclass")
        arcpy.CopyFeatures_management(in_points, tmpPt)
        lyrSites = "lyrSites"
        arcpy.MakeFeatureLayer_management(tmpPt, lyrSites)

        # populate samplePop dictionary
        # elements are: j,i,chosen,cellno
        samplePop = {}
        sRows = arcpy.SearchCursor(lyrSites)
        OIDField = arcpy.Describe(lyrSites).OIDFieldName
        CHOSEN = 0
        CELLNO = 0
        for sRow in sRows:
            pnt = sRow.shape.getPart()
            try: jFloat,iFloat = XYtoJI(pnt.X, pnt.Y, procExt, proc_cell)
            except:  # x,y invalid for extent
                continue
            try:
                j, i = int(jFloat), int(iFloat)
                if R[j,i] == numpyNoData: continue   # get rid of noData
                pid = sRow.getValue(OIDField)
                samplePop[pid] = [jFloat,iFloat,CHOSEN,CELLNO]
            except: pass
        del sRow, sRows

        GPMsg(fmtI.format("  Available sample points:",
                          len(samplePop)))

        # Apply iterating grids method
        dataCells = (R == numpyData).sum()
        totCells = J*I
        inflate = Np * (dataCells / totCells)
        estimatedCellNum = Np * totCells / dataCells

        GPMsg("Making initial guess at grid characteristics")

        convergence = False
        for q in range(iter_count):
            GPMsg("  Beginning iteration #{0}".format(q+1))
            GPMsg("    Defining cell characteristics")

            # calculate grid parameters
            length = int( 0.5 + (J*I/estimatedCellNum)**0.5)
            offsetJ = int(-random.choice(numpy.arange(length-1)))
            offsetI = int(-random.choice(numpy.arange(length-1)))

            # assign CELLNO for this iteration (zero-based counting)
            cols = 1 + int(I - offsetI - 0.001)/length      # number of columns
            maxCellIndex = ((J-1)-offsetJ)/length*cols + ((I-1)-offsetI)/length
            GPMsg(("    Constructing a grid with {0} cells"
                   ).format(maxCellIndex+1))
            subCellLists = {}
            for num in range(maxCellIndex+1): subCellLists[num] = []
            countsByCELLNO = [0 for num in range(maxCellIndex+1)]
            candidates = []
            for key in samplePop.keys():
                j, i = int(samplePop[key][0]), int(samplePop[key][1])
                cellno = (j - offsetJ)/length*cols + (i - offsetI)/length
                samplePop[key][3] = cellno
                countsByCELLNO[cellno] = countsByCELLNO[cellno] + 1
                subCellLists[cellno].append(key)
                candidates.append(key)

            # prepare information about the cells' contents
            Nct, Na = 0, 0
            sufficientCellsList = []
            notSufficientCellsList = []
            for num in range(maxCellIndex+1):
                if countsByCELLNO[num] < Nr:
                    Na = Na + countsByCELLNO[num]
                    notSufficientCellsList.append(num)
                else:
                    Nct += 1
                    sufficientCellsList.append(num)

            # assess convergence and exiting criteria
            GPMsg("    Checking for convergence")
            ave = int(sum(countsByCELLNO) /
                      sum([True for e in countsByCELLNO if e <> 0]))
            condition1 = bool( iter_frac <= Nct/float(Np) <= 1.0 )
            condition2 = bool( Na >= Nr*(Np - Nct) )

            ## Debug code - display convergence criteria on each iteration
            ##
            ## GPMsg("(1) Fa <= (Nct / Np ) <= 1     %s <= ( %s / %s ) <= 1  [%.2f : %s]" % (iter_frac, Nct, Np, Nct/float(Np),condition1))
            ## GPMsg("(2) Na >= Nr(Np - Nct)         %s >= %s(%s - %s) [%.2f : %s]" % (Na, Nr, Np, Nct, Nr*(Np - Nct), condition2))

            ave = int(sum(countsByCELLNO) /
                      sum([True for e in countsByCELLNO if e <> 0]))

            if condition1 and condition2:
                convergence = True
                break
            GPMsg("w",
                  '    Iteration failed to converge using {0} cells'.format(
                      maxCellIndex + 1))
            GPMsg("w",('    Populated cells contain an average of {0} sites.'
                       ).format(ave))
            GPMsg("w",('    {0} cells contained at least {0} sites'
                       ).format(len(sufficientCellsList), Nr))

            if q + 1 == iter_count:
                GPMsg("e",
                      ('Maximum iterations of {0} reached without convergence'
                       ).format(iter_count))
                break

            # calculate new estimates of parameters for next iteration
            if Nct == 0:
                GPMsg("w",'    Solution was overshot, '
                      'may need smaller damping factor')
                estimatedCellNumber = inflate
            else:
                if 0.6 < Nct/float(Np) < 1.4: damping = iter_damp
                else: damping = 1.0
                ##  print Nct/float(Np), damping
                estimatedCellNum = estimatedCellNum * \
                    ( damping*(Np - Nct)/Nct + 1)
            GPMsg("w","    Convergence failed on iteration #{0}".format(q + 1))
            GPMsg("w",'    Next iteration will use at least {0} cells'.format(
                int(estimatedCellNum)))

        if convergence:
            GPMsg(('  Conversion succeeded on iteration #{0} using {1} cells'
                   ).format(q+1, maxCellIndex+1))
            GPMsg(('  Populated cells contain an average of {0} sites.'
                   ).format(ave))
            GPMsg(('  {0} cells contained at least {1} sites.'
                   ).format(len(sufficientCellsList), Nr))
            ## GPMsg(sufficientCellsList)

            kpick = 1 # tally of rounds

            GPMsg(("Selecting from {0} available points..."
                   ).format(len(samplePop)))

            # Stage I selection from cells with sufficient sites
            GPMsg('Selecting from well populated cells...')
            for cellno in sufficientCellsList:
                GPMsg('  ({0}) Selecting site #{1}...'.format(kpick, cellno+1))
                numSamples = len(subCellLists[cellno])
                if numSamples < Nr:
                    sampleSize = numSamples
                    GPMsg("w", ("Selection quantity requested ({0}) "
                                "exceeds available population ({1})"
                                ).format(Nr, numSamples))
                else: sampleSize = Nr
                selection = random.sample(subCellLists[cellno],sampleSize)

                # assign CHOSEN value to samplePop, and remove any picked sites
                for n, pt in enumerate(selection):
                    # assign CHOSEN value
                    # 1: primary site;  2,3...: alternate sites
                    samplePop[pt][2] = n + 1
                    # remove from cell/site and global candidate lists
                    try:
                        subCellLists[cellno].remove(pt)
                        candidates.remove(pt)
                    except:
                        pass

                # remove distance excluded sites
                if exdist > 0.0:
                    exlist = []
                    for pt in selection:
                        p1 = [samplePop[pt][0], samplePop[pt][1]]
                        for other in candidates:
                            p2 = [samplePop[other][0], samplePop[other][1]]
                            if EucDist(p2, p1) < jiRadius:
                                exlist.append(other)
                    # unique the list of sites exlist excluded
                    exlist = list(set(exlist))
                    # tag the sites excluded
                    for expt in exlist:
                        # assign CHOSEN value to -1 (excluded)
                        samplePop[expt][2] = -1
                        # remove from cell list, and candidate list
                        try:
                            subCellLists[samplePop[expt][3]].remove(expt)
                            candidates.remove(expt)
                        except:
                            pass

                # remove all remaining sites in this cell
                # from general candidates list (for performance reasons)
                for pt in subCellLists[cellno]:
                    try: candidates.remove(pt)
                    except: pass
                kpick += 1

            # Stage II selection: from less populated cells
            aggregateCellNum = maxCellIndex
            for key in candidates: samplePop[key][3] = 0
            sitesNeeded = Np - len(sufficientCellsList)
            if sitesNeeded:
                GPMsg('Selecting from aggregate of less-populated cells...')
                for pick_round in range(sitesNeeded):
                    siteno = aggregateCellNum * 100 + pick_round
                    GPMsg(('  ({0}) Selecting site #{1}...'
                           ).format(kpick, siteno+1))
                    numSamples = len(candidates)
                    if numSamples < Nr:
                        sampleSize = numSamples
                        GPMsg("w", ("Selection quantity requested ({0}) "
                                    "exceeds available population ({1})"
                                    ).format(Nr, numSamples))
                    else:
                        sampleSize = Nr
                    selection = random.sample(candidates, sampleSize)

                    # assign CHOSEN in samplePop, and remove picked sites
                    for n, pt in enumerate(selection):
                        samplePop[pt][2] = n + 1                # set CHOSEN
                        samplePop[pt][3] = siteno
                        try:
                            candidates.remove(pt)
                        except:
                            pass

                    # remove distance excluded sites
                    if exdist > 0.0:
                        exlist = []
                        for pt in selection:
                            p1 = [samplePop[pt][0], samplePop[pt][1]]
                            for other in candidates:
                                p2 = [samplePop[other][0], samplePop[other][1]]
                                if EucDist(p2, p1) < jiRadius:
                                    exlist.append(other)
                        # unique the list of sites exlist excluded
                        exlist = list(set(exlist))
                        # tag the sites excluded
                        for expt in exlist:
                            # assign CHOSEN value to -1 (excluded)
                            samplePop[expt][2] = -1
                            # remove from candidate list
                            try:
                                candidates.remove(expt)
                            except:
                                pass
                    kpick += 1
                # assign remaining candidates
                # (these are all in insuff cells) to agg cell
                for other in candidates:
                    samplePop[other][3] = aggregateCellNum - 1

            # generate the numpy areaMap, add 1 for output for 1-based counting
            cellMap = numpy.zeros((J,I), dtype=int)
            for j in range(J):
                for i in range(I):
                    cellno = (j - offsetJ)/length*cols + (i - offsetI)/length
                    if cellno in sufficientCellsList: cellMap[j,i] = cellno + 1
                    else: cellMap[j,i] = aggregateCellNum

            # Create polygon feature class of sampling cells
            GPMsg('Creating sampling-cell polygons...')

            # do raster processing in a folder (make a grid)
            arcpy.env.workspace = arcpy.env.scratchWorkspace = ScratchFolder()
            tmpRas = arcpy.NumPyArrayToRaster(
                cellMap, arcpy.Point(left, bottom),
                proc_cell, proc_cell, numpyNoData)
            arcpy.env.workspace = arcpy.env.scratchWorkspace = outWS

            # Convert raster to cell polygons
            tmpCells = arcpy.CreateScratchName("xxcells", "", "featureclass")
            arcpy.RasterToPolygon_conversion(tmpRas, tmpCells,
                                             "NO_SIMPLIFY", "VALUE")
            tmpRas = tmpRas.catalogPath # easier cleanup
            lyrCells = "lyrCells"
            arcpy.MakeFeatureLayer_management(tmpCells, lyrCells)
            if arcpy.ListFields(lyrCells, "GRIDCODE"): valField = "GRIDCODE"
            else: valField = "GRID_CODE"
            arcpy.AddField_management(lyrCells, "CELLNO", "LONG")
            expr = "!{0}!".format(valField)
            arcpy.CalculateField_management(lyrCells, "CELLNO",
                                            expr, "PYTHON_9.3")
            arcpy.DeleteField_management(lyrCells, ["Id", valField])

            # Record sample results
            GPMsg('Recording sample results...')
            arcpy.AddField_management(lyrSites, "CHOSEN", "SHORT")
            arcpy.AddField_management(lyrSites, "CELLNO", "SHORT")

            iRows = arcpy.UpdateCursor(lyrSites, "", "",
                                       OIDField + ";CHOSEN;CELLNO")
            offset = 0
            for iRow in iRows:
                OID = iRow.getValue(OIDField)
                try:
                    iRow.chosen = samplePop[OID][2]
                    iRow.cellno = samplePop[OID][3] + 1
                    iRows.updateRow(iRow)
                except:
                    pass
            del iRow, iRows
            iRow, iRows = None, None # for cleanup routine

            # remove points marked <null> --
            # these are sites that did not discretize into
            # data cells in the areaMap
            arcpy.SelectLayerByAttribute_management(lyrSites, "",
                                                    "CHOSEN IS NOT NULL")

            arcpy.CopyFeatures_management(lyrSites, out_points)
            arcpy.CopyFeatures_management(lyrCells, out_poly)

    except arcpy.ExecuteError:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        arcpy.AddError(tbinfo.strip())
        GPMsg(arcpy.GetMessages(0))
    except Exception, xmsg:
        tbinfo = traceback.format_tb(sys.exc_info()[2])[0]
        arcpy.AddError(tbinfo + str(xmsg))
    finally:
        for var in [Row,Rows,iRow,iRows]:
            try: del var
            except: pass
        for f in [lyrCells, lyrSites,
                  tmpPt, tmpCells, procRas, tmpRas, rasWK]:
            if f:
                try:
                    arcpy.Delete_management(f)
                except:
                    pass

    return True


if __name__ == '__main__':
    RunSample()

