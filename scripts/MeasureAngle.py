# calculate an azimuth angle from a interactively entered
# line (feature set)
#
# Curtis Price, cprice@usgs.gov,  9/18/2013 11:51:10 AM

# MEASURE AN ANGLE

import math
import arcpy



def get_angle(xy1, xy2):
  """Calculate azimuth angle from two points. (Zero is north.)"""
  import math
  try:
    # ArcPy point objects
    x1, y1, x2, y2 = xy1.X, xy1.Y, xy2.X, xy2.Y
  except:
    # xy strings, e.g. "0 0"
    x1, y1 = [float(x) for x in xy1.split()]
    x2, y2 = [float(x) for x in xy2.split()]
  dx, dy = (x2 - x1, y2 - y1)
  return 90 - math.degrees(math.atan2(dy, dx))

# script tool interface
if __name__ == '__main__':
  try:
    # read feature set
    line = arcpy.GetParameterAsText(0)
    # get first and last point of a line
    SHAPE = arcpy.Describe(line).shapeFieldName
    Rows = arcpy.SearchCursor(line,"","",SHAPE)
    feat = Rows.next().getValue(SHAPE)
    pt1 = feat.firstPoint
    pt2 = feat.lastPoint
    angle = get_angle(pt1, pt2)
    msg1 = "  First point: {0:.1f}, {0:.1f}".format(pt1.X, pt1.Y)
    msg2 = "  Last point:  {0:.1f}, {0:.1f}".format(pt2.X, pt2.Y)
    msg3 = "  Azimuth angle (in degrees):     {0:.1f}".format(angle)
    msg4 = "  Rotation angle for Sample tool: {0:.1f}".format(angle * -1)
    arcpy.AddMessage("{0}\n{1}\n{2}\n{3}".format(msg1, msg2, msg3, msg4))
    arcpy.SetParameterAsText(1, angle)
  except Exception, msg:
    arcpy.AddError(str(msg))
    arcpy.AddIDMessage("Error", 152)