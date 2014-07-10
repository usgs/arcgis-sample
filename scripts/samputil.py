"""Python utilities for sample scripts

"Although this program has been used by the USGS, no warranty, expressed or
implied, is made by the USGS or the United States Government as to the
accuracy and functioning of the program and related program material nor shall
the fact of distribution constitute any such warranty, and no responsibility
is assumed by the USGS in connection therewith."

Author: Curtis Price, cprice@usgs.gov
"""
import sys
import os
import arcpy

# Set up timer function for GPMsg()
# Matthew Collier, USGS, mcollier@usgs.gov
import time
global dt
t0 = time.clock()      # for execution time reporting at end of script
dt = [0,t0]            # set up a time object to report on differences
def timeDiff(T=None):  # define function to juggle and print time differences
    """Timer function for use in GPMsg()"""
    if not T: T=dt
    T[0] = T[1]
    T[1] = time.clock()
    return T[1] - T[0]


def GetCount(input):
    """Count rows in  input, returns an integer

    (modified from ESRI HelperFunctions.py)
    """
    try:
        Result = arcpy.GetCount_management(input)
        return int(Result.getOutput(0))
    except Exception, msg:
        raise Exception, str(msg) + "\nCould not get count from %s" % input


class MsgError(Exception):
    """Exception class for reporting a user error from a script

    Example

      raise MsgError, "Houston, we have a problem...."

    """


class GPModeThing:
    """A little Python class to keep track of print mode for GPMsg()

    See the help for the GPMsg() function for details.
      """
    def __init__(self):
        self.data = "gp"  # set default: "gp"

    def __call__(self,strMode=None):
        #  "" strMode returns current value"""
        if strMode:
            # check argument to make sure it is valid
            strMode = strMode.lower()
            if strMode not in ["gp","print","both","off"]:
                print 'Valid values are: "gp","print","both", or "off"'
            else:
                self.data = strMode
        return self.data

# initialize it
GPMode = GPModeThing()


def GPMsg(sev=None,msg=None,msgType=None):
    """Send a message to the geoprocessor,"python-print", or both

    Geoprocessing messages displayed with methods like "arcpy.AddMessage"
    may be visible in ArcGIS or in tty python output, but when
    running the script within IDE environments like IDLE or Wing,
    invisible, or display garbled. This method allows you
    to specify to send a message to python-print, just the geoprocessor, or
    to both places. The syntax also allows for a little less typing
    than the gp messaging methods it calls.

    dependencies

    GPModeThing

    arguments

    sev - severity code / message option

      "Message","Warning","Error","Return","Time","Gpmessages"
      (only the first letter is required)

     msg - text message to display

      A special syntax for msg is used to support arcpy.AddIDMessage().
      (Spaces are required between each argument!)

        ID <MessageID> {AddArgument1} {AddArgument2}

      For example, to do this:
        arcpy.AddIDMessage("Error", 12, outFeatureClass)
      You can use this syntax with GPMsg:
        GPMsg("Error","ID %s %s" % (12,outFeatureClass))
      If a message argument contains a space, you can use | to separate
        the second argument so it will parse correctly:
         GPMsg("Error","ID %s %s|%s" % (328,"Input dataset",outFeatureClass))
      (Please only use error numbers documented in the ArcGIS help!)

     msgType - Where to send the message. If this argument is given
        the destination will stay the same until GPMode is used or
        the argument is given again.

        "gp"      Geoprocessor (default) (arcpy or 9.x gp object)
        "print"   Python print
        "both"    Both places
        "off"     Nothing prints anywhere (use with care)
        None      Use current value of GPMode()

    examples

       GPMode("print")  # change default output to python-print
       GPMsg("This is a message") # default output, print ONLY
       GPMsg("t","The time is now","gp") # output to ArcGIS ONLY
       GPMsg() # print arcpy.AddMessages(0) GP messages to ARCGIS only
       GPMsg("w","ID 345","both") # use arcpy.AddIDMessage, output to gp,print
       GPMode("off") # no messages printed

       Output:

       This is a message
       10:40:05 The time is now
       Executing: CopyFeatures poly_inout E:\work\poly_inout_CopyFeatures.shp # 0 0 0
       Start Time: Wed Apr 07 11:11:58 2010
       Executed (CopyFeatures) successfully.
       End Time: Wed Apr 07 11:11:58 2010 (Elapsed Time: 0.00 seconds)
       WARNING 000345: Input must be a workspace
    """


    # support shorthand usage: GPMsg("message") and GPMsg()
    if sev != None and msg == None:
        # GPMsg("message") ->  GPMsg("","message")
        sev,msg = None,sev
    elif sev == None and msg == None:
        # GPMsg() -> GPMsg("Message",arcpy.GetMessages(0))
        sev,msg  = "g",None

    if not msgType:
        msgType = GPMode()  # use current value of GPMode
    else:
        msgType = GPMode(msgType) # set GPMode (and remember for next GPMsg)

    if msgType == "off":
        # Do not print anything! (like AML "&messages &off &all")
        return

    # decode severity to a code 0 thru 5:
    # sev  isev description
    #        0  message
    # "w"    1  warning
    # "e"    2  error
    # "r"    3  returnmessage
    # "t"    4  message with clock time & sec since last check
    # "g"    5  return arcpy.GetMessages(0)
    dictSev = {"w":1, "e":2, "r":3, "t":4,"g":5 }
    try:
        sev = str(sev).lower()[:1] # "Warning" -> "w"
        isev = dictSev[sev]
    except:
        isev = 0

    # support arcpy.AddIDMessage
    IDMessage = False # assume this isn't an id message
    try:
        # Usage: "ID <msgID> {|} {arg1} {|} {arg2}"
        lstMsg = msg.split()
        if lstMsg[0].lower() != "id": raise
        try:
            MessageID = int(lstMsg[1])
            if MessageID <= 0 or MessageID > 99999: raise
        except:
            GPMsg("w","GPMsg: Invalid message ID: %s" % MessageID)
            raise
        xmsg = " ".join(lstMsg[2:])
        if xmsg.find("|") > -1:
            lstMsg = xmsg.split("|")
        else:
            lstMsg = xmsg.split()

        IDMessage = True
        # capture AddArgument1, AddArgument2
        try:
            IDArg1, IDArg2 = "", ""
            IDArg1 = lstMsg[0]
            IDArg2 = lstMsg[1]
        except:
            pass
    except:
        pass

    # send our message

    if msgType.lower() in ["gp","both"]:
        # send message to geoprocessor
        if isev == 0:
            arcpy.AddMessage(msg)
        elif isev == 1:
            if not IDMessage:
                arcpy.AddWarning(msg)
            else:
                arcpy.AddIDMessage("Warning",MessageID,IDArg1,IDArg2)
        elif isev == 2:
            if not IDMessage:
                arcpy.AddError(msg)
            else:
                arcpy.AddIDMessage("Error",MessageID,IDArg1,IDArg2)
        elif isev == 3:
            arcpy.AddReturnMessage(int(msg))
        elif isev == 4:
            curTime = time.strftime("%H:%M:%S ", time.localtime())
            elpTime = timeDiff()
            arcpy.AddMessage("%s %.2f %s" % (curTime,elpTime,msg))
        elif isev == 5:
            arcpy.AddMessage(arcpy.GetMessages(0))
    if msgType.lower() in ["print","both"]:
        # python-print messages
        SevLabel = ["","WARNING","ERROR"]
        if isev == 0:
            print msg
        elif isev in [1,2]:
            if not IDMessage:
                print("%s: %s" % (SevLabel[isev],msg))
            else:
                print("%s %06d: %s" % (SevLabel[isev],MessageID,lstMsg))
        elif isev == 3:
            print(arcpy.GetMessage(int(msg)))
        elif isev == 4:
            curTime = time.strftime("%H:%M:%S ", time.localtime())
            elpTime = timeDiff()
            print("%s %.2f %s" % (curTime,elpTime,msg))
        elif isev == 5:
            msg = arcpy.GetMessages(0)
            if len(msg) > 0: print(msg)


def ScratchFolder():
    """Create a scratch folder

    1) if current workspace is a folder, in the current workspace
    2) if a .mdb or .gdb, parallel to it
    3) in TEMP
    """
    import arcgisscripting
    gp = arcgisscripting.create()
    try:
        # this works at 10.1
        sw = gp.scratchFolder
    except:
        try:
            sw = gp.scratchWorkspace
            swType = gp.Describe(sw).dataType
            if swType == "Folder":
                sw = os.path.join(sw,"scratch")
            elif swType == "Workspace":
                pth = os.path.dirname(sw)
                if not gp.Exists(pth): raise
                sw = os.path.join(pth,"scratch")
        except:
            # put it in TEMP
            sw = os.path.join(os.environ["TEMP"],"scratch")
        finally:
            if not gp.Exists(sw):
                os.mkdir(sw)
                print "created " + sw
    return sw

