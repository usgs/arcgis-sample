------------------------------------------------------------------------------
Getting Started:
This package is supported under ArcGIS Desktop 10.0/10.1/10.2 software. No support exists for previous versions of ArcGIS.

This folder contains the software, associated files, and test datasets. You should get started by reviewing pages 1-14 of Scott(1990). It will explain the original methodology. Every reasonable effort has been made to closely follow this methodology. However, the description of how to run the software from page 15 onwards may be ignored as it describes the original implementation using ArcInfo Workstation software (now deprecated by Esri). This new implementation will be accessed entirely through an ArcGIS 10.x Python script tool.

Once you are comfortable with the methodology, open up the sampleTest MXD file. You'll find two data frames for testing purposes. The first is a simplified study area illustrated on page 21 of Scott(1990) and used in his example from pages 30 to 50. While various layers may be present, they may not display if not yet created. So open ArcToolbox, find the 'Sample' toolbox, and run the 'Sample' script. Built-in help will appear in the help pane, and selections will automatically grey out if not active in a particular context. Once the layers are built, you can select them in the table of contents and view them. Take some time to try out the different selection methods, and vary the parameters as you see fit. For instance, you can generate a grid of test points as described in the methodology, or you could use the 'testarea_point' feature class provided with the geodatabase. Symbology has been provided as a suggestion, but should be tailored for your individual needs.

Finally, move on to the more complex test study area in the other data frame when you are ready (Or, your own data!). This complex area is based on the example provided on pages 50 through 53 in Scott 1990. No test points are provided with this data frame, but you can easily generate test grids as above.

References:

Squillace, Paul J., and Price, Curtis V., 1996, Urban land-use study plan for the National Water-Quality Assessment Program, USGS Open-File Report 96-217, 19 p. http://pubs.er.usgs.gov/publication/ofr96217

Scott, Jonathan C., 1990, Computerized stratified random site-selection approaches for design of a ground-water-quality sampling network, USGS Water-Resources Investigations Report 90-4101, 109 p.  http://pubs.er.usgs.gov/publication/wri904101

------------------------------------------------------------------------------
Description of files and subfolders:

00README.txt
    This file.
LICENSE.txt
    Software license information
sampleData.gdb
    A collection of test data for the sample tools in Esri
    ArcGIS 10.0 geodatabase format.
scratch
    A folder to receive intermediary geoprocessing files. You
    may wish to occasionally check this folder to clear it of
    unnecessary files.
scripts
    A folder with supporting Python code and ArcGIS layer files.
sampleTest.mxd
    An ArcGIS map document arranged to quickly illustrate the package.
Sample Tools.tbx
    ArcGIS toolbox used to access the tools from ArcGIS software

------------------------------------------------------------------------------