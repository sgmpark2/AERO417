# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2022.HF8 replay file
# Internal Version: 2023_03_27-08.24.44 RELr424 176957
# Run by sglallpr on Sat Nov  2 14:13:30 2024
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=177.884902954102, 
    height=217.677780151367)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
openMdb('Dyanamictruss2.cae')
#: The model database "M:\Documents\AERO417\AERO417\Coursework1\Abaqus\Dynamic\Dyanamictruss2.cae" has been opened.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p = mdb.models['Model-1'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
a = mdb.models['Model-1'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
o3 = session.openOdb(
    name='M:/Documents/AERO417/AERO417/Coursework1/Abaqus/Dynamic/Job-2.odb')
#: Model: M:/Documents/AERO417/AERO417/Coursework1/Abaqus/Dynamic/Job-2.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       2
#: Number of Node Sets:          9
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
odb = session.odbs['M:/Documents/AERO417/AERO417/Coursework1/Abaqus/Dynamic/Job-2.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Spatial displacement: U1 at Node 2 in NSET Node  2', 
    steps=('transient', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.XYPlot('XYPlot-1')
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
odb = session.odbs['M:/Documents/AERO417/AERO417/Coursework1/Abaqus/Dynamic/Job-2.odb']
xy1 = xyPlot.XYDataFromHistory(odb=odb, 
    outputVariableName='Spatial displacement: U2 at Node 2 in NSET Node  2', 
    steps=('transient', ), suppressQuery=True, __linkedVpName__='Viewport: 1')
c1 = session.Curve(xyData=xy1)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
chart.setValues(curvesToPlot=(c1, ), )
session.charts[chartName].autoColor(lines=True, symbols=True)
