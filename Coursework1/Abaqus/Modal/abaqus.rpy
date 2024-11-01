# -*- coding: mbcs -*-
#
# Abaqus/CAE Release 2022.HF8 replay file
# Internal Version: 2023_03_27-08.24.44 RELr424 176957
# Run by sglallpr on Fri Nov  1 14:24:09 2024
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
openMdb('Modaltruss.cae')
#: The model database "M:\Documents\AERO417\AERO417\Coursework1\Abaqus\Modal\Modaltruss.cae" has been opened.
session.viewports['Viewport: 1'].setValues(displayedObject=None)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=ON)
p = mdb.models['Truss'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
a = mdb.models['Truss'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
o3 = session.openOdb(
    name='M:/Documents/AERO417/AERO417/Coursework1/Abaqus/Modal/Job-1.odb')
#: Model: M:/Documents/AERO417/AERO417/Coursework1/Abaqus/Modal/Job-1.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       3
#: Number of Node Sets:          6
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
a = mdb.models['Truss'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(step='modal')
mdb.models['Truss'].steps['modal'].setValues(numEigen=5)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=OFF)
mdb.Job(name='Job-2', model='Truss', description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
    multiprocessingMode=DEFAULT, numCpus=1, numGPUs=0)
mdb.jobs['Job-2'].submit(consistencyChecking=OFF)
#: The job input file "Job-2.inp" has been submitted for analysis.
#: Job Job-2: Analysis Input File Processor completed successfully.
#: Job Job-2: Abaqus/Standard completed successfully.
#: Job Job-2 completed successfully. 
session.viewports['Viewport: 1'].setValues(
    displayedObject=session.odbs['M:/Documents/AERO417/AERO417/Coursework1/Abaqus/Modal/Job-1.odb'])
o3 = session.openOdb(
    name='M:/Documents/AERO417/AERO417/Coursework1/Abaqus/Modal/Job-2.odb')
#: Model: M:/Documents/AERO417/AERO417/Coursework1/Abaqus/Modal/Job-2.odb
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       3
#: Number of Node Sets:          6
#: Number of Steps:              1
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].partDisplay.setValues(sectionAssignments=ON, 
    engineeringFeatures=ON)
session.viewports['Viewport: 1'].partDisplay.geometryOptions.setValues(
    referenceRepresentation=OFF)
p1 = mdb.models['Truss'].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
a = mdb.models['Truss'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].setValues(
    displayedObject=session.odbs['M:/Documents/AERO417/AERO417/Coursework1/Abaqus/Modal/Job-2.odb'])
o3 = session.openOdb(
    name='M:/Documents/AERO417/AERO417/Coursework1/Abaqus/Modal/Job-2.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o3)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    DEFORMED, ))
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    UNDEFORMED, ))
a = mdb.models['Truss'].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
    adaptiveMeshConstraints=ON)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(loads=ON, bcs=ON, 
    predefinedFields=ON, connectors=ON, adaptiveMeshConstraints=OFF)
