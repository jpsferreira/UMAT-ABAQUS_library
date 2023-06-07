# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 6.14-3 replay file
# Internal Version: 2015_02_02-21.14.46 134785
# Run by jferreira on Fri Dec 18 14:10:08 2015
#

# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
import os
from abaqus import *
from abaqusConstants import *

session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=243.416656494141, 
    height=165.277770996094)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()

odb_file = os.path.join(os.getcwd(), 'cube_gho_nld.odb')

o1 = session.openOdb(name=odb_file)
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
odb = session.odbs[odb_file]
session.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=(('U', 
    NODAL, ((COMPONENT, 'U2'), )), ('S', INTEGRATION_POINT, ((INVARIANT, 
    'Max. Principal'), )), ('SDV_DET', INTEGRATION_POINT), ), nodePick=((
    'PART-1-1', 1, ('[#1 ]', )), ), )
xy1 = session.xyDataObjects['U:U2 PI: PART-1-1 N: 1']
xy2 = session.xyDataObjects['S:Max Principal (Avg: 75%) PI: PART-1-1 N: 1']
xy3 = combine(xy1+1, xy2)
xy3.setValues(
    sourceDescription='combine ( "U:U2 PI: PART-1-1 N: 1"+1,"S:Max Principal (Avg: 75%) PI: PART-1-1 N: 1" )')
tmpName = xy3.name
session.xyDataObjects.changeKey(tmpName, 'Displacement vs SMax_Principal')
x0 = session.xyDataObjects['Displacement vs SMax_Principal']
session.writeXYReport(fileName='output.txt', xyData=(x0, ))

