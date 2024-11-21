#!/usr/bin/python


import sys
import os
from bishopUtils import *
from ImplicitShapes import *


beginJob()

#
#=========== pyro construction ============================
#

allparms =  AllParameters(PyroSphereParameters())

CmdLineHelp( "-h" )


allparms["-pyrogamma"] = 0.75
allparms["-amp"] = 1.5
allparms["-octaves"] = 6.0
allparms["-roughness"] = 0.45
allparms["-translate"] = [ 0, 0, 0.34 ]
allparms["-name"] = "host_pyro.exr"
allparms["-scatter"] = 10.0
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-fov"] = 70.0
allparms["-ds"] = 0.01
allparms["-near"] = 4.0
allparms["-radius"] = 2.5
allparms["-amp"] = 1.0
allparms["-maxpathlength"] = 10.0
allparms["-dsmcell"] = [ 0.04, 0.04, 0.04 ]
allparms["-dsmpartition"] = 10
allparms["-eye"] = [ 0, 0, 9 ]
allparms["-view"] = [ 0, 0, -1 ]
allparms["-subsamples"] = 5
#allparms["-size"] = [ 960, 540 ]
allparms["-clamp"] = 0.2

allparms["-nbthreads"] = 3
allparms["-color"] = [1, 1, 1]

allparms["-lightstyle"] = "keyrimfill"

allparms = CmdLineFindValues( allparms )
PrintParameters( allparms )

print "Building Pyro" 
volumeObject = PyroSphere(allparms)

clampv = float( allparms["-clamp"] )

llc = Vector( float(allparms["-dsmllc"][0]), float(allparms["-dsmllc"][1]), float(allparms["-dsmllc"][2]) )
urc = llc + Vector( float(allparms["-dsmlength"][0]), float(allparms["-dsmlength"][1]), float(allparms["-dsmlength"][2]) )
res = Vector( float(allparms["-dsmcell"][0]), float(allparms["-dsmcell"][1]), float(allparms["-dsmcell"][2]) )
gridb = makeGridBox( llc, urc, res )
densityGrid = makeGrid( gridb, 0.0, int(allparms["-dsmpartition"]) )

print( str(nx(gridb)) + " " + str(ny(gridb)) + " " + str(nz(gridb)) )

stamp( densityGrid, volumeObject, 1 )


ambient = allparms["-ambient"]
ambientColor = constant( Color( float(ambient[0]), float( ambient[1] ), float( ambient[2] ), 1.0  ) ) 
colr = allparms["-color"]
litColor = constant( Color( float(colr[0]), float( colr[1]  ), float( colr[2] ), 1.0  ) )

volumeObject = gridded(densityGrid)
densityField = clamp( volumeObject/constant(clampv), 0, 1.0 )

#
#========================= End of scene definition =========================
#


StandardRender( densityField, allparms )

endJob()
