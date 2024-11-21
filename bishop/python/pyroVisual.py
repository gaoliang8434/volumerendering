#!/usr/bin/python


import sys
import os
from bishoputils import *
from ImplicitShapes import *


beginJob()

#
#=========== pyro construction ============================
#

allparms =  AllParameters(PyroSphereParameters())

CmdLineHelp( "-h" )


allparms["-pyrogamma"] = 0.2
allparms["-amp"] = 1.5
allparms["-octaves"] = 6.0
allparms["-roughness"] = 0.5
allparms["-name"] = "pyro.exr"
allparms["-scatter"] = 2.5
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-size"] = [ 960, 540 ]
allparms["-dsmNxNyNz"] = [ 100, 100, 100 ]
allparms["-dsmllc"] = [ -4,-4,-4 ]
allparms["-dsmurc"] = [ 4,4,4 ]
allparms["-dsmdxdydz"] = [ 0.01,0.01,0.01 ]
allparms["-fov"] = 80.0
allparms["-ds"] = 0.025
allparms["-near"] = 0.1 
allparms["-maxpathlength"] = 100.0
allparms["-renderllc"] = [-4,-4,-4]
allparms["-renderurc"] = [4,4,4]


allparms = CmdLineFindValues( allparms )
PrintParameters( allparms )

print "Building Pyro" 
volumeObject = PyroSphere(allparms)


llc = Vector(-4,-4,-4)
urc = Vector(4,4,4)
cellSize = float(allparms["-ds"])
vdbgb = makeGridBox( llc, urc, Vector( cellSize, cellSize, cellSize) )
writeOpenVDB( "pyroclasticSphere.vdb", -volumeObject, vdbgb )



#
#========================= End of scene definition =========================
#



endJob()
