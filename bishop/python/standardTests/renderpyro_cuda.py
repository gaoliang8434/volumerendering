#!/usr/bin/python


import sys
import os
from bishopUtils import *
from ImplicitShapes import *
from cuBishopUtils import *


beginJob()

#
#=========== pyro construction ============================
#

allparms =  AllParameters(PyroSphereParameters())
CmdLineHelp( "-h" )


allparms["-pyrogamma"] = 0.8
allparms["-octaves"] = 10.0
allparms["-roughness"] = 0.5
allparms["-translate"] = [ 0, 0, 0 ]
allparms["-radius"] = 2.25
allparms["-amp"] = 1.5
allparms["-freq"] = 1.5

allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-color"] = [1, 1, 1]
allparms["-fov"] = 70.0
allparms["-near"] = 4.0

allparms["-maxpathlength"] = 10.0
allparms["-eye"] = [ 0, 0, 15 ]
allparms["-view"] = [ 0, 0, -1 ]
allparms["-clamp"] = 0.2

allparms["-dsmllc"] = [-5, -5, -5]
allparms["-dsmlength"] = [10, 10, 10]
allparms["-dsmcell"] = [ 0.025, 0.025, 0.025 ]
allparms["-dsmpartition"] = 10

allparms["-renderllc"] = [-5, -5, -5]
allparms["-renderurc"] = [5, 5, 5]

allparms["-lightP"] = [[15, 90, 15], [0, -90, 0], [0, 0, -90]]
allparms["-lightCd"] = [[0, 0, 1], [1, 0, 0], [0, 1, 0]]
allparms["-lightstyle"] = "custom"

allparms["-scatter"] = 9.0
allparms["-subsamples"] = 2
allparms["-nbthreads"] = 3
allparms["-nbpatches"] = 100
allparms["-ds"] = 0.005


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

density = clamp( volumeObject, 0.0, 10000.0 )
stamp( densityGrid, density, 1 )


ambient = allparms["-ambient"]
ambientColor = constant( Color( float(ambient[0]), float( ambient[1] ), float( ambient[2] ), 1.0  ) ) 
colr = allparms["-color"]
litColor = constant( Color( float(colr[0]), float( colr[1]  ), float( colr[2] ), 1.0  ) )

volumeObject = gridded(densityGrid)
densityField = clamp( volumeObject/constant(clampv), 0, 1.0 )

#
#========================= End of scene definition =========================
#

allparms["-name"] = "cuda_pyro_AABB.exr"
print("Rendering with CUDA...")
initCUDA(densityGrid)
cuStandardRender( densityGrid, allparms )

# allparms["-name"] = "cpu_pyro.exr"
# print("Rendering with CPU...")
# StandardRender( densityField, allparms )

endJob()
