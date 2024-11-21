#!/usr/bin/python



from bishopUtils import *
from cuBishopUtils import *


beginJob()

csgparameters = {"-radius":2.5, "-center":[0.0, 0.0, 0.0],
		 "-clamp":0.1,
		 "-doturntable":0
		 }



parms = AllParameters( csgparameters )
CmdLineHelp("-h")
parms["-ds"] = 0.1
parms["-near"] = 4.0
parms["-maxpathlength"] = 10.0
parms["-fov"] = 65.0
parms["-scatter"] = 1.0
parms["-nbthreads"] = 3
parms["-dsmlength"] = [ 10, 10, 10 ]
parms["-dsmcell"] = [ 0.05, 0.05, 0.05 ]
parms["-dsmllc"] = [ -5, -5, -5 ]
parms["-dsmpartition"] = 10
parms["-nodsm"] = 0
parms["-color"] = [1, 1, 1]
parms["-ambient"] = [1, 1, 1]

parms["-lightstyle"] = "keyrimfill"

parms = CmdLineFindValues( parms )
PrintParameters( parms )


##### Sphere
sphereRadius = float(parms["-radius"])
sphereCenter = Vector( float(parms["-center"][0]),  float(parms["-center"][1]),  float(parms["-center"][2]) )
v1 = Sphere( sphereCenter, sphereRadius )

#parms["-case"] = "Mask" + str(parms["-ds"])
#v1 = mask(v1)

parms["-case"] = "Clamp" + str(parms["-ds"]) + "Scatter" + str(parms["-scatter"]) + "GriddedImplicit"
clampv = float( parms["-clamp"] )

llc = Vector( float(parms["-dsmllc"][0]), float(parms["-dsmllc"][1]), float(parms["-dsmllc"][2]) )
NND = Vector( float(parms["-dsmlength"][0]), float(parms["-dsmlength"][1]), float(parms["-dsmlength"][2]) )
res = Vector( float(parms["-dsmcell"][0]), float(parms["-dsmcell"][1]), float(parms["-dsmcell"][2]) )
gb = makeGridBox(llc, llc+NND, res)
density = makeGrid(gb, 0.0, int(parms["-dsmpartition"]))

print( str(nx(gb)) + " " + str(ny(gb)) + " " + str(nz(gb)) )
print(type(density))
stamp( density, v1, 1 )

print("### Begin CUDA render ###")
parms["-name"] = "cuda_test.exr"
initCUDA(density)
cuStandardRender( density, parms )
print("### End CUDA render ###")

print("### Begin CPU render ###")
parms["-name"] = "cpu_test.exr"
v1 = gridded(density)
v1 = clamp( v1/constant(clampv), 0.0, 1.0 )
StandardRender( v1, parms )
print("### End CPU render ###")

endJob()
