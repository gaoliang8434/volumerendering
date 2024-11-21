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
allparms["-freq"] = 2.124
allparms["-octaves"] = 1.0
allparms["-roughness"] = 0.5
allparms["-amp"] = 1.0
allparms["-pyrogamma"] = 1.0

allparms = CmdLineFindValues( allparms )
PrintParameters( allparms )

noiseparms = Noise_t()
noiseparms.frequency = float(allparms["-freq"])
noiseparms.roughness = float(allparms["-roughness"] )
noiseparms.amplitude = float(allparms["-amp"])
noiseparms.octaves = float( allparms["-octaves"])

llc = Vector(-4,-4,-4)
urc = Vector(4,4,4)
cellSize = float(allparms["-ds"])
vdbgb = makeGridBox( llc, urc, Vector( cellSize, cellSize, cellSize) )

volumeObject = CsgBox( Vector(0,0,0), 1.5, 5.0 )
#writeOpenVDB( "csgbox.vdb", -volumeObject, vdbgb )
noise = perlin( noiseparms )
Y = ImplicitSurfacePoint( volumeObject, 0.2, 20 )
von =  volumeObject + warp( pow( abs(SFNoise(noise)), float(allparms["-pyrogamma"]))*constant(1000.0), Y )
#writeOpenVDB( "pyrocsgbox.vdb", -von, vdbgb )





noiseparms.frequency = float(allparms["-freq"]) * 0.37
volumeObject = Torus( Vector(0,0,0), Vector(0,0,1), 3.0, 0.4 )
#writeOpenVDB( "torus.vdb", -volumeObject, vdbgb )
noise = perlin( noiseparms )
Y = ImplicitSurfacePoint( volumeObject, 0.2, 20 )
von =  volumeObject + warp( pow( abs(SFNoise(noise)), float(allparms["-pyrogamma"]))*constant(40.0), Y )
#writeOpenVDB( "pyrotorus.vdb", -von, vdbgb )




llc = Vector(-7,-7,-7)
urc = Vector(7,7,7)
cellSize = float(allparms["-ds"])
vdbgb = makeGridBox( llc, urc, Vector( cellSize, cellSize, cellSize) )
noiseparms.frequency = 1.0
volumeObject = Icosahedron()
#writeOpenVDB( "icosahedron.vdb", -volumeObject, vdbgb )
noise = perlin( noiseparms )
Y = ImplicitSurfacePoint( volumeObject, 0.001, 100 )
von =  volumeObject + warp( pow( abs(SFNoise(noise)), float(allparms["-pyrogamma"]))*constant(1.0), Y )
writeOpenVDB( "pyroicosahedron.vdb", -von, vdbgb )

#
#========================= End of scene definition =========================
#



endJob()
