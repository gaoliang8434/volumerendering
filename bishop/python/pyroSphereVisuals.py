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
allparms["-freq"] = 1.124
allparms["-octaves"] = 1.0
allparms["-roughness"] = 0.5
allparms["-amp"] = 1.0
allparms["-pyrogamma"] = 0.5

allparms = CmdLineFindValues( allparms )
PrintParameters( allparms )

noiseparms = Noise_t()
noiseparms.frequency = float(allparms["-freq"])
noiseparms.roughness = float(allparms["-roughness"] )
noiseparms.amplitude = float(allparms["-amp"])
noiseparms.octaves = float( allparms["-octaves"])

llc = Vector(-3,-3,-3)
urc = Vector(3,3,3)
cellSize = float(allparms["-ds"])
vdbgb = makeGridBox( llc, urc, Vector( cellSize, cellSize, cellSize) )

volumeObject = Sphere( Vector(0,0,0), 2.0 )
#writeOpenVDB( "csgbox.vdb", -volumeObject, vdbgb )
noise = perlin( noiseparms )
Y = ImplicitSurfacePoint( volumeObject, 0.2, 15 )
von =  volumeObject + warp( abs(SFNoise(noise)), Y )
#writeOpenVDB( "pyrospherepos.vdb", -von, vdbgb )





von =  volumeObject - warp( abs(SFNoise(noise)), Y )
#writeOpenVDB( "pyrosphereneg.vdb", -von, vdbgb )




von =  volumeObject + warp( SFNoise(noise), Y )
#writeOpenVDB( "pyrosphereposneg.vdb", -von, vdbgb )













allparms["-pyrogamma"] = 0.2
volumeObject = PyroSphere(allparms)
writeOpenVDB( "pyrospheregama0.2.vdb", -volumeObject, vdbgb )


allparms["-pyrogamma"] = 0.5
volumeObject = PyroSphere(allparms)
writeOpenVDB( "pyrospheregama0.5.vdb", -volumeObject, vdbgb )


allparms["-pyrogamma"] = 1.0
volumeObject = PyroSphere(allparms)
writeOpenVDB( "pyrospheregama1.0.vdb", -volumeObject, vdbgb )


allparms["-pyrogamma"] = 1.5
volumeObject = PyroSphere(allparms)
writeOpenVDB( "pyrospheregama1.5.vdb", -volumeObject, vdbgb )

allparms["-pyrogamma"] = 2.0
volumeObject = PyroSphere(allparms)
writeOpenVDB( "pyrospheregama2.0.vdb", -volumeObject, vdbgb )











#
#========================= End of scene definition =========================
#



endJob()
