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
allparms["-dsmNxNyNz"] = [ 300, 300, 300 ]
allparms["-dsmllc"] = [ -4,-4,-4 ]
allparms["-dsmurc"] = [ 4,4,4 ]
allparms["-dsmlength"] = [ 8,8,8 ]
allparms["-dsmcell"] = [ 0.01,0.01,0.01 ]
allparms["-dsmdxdydz"] = [ 0.02,0.02,0.02 ]
allparms["-fov"] = 80.0
allparms["-ds"] = 0.015
allparms["-near"] = 0.1 
allparms["-maxpathlength"] = 100.0
allparms["-freq"] = 1.124*0.3
allparms["-octaves"] = 3.0
allparms["-roughness"] = 0.75
allparms["-amp"] = 1.0
allparms["-pyrogamma"] = 0.75

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

allparms["-renderllc"] = [ llc.X(), llc.Y(), llc.Z()   ]
allparms["-renderurc"] = [ urc.X(), urc.Y(), urc.Z()   ]





radius = 2.0
volumeObject = Sphere( Vector(0,0,0), radius )
noise = perlin( noiseparms )












#writeOpenVDB( "sphere.vdb", -volumeObject, vdbgb )

allparms["-name"] = "original_sphere.exr"
densityField = clamp( volumeObject/constant(0.1) , 0, 1 )
#StandardRender( densityField, allparms )

iX = unitvector( grad(volumeObject) ) * constant(radius)

noise_disp = pow( abs( warp(SFNoise(noise), iX) ), float(allparms["-pyrogamma"])   )*grad(volumeObject)*constant(0.75)

X = identity() + noise_disp

sphere_disp = warp( volumeObject, X )
#writeOpenVDB( "sphere_disp.vdb", -sphere_disp, vdbgb )

allparms["-name"] = "displaced_sphere.exr"
densityField = clamp( sphere_disp/constant(0.1) , 0, 1 )
#StandardRender( densityField, allparms )








volumeObject = Torus( Vector(0,0,0), Vector( 1,1,1), radius, 0.2*radius )

allparms["-name"] = "original_torus.exr"
densityField = clamp( volumeObject/constant(0.1) , 0, 1 )
#StandardRender( densityField, allparms )

iX = unitvector( grad(volumeObject) ) * constant(radius)

noise_disp = abs( warp(SFNoise(noise), iX)) * unitvector(grad(volumeObject))*constant(0.75)

X = identity() + noise_disp

sphere_disp = warp( volumeObject, X )
#writeOpenVDB( "torus_disp.vdb", -sphere_disp, vdbgb )

allparms["-name"] = "displaced_torus.exr"
densityField = clamp( sphere_disp/constant(0.1) , 0, 1 )
#StandardRender( densityField, allparms )













volumeObject = Torus( Vector(0,0,0), Vector( 1,1,1), radius, 0.2*radius )

allparms["-name"] = "original_torus.exr"
densityField = clamp( volumeObject/constant(0.1) , 0, 1 )
#StandardRender( densityField, allparms )

iX = ImplicitSurfacePoint( volumeObject, 0.1*radius, 10 )

noise_disp = abs( warp(SFNoise(noise), iX)) * unitvector(grad(volumeObject))*constant(0.75)

X = identity() + noise_disp

sphere_disp = warp( volumeObject, X )
writeOpenVDB( "torus_IPT_disp.vdb", -sphere_disp, vdbgb )

allparms["-name"] = "displaced_IPT_torus.exr"
densityField = clamp( sphere_disp/constant(0.1) , 0, 1 )
#StandardRender( densityField, allparms )





#
#========================= End of scene definition =========================
#



endJob()
