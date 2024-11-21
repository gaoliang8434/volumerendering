#!/usr/bin/python


import sys
import os
from bishoputils import *

beginJob()

# From 
# Back and forth error compensation and correction mehtods for semi-lagrangian schemes with application to level set interface computations
# TODD F. DUPONT AND YINGJIE LIU
def advectBFECC( field, velocity, timestep ):
	step1 = advect( field, velocity, timestep )
	step2 = advect( step1, velocity, -timestep )
	step3 = field * constant(1.5) - step2 * constant(0.5)
	step4 = advect( step3, velocity, timestep )
	return step4


# From 
# An Unconditionally Stable MacCormack Method
# Andrew Selle	Ronald Fedkiw	ByungMoon Kim Yingjie Liu	Jarek Rossignac
def advectModifiedMacCormack( field, velocity, timestep ):
	step1 = advect( field, velocity, timestep )
	step2 = advect( step1, velocity, -timestep )
	step3 = step1 + (field - step2) * constant(0.5)
	return step3

def PutInGrid( field, gb, defvalue ):
	grid = makeGrid( gb, defvalue )
	stamp( grid, field, 1 )
	return gridded( grid )

allparms = AllParameters( { "-nbadvectioniterations":0 } )

CmdLineHelp( "-h" )

allparms["-name"] = "gscfd.exr"
allparms["-scatter"] = 5.0
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-fov"] = 20.0
allparms["-eye"] = [0, 0, 20]
allparms["-view"] = [0, 0, 0]
allparms["-up"] = [0, 1, 0]
allparms["-ds"] = 0.001
allparms["-near"] = 10.0
allparms["-maxpathlength"] = 15.0
allparms["-dsmcell"] = [ 0.02, 0.02, 0.02 ] 
allparms["-dsmlength"] = [ 15, 15, 15 ]
allparms["-dsmllc"] = [ -7.5, -7.5, -7.5 ]
allparms["-case"] = "CompositionTestA"
allparms["-version"] = 1
allparms["-subsamples"] = 1
allparms["-size"] = [1920, 1080]
allparms["-aspect"] = 1920.0/1080.0
allparms["-nbthreads"] = 4

allparms = CmdLineFindValues( allparms )

noiseParms1 = Noise_t()
noiseParms1.octaves = 1.0
noiseParms1.freq = 0.331839
noiseParms1.translate = Vector(0,0,0)
perlin1 = perlin( noiseParms1 )
noisefield1 = SFNoise( perlin1, 0.01 )


noiseParms2 = Noise_t()
noiseParms2.octaves = 1.0
noiseParms2.freq = 0.25865758
noiseParms2.translate = Vector(-0.3,1.5,0.6)
perlin2 = perlin( noiseParms2 )
noisefield2 = SFNoise( perlin2, 0.01 )


nrange = 4
dilationScale = math.pow(2, nrange)

velocityField = cross( grad(noisefield1), grad(noisefield2) ) 




baseT = (1.0/24.0)*0.5*dilationScale

L = 3.0 
dL = 0.01
gridbox = makeGridBox( Vector(-L,-L,-L), Vector(L,L,L), Vector( dL, dL, dL ) )
velocityField = FFTDivFree( gridbox, velocityField )

baseVolume =  clamp( Sphere( Vector(0,0,0), 1.0 )/constant(0.1), 0.0, 1.0 )
cfdmap = identity()

T = baseT
cfdmap = PutInGrid(  gradientStretchCM( velocityField, T, 10 ) - identity(), gridbox, Vector(0.0,0.0,0.0) ) + identity()
density = warp(baseVolume, cfdmap )
StandardRender( density, allparms )


allparms["-case"] = "CompositionTestB"


nrange = 8
dilationScale = math.pow(2, nrange)
T = baseT/dilationScale
cfdmap = PutInGrid(  gradientStretchCM( velocityField, T, 1 ) - identity(), gridbox, Vector(0.0,0.0,0.0) ) + identity()
for i in range(0,nrange):
	cfdmap = PutInGrid( warp(cfdmap, cfdmap ) - identity(), gridbox, Vector(0,0,0) ) + identity()
density = warp(baseVolume, cfdmap )
StandardRender( density, allparms )








endJob()
