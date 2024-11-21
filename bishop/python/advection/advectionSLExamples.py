#!/usr/bin/python


import sys
import os
from bishoputils import *



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
	print "Generating gridded version of field"
	grid = makeGrid( gb, defvalue )
	stamp( grid, field, 1 )
	print "DONE Generating gridded version of field"
	return gridded( grid )

allparms = AllParameters( { "-nbadvectioniterations":0, "-advectiontime":0.15 } )

CmdLineHelp( "-h" )

allparms["-name"] = "advectionNoise.exr"
allparms["-scatter"] = 5.0
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-fov"] = 30.0
allparms["-ds"] = 0.003
allparms["-near"] = 3.0
allparms["-maxpathlength"] = 14.0
allparms["-dsmcell"] = [ 0.01, 0.01, 0.01 ] 
allparms["-dsmlength"] = [ 14, 14, 14 ]
allparms["-dsmllc"] = [ -7, -7, -7 ]
allparms["-case"] = ""
allparms["-version"] = 6
allparms["-subsamples"] = 10
allparms["-size"] = [1920, 1080]
allparms["-aspect"] = 1920.0/1080.0
allparms["-nbthreads"] = 4

allparms = CmdLineFindValues( allparms )

noiseParms1 = Noise_t()
noiseParms1.octaves = 2.0
noiseParms1.freq = 2.0
noiseParms1.translate = Vector(0,0,0)
perlin1 = perlin( noiseParms1 )
noisefield1 = SFNoise( perlin1 )


noiseParms2 = Noise_t()
noiseParms2.octaves = 2.0
noiseParms2.freq = 1.5
noiseParms2.translate = Vector(-0.3,1.5,0.6)
perlin2 = perlin( noiseParms2 )
noisefield2 = SFNoise( perlin2 )

velocityField = cross( grad(noisefield1), grad(noisefield2) )


L = 7.0 
dL = 0.05
gridbox = makeGridBox( Vector(-L,-L,-L), Vector(L,L,L), Vector( dL, dL, dL ) )


velocityField = FFTDivFree( gridbox, velocityField )


baseVolume =  clamp( cutout( Sphere( Vector(0,0,0), 1.0 ),   Sphere( Vector(1,0,0), 0.7 )    )/constant(0.1), 0.0, 1.0 )
baseMap = identity()
#fields = [ PutInGrid(baseVolume,gridbox,0.0), PutInGrid(baseVolume,gridbox,0.0), PutInGrid(baseVolume,gridbox,0.0) ]

hdisplacement = Vector( 3.0, 0.0, 0.0 )

T = float(allparms["-advectiontime"])
dL = 0.02
gridbox = makeGridBox( Vector(-L,-L,-L), Vector(L,L,L), Vector( dL, dL, dL ) )
advmap = PutInGrid(  gradientStretchCM( velocityField, T, 5 ) - identity(), gridbox, Vector(0.0,0.0,0.0) ) + identity()

density = warp( baseVolume, advmap ) 
allparms["-case"] = "GS"
StandardRender( density, allparms )


