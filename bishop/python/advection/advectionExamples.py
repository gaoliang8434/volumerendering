#!/usr/bin/python


import sys
import os
from vr import *
from vrutils import *



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

allparms = AllParameters( { "-nbadvectioniterations":1 } )

CmdLineHelp( "-h" )

allparms["-name"] = "advection.exr"
allparms["-scatter"] = 2.0
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-fov"] = 60.0
allparms["-ds"] = 0.01
allparms["-near"] = 3.0
allparms["-maxpathlength"] = 14.0
allparms["-dsmcell"] = [ 0.02, 0.02, 0.02 ] 
allparms["-dsmlength"] = [ 10, 10, 10 ]
allparms["-dsmllc"] = [ -5, -5, -5 ]
allparms["-case"] = "Comparison"
allparms["-version"] = 1
allparms["-subsamples"] = 10
allparms["-size"] = [1920, 540]
allparms["-aspect"] = 1920.0/540.0

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

L = 3.0 
dL = 0.05
gridbox = makeGridBox( Vector(-L,-L,-L), Vector(L,L,L), Vector( dL, dL, dL ) )

baseVolume =  clamp( Sphere( Vector(0,0,0), 1.0 )/constant(0.1), 0.0, 1.0 )
baseMap = identity()
fields = [ PutInGrid(baseVolume,gridbox,0.0), PutInGrid(baseVolume,gridbox,0.0), PutInGrid(baseVolume,gridbox,0.0) ]
maps = [ identity(), identity(), identity() ]


T = 0.2
nbIterations = int(allparms["-nbadvectioniterations"])
dt = T/float(nbIterations)


for i in range(0,nbIterations):
	print "Iteration " + str(i)
	fields[0] = PutInGrid( advect( fields[0], velocityField, dt ), gridbox, 0.0 )
	fields[1] = PutInGrid( advectBFECC( fields[1], velocityField, dt ), gridbox, 0.0 )
	fields[2] = PutInGrid( advectModifiedMacCormack( fields[2], velocityField, dt ), gridbox, 0.0 )
	maps[0] = PutInGrid( advect( maps[0], velocityField, dt ) - identity(), gridbox, Vector(0.0,0.0,0.0) ) + identity()
	maps[1] = PutInGrid( advectBFECC( maps[1], velocityField, dt ) - identity(), gridbox, Vector(0.0,0.0,0.0) ) + identity()
	maps[2] = PutInGrid( advectModifiedMacCormack( maps[2], velocityField, dt ) - identity(), gridbox, Vector(0.0,0.0,0.0) ) + identity()




# Displace and Render

hdisplacement = Vector( 3.0, 0.0, 0.0 )

density = constant(0.0)
for i in range(0, len(fields) ):
	density = density + translate( fields[i], hdisplacement * float( i - int(len(fields)/2.0) )  )
allparms["-case"] = "GriddedFieldsSteps" + str(nbIterations)
#StandardRender( density, allparms )

allparms["-ds"] = 0.0015

density = constant(0.0)
for i in range(0, len(fields) ):
	density = density + translate( warp( baseVolume, maps[i] ), hdisplacement * float( i - int(len(fields)/2.0) ) ) 
allparms["-case"] = "GriddedMapsSteps" + str(nbIterations)
StandardRender( density, allparms )


