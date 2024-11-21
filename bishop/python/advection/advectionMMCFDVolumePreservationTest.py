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

allparms["-name"] = "bfecccfd.exr"
allparms["-scatter"] = 5.0
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-fov"] = 20.0
allparms["-eye"] = [0, 0, 20]
allparms["-view"] = [0, 0, 0]
allparms["-up"] = [0, 1, 0]
allparms["-ds"] = 0.05
allparms["-near"] = 10.0
allparms["-maxpathlength"] = 15.0
allparms["-dsmcell"] = [ 0.02, 0.02, 0.02 ] 
allparms["-dsmlength"] = [ 15, 15, 15 ]
allparms["-dsmllc"] = [ -7.5, -7.5, -7.5 ]
allparms["-case"] = "VolumePreservationTest"
allparms["-version"] = 1
allparms["-subsamples"] = 1
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

nrange = 4 
dilationScale = math.pow(2, nrange)

velocityField = cross( grad(noisefield1), grad(noisefield2) ) 


velocityLeVeque = component(  constant(2.0)*pow( sin( constant(3.14159265)*xIdentity())  , 2.0 ) * sin( constant(2.0*3.14159265)*yIdentity() ) * sin( constant(2.0*3.14159265)*zIdentity() ),    - pow( sin( constant(3.14159265)*yIdentity())  , 2.0 ) * sin( constant(2.0*3.14159265)*xIdentity() ) * sin( constant(2.0*3.14159265)*zIdentity() ),     - pow( sin( constant(3.14159265)*zIdentity())  , 2.0 ) * sin( constant(2.0*3.14159265)*xIdentity() ) * sin( constant(2.0*3.14159265)*yIdentity() )        )

baseT = (1.0/24.0)*0.5

L = 1.0 
dL = 0.01
gridbox = makeGridBox( Vector(-L,-L,-L), Vector(L,L,L), Vector( dL, dL, dL ) )
#velocityField = FFTDivFree( gridbox, velocityField )

baseVolume =  clamp( Sphere( Vector(0,0,0), 1.0 )/constant(0.1), 0.0, 1.0 )

#T = baseT*dilationScale
#cfdmap = gradientStretchCM( velocityField, T, 500 )



T = baseT*dilationScale
nrange = 4
DT = T / math.pow(2.0,nrange)
print "# Modified MacCormack T = " + str(T) + "   Number of folds = " + str(nrange) + "    DT = " + str(DT)
cfdmap = advectModifiedMacCormack( identity(),  velocityLeVeque, DT )
#cfdmap = advectModifiedMacCormack( identity(),  velocityField, DT )
for n in range(0,nrange):
	cfdmap = warp(cfdmap,cfdmap)

volumepres = det( grad(cfdmap) )



x = -L
y = 0.1
z = 0.1
while x <= L:
	value = evaluate( volumepres, Vector( x, y, z ) )
	print str(x) + " " + str(value)
	x = x + dL

#density = warp(baseVolume, cfdmap )
#StandardRender( volumepres, allparms )
#StandardTwoPartRender( volumepres, constant(0.0), constant(Color(1,1,1,1)), constant( Color(0,0,0,0)), allparms )

#allparms["-case"] = "CompositionTestB"


#T = baseT
#cfdmap = PutInGrid(  advect( identity(), velocityField, T ) - identity(), gridbox, Vector(0.0,0.0,0.0) ) + identity()
#for i in range(0,nrange):
#	cfdmap = PutInGrid( warp(cfdmap, cfdmap ) - identity(), gridbox, Vector(0,0,0) ) + identity()
#density = warp(baseVolume, cfdmap )
#StandardRender( density, allparms )








endJob()
