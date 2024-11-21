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

allparms = AllParameters( { "-nbadvectioniterations":0, "-gsnbintervals":10, "-gsgridres":0.01, "-interpolations":1 } )

CmdLineHelp( "-h" )

allparms["-name"] = "sl_interp.exr"
allparms["-scatter"] = 20.0
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-fov"] = 20.0
allparms["-eye"] = [0.35, 0.35, -10]
allparms["-view"] = [0.35, 0.35, 0.35]
allparms["-up"] = [0, 1, 0]
allparms["-ds"] = 0.005
allparms["-near"] = 0.1
allparms["-maxpathlength"] = 100.0
allparms["-dsmcell"] = [ 0.005, 0.005, 0.005 ] 
allparms["-dsmlength"] = [ 1, 1, 1 ]
allparms["-dsmllc"] = [ 0,0,0 ]
allparms["-renderllc"] = [ 0,0,0 ]
allparms["-renderurc"] = [ 1,1,1 ]
allparms["-case"] = "ShearTest"
allparms["-version"] = 3
allparms["-subsamples"] = 1
allparms["-size"] = [1920, 1080]
allparms["-aspect"] = 1920.0/1080.0
allparms["-nbthreads"] = 4
allparms["-gsnbintervals"] = 1 
allparms["-lightstyle"] = "keyrimfill" 

allparms = CmdLineFindValues( allparms )


# This test as described in
#
#a"A 3D unsplit-advection volume tracking algorithm with planarity-preserving interface reconstruction", Liovic, Rudman, Liow, Lakehal, and Kothe, Computers and Fluids, 35, (2006), 1011-1032.
#  doi:10.1016/j.compfluid.2005.09.003
#

interpolations = int( allparms["-interpolations"])
allparms["-name"] = "sl_interp_" + str(interpolations) + ".exr"

umax = 1.0
R = 0.5
r = pow( pow( xIdentity()-constant(0.5)   , 2 ) + pow( yIdentity() - constant(0.5), 2 ) , 0.5 ) 
velocityShear = component(  
                              pow( sin( constant(3.14159265)*xIdentity())  , 2.0 ) * sin( constant(2.0*3.14159265)*yIdentity() ), 
			                  - pow( sin( constant(3.14159265)*yIdentity())  , 2.0 ) * sin( constant(2.0*3.14159265)*xIdentity() ),
			                  constant(umax)*pow( constant(1.0) - r/constant(R), 2.0 )        )



L = 1.0 
dL = float( allparms["-gsgridres"] )

baseT = (3.0/150.0)*3.0

gridbox = makeGridBox( Vector(0,0,0), Vector(L,L,L), Vector( dL, dL, dL ) )
Xgridbox = makeGridBox( Vector(0,0,0), Vector(L,L,L), Vector( dL, dL, dL )*0.5 )

setInterpOrder( gridbox, interpolations )
setInterpOrder( Xgridbox, interpolations )


startsphere =  Sphere( Vector(0.5,0.75,0.25), 0.15 )
baseVolume =  clamp( startsphere/constant(0.01), 0.0, 1.0 )
allparms["-renderllc"] = [0,0,0]
allparms["-renderurc"] = [L,L,L]

Tperiod = 3.0
T = 0
nrange = 0
DT = baseT / math.pow(2.0,nrange)
gsSteps = int( allparms["-gsnbintervals"] )
print "# Gradient Stretch T = " + str(T) + "   Number of folds = " + str(nrange) + "    DT = " + str(DT) + "   GS Steps = " + str(gsSteps)

cfdmap = identity()


nbframes = int( Tperiod / DT ) + 1


for frame in range(1,nbframes):
	allparms["-turntableframe"] = frame
	density = warp(baseVolume, cfdmap )
	if frame==1:
		PrintParameters(allparms)
		StandardRender( density, allparms )
	if frame==int(nbframes/2):
		PrintParameters(allparms)
		StandardRender( density, allparms )
	if frame==nbframes-1:
		PrintParameters(allparms)
		StandardRender( density, allparms )
	T = T + DT
	velocity = velocityShear * constant( math.cos( 3.14159265*T/Tperiod ) )
	#DDT = DT / math.pow(2,8)
	stepmap = advect( identity(), velocity, DT )
	#for n in range(0,8):
	#	stepmap = warp(stepmap,stepmap)
	cfdmap = PutInGrid( warp(cfdmap,stepmap) - identity(), Xgridbox, Vector(0,0,0) ) + identity()
	#volumepres = det( grad(stepmap) )
	#x = -L
	#y = 0.1
	#z = 0.1
	#while x <= L:
	#	value = evaluate( volumepres, Vector( x, y, z ) )
	#	print str(x) + " " + str(value)
	#	x = x + dL


#finalsphere = warp( startsphere, cfdmap )
#error = Union( cutout( startsphere, finalsphere ) , cutout( finalsphere, startsphere ) )
#errordensity = mask(error)
#allparms["-case"] = "ShearTestError"
#StandardRender( errordensity, allparms )




endJob()
