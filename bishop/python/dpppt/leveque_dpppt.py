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

allparms["-name"] = "gsdpppt_leveque.exr"
allparms["-scatter"] = 30.0
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-fov"] = 45.0
allparms["-eye"] = [0.5, -2, 0.5]
allparms["-view"] = [0.5, 0.5, 0.5]
allparms["-up"] = [0, 0, 1]
allparms["-ds"] = 1.0/2000.0
allparms["-near"] = 0.1
allparms["-maxpathlength"] = 100.0
allparms["-dsmcell"] = [ 0.0025, 0.0025, 0.0025 ] 
allparms["-dsmlength"] = [ 1, 1, 1 ]
allparms["-dsmllc"] = [ 0,0,0 ]
allparms["-renderllc"] = [ 0,0,0 ]
allparms["-renderurc"] = [ 1,1,1 ]
allparms["-case"] = "LeVequeTest"
allparms["-version"] = 7
allparms["-subsamples"] = 1
allparms["-size"] = [1920, 1080]
allparms["-aspect"] = 1920.0/1080.0
allparms["-nbthreads"] = 8
allparms["-lightCd"] = [ [1,0,0], [0,1,0], [0,0,1] ]
allparms["-lightP"] = [ [10,0,0], [0,-10,10], [0,10,1] ]
allparms["-lightstyle"] = "custom"

allparms = CmdLineFindValues( allparms )


velocityLeVeque = component(  
                              constant(2.0)*pow( sin( constant(3.14159265)*xIdentity() )  , 2.0 ) * sin( constant(2.0*3.14159265)*yIdentity() ) * sin( constant(2.0*3.14159265)*zIdentity() ), 
			                  - pow( sin( constant(3.14159265)*yIdentity() )  , 2.0 ) * sin( constant(2.0*3.14159265)*xIdentity() ) * sin( constant(2.0*3.14159265)*zIdentity() ),
			                  - pow( sin( constant(3.14159265)*zIdentity() )  , 2.0 ) * sin( constant(2.0*3.14159265)*xIdentity() ) * sin( constant(2.0*3.14159265)*yIdentity() )        )



baseT = 3.0/150.0

L = 1.0
dL = 0.01
gridbox = makeGridBox( Vector(0,0,0), Vector(L,L,L), Vector( dL, dL, dL ) )



baseVolume =  clamp( Sphere( Vector(0.35,0.35,0.35), 0.15 )/constant(0.01), 0.0, 1.0 )


T = 0
nrange = 0
DT = baseT / math.pow(2.0,nrange)
gsSteps = 6 
print "# Gradient Stretch T = " + str(T) + "   Number of folds = " + str(nrange) + "    DT = " + str(DT) + "   GS Steps = " + str(gsSteps)


Omega = constant(Vector(0,0,0))

boundingSphere = -Sphere( Vector(L/2.0,L/2.0,L/2.0), 1.0 )






cfdmap = identity()
Tperiod = 3.0
nframes = int( Tperiod/DT ) + 2
DT = Tperiod / float(nframes-1)
print "total number of frames = " + str(nframes-1)
for frame in range(1,nframes):
	allparms["-turntableframe"] = frame
	PrintParameters(allparms)
	density = warp(baseVolume, cfdmap )
	StandardRender( density, allparms )
	T = T + DT
	print "Frame " + str(frame) + "  T " + str(T)
	velocityLeVequeT = velocityLeVeque * constant( math.cos( 3.14159265 * T / Tperiod ) )
	DDT = DT/math.pow(2,4)
	stepmap = gradientStretchCM( velocityLeVequeT, DT, gsSteps )
	#for i in range(0,4):
	#	stepmap = warp(stepmap, stepmap )
	stepmap = PutInGrid( stepmap-identity(), gridbox, Vector(0,0,0) ) + identity()
        vorticity = curl(velocityLeVequeT)
        Omega = warp(Omega,stepmap) + vorticity * constant(DT*0.5)
        Omega = PutInGrid( Omega, gridbox, Vector(0,0,0) )
        gradX = rotation(Omega)
	XX = gradientDisplacement( gradX, boundingSphere, dL/10.0 ) + identity()
	cfdmap = PutInGrid( warp(cfdmap,XX) - identity(), gridbox, Vector(0,0,0) ) + identity()






endJob()
