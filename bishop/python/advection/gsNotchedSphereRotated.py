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

allparms["-name"] = "gsRotated.exr"
allparms["-scatter"] = 5.0
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-fov"] = 50.0
allparms["-eye"] = [0, 0, 200]
allparms["-view"] = [0, 0, 0]
allparms["-up"] = [0, 1, 0]
allparms["-ds"] = 0.05
allparms["-near"] = 100.0
allparms["-maxpathlength"] = 200.0
allparms["-dsmcell"] = [ 0.05, 0.05, 0.05 ] 
allparms["-dsmlength"] = [ 100, 100, 100 ]
allparms["-dsmllc"] = [ -50,-50,-50 ]
allparms["-case"] = "ZalesakSphere"
allparms["-version"] = 1
allparms["-subsamples"] = 1
allparms["-size"] = [1920, 1080]
allparms["-aspect"] = 1920.0/1080.0
allparms["-nbthreads"] = 4

allparms = CmdLineFindValues( allparms )


nrange = 0
dilationScale = math.pow(2, nrange)

velocityField = component( constant(1.0/10.0)*( constant(00.0)-yIdentity() ) , constant(1.0/10.0)*( xIdentity() - constant(00.0) ) , constant(0.0) )



baseT = 2.0*3.14159265

L = 100 
dL = 0.5
gridbox = makeGridBox( Vector(-L,-L,-L), Vector(L,L,L), Vector( dL, dL, dL ) )
periodicgridbox = makeGridBox( Vector(-L,-L,-L), Vector(L,L,L), Vector( dL, dL, dL ) )
setPeriodic(periodicgridbox)

#velocityField = PutInGrid( velocityField, periodicgridbox, Vector(0,0,0) )

baseVolume =  Sphere( Vector(0,25,0), 15.0 )
notch = HardBox( Vector( -2.5,10,-L), Vector( 2.5, 22.5, L)  )
baseVolume = cutout( baseVolume, notch )
density0 = clamp( baseVolume/constant(1.0), 0, 1.0 )


#for j in range(0,400):
#	x = Vector( -L + dL*j, 0,0 )
#	vf = evaluate( velocityField, x )
#	gvf = evaluate( grad(velocityField), x )
#	print str( x.X() ) + "  " + str(vf) + "      [ " + str(gvf.Get(0,0)) + ", " + str(gvf.Get(1,0)) + ", " + str(gvf.Get(2,0)) + "         " + str(gvf.Get(0,1)) + ", " + str(gvf.Get(1,1)) + ", " + str(gvf.Get(2,1)) + "         " + str(gvf.Get(0,2)) + ", " + str(gvf.Get(1,2)) + ", " + str(gvf.Get(2,2)) + " ]        "



cfdmap = identity()
for  i in range(1,41+1,1):
	allparms["-turntableframe"] = i
	T = baseT * (i-1)
	stepmap =  PutInGrid( gradientStretchCM( velocityField, T, 1 ) - identity(), gridbox, Vector(0,0,0) ) + identity()
	density = warp(density0, stepmap )
	PrintParameters( allparms )
	StandardRender( density, allparms )



endJob()
