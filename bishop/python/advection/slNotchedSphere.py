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

allparms["-name"] = "sl.exr"
allparms["-scatter"] = 5.0
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-fov"] = 50.0
allparms["-eye"] = [0, 0, 200]
allparms["-view"] = [0, 0, 0]
allparms["-up"] = [0, 1, 0]
allparms["-ds"] =0.1 
allparms["-near"] = 100.0
allparms["-maxpathlength"] = 200.0
allparms["-dsmcell"] = [ 0.1,0.1,0.1 ] 
allparms["-dsmlength"] = [ 200, 200, 200 ]
allparms["-dsmllc"] = [ -100,-100,-100 ]
allparms["-case"] = "ZalesakSphere"
allparms["-version"] = 1
allparms["-subsamples"] = 1
allparms["-size"] = [1920, 1080]
allparms["-aspect"] = 1920.0/1080.0
allparms["-nbthreads"] = 4
allparms["-lightstyle"] = "keyrimfill"

allparms = CmdLineFindValues( allparms )


nrange = 0
dilationScale = math.pow(2, nrange)

velocityField = component( constant(3.14159265/314.0)*( constant(00.0)-yIdentity() ) , constant(3.14159265/314.0)*( xIdentity() - constant(00.0) ) , constant(0.0) )



baseT = 6.28

L = 100 
dL = 1.0
gridbox = makeGridBox( Vector(-L,-L,-L), Vector(L,L,L), Vector( dL, dL, dL ) )
allparms["-renderllc"] = [-L,-L,-L]
allparms["-renderurc"] = [L,L,L]

#velocityField = PutInGrid( velocityField, gridbox, Vector(0,0,0) )

baseVolume =  Sphere( Vector(0,25,0), 15.0 )
notch = HardBox( Vector( -2.5,10,-L), Vector( 2.5, 22.5, L)  )
baseVolume = cutout( baseVolume, notch )
density0 = clamp( baseVolume/constant(1.0), 0, 1.0 )






cfdmap = identity()
T = baseT
stepmap = PutInGrid(  advect( identity(), velocityField, T ) - identity(), gridbox, Vector(0.0,0.0,0.0) ) + identity()
for  i in range(1,101+1):
	print "Frame " + str(i)
	allparms["-turntableframe"] = i
	density = warp(density0, cfdmap )
	if i==101:
		StandardRender( density, allparms )
	cfdmap = PutInGrid( warp(cfdmap, stepmap) -identity(), gridbox, Vector(0,0,0) ) + identity()



endJob()
