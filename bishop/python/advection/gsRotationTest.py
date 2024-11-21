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

allparms["-name"] = "gs.exr"
allparms["-scatter"] = 5.0
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-fov"] = 50.0
allparms["-eye"] = [0, 0, 200]
allparms["-view"] = [0, 0, 0]
allparms["-up"] = [0, 1, 0]
allparms["-ds"] = 0.5
allparms["-near"] = 100.0
allparms["-maxpathlength"] = 200.0
allparms["-dsmcell"] = [ 0.2, 0.2, 0.2 ] 
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

velocityField = component( constant(3.14159265/314.0)*( constant(00.0)-yIdentity() ) , constant(3.14159265/314.0)*( xIdentity() - constant(00.0) ) , constant(0.0) )

procGrad = grad(velocityField)


baseT = 6.28*10

L = 150 
dL = 1.0
gridbox = makeGridBox( Vector(-L,-L,-L), Vector(L,L,L), Vector( dL, dL, dL ) )



velocityField = PutInGrid( velocityField, gridbox, Vector(0,0,0) )

gridGrad = grad(velocityField)

for  x in range(-L,L,int(dL)):
	pos = Vector( x, 0, 0 )
	procvalue = evaluate( procGrad, pos )
	gridvalue = evaluate( gridGrad, pos )
	print "pos " + str(pos)
	print "proc  " + str(procvalue.Get(0,0) ) + " " + str(procvalue.Get(0,1) ) + " " + str(procvalue.Get(0,2) )
	print "      " + str(procvalue.Get(1,0) ) + " " + str(procvalue.Get(1,1) ) + " " + str(procvalue.Get(1,2) )
	print "      " + str(procvalue.Get(2,0) ) + " " + str(procvalue.Get(2,1) ) + " " + str(procvalue.Get(2,2) )
	print "grid  " + str(gridvalue.Get(0,0) ) + " " + str(gridvalue.Get(0,1) ) + " " + str(gridvalue.Get(0,2) )
	print "      " + str(gridvalue.Get(1,0) ) + " " + str(gridvalue.Get(1,1) ) + " " + str(gridvalue.Get(1,2) )
	print "      " + str(gridvalue.Get(2,0) ) + " " + str(gridvalue.Get(2,1) ) + " " + str(gridvalue.Get(2,2) )
	print " "	



baseVolume =  Sphere( Vector(0,25,0), 15.0 )
notch = HardBox( Vector( -2.5,10,-L), Vector( 2.5, 22.5, L)  )
baseVolume = cutout( baseVolume, notch )
density0 = clamp( baseVolume/constant(1.0), 0, 1.0 )



T = baseT
XField = identity() - velocityField * constant(T) * sinch( grad(velocityField) * constant(T) )


cfdmap = identity()
stepmap10 = gradientStretchCM( velocityField, T, 5000 ) 
stepmap01 = gradientStretchCM( velocityField, T, 1 ) 
#for  x in range(-L,L,int(dL)):
#	pos = Vector( x, 0, 0 )
#	print "pos " + str(pos)
#	ev01 = evaluate( stepmap01, pos )
#	ev10 = evaluate( stepmap10, pos )
#	diff = (ev01-ev10).magnitude()
#	print str(x) + "   " + str(ev01) + " " + str(ev10) + "   " + str(diff)
#	gu = evaluate( constant(T) * sinch( grad(velocityField)*constant(T)) , pos )
#	print "gu  " + str(gu.Get(0,0) ) + " " + str(gu.Get(0,1) ) + " " + str(gu.Get(0,2) )
#	print "    " + str(gu.Get(1,0) ) + " " + str(gu.Get(1,1) ) + " " + str(gu.Get(1,2) )
#	print "    " + str(gu.Get(2,0) ) + " " + str(gu.Get(2,1) ) + " " + str(gu.Get(2,2) )
#	vf = evaluate( velocityField, pos )
#	print "V " + str(vf)
#	gu = evaluate( grad(velocityField) , pos )
#	print "gV  " + str(gu.Get(0,0) ) + " " + str(gu.Get(0,1) ) + " " + str(gu.Get(0,2) )
#	print "    " + str(gu.Get(1,0) ) + " " + str(gu.Get(1,1) ) + " " + str(gu.Get(1,2) )
#	print "    " + str(gu.Get(2,0) ) + " " + str(gu.Get(2,1) ) + " " + str(gu.Get(2,2) )
#	gu = evaluate( sinch( grad(velocityField) * constant(T) ) , pos )
#	print "sg  " + str(gu.Get(0,0) ) + " " + str(gu.Get(0,1) ) + " " + str(gu.Get(0,2) )
#	print "    " + str(gu.Get(1,0) ) + " " + str(gu.Get(1,1) ) + " " + str(gu.Get(1,2) )
#	print "    " + str(gu.Get(2,0) ) + " " + str(gu.Get(2,1) ) + " " + str(gu.Get(2,2) )
#	X = evaluate( XField, pos )
#	print "X  " + str(X)
#	print "Magnitudes:  |pos| = " + str(pos.magnitude()) + "    |ev01| = " + str(ev01.magnitude()) + "    |ev10| = " + str(ev10.magnitude()) + "      |X| = " + str(X.magnitude())
#


endJob()
