#!/usr/bin/python

from bishoputils import *

parms = AllParameters( {} )
parms = CmdLineFindValues( parms )
PrintParameters( parms )

X = identity()

gb = makeGridBox( Vector(-20,-20,-20), Vector(20,20,20), Vector(0.5,0.5,0.5) )

grid = makeGrid( gb, Vector(0,0,0) );

print "Starting surface convolution"
GreenSurfaceConvolve( X, grid );

Xconvolve = gridded(grid)

MM = constant( unitMatrix() )
mgrid = makeGrid( gb, Vector(0,0,0) )

print "Starting volume convolution"
GreenConvolve( MM, mgrid )

XXconvolve = gridded( mgrid )

Xtotal = Xconvolve + XXconvolve

gradMap = grad(Xtotal)
detGradMap = det(gradMap)

x = -10.0
y = -10.0
z = -10.0
dx = 0.2

while x <= 10.0:
	P = Vector(x,y,z)
	value = evaluate( Xconvolve, P )
	mvalue = evaluate( XXconvolve, P )
	tvalue = evaluate( Xtotal, P )
	dvalue = evaluate( detGradMap, P )
	print str(P.X()) + " " + str(P.Y()) + " " + str(P.Z()) + "      " + str(value.X()) + " " + str(value.Y()) + " " + str(value.Z())  + "      " + str(mvalue.X()) + " " + str(mvalue.Y()) + " " + str(mvalue.Z())   + "      " + str(tvalue.X()) + " " + str(tvalue.Y()) + " " + str(tvalue.Z()) + "       " + str(dvalue)

	x += dx
	y += dx
	z += dx
