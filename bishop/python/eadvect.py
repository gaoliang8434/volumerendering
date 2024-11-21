#!/usr/bin/python

from bishoputils import *

parms = AllParameters({})


def XDistribute( X, gbox ):
	delx = dx(gbox)
	dely = dy(gbox)
	delz = dz(gbox)
	cdx = constant(Vector(2.0*delx,0,0))
	cdy = constant(Vector(0,2.0*dely,0))
	cdz = constant(Vector(0,0,2.0*delz))
	result = warp( X, identity() + cdx )
	result += warp( X, identity() + cdy )
	result += warp( X, identity() + cdz )
	result += warp( X, identity() - cdx )
	result += warp( X, identity() - cdy )
	result += warp( X, identity() - cdz )
	result *= constant(1.0/6.0)
	return result

def MDistribute( X, M, gbox ):
	delx = dx(gbox)
	dely = dy(gbox)
	delz = dz(gbox)
	cdx = constant(Vector(delx,0,0))
	cdy = constant(Vector(0,dely,0))
	cdz = constant(Vector(0,0,delz))
	result = cdx*warp(M,warp(X,identity()-cdx))
	result += cdy*warp(M,warp(X,identity()-cdy))
	result += cdz*warp(M,warp(X,identity()-cdz))
	result -= cdx*warp(M,warp(X,identity()+cdx))
	result -= cdy*warp(M,warp(X,identity()+cdy))
	result -= cdz*warp(M,warp(X,identity()+cdz))
	result *= constant(1.0/3.0)
	return result


def XSolve( U, delt, gbox, nbiter ):
	gU = grad(U) * constant(delt)
	M = exp(gU)
	dM = det(M)
	ShowStatistics(dM)
	X = identity() - U*constant(delt)
	#X = identity()
	for i in range(0,nbiter):
		#X = MDistribute( identity(), M, gbox ) + XDistribute( X, gbox )
		X = MDistribute( X, M, gbox ) + XDistribute( X, gbox )
		grid = makeGrid( gbox, Vector(0,0,0) )
		stamp( grid, X-identity() )
		X = identity() + gridded(grid)
		ShowStatistics( det(grad(X)) )
	return X





def LeVequeVelocity():
	X = xComponent( identity() )
	Y = yComponent( identity() )
	Z = zComponent( identity() )
	Pi = constant( math.pi )
	TwoPi = constant( 2.0*math.pi )
	U = constant(2.0) * sin( X*Pi )*sin( X*Pi ) * sin( TwoPi*Y ) * sin( TwoPi*Z )
	V = constant(-1.0)* sin( X*TwoPi ) * sin( Y*Pi )*sin( Y*Pi ) * sin( TwoPi*Z )
	W = constant(-1.0)* sin( X*TwoPi ) * sin( Y*TwoPi ) * sin( Z*Pi )*sin( Z*Pi )
	return XYZ( U,V,W )


def ShowPreservation( Y ):
	dY = grad(Y)
	detdY = det(dY)
	for j in range(0,100):
		P = Vector( 0.01*j, 0.01*j, 0.01*j )
		dmvalue = evaluate( detdY, P )
		print str(j) + "  " + str(dmvalue)

def ShowStatistics( detdY ):
	mean = 0.0
	stddev = 0.0
	nbterms = 0
	maxvalue = -1000.0
	minvalue = 1000.0
	for k in range(0,100):
		for j in range(0,100):
			for i in range(0,100):
				P = Vector( 0.01*i, 0.01*j, 0.01*k )
				value = evaluate(detdY, P )
				mean += value
				stddev += value*value
				nbterms += 1
				maxvalue = max( maxvalue, value )
				minvalue = min( minvalue, value )
	mean /= float(nbterms)
	stddev /= float(nbterms)
	stddev -= mean*mean
	if stddev >= 0:
		stddev = math.sqrt(stddev)
	else:
		stddev = -math.sqrt(-stddev)
	print "max " + str(maxvalue) + "  min " + str(minvalue) + "  mean " + str(mean) + "   stddev " + str(stddev)


sphere = constant(0.15*0.15) - (identity()-constant(Vector(0.35,0.35,0.35)))*(identity()-constant(Vector(0.35,0.35,0.35)))

velocity = LeVequeVelocity()
dt = 0.1
box = makeGridBox( Vector(0,0,0), Vector(1,1,1), Vector(1.0/128.0, 1.0/128.0, 1.0/128.0) )

Y = FFTVolumePreserve( box, velocity*constant(dt) )
#Y = XSolve( velocity, dt, box, 10 )
ShowPreservation( Y )
ShowStatistics( det(grad(Y)) )

parms["-fov"] = 10.0
parms["-view"] = [0.35,0.35,0.35]
parms["-ds"] = 0.02
parms["-near"] = 8.5
parms["-maxpathlength"] = 3.5
parms["-eye"][0] += parms["-view"][0]
parms["-eye"][1] += parms["-view"][1]
parms["-eye"][2] += parms["-view"][2]
parms["-name"] = "leveque.exr"
parms["-case"] = "SlvI"



X = identity()


for i in range(1,11):
	print "Frame " + str(i)
	parms["-turntableframe"] = int(i)
	grid = makeGrid( box, Vector(0,0,0) )
	stamp( grid, X-identity() )
	X = identity() + gridded(grid)
	density = mask(warp(sphere,X))
	PrintParameters(parms)
	StandardRender( density, parms )
	#X = advect( X, velocity, dt )
	X =  warp(X,Y)
	#ShowPreservation( X )
	ShowStatistics( det(grad(X)) )




beginJob()












endJob()
