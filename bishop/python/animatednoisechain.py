#!/usr/bin/python


import sys
import os
from bishoputils import *
from ImplicitShapes import *



def TravelingTrail( trailage, nbTrailPoints, chain  ):
	print "TravelingTrail"
	dx = trailage / (nbTrailPoints-1 )
	trail = AnchorChain()
	for i in range(0,nbTrailPoints):
		x = i * dx
		xtrailhead = x + (1.0 - trailage)
		anchor = evaluateAnchorChain( chain, xtrailhead )
		localanchor = evaluateAnchorChain( chain, x )
		anchor.P = localanchor.P 
		anchor.v = localanchor.v 
		anchor.tanget = localanchor.tangent
		anchor.normal = localanchor.normal
		anchor.binormal = localanchor.binormal
		trail.push_back( anchor )
	return trail


def StaticTrail( trailage, nbTrailPoints, chain  ):
	print "StaticTrail"
	dx = trailage / (nbTrailPoints-1 )
	trail = AnchorChain()
	for i in range(0,nbTrailPoints):
		x = i * dx
		anchor = evaluateAnchorChain( chain, x )
		trail.push_back( anchor )
	return trail



def CoiledSpring( nbTurns, nbVerts, springRadius, springGap, ybase ):
	translate = Vector(0,0,0)
	chain = AnchorChain()
	for j in range(0,nbTurns):
		for i in range(0,nbVerts):
			theta = i * du
			ct = math.cos(theta)*springRadius
			st = math.sin(theta)*springRadius
			y = (i + j*nbVerts) * springGap - ybase
			noiseparms.P = Vector(ct, y, st)
			ct = math.cos(theta+du)*springRadius
			st = math.sin(theta+du)*springRadius
			y = (i + j*nbVerts) * springGap - ybase + springGap
			noiseparms.tangent = (Vector(ct, y, st) - noiseparms.P).unitvector()
			ct = math.cos(theta-du)*springRadius
			st = math.sin(theta-du)*springRadius
			y = (i + j*nbVerts) * springGap - ybase - springGap
			T0 = (noiseparms.P - Vector(ct, y, st)).unitvector()
			noiseparms.normal = noiseparms.tangent - T0
			if i == 0 and j == 0:
				noiseparms.normal = Vector( -noiseparms.tangent.Z(), 0.0, noiseparms.tangent.X() )
			noiseparms.normal = noiseparms.normal - noiseparms.tangent*( noiseparms.tangent*noiseparms.normal )
			noiseparms.normal = noiseparms.normal.unitvector()
			noiseparms.binormal = cross_product(noiseparms.tangent, noiseparms.normal).unitvector()
			if i == 0 and j == 0:
				P0 = noiseparms.P
			noiseparms.radius = thickness * 2.0
			noiseparms.translate = P0
			noiseparms.gamma = 0.3
			noiseparms.amplitude = 0.2
			noiseparms.frequency = 2.33428293
			noiseparms.octaves = 4.2
			noiseparms.roughness = 0.65
			noiseparms.pscale = 1.0 
			noiseparms.falloff=  1.0
			dP = noiseparms.P - P0
			noiseparms.translate = noiseparms.translate + noiseparms.tangent * dP.magnitude()
			chain.push_back( noiseparms )
			P0 = noiseparms.P;
	return chain



beginJob()

#
#=========== trail construction ============================
#

allparms = AllParameters( {"-usenoise":"yes", "-static":"no"} )

CmdLineHelp( "-h" )

allparms["-pyrogamma"] = 0.75
allparms["-amp"] = 1.5
allparms["-octaves"] = 6.0
allparms["-roughness"] = 0.45
allparms["-translate"] = [ 0, 0, 0.34 ]
allparms["-name"] = "animatedchain.exr"
allparms["-scatter"] = 5.0
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-fov"] = 70.0
allparms["-ds"] = 0.025
allparms["-near"] = 7.0
allparms["-radius"] = 2.5
allparms["-amp"] = 1.0
allparms["-maxpathlength"] = 10.0
allparms["-dsmcell"] = [ 0.03, 0.03, 0.03 ] 
allparms["-dsmlength"] = [ 6, 6, 6 ]
allparms["-dsmllc"] = [ -3, -3, -3 ]
allparms["-case"] = "Noise"

allparms = CmdLineFindValues( allparms )


center = [0, 0, 0]
radii = []
centers = []
amp = allparms["-amp"]
colors = []
noiseparms = Noise_t()
noiseparms.radius = 0.25

chain = AnchorChain()
thickness = 0.15
springradius = 2.0
dz = 7.5*thickness
nbVerts = 36
nbTurns = 4

ddz = dz/(nbVerts-1.0)
z = -2.0
P0 = Vector(0,0,0)
T0 = Vector(1,1,1)

du = 2.0*3.14159265/nbVerts
L = ddz/du


translate = Vector(0,0,0)
for j in range(0,nbTurns):
	for i in range(0,nbVerts):
		theta = i * du
		ct = math.cos(theta)*springradius
		st = math.sin(theta)*springradius
		z = (i + j*nbVerts) * ddz - 2.0
		noiseparms.P = Vector(ct, z, st)
		ct = math.cos(theta+du)*springradius
		st = math.sin(theta+du)*springradius
		z = (i + j*nbVerts) * ddz - 2.0 + ddz
		noiseparms.tangent = (Vector(ct, z, st) - noiseparms.P).unitvector()
		ct = math.cos(theta-du)*springradius
		st = math.sin(theta-du)*springradius
		z = (i + j*nbVerts) * ddz - 2.0 - ddz
		T0 = (noiseparms.P - Vector(ct, z, st)).unitvector()
		noiseparms.normal = noiseparms.tangent - T0
		if i == 0 and j == 0:
			noiseparms.normal = Vector( -noiseparms.tangent.Z(), 0.0, noiseparms.tangent.X() )
		noiseparms.normal = noiseparms.normal - noiseparms.tangent*( noiseparms.tangent*noiseparms.normal )
		noiseparms.normal = noiseparms.normal.unitvector()
		noiseparms.binormal = cross_product(noiseparms.tangent, noiseparms.normal).unitvector()
		if i == 0 and j == 0:
			P0 = noiseparms.P
		noiseparms.radius = thickness
		noiseparms.translate = P0
		noiseparms.gamma = 0.3
		noiseparms.amplitude = 0.2
		noiseparms.frequency = 2.33428293
		noiseparms.octaves = 2.2
		noiseparms.roughness = 0.65
		dP = noiseparms.P - P0
		noiseparms.translate = noiseparms.translate + noiseparms.tangent * dP.magnitude()
		chain.push_back( noiseparms )
		P0 = noiseparms.P;



PrintParameters( allparms )


llc = Vector( float(allparms["-dsmllc"][0]), float(allparms["-dsmllc"][1]), float(allparms["-dsmllc"][2]) )
urc = llc + Vector( float(allparms["-dsmlength"][0]), float(allparms["-dsmlength"][1]), float(allparms["-dsmlength"][2]) )
res = Vector( 0.025, 0.025, 0.025 )
gridb = makeGridBox( llc, urc, res )

nbFrames = 300
for frame in range(299,nbFrames):
	print "FRAME " + str(frame+1)
	allparms["-turntableframe"] = int(frame+1)
	trailage = float(frame)/float(nbFrames-1)
	nbTrailPoints = int( float(trailage) * float(nbVerts*nbTurns) )
	print "Number of points: " + str(nbTrailPoints)
	if nbTrailPoints < 2:
		nbTrailPoints = 2
	trail = AnchorChain
	if allparms["-static"] != "yes":
		trail = TravelingTrail( trailage, nbTrailPoints, chain  )
	if allparms["-static"] == "yes":
		trail = StaticTrail( trailage, nbTrailPoints, chain  )
	volumeObject = PiecewiseCurveField( trail )
	if allparms["-usenoise"] == "yes":
		volumeObject = PiecewiseNoiseCurveField( trail )
		print "Using PiecewiseNoiseCurveField"
	densityGrid = makeGrid( gridb, 0.0 )
	stamp( densityGrid, volumeObject )
	volumeObject = gridded(densityGrid)
	densityField = clamp( volumeObject, 0, 1.0 )
	StandardRender( densityField, allparms )

endJob()


