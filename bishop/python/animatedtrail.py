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
	noiseparms = Noise_t()
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
			noiseparms.radius = thickness
			noiseparms.translate = P0
			noiseparms.gamma = 0.3
			noiseparms.amplitude = 0.2
			noiseparms.frequency = 2.33428293
			noiseparms.octaves = 4.2
			noiseparms.roughness = 0.65
			dP = noiseparms.P - P0
			noiseparms.translate = noiseparms.translate + noiseparms.tangent * dP.magnitude()
			chain.push_back( noiseparms )
			P0 = noiseparms.P;
	return chain





def TrailPath( nbAnchors, startPosition, endPosition, thickness ):
	translate = Vector(0,0,0)
	chain = AnchorChain()
	for i in range(0,nbAnchors):
		noiseparms = Noise_t()
		t = float(i)/float(nbAnchors-1)
		noiseparms.P = startPosition * (1.0-t) + endPosition * t
		noiseparms.tangent = (endPosition-startPosition).unitvector()
		noiseparms.normal = Vector(noiseparms.tangent.Z(), 0.0, -noiseparms.tangent.X() ).unitvector()
		noiseparms.binormal = cross_product(noiseparms.tangent, noiseparms.normal).unitvector()
		noiseparms.radius = thickness
		noiseparms.translate = noiseparms.P - startPosition 
		noiseparms.gamma = 0.6
		noiseparms.amplitude = 0.1
		noiseparms.frequency = 2.33428293
		noiseparms.octaves = 4.2
		noiseparms.roughness = 0.65
		chain.push_back( noiseparms )
	return chain

def FinalTrailPath( nbAnchors, startPosition, endPosition, thickness ):
	translate = Vector(0,0,0)
	chain = AnchorChain()
	for i in range(0,nbAnchors):
		noiseparms = Noise_t()
		t = float(i)/float(nbAnchors-1)
		noiseparms.P = startPosition * (1.0-t) + endPosition * t
		noiseparms.tangent = (endPosition-startPosition).unitvector()
		noiseparms.normal = Vector(noiseparms.tangent.Z(), 0.0, -noiseparms.tangent.X() ).unitvector()
		noiseparms.binormal = cross_product(noiseparms.tangent, noiseparms.normal).unitvector()
		noiseparms.radius = thickness*4.0
		noiseparms.translate = noiseparms.P - startPosition 
		noiseparms.gamma = 0.3
		noiseparms.amplitude = 0.5
		noiseparms.frequency = 2.33428293
		noiseparms.octaves = 3.2
		noiseparms.roughness = 0.35
		chain.push_back( noiseparms )
	return chain




def AgeTrain( chain1, chain2, earlyAge, lateAge ):
	chain = AnchorChain()
	for i in range(0, chain1.size() ):
		a1 = chain1[i]
		a2 = chain2[i]
		age = earlyAge + i * ( lateAge-earlyAge )/( chain1.size()-1 )
		chain.push_back( interpolateAnchors( a1, a2, age )  )
	return chain

beginJob()

#
#=========== pyro construction ============================
#

allparms = AllParameters( {"-usepyro":"yes", "-static":"yes"} )

CmdLineHelp( "-h" )

allparms["-pyrogamma"] = 0.75
allparms["-amp"] = 1.5
allparms["-octaves"] = 6.0
allparms["-roughness"] = 0.45
allparms["-translate"] = [ 0, 0, 0.34 ]
allparms["-name"] = "animatedTrail.exr"
allparms["-scatter"] = 10.0
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-fov"] = 70.0
allparms["-ds"] = 0.028
allparms["-near"] = 0.0
allparms["-radius"] = 2.5
allparms["-amp"] = 1.0
allparms["-maxpathlength"] = 50.0
allparms["-dsmcell"] = [ 0.05, 0.05, 0.05 ] 
allparms["-dsmlength"] = [ 20, 20, 30 ]
allparms["-dsmllc"] = [ -10, -10, -30 ]
allparms["-case"] = "PyroHQ"
allparms["-subsamples"] = 5
#allparms["-size"] = [ 960, 540 ]

allparms = CmdLineFindValues( allparms )


center = [0, 0, 0]
radii = []
centers = []
amp = allparms["-amp"]
colors = []

chain = AnchorChain()
thickness =  0.1

nbAnchors = 60
startPosition = Vector( 10, 10, 10 )
endPosition = Vector( 0, -7, -30)
earlychain = TrailPath( nbAnchors, startPosition, endPosition, thickness )
latechain = FinalTrailPath( nbAnchors, startPosition, endPosition, thickness )

PrintParameters( allparms )

llc = Vector( float(allparms["-dsmllc"][0]), float(allparms["-dsmllc"][1]), float(allparms["-dsmllc"][2]) )
urc = llc + Vector( float(allparms["-dsmlength"][0]), float(allparms["-dsmlength"][1]), float(allparms["-dsmlength"][2]) )
res = Vector( float(allparms["-dsmcell"][0]), float(allparms["-dsmcell"][1]), float(allparms["-dsmcell"][2]) )
gridb = makeGridBox( llc, urc, res )

nbFrames = 100
for frame in range(99,nbFrames):
	print "FRAME " + str(frame+1)
	allparms["-turntableframe"] = int(frame+1)
	trailage =  0.3 + 0.7*float(frame)/float(nbFrames-1)
	nbTrailPoints = int( float(trailage) * float(nbAnchors) )
	print "Number of points: " + str(nbTrailPoints)
	if nbTrailPoints < 2:
		nbTrailPoints = 2
	trail = AnchorChain
	if allparms["-static"] != "yes":
		early = TravelingTrail( trailage, nbTrailPoints, earlychain )
		late = TravelingTrail( trailage, nbTrailPoints, latechain )
		trail = AgeTrain( early, late, trailage, 0.0  )
	if allparms["-static"] == "yes":
		early = StaticTrail( trailage, nbTrailPoints, earlychain )
		late = StaticTrail( trailage, nbTrailPoints, latechain )
		trail = AgeTrain( early, late, trailage, 0.0  )
	volumeObject = PiecewiseCurveField( trail )
	if allparms["-usepyro"] == "yes":
		volumeObject = PiecewisePyroCurveField( trail )
		print "Using PiecewisePyroCurveField"
	densityGrid = makeGrid( gridb, 0.0 )
	stamp( densityGrid, volumeObject )
	volumeObject = gridded(densityGrid)
	densityField = clamp( volumeObject/constant(0.05), 0, 1.0 )
	StandardRender( densityField, allparms )

endJob()


