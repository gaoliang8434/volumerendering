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
			noiseparms = Noise_t()
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


def TornadoChain( nbVerts, bottom, top, bottomradius, topradius, parameters  ):
	print "Tornado"
	dx = 1.0 / (nbVerts-1 )
	trail = AnchorChain()
	tangent = (top-bottom).unitvector()
	normal = Vector(1,0,0)
	normal = normal - tangent*dot_product(tangent,normal)
	normal = normal.unitvector()
	binormal = cross_product(normal,tangent).unitvector()
	print "tangent  " + str(tangent)
	print "normal   " + str(normal)
	print "binormal " + str(binormal)
	P0 = Vector(0,0,0)
	for i in range(0,nbVerts):
		x = i * dx
		noiseparms = Noise_t()
		noiseparms.P = top * (1.0-x) + bottom * x
		noiseparms.radius = topradius*(1.0-x) + bottomradius * x
		noiseparms.tangent = tangent
		noiseparms.normal = normal
		noiseparms.binormal = binormal
		noiseparms.translate = top - noiseparms.P
		noiseparms.gamma = float(parameters["-pyrogamma"])
		noiseparms.amplitude = float(parameters["-amp"])
		noiseparms.frequency = 2.33428293
		noiseparms.octaves = float(parameters["-octaves"])
		noiseparms.roughness = float(parameters["-roughness"])
		trail.push_back( noiseparms )
	return trail



def ProfileTornadoChain( nbVerts, bottom, top, bottomradius, topradius, rate, parameters  ):
	print "Tornado"
	dx = 1.0 / (nbVerts-1 )
	trail = AnchorChain()
	tangent = (top-bottom).unitvector()
	normal = Vector(1,0,0)
	normal = normal - tangent*dot_product(tangent,normal)
	normal = normal.unitvector()
	binormal = cross_product(normal,tangent).unitvector()
	P0 = Vector(0,0,0)
	for i in range(0,nbVerts):
		x = i * dx
		xx =  math.exp( - rate*x )
		noiseparms = Noise_t()
		noiseparms.P = top * (1.0-x) + bottom * x
		noiseparms.radius = topradius*xx 
		if noiseparms.radius < bottomradius:
			noiseparms.radius = bottomradius
		noiseparms.capradius = 0.0*noiseparms.radius
		noiseparms.tangent = tangent
		noiseparms.normal = normal
		#noiseparms.normal = normal * (1.0-xx) + tangent * xx
		noiseparms.normal = noiseparms.normal.unitvector()
		noiseparms.binormal = cross_product( noiseparms.normal, noiseparms.tangent ).unitvector()
		noiseparms.translate = top - noiseparms.P
		noiseparms.gamma = float(parameters["-pyrogamma"])
		noiseparms.amplitude = float(parameters["-amp"])
		#noiseparms.frequency = 0.2 * 2.33428293 * topradius / noiseparms.radius
		noiseparms.frequency = 2.33428293
		noiseparms.octaves = float(parameters["-octaves"])
		noiseparms.roughness = float(parameters["-roughness"])
		trail.push_back( noiseparms )
	return trail


		




def RotateTrail( trail, angle ):
	rotatedTrail = AnchorChain()
	for node in trail:
		theta = angle
		cth = math.cos( theta )
		sth = math.sin( theta )
		normal = node.normal * cth + node.tangent * ( node.tangent*node.normal ) * (1.0-cth) + cross_product( node.tangent, node.normal ) * sth
		binormal = node.binormal * cth + node.tangent * ( node.tangent*node.binormal ) * (1.0-cth) + cross_product( node.tangent, node.binormal ) * sth
		node.normal = normal
		node.binormal = binormal
		rotatedTrail.push_back(node)
	return rotatedTrail

def RadiusRotateTrail( trail, angle, maxradius, minradius, bias ):
	rotatedTrail = AnchorChain()
	for node in trail:
		rotationfactor = (1.0 - (node.radius-minradius)/(maxradius-minradius))
		if rotationfactor < bias:
			rotationfactor = bias
		theta = angle * rotationfactor
		cth = math.cos( theta )
		sth = math.sin( theta )
		normal = node.normal * cth + node.tangent * ( node.tangent*node.normal ) * (1.0-cth) + cross_product( node.tangent, node.normal ) * sth
		binormal = node.binormal * cth + node.tangent * ( node.tangent*node.binormal ) * (1.0-cth) + cross_product( node.tangent, node.binormal ) * sth
		node.normal = normal
		node.binormal = binormal
		rotatedTrail.push_back(node)
	return rotatedTrail

def  TranslateNoise( trail, displacement ):
	trail2 = AnchorChain()
	for node in trail:
		node.translate = node.translate + displacement
		trail2.push_back(node)
	return trail2




beginJob()

#
#=========== pyro construction ============================
#

allparms = AllParameters( {"-usepyro":"yes", "-static":"yes"} )

CmdLineHelp( "-h" )

allparms["-pyrogamma"] = 0.3
allparms["-amp"] = 0.4
allparms["-octaves"] = 4.2
allparms["-roughness"] = 0.6
allparms["-translate"] = [ 0, 0, 0.34 ]
allparms["-name"] = "tornado.exr"
allparms["-scatter"] = 2.0
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-fov"] = 70.0
allparms["-ds"] = 0.02
allparms["-near"] = 7.0
allparms["-radius"] = 2.5
allparms["-maxpathlength"] = 10.0
allparms["-dsmcell"] = [ 0.025, 0.025, 0.025 ] 
allparms["-dsmlength"] = [ 8, 8, 8 ]
allparms["-dsmllc"] = [ -4, -4, -4 ]
allparms["-case"] = "ProfilePyroUpwardFlow"
allparms["-version"] = 1
allparms["-subsamples"] = 1
allparms["-lightstyle"] = "custom"
allparms["-lightP"] = [[10, 0, 15], [-5, -40, 10] ]
allparms["-lightCd"] = [[0.5,0.5,0.5], [1,1,1]]
allparms["-backgroundcolor"] = [ 0.5, 0.5, 0.6, 0.0 ]
allparms = CmdLineFindValues( allparms )


center = [0, 0, 0]
radii = []
centers = []
amp = allparms["-amp"]
colors = []
noiseparms = Noise_t()
noiseparms.radius = 0.25

tornadotopradius = 2.0
tornadobottomradius = 0.75
rotationbias = 0.1
tornadotop = Vector( 0, 3, 0 )
tornadobottom = Vector( 0, -3, 0 )
nbVerts = 36
#chain = TornadoChain( nbVerts, tornadobottom, tornadotop, tornadobottomradius, tornadotopradius, allparms )

deviationpoint = 0.5
deviationwidth = 0.005
deviationdirection = Vector( 1, 0, 0 )


tornadotopradius = 3.0 
tornadobottomradius = 0.2
profilerate = 5.0 
chain = ProfileTornadoChain( nbVerts, tornadobottom, tornadotop, tornadobottomradius, tornadotopradius, profilerate, allparms )
#DeviateTornadoChain( chain, deviationpoint, deviationwidth, deviationdirection )




llc = Vector( float(allparms["-dsmllc"][0]), float(allparms["-dsmllc"][1]), float(allparms["-dsmllc"][2]) )
urc = llc + Vector( float(allparms["-dsmlength"][0]), float(allparms["-dsmlength"][1]), float(allparms["-dsmlength"][2]) )
res = Vector( float(allparms["-dsmcell"][0]), float(allparms["-dsmcell"][1]), float(allparms["-dsmcell"][2]) )
#res = Vector( 0.02, 0.02, 0.02 )
gridb = makeGridBox( llc, urc, res )

PrintParameters( allparms )

nbFrames = 300
dtheta = float(5.0)
for frame in range(10,nbFrames+1):
	print "FRAME " + str(frame)
	allparms["-turntableframe"] = int(frame)
	trailage = float(frame)/float(nbFrames-1)
	nbTrailPoints = int( float(trailage) * float(nbVerts) )
	print "Number of points: " + str(nbTrailPoints)
	if nbTrailPoints < 2:
		nbTrailPoints = 2
	trail = AnchorChain()
	if allparms["-static"] != "yes":
		trail = TravelingTrail( trailage, nbTrailPoints, chain  )
	if allparms["-static"] == "yes":
		trail = StaticTrail( trailage, nbTrailPoints, chain  )
	trail = RotateTrail( trail, dtheta * float(frame-1) * (3.14159265/180.0) )
	displacement = Vector( 0.0, 1.0, 0.0 ) * 0.04 * float(frame)
	trail = TranslateNoise( trail, displacement )
	for node in trail:
		print str(node.P) + " "  + str(node.radius)
	volumeObject = PiecewiseCurveField( trail )
	if allparms["-usepyro"] == "yes":
		volumeObject = PiecewisePyroCurveField( trail )
		print "Using PiecewisePyroCurveField"
	densityGrid = makeGrid( gridb, 0.0 )
	stamp( densityGrid, volumeObject, 1 )
	volumeObject = gridded(densityGrid)
	densityField = clamp( volumeObject/constant(0.1), 0, 1.0 )
	StandardRender( densityField, allparms )

endJob()


