#!/usr/bin/python


import sys
import os
from bishoputils import *
from ImplicitShapes import *



def evaluateParticleGroup( anchorList, x ):
	anchorsize = int(anchorList.size())
	anchor = int(x * anchorsize)
	w = x * float(anchorsize) - float(anchor)
   	if anchor < 0:
		anchor = 0
		w = 0
	if anchor >= anchorsize-1:
		anchor = anchorsize-2
		w = 1
	value = Particle()
	Interpolate( anchorList[anchor], anchorList[anchor+1], w, value )
	return value


def TravelingTrail( trailage, nbTrailPoints, chain  ):
	print "TravelingTrail"
	dx = trailage / (nbTrailPoints-1 )
	trail = ParticleGroup()
	id = chain[chain.size()-1]._id - nbTrailPoints
	for i in range(0,nbTrailPoints):
		x = i * dx
		xtrailhead = x + (1.0 - trailage)
		anchor = evaluateParticleGroup( chain, xtrailhead )
		localanchor = evaluateParticleGroup( chain, x )
		anchor._P = localanchor.P() 
		anchor._v = localanchor.v() 
		anchor._up = localanchor.up()
		anchor._normal = localanchor.normal()
		anchor._right = localanchor.right()
		anchor._id = id + i
		trail.push_back( anchor )
	if trail.size() > 0:
		trail[trail.size()-1]._opacity = 0
	return trail


def TravelingTrail2( trailage, nbTrailPoints, chain  ):
	print "TravelingTrail"
	dx = trailage / (nbTrailPoints-1 )
	trail = ParticleGroup()
	id = chain[chain.size()-1]._id - nbTrailPoints
	for i in range(0,nbTrailPoints):
		x = i * dx
		anchor = chain[ chain.size() + i - nbTrailPoints ]
		localanchor = chain[i]
		#localanchor = evaluateParticleGroup( chain, x )
		anchor._P = localanchor.P() 
		anchor._v = localanchor.v() 
		anchor._up = localanchor.up()
		anchor._normal = localanchor.normal()
		anchor._right = localanchor.right()
		anchor._translate = localanchor.translate()
		anchor._wisp_translate = localanchor.wispTranslate()
		trail.push_back( anchor )
	trail.push_back(anchor)
	return trail




def StaticTrail( trailage, nbTrailPoints, chain  ):
	print "StaticTrail"
	dx = trailage / (nbTrailPoints-1 )
	trail = ParticleGroup()
	for i in range(0,nbTrailPoints):
		x = i * dx
		anchor = evaluateParticleGroup( chain, x )
		anchor._id = i
		trail.push_back( anchor )
	if trail.size() > 0:
		trail[trail.size()-1]._opacity = 0
	return trail


def StaticTrail2( trailage, nbTrailPoints, chain  ):
	print "StaticTrail"
	dx = trailage / (nbTrailPoints-1 )
	trail = ParticleGroup()
	for i in range(0,nbTrailPoints-1):
		anchor = chain[i]
		trail.push_back( anchor )
	anchor = evaluateParticleGroup( chain, trailage )
	trail.push_back( anchor )
	return trail





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
allparms["-name"] = "coil.exr"
allparms["-scatter"] = 5
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-fov"] = 70.0
allparms["-ds"] = 0.02
allparms["-near"] = 7.0
allparms["-radius"] = 2.5
allparms["-capradius"] = 1.0
allparms["-amp"] = 1.0
allparms["-maxpathlength"] = 10.0
allparms["-dsmcell"] = [ 0.02, 0.02, 0.02 ] 
allparms["-dsmlength"] = [ 6, 6, 6 ]
allparms["-dsmllc"] = [ -3, -3, -3 ]
allparms["-dsmsamples"] = 1
allparms["-case"] = "StaticWispChain"
allparms["-version"] = 1
allparms["-subsamples"] = 1
allparms["-lightP"] = [ [0, 20, 0] ]
allparms["-lightCd"] = [ [1, 1, 1] ]
allparms["-nbthreads"] = 16
allparms["-size"] = [ 960, 540 ]
allparms["-eye"] = [ 0, 3, 10 ]
allparms = CmdLineFindValues( allparms )


center = [0, 0, 0]
radii = []
centers = []
amp = allparms["-amp"]
colors = []

thickness = 0.15
capthickness = 0.15
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

topColor = Color( 1.0, 0.2, 1.0, 1.0 )
bottomColor = Color( 0.2, 1.0, 0.2, 1.0 )

chain = ParticleGroup()
translate = Vector(0,0,0)
for j in range(0,nbTurns):
	for i in range(0,nbVerts):
		noiseparms = Particle()
		theta = i * du
		ct = math.cos(theta)*springradius
		st = math.sin(theta)*springradius
		z = (i + j*nbVerts) * ddz - 2.0
		noiseparms._P = Vector(ct, z, st)
		ct = math.cos(theta+du)*springradius
		st = math.sin(theta+du)*springradius
		z = (i + j*nbVerts) * ddz - 2.0 + ddz
		noiseparms._normal = (Vector(ct, z, st) - noiseparms.P()).unitvector()
		ct = math.cos(theta-du)*springradius
		st = math.sin(theta-du)*springradius
		z = (i + j*nbVerts) * ddz - 2.0 - ddz
		T0 = (noiseparms.P() - Vector(ct, z, st)).unitvector()
		noiseparms._right = noiseparms.normal() - T0
		if i == 0 and j == 0:
			noiseparms._right = Vector( -noiseparms.normal().Z(), 0.0, noiseparms.normal().X() )
		noiseparms._right = noiseparms.right() - noiseparms.normal()*( noiseparms.normal()*noiseparms.right() )
		noiseparms._right = noiseparms.right().unitvector()
		noiseparms._up = cross_product(noiseparms.normal(), noiseparms.right()).unitvector()
		if i == 0 and j == 0:
			P0 = noiseparms.P()
		noiseparms._pscale = thickness
		noiseparms._translate = noiseparms.P()
		noiseparms._wisp_translate = noiseparms.P() 
		noiseparms._freq = 1.33428293
		noiseparms._wisp_freq = 1.33428293
		noiseparms._wisp_octaves = 2
		noiseparms._octaves = 2
		noiseparms._roughness = 0.65
		noiseparms._wisp_displacementScale= 1.5
		noiseparms._nbWisps = 1000
		noiseparms._id = (i + j*nbVerts) 
		dP = noiseparms.P() - P0
		#noiseparms._normal = noiseparms.normal() + noiseparms.right() * dP.magnitude()
		epsilon = float(i + j*nbVerts)/float( nbVerts*nbTurns -1 )
		noiseparms._Cd = topColor * epsilon + bottomColor * (1.0-epsilon)
		chain.push_back( noiseparms )
		P0 = noiseparms.P();



PrintParameters( allparms )


llc = Vector( float(allparms["-dsmllc"][0]), float(allparms["-dsmllc"][1]), float(allparms["-dsmllc"][2]) )
urc = llc + Vector( float(allparms["-dsmlength"][0]), float(allparms["-dsmlength"][1]), float(allparms["-dsmlength"][2]) )
#res = Vector( float(allparms["-dsmcell"][0]), float(allparms["-dsmcell"][1]), float(allparms["-dsmcell"][2]) )
res = Vector( 0.025, 0.025, 0.025 )
gridb = makeGridBox( llc, urc, res )

nbFrames = 300
for frame in range(1,nbFrames):
	print "FRAME " + str(frame)
	allparms["-turntableframe"] = int(frame)
	trailage = float(frame)/float(nbFrames-1)
	nbTrailPoints = int( float(trailage) * float(nbVerts*nbTurns) )
	print "Number of points: " + str(nbTrailPoints)
	if nbTrailPoints < 2:
		nbTrailPoints = 2
	trail = ParticleGroup
	if allparms["-static"] != "yes":
		trail = TravelingTrail2( trailage, nbTrailPoints, chain  )
	if allparms["-static"] == "yes":
		trail = StaticTrail2( trailage, nbTrailPoints, chain  )
	densityGrid = makeGrid( gridb, 0.0 )
	colorGrid = makeGrid( gridb, Color(0,0,0,0) )
	StampSplineWisps( densityGrid, colorGrid, trail )
	densityField = gridded(densityGrid)
	colorField = gridded(colorGrid) / ( densityField + constant(0.000001) )
	StandardTwoPartRender( constant(0.0), densityField, constant(Color(0,0,0,0)), colorField, allparms )

endJob()


