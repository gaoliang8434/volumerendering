#!/usr/bin/python


import sys
import os
from bishoputils import *
from ImplicitShapes import *


def ObjParameters():
	parms = {  "-objfilename":"", "-objscale":1.0, "-objtranslate":[0,0,0], "-levelsetNxNyNz":[100, 100, 100], 
		           "-lsthresh":0.1, "-clamp":10.0, "-saveLevelset":"", "-readLevelset":"" }
	return parms



BlueRidge     = Color( 0.227, 0.286, 0.345, 1)
BowmanField   = Color( 0.337, 0.380, 0.153, 1)
HowardsRock   = Color( 0.408, 0.361, 0.325, 1)
ClemsonOrange = Color( 0.965, 0.404, 0.2  , 1)
Regalia       = Color(   82.0/255.0, 45.0/255.0, 128.0/255.0, 1)

ClemsonColors = [ BlueRidge, BowmanField, HowardsRock, ClemsonOrange, Regalia  ]

random.seed( 48475 )

def GenerateBunnies( nbBunnies, chain, nbFrames ):
	collection = []
	for i in range(0,nbBunnies):
		P = chain[chain.size()-1].P()*1.0
		startTime = int( random.random() * float(nbFrames) )
		orientation = Vector( random.random()-0.5, random.random()-0.5, random.random()-0.5 ).unitvector()
		orientationAngle = random.random() * 2.0 * 3.14159265
		rotationAxis = Vector( random.random()-0.5, random.random()-0.5, random.random()-0.5 ).unitvector()
		rotationAngle = 0.0
		color = ClemsonColors[ int( random.random() * 3848404048.0 ) % len(ClemsonColors) ]
		scale = random.random() * 0.3
		if scale < 0.1:
			scale = 0.1
		bunny = [ P, startTime, orientation, orientationAngle, rotationAxis, rotationAngle, color, scale, 0.0 ]
		print "Bunny location " + str(P) + " startTime " + str(startTime) + " orientation " + str(orientation) + " color " + str(color) 
		collection.append( bunny )
	return collection
		


def TornadoChain( nbVerts, bottom, top, bottomradius, topradius, parameters  ):
	print "Tornado"
	dx = 1.0 / (nbVerts-1 )
	trail = ParticleGroup()
	normal = (top-bottom).unitvector()
	right = Vector(1,0,0)
	right = right - normal*dot_product(right,normal)
	right = right.unitvector()
	up = cross_product(right,normal).unitvector()
	print "normal   " + str(normal)
	print "right    " + str(right)
	print "up       " + str(up)
	print "radii " + str(bottomradius) + "  " + str(topradius)
	P0 = Vector(0,0,0)
	topcolor = Color( 0.227,0.286,0.345,1)
	bottomcolor = Color( 0.965,0.404,0.2,1)
	rate = 5.0
	for i in range(0,nbVerts):
		x = i * dx
		noiseparms = Particle()
		noiseparms._P = top * (1.0-x) + bottom * x
		xx =  math.exp( - rate*x )
		radius = topradius * xx
		if radius < bottomradius:
			radius = bottomradius
		noiseparms._pscale = radius
		noiseparms._normal = normal
		noiseparms._right = right
		noiseparms._up = up
		noiseparms._translate = top - noiseparms.P()
		noiseparms._wisp_translate = top - noiseparms.P()
		noiseparms._freq = 2.33428293
		noiseparms._wisp_freq = 2.33428293
		noiseparms._octaves = float(parameters["-octaves"])
		noiseparms._roughness = float(parameters["-roughness"])
		noiseparms._nbWisps = int(50000 * math.pow( radius/topradius,2));
		noiseparms._opacity = 1.0
		noiseparms._id = i
		xxx = math.tanh( (x-0.65)*rate )
		noiseparms._Cd = bottomcolor * xxx + topcolor * (1.0-xxx)
		noiseparms._wisp_displacementScale = 0.2
		trail.push_back( noiseparms )
	return trail


def CurvedTornadoChain( nbVerts, bottom, top, bottomradius, topradius, deflection, parameters  ):
	print "Tornado"
	dx = 1.0 / (nbVerts-1 )
	trail = ParticleGroup()
	normal = (top-bottom).unitvector()
	P0 = top
	topcolor = BlueRidge
	bottomcolor = ClemsonOrange/2.0
	rate = 5.0
	deflrate = 10.0
	translate = Vector(0,0,0)
	deflectionPoint = 0.75
	for i in range(0,nbVerts):
		x = i * dx
		xxx = (math.tanh( (x-deflectionPoint)*rate ) + 1.0)/2.0
		xx =  math.exp( - rate*x )
		xxxx = (math.tanh( (x-deflectionPoint)*deflrate ) + 1.0)/2.0
		xxxxp = ( 1.0 - math.tanh( (x-deflectionPoint)*deflrate ) *  math.tanh( (x-deflectionPoint)*deflrate ) )/2.0
		radius = topradius * xx
		if radius < bottomradius:
			radius = bottomradius
		noiseparms = Particle()
		noiseparms._P = top * (1.0-x) + bottom * x + Vector( deflection*xxx, 0, 0 )
		print "Particle " + str(i) + " location " + str( noiseparms.P() )
		noiseparms._pscale = radius
		noiseparms._normal = (bottom-top + Vector( deflection*deflrate*xxxxp, 0, 0 )).unitvector()
		noiseparms._right =  ( Vector( (top-bottom).magnitude(), 0, 0 ) - (bottom-top)*xxxxp*deflection*deflrate ).unitvector()
		noiseparms._up = Vector(0,0,1)
		translate = translate + P0-noiseparms.P()
		P0 = noiseparms.P()
		noiseparms._translate = translate
		noiseparms._wisp_translate = translate
		noiseparms._freq = 2.33428293
		noiseparms._wisp_freq = 2.33428293
		noiseparms._octaves = float(parameters["-octaves"])
		noiseparms._roughness = float(parameters["-roughness"]) * (1.0-xxxx) + 0.1 * xxxx
		noiseparms._nbWisps = int(5000000 * math.pow( radius/topradius,2));
		noiseparms._opacity = 1.0
		noiseparms._id = i
		noiseparms._Cd = bottomcolor * xxxx  + topcolor * (1.0-xxxx)
		noiseparms._wisp_displacementScale = 0.2
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
	rotatedTrail = ParticleGroup()
	for node in trail:
		theta = angle
		cth = math.cos( theta )
		sth = math.sin( theta )
		right = node.right() * cth + node.normal() * ( node.right()*node.normal() ) * (1.0-cth) + cross_product( node.right(), node.normal() ) * sth
		up = node.up() * cth + node.normal() * ( node.up()*node.normal() ) * (1.0-cth) + cross_product( node.up(), node.normal() ) * sth
		node._right = right
		node._up = up
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
	trail2 = ParticleGroup()
	previousNode = "None"
	for node in trail:
		if previousNode != "None":
			shift = previousNode.P() - node.P()
			shift = shift * displacement.magnitude()
			node._translate = node.translate() + shift
			node._wisp_translate = node.wispTranslate() + shift
		else:
			node._translate = node.translate() - node.normal() * displacement.magnitude()
			node._wisp_translate = node.wispTranslate() - node.normal() * displacement.magnitude()
		previousNode = node
		trail2.push_back(node)
	return trail2




# Teapot creation
def GetLevelset( filename, parameters ):
	p = ObjParser()
	p.ParseFile( filename )
	g = TriangleGeometry()
	g.setScaling( allparms["-objscale"] )
	fillresult = p.Fill(g)
	if fillresult == 0:
		print "Could not read geometry from file " + allparms["-objfilename"]
		exit()
	print "Number of vertices = " + str( g.nbVertices() ) + "   number of faces = " + str( g.nbFaces() )
	dims = g.URC() - g.LLC()
	center = (g.URC() + g.LLC())*0.5;
	dims *= 1.2;
	llc = center - dims*0.5;
	urc = llc + dims;
	print "Obj BB:   ",
	print str(llc.X()) + " " +  str(llc.Y()) + " " +  str(llc.Z()) + "  X  " +  str(urc.X()) + " " +  str(urc.Y()) + " " +  str(urc.Z())
	print "Obj Size: " + str(dims.X()) + " X " + str(dims.Y()) + " X " + str(dims.Z()) 
	cell = Vector( float(allparms["-dsmcell"][0]),  float(allparms["-dsmcell"][1]),  float(allparms["-dsmcell"][2])  )*0.5
	objgb = makeGridBox( llc, urc, cell )
	sgrid = makeGrid( objgb, 0.0 )
	RayMarchLevelSet( g, sgrid )
	lssdf = gridded(sgrid)
	return [ lssdf, objgb ]




beginJob()

#
#=========== pyro construction ============================
#

allparms = AllParameters( ObjParameters() )

CmdLineHelp( "-h" )

allparms["-objfilename"] = "../819/models/cleanteapot.obj"
allparms["-pyrogamma"] = 0.3
allparms["-amp"] = 0.4
allparms["-octaves"] = 2.2
allparms["-roughness"] = 0.6
allparms["-translate"] = [ 0, 0, 0.34 ]
allparms["-name"] = "wispTornado.exr"
allparms["-scatter"] = 2.0
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-fov"] = 90.0
allparms["-ds"] = 0.01
allparms["-near"] = 6.0
allparms["-radius"] = 2.5
allparms["-maxpathlength"] = 10.0
allparms["-dsmcell"] = [ 0.04, 0.04, 0.04 ] 
allparms["-dsmlength"] = [ 8, 8, 8 ]
allparms["-dsmllc"] = [ -4, -4, -4 ]
allparms["-case"] = "TeapotBunny"
allparms["-version"] = 2
allparms["-subsamples"] = 1
allparms["-lightstyle"] = "custom"
allparms["-lightP"] = [[10, 0, 15], [-5, -40, 10] ]
allparms["-lightCd"] = [[0.5,0.5,0.5], [1,1,1]]
#allparms["-size"] = [960,540]

allparms = CmdLineFindValues( allparms )




# Teapot creation
teapotData = GetLevelset( "../819/models/cleanteapot.obj", allparms )
teapotsdf = scale( teapotData[0], Vector(1.0,1.0,1.0)*2.0 )
# first attempt at putting the spout in the correct place
shift = llc(teapotData[1]) - Vector(-0.5,3.25,0)
teapotsdf = translate(teapotsdf, shift)
teapotdensity = clamp( teapotsdf*constant(10.0), 0.0, 1.0 )*constant(10.0)
teapotcolor = constant(BowmanField) * mask(teapotsdf)


# bunny creation
bunnyData = GetLevelset( "../819/models/cleanbunny.obj", allparms )





# Build tornado
center = [0, 0, 0]
radii = []
centers = []
amp = allparms["-amp"]
colors = []
noiseparms = Noise_t()
noiseparms.radius = 0.5

tornadotopradius = 3.0
tornadobottomradius = 0.25
rotationbias = 0.1
tornadotop = Vector( 0, 3, 0 )
tornadobottom = Vector( 0, -3, 0 )
nbVerts = 36
defl = 1.5
chain = CurvedTornadoChain( nbVerts, tornadobottom, tornadotop, tornadobottomradius, tornadotopradius, defl, allparms )
#chain = TornadoChain( nbVerts, tornadobottom, tornadotop, tornadobottomradius, tornadotopradius, allparms )


llc = Vector( float(allparms["-dsmllc"][0]), float(allparms["-dsmllc"][1]), float(allparms["-dsmllc"][2]) )
urc = llc + Vector( float(allparms["-dsmlength"][0]), float(allparms["-dsmlength"][1]), float(allparms["-dsmlength"][2]) )
#res = Vector( float(allparms["-dsmcell"][0]), float(allparms["-dsmcell"][1]), float(allparms["-dsmcell"][2]) )
res = Vector( 0.01, 0.01, 0.01 )
gridb = makeGridBox( llc, urc, res )


allparms["-dsmllc"][1] = float(allparms["-dsmllc"][1]) - 4.0
allparms["-dsmlength"][1] = float(allparms["-dsmlength"][1]) + 4.0
allparms["-dsmcell"] = [ 0.01, 0.01, 0.01 ]

nbFrames = 300


# set up bunnies
nbBunnies = 100
bunnies = GenerateBunnies( nbBunnies, chain, nbFrames ) 



renderStartFrame = 1 
PrintParameters( allparms )
dt = 0.1 * 300.0/float(nbFrames)
dtheta = float(5.0)
displacement = Vector( 0.0, 1.0, 0.0 ) * 0.04 
for frame in range(1,nbFrames+1):
	print "\n\nFRAME " + str(frame)
	allparms["-turntableframe"] = int(frame)
	chain = TranslateNoise( chain, displacement )
	trail = RotateTrail( chain, dtheta * float(frame-1) * (3.14159265/180.0) )
	for bunny in bunnies:
		if frame > bunny[1]:
			bunny[0] += Vector( random.random() - 0.5, random.random(), random.random()-0.5 ) * dt
			bunny[5] += 15.0 * random.random() * 300.0/float(nbFrames)
			if bunny[0].Y() > 3.5:
				bunny[8] -= 0.1
			else:
				bunny[8] += 0.1
			if bunny[8] < 0:
				bunny[8] = 0.0
			if bunny[8] > 1.0:
				bunny[8] = 1.0
		print "Bunny P " + str(bunny[0])
	if frame >= renderStartFrame:
		densityGrid = makeGrid( gridb, 0.0 )
		colorGrid = makeGrid( gridb, Color(0,0,0,0) )
		StampSplineWisps( densityGrid, colorGrid, trail )
		density = gridded(densityGrid)  + teapotdensity
		color = (gridded(colorGrid) / ( density + constant(0.0000001) ))  + teapotcolor
		#density = constant(0)
		#color = constant(Color(0,0,0,0))
		nb = 0
		bunnydensity = constant(0.0)
		bunnycolor = constant(Color(0,0,0,0))
		for bunny in bunnies:
			if frame >= bunny[1]:
				if bunny[0].Y() < 4.0:
					bunnySDF = scale( bunnyData[0], Vector(1,1,1)*bunny[7] )
					bunnySDF = rotate( bunnySDF, bunny[2]*bunny[3] )
					bunnySDF = rotate( bunnySDF, bunny[4]*bunny[5] )
					bunnySDF = translate( bunnySDF, bunny[0] )
					bunnydensity = bunnydensity + clamp( bunnySDF/constant(0.01), 0.0, 1.0 )*constant(10.0 * bunny[8])
					bunnycolor = bunnycolor + constant(bunny[6])*mask(bunnySDF)
					nb = nb + 1
		bundengrid = makeGrid(gridb, 0.0 )
		buncolgrid = makeGrid(gridb, Color(0,0,0,0))
		stamp( bundengrid, bunnydensity, 1 )
		bunnydensity = gridded( bundengrid )
		stamp( buncolgrid, bunnycolor, 1 )
		bunnycolor = gridded( buncolgrid )
		density = density + bunnydensity
		color = color + bunnycolor
		StandardTwoPartRender( constant(0.0), density, constant(Color(0,0,0,0)), color, allparms )

endJob()


