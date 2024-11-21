#!/usr/bin/python


import sys
import os
from bishoputils import *
from ImplicitShapes import *



beginJob()

#
#=========== pyro construction ============================
#

allparms = AllParameters( {} )

CmdLineHelp( "-h" )

allparms["-pyrogamma"] = 0.75
allparms["-amp"] = 1.5
allparms["-octaves"] = 6.0
allparms["-roughness"] = 0.45
allparms["-translate"] = [ 0, 0, 0.34 ]
allparms["-name"] = "piecewise.exr"
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
allparms["-backgroundcolor"] = [ 0.5, 0.5, 0.5, 0.0 ]
allparms = CmdLineFindValues( allparms )


center = [0, 0, 0]
radii = []
centers = []
amp = allparms["-amp"]
colors = []
noiseparms = Noise_t()
noiseparms.radius = 0.25

chain = AnchorChain()
thickness = 0.25
springradius = 2.0
dz = 3.0*thickness
nbVerts = 36
nbTurns = 4

ddz = dz/(nbVerts-1.0)
z = -2.0
for j in range(0,nbTurns):
	for i in range(0,nbVerts):
		theta = i * 2.0 * 3.14159265/nbVerts
		ct = math.cos(theta)*springradius
		st = math.sin(theta)*springradius
		noiseparms.P = Vector(ct, z, st)
		noiseparms.radius = thickness
		chain.push_back( noiseparms )
		z = z + ddz


allparms = CmdLineFindValues( allparms )
PrintParameters( allparms )


llc = Vector( float(allparms["-dsmllc"][0]), float(allparms["-dsmllc"][1]), float(allparms["-dsmllc"][2]) )
urc = llc + Vector( float(allparms["-dsmlength"][0]), float(allparms["-dsmlength"][1]), float(allparms["-dsmlength"][2]) )
res = Vector( 0.01,0.01,0.01 )
gridb = makeGridBox( llc, urc, res )



volumeObject = PiecewiseCurveField( chain )
densityGrid = makeGrid( gridb, 0.0 )
stamp( densityGrid, volumeObject, 1 )
volumeObject = gridded(densityGrid)
densityField = clamp( volumeObject/constant(0.1), 0, 1.0 )
StandardRender( densityField, allparms )

endJob()


