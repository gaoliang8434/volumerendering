#!/usr/bin/python


import sys
import os
from bishoputils import *
from ImplicitShapes import *

beginJob()


allparms =  AllParameters(PyroSphereParameters())

CmdLineHelp( "-h" )



allparms["-name"] = "dpppt_displ.exr"
allparms["-scatter"] = 1.5
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-size"] = [ 960, 540 ]
allparms["-dsmNxNyNz"] = [ 300, 300, 300 ]
allparms["-dsmllc"] = [ -4,-4,-4 ]
allparms["-dsmurc"] = [ 4,4,4 ]
allparms["-dsmlength"] = [ 8,8,8 ]
allparms["-dsmcell"] = [ 0.01,0.01,0.01 ]
allparms["-dsmsamples"] = 1
allparms["-dsmdxdydz"] = [ 0.02,0.02,0.02 ]
allparms["-fov"] = 60.0
allparms["-ds"] = 0.0025
allparms["-nbthreads"] = 4
allparms["-subsamples"] = 4
allparms["-near"] = 0.1
allparms["-maxpathlength"] = 100.0
allparms["-freq"] = 1.124*0.5
allparms["-octaves"] = 3.0
allparms["-roughness"] = 0.75
allparms["-amp"] = 1.0

allparms = CmdLineFindValues( allparms )

noiseparms = Noise_t()
noiseparms.frequency = float(allparms["-freq"])
noiseparms.roughness = float(allparms["-roughness"] )
noiseparms.amplitude = float(allparms["-amp"])
noiseparms.octaves = float( allparms["-octaves"])

llc = Vector(-4,-4,-4)
urc = Vector(4,4,4)
cellSize = 0.025
vdbgb = makeGridBox( llc, urc, Vector( cellSize, cellSize, cellSize) )
outcellSize = 0.025
vdbgbout = makeGridBox( llc, urc, Vector( outcellSize, outcellSize, outcellSize) )

allparms["-renderllc"] = [ llc.X(), llc.Y(), llc.Z()   ]
allparms["-renderurc"] = [ urc.X(), urc.Y(), urc.Z()   ]
PrintParameters( allparms )


radius = 1.5

vol = Sphere( Vector(0,0,0), radius )
boundingSphere = Sphere( Vector(0,0,0), 0.2 )

#for f in range(2,101):
#    allparms["-turntableframe"] = f
frame = int(allparms["-turntableframe"])
nbfolds = 10
foldfactor = math.pow( 2.0, nbfolds )

noiseparms.translate = Vector(0, 0, -frame*0.1 )
noise = perlin( noiseparms )
om1 = SFNoise( noise )*constant(1.9/foldfactor)
om2 = translate( om1, Vector( 0.7563, -0.2029, 0.3 ) ) 
om3 = translate( om1, Vector( -0.48484, 0.1208, 0.75757 ) ) 
omega = component( om1, om2, om3 )
noiseM = rotation( omega )
XX = gradientDisplacement( noiseM, boundingSphere, 0.05 ) + identity()
Xgrid = makeGrid( vdbgb, Vector(0,0,0) )
stamp( Xgrid, XX - identity(), 1 )
XX = gridded( Xgrid ) + identity()
for f in range(0,nbfolds):
    XX = warp( XX, XX )
    Xgrid = makeGrid( vdbgb, Vector(0,0,0) )
    stamp( Xgrid, XX - identity(), 1 )
    XX = gridded( Xgrid ) + identity()
displ = warp( vol, XX )
fname = "test." + formattedFrame(frame) + ".vdb"
writeOpenVDB( fname, -displ, vdbgbout )
#density = clamp( displ/constant(cellSize), 0.0, 1.0 )
#StandardRender( density, allparms )


endJob()
