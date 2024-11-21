#!/usr/bin/python

from bishoputils import *


from ImplicitShapes import *


parms = AllParameters( {} )

parms["-name"] = "shadow.exr"
parms = CmdLineFindValues(parms)
CmdLineHelp( "-h" )
PrintParameters( parms )



#
# initialize velocity noise
#
noiseparm = Noise_t()
noiseparm.octaves = 3.0
noiseparm.roughness = 0.9

nmx = perlin( noiseparm )

noiseparm.translate = Vector(0.1, -0.4, 0.222 )

nmy = perlin( noiseparm )

noiseparm.translate = Vector(0.02343, 0.5873, -0.444 )

nmz = perlin( noiseparm )

pfnxf = SFNoise( nmx )
pfnyf = SFNoise( nmy )
pfnzf = SFNoise( nmz )
velocitynoise = component( pfnxf, pfnyf , pfnzf ) 



shadow = CappedCylinder( Vector(0,0,0), Vector(0,1,0), 8.0, 2.0 )
shadowDensity = clamp( shadow/constant(0.1), 0, 1 )



X = identity()
gb = makeGridBox( Vector(-5,-5,-5), Vector(5,5,5), Vector(0.5, 0.5, 0.5) )



for i in range(0,10):
	X = advect( X, velocitynoise, 1.0 )
	Xgrid = makeGrid( gb, Vector(0,0,0) )
	stamp( Xgrid, X-identity(), 1 )
	X = identity() + gridded( Xgrid )
shadowDensity = warp( shadowDensity, X )



print shadowDensity

StandardRender( shadowDensity, parms )

