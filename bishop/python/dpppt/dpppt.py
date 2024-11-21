#!/usr/bin/python


import sys
import os
from bishoputils import *
from ImplicitShapes import *




allparms =  AllParameters(PyroSphereParameters())

CmdLineHelp( "-h" )



allparms["-pyrogamma"] = 0.2
allparms["-amp"] = 1.0
allparms["-octaves"] = 6.0
allparms["-roughness"] = 0.5
allparms["-name"] = "pyro.exr"
allparms["-scatter"] = 2.5
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-size"] = [ 960, 540 ]
allparms["-dsmNxNyNz"] = [ 300, 300, 300 ]
allparms["-dsmllc"] = [ -4,-4,-4 ]
allparms["-dsmurc"] = [ 4,4,4 ]
allparms["-dsmlength"] = [ 8,8,8 ]
allparms["-dsmcell"] = [ 0.01,0.01,0.01 ]
allparms["-dsmdxdydz"] = [ 0.02,0.02,0.02 ]
allparms["-fov"] = 80.0
allparms["-ds"] = 0.05
allparms["-near"] = 0.1
allparms["-maxpathlength"] = 100.0
allparms["-freq"] = 1.124*0.5
allparms["-octaves"] = 3.0
allparms["-roughness"] = 0.75
allparms["-amp"] = 1.0
allparms["-pyrogamma"] = 0.75

allparms = CmdLineFindValues( allparms )

noiseparms = Noise_t()
noiseparms.frequency = float(allparms["-freq"])
noiseparms.roughness = float(allparms["-roughness"] )
noiseparms.amplitude = float(allparms["-amp"])
noiseparms.octaves = float( allparms["-octaves"])

llc = Vector(-3,-3,-3)
urc = Vector(3,3,3)
cellSize = float(allparms["-ds"])
vdbgb = makeGridBox( llc, urc, Vector( cellSize, cellSize, cellSize) )

allparms["-renderllc"] = [ llc.X(), llc.Y(), llc.Z()   ]
allparms["-renderurc"] = [ urc.X(), urc.Y(), urc.Z()   ]
PrintParameters( allparms )


radius = 1.5

vol = Sphere( Vector(0,0,0), radius )

noise = perlin( noiseparms )
om1 = SFNoise( noise )*constant(0.6)
om2 = translate( om1, Vector( 0.7563, -0.2029, 0.3 ) ) 
om3 = translate( om1, Vector( -0.48484, 0.1208, 0.75757 ) ) 

omega = component( om1, om2, om3 )
#omega = omega * clamp( (vol + constant(radius/3.0))/constant(5.0*cellSize), 0.000001, 1.0 )




noiseM = rotation( omega )


boundingSphere = Sphere( Vector(0,0,0), 0.2 )

XX = gradientDisplacement( noiseM, boundingSphere, 0.05 )
Xgrid = makeGrid( vdbgb, Vector(0,0,0) )
stamp( Xgrid, XX, 1 )
XX = gridded( Xgrid ) + identity()



error = det( grad(XX) ) - det(noiseM)

x = -3.0
for i in range(0,101):
    xv = Vector( x,0,0 )
    value = evaluate( error, xv )
    displacement = evaluate( XX-identity(), xv )
    print str(x) + "  " + str(value) + "      " + str(displacement)
    x = x + 6.0/100.0



#writeOpenVDB( "dpppt_sphere_error.vdb", -error, vdbgb )

displ = warp( vol, XX )
writeOpenVDB( "dpppt_sphere_disp.vdb", -displ, vdbgb )

dispscale = 0.43445
omega = scale( omega, Vector(dispscale,dispscale,dispscale) )
omega = rotate( omega, Vector(0,1.7,0) )
noiseM = rotation( omega )
XX2 = gradientDisplacement( noiseM, boundingSphere, 0.02 )
Xgrid = makeGrid( vdbgb, Vector(0,0,0) )
stamp( Xgrid, XX2, 1 )
XX2 = gridded( Xgrid ) + identity()


displ2 = warp( warp( vol, XX2 ), XX )
writeOpenVDB( "dpppt_sphere_disp2.vdb", -displ2, vdbgb )

#density = clamp( displ/constant(0.05), 0.0, 1.0 )
#allparms["-name"] = "dpppt_displ.exr"
#StandardRender( density, allparms )

