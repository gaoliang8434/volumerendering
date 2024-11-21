#!/usr/bin/python


import sys
import os
from bishoputils import *
from ImplicitShapes import *


beginJob()

#
#=========== pyro construction ============================
#

allparms =  AllParameters(PyroSphereParameters())

CmdLineHelp( "-h" )


allparms["-pyrogamma"] = 0.5
allparms["-amp"] = 1.5
allparms["-octaves"] = 4.0
allparms["-roughness"] = 0.7
allparms["-name"] = "testpyrodisplacement.exr"
allparms["-scatter"] = 10.0
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-size"] = [ 960, 540 ]
allparms["-dsmNxNyNz"] = [ 300, 300, 300 ]
allparms["-dsmlength"] = [ 7,7,7 ]
allparms["-dsmllc"] = [ -3.5, -3.5, -3.5 ]
allparms["-fov"] = 70.0
allparms["-ds"] = 0.01
allparms["-ambient"] = [ 0.0, 0.0, 0.0 ]
allparms["-near"] = 4.0
allparms["-maxpathlength"] = 10.0
allparms["-size"] = [ 960, 540 ]



allparms = CmdLineFindValues( allparms )

ambient = allparms["-ambient"]
ambientColor = constant( Color( float(ambient[0]), float( ambient[1] ), float( ambient[2] ), 1.0  ) ) 
colr = allparms["-color"]
litColor = constant( Color( float(colr[0]), float( colr[1]  ), float( colr[2] ), 1.0  ) )
axis = Vector(1,1,1).unitvector()
basesdf = Torus( Vector(0,0,0), axis, 3.0, 0.25 )


angle = float(int(allparms["-turntableframe"])-1.0)*3.14159265*2.0/100.0
print "Angle " + str(angle)
testrange = math.sin( angle )*1.0/3.0
baseCenter = axis * testrange
displacementNoise = PerlinFractalSum()
noiseparms = Noise_t()
displacementNoise.getParameters( noiseparms )
noiseparms.wavelength = 3.0
noiseparms.octaves = 4.2
noiseparms.translate = axis*testrange
displacementNoise.setParameters( noiseparms )
displacementNoise.getParameters( noiseparms )
print "Noise translate: " + str( noiseparms.translate.X() ) + " " + str( noiseparms.translate.Y() ) + " " + str( noiseparms.translate.Z() )
dispN = pow(abs( SFNoise( displacementNoise ) ),0.5)





ix = identity() - constant(baseCenter)
surfaceOfSphere = unitvector(ix)
displacement = warp( dispN, surfaceOfSphere )*constant(10.1)





sdf = basesdf + displacement
densityField = clamp( sdf/constant(0.1), 0, 1.0 )
#
#========================= End of scene definition =========================
#
PrintParameters( allparms )
StandardRender( densityField, allparms )


