#!/usr/bin/python


import sys
import os
from bishoputils import *
from ImplicitShapes import *

parms = CmdLineFindValues(Turntable(AllParameters(ImplicitJackParameters())))
CmdLineHelp( "-h" )



p1 = Vector( 0,0,0)
p2 = Vector( 1,0,0 )
p3 = Vector( 0,0,1 )
p4 = Vector( 0,1,0 )

tri1 = Triangle( p1, p2, p3 )
tri2 = Triangle( p1, p4, p2 )
tri3 = Triangle( p1, p3, p4 )
tri4 = Triangle( p2, p4, p3 )



sdf = SignedDistance( tri1 )
sdf.setMaxDistance( 123456789.0 )
sdf.addTriangle( tri2 )
sdf.addTriangle( tri3 )
sdf.addTriangle( tri4 )

clampvalue = 2.0
vvol = ConstantVolume( 16.3 );
scaledSDF = KeepVolume(MultiplyVolume( sdf, float(clampvalue) ))
objdensity = KeepVolume(ClampVolume( sdf, 0.0, 1.0 ))

center = (p1 + p2 + p3 + p4)/4.0

y = 0.15
z = -10.0 
x = 0.25
dz = 0.01
while z < 10.0:
	P = center + Vector( z,y,z )
	ssdf = scaledSDF.eval(P)
	den = objdensity.eval(P)
	s = sdf.eval(P)
	print "profile " + str(z) + " " + str(s) + " "  + str(den) + " " + str(ssdf)
	z = z + dz

sys.exit()

ambientColor = ConstantColor( Color( float(parms["-ambient"][0]), float( parms["-ambient"][1] ), float( parms["-ambient"][2] ), 1.0  ) )
objColor = ConstantColor( Color( float(parms["-color"][0]), float( parms["-color"][1]  ), float( parms["-color"][2] ), 1.0  ) )


#
#========================= End of scene definition =========================
#

PrintParameters( parms )


print "Building render object" 
camera = MakeCamera( parms )
image = MakeImage( parms )
renderData = MakeRenderData( objdensity, objdensity, objColor, ambientColor, parms )
#renderData = MakeRenderData( objdensity, objdensity, objColor, ambientColor, parms )


dsms = MakeDSMs( renderData, parms )
dsmVolumes = PackageDSMs( dsms, renderData )
RenderLoop( camera, image, renderData, parms )

OutputImage( image, parms )

LogIt("CLOSE")
