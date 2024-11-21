#!/usr/bin/python


import sys
import os
from bishoputils import *
from ImplicitShapes import *

print versionString()


cparms = CmdLineFindValues( CameraParameters() )
iparms = CmdLineFindValues( ImageParameters() )
rparms = CmdLineFindValues( RenderDataParameters() )
loopparms = CmdLineFindValues( RenderLoopParameters() )
dsmparms = CmdLineFindValues( DSMParameters() )




#
#=========== ellipse construction ============================
#

parms = CmdLineFindValues( EllipseParameters() )
CmdLineHelp( "-h" )

print "Building Ellipse" 
volumeObject = Ellipse(parms)
ambient = rparms["-ambient"]
ambientColor = ConstantColor( Color( float(ambient[0]), float( ambient[1] ), float( ambient[2] ), 1.0  ) ) 
colr = rparms["-color"]
litColor = ConstantColor( Color( float(colr[0]), float( colr[1]  ), float( colr[2] ), 1.0  ) )


#
#========================= End of scene definition =========================
#


PrintParameters( cparms )
PrintParameters( iparms )
PrintParameters( rparms )
PrintParameters( loopparms )
PrintParameters( dsmparms )
PrintParameters( parms )

camera = MakeCamera( cparms )
image = MakeImage( iparms )
renderData = MakeRenderData( volumeObject, litColor, ambientColor, rparms )


dsms = MakeDSMs( renderData, dsmparms )
dsmVolumes = PackageDSMs( dsms, renderData )
RenderLoop( camera, image, renderData, loopparms )

OutputImage( image, iparms )
