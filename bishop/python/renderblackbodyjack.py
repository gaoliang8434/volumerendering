#!/usr/bin/python


import sys
import os
from bishoputils import *
from ImplicitShapes import *



cparms = CmdLineFindValues( CameraParameters() )
iparms = CmdLineFindValues( ImageParameters() )
rparms = CmdLineFindValues( RenderDataParameters() )
loopparms = CmdLineFindValues( RenderLoopParameters() )
dsmparms = CmdLineFindValues( DSMParameters() )




#
#=========== Jack construction ============================
#
parms = CmdLineFindValues(ImplicitJackParameters())
CmdLineHelp( "-h" )

print "Building jack" 
jack = ImplicitJack(parms)
scatterdensity = ConstantVolume( 0.0 )
bb = BlackBodyVolume( KeepVolume( MultiplyVolume( KeepVolume(MultiplyVolume(jack,jack)), 3000.0 ) ) );
ambientColor = KeepVolume(bb);
jackColor = ConstantColor( Color( 0,0,0, 1.0  ) )


#
#========================= End of scene definition =========================
#

display = dict( cparms.items() + iparms.items() + rparms.items() + loopparms.items() + dsmparms.items() + parms.items() )
PrintParameters( display )


print "Building render objects" 
camera = MakeCamera( cparms )
image = MakeImage( iparms )
renderData = MakeRenderData( scatterdensity, jack, jackColor, ambientColor, rparms )


dsms = MakeDSMs( renderData, dsmparms )
dsmVolumes = PackageDSMs( dsms, renderData )
RenderLoop( camera, image, renderData, loopparms )

OutputImage( image, iparms )

LogIt("CLOSE")
