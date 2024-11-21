#!/usr/bin/python


import sys
import os
from bishoputils import *
from ImplicitShapes import *

beginJob()


#
#=========== Jack construction ============================
#
parms = CmdLineFindValues(AllParameters(ImplicitJackParameters()))
parms = Turntable( parms )
parms["-backgroundcolor"] = [ 0.5,0.5,0.6,0.0]
parms["-ds"] = 0.03
parms["-scatter"] = 2.0
parms["-dsmcell"] = [ 0.03, 0.03, 0.03 ]
CmdLineHelp( "-h" )

print "Building jack" 
jack = ImplicitJack(parms)
ambientColor = constant( Color( float(parms["-ambient"][0]), float( parms["-ambient"][1] ), float( parms["-ambient"][2] ), 1.0  ) )
jackColor = constant( Color( float(parms["-color"][0]), float( parms["-color"][1]  ), float( parms["-color"][2] ), 1.0  ) )


#
#========================= End of scene definition =========================
#

PrintParameters( parms )

print jack

print "Building render objects" 
renderData = MakeRenderData( jack, jack, jackColor, ambientColor, parms )
dsms = MakeDSMs( renderData, parms )
dsmVolumes = PackageDSMs( dsms, renderData )
image = RenderVolume( renderData, parms )

OutputImage( image, parms )

endJob()
