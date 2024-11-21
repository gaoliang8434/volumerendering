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


allparms["-pyrogamma"] = 0.75
allparms["-amp"] = 1.5
allparms["-octaves"] = 6.0
allparms["-roughness"] = 0.45
allparms["-translate"] = [ 0, 0, 0.34 ]
allparms["-name"] = "pyro.exr"
allparms["-scatter"] = 5.0
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-fov"] = 70.0
allparms["-ds"] = 0.035
allparms["-near"] = 4.0
allparms["-radius"] = 2.5
allparms["-amp"] = 1.0
allparms["-maxpathlength"] = 10.0
allparms["-dsmcell"] = [ 0.025, 0.025, 0.025 ] 
allparms["-eye"] = [ 1, 0, 9 ]
allparms["-subsamples"] = 15
allparms["-size"] = [ 960, 540 ]


allparms = CmdLineFindValues( allparms )
PrintParameters( allparms )

ambient = allparms["-ambient"]
ambientColor = constant( Color( float(ambient[0]), float( ambient[1] ), float( ambient[2] ), 1.0  ) ) 
colr = allparms["-color"]
litColor = constant( Color( float(colr[0]), float( colr[1]  ), float( colr[2] ), 1.0  ) )

dx = 0.1

for i in range(1,101):
	print "Frame " + str(i) + "  Building Pyro" 
	allparms["-turntableframe"] = i
	volumeObject = PyroSphere(allparms)
	densityField = clamp( volumeObject/constant(0.1), 0, 1.0 )
	StandardRender( densityField, allparms )
	if i < 50:
		allparms["-translate"][0] = allparms["-translate"][0] - dx
	if i >= 50:
		allparms["-translate"][2] = allparms["-translate"][2] - dx

cmd = "dpaffmpeg -base " + GenerateStandardBasename(allparms) + " -usemask -useslate"
os.system(cmd)
endJob()


