#!/usr/bin/python

import math
from bishoputils import *
from ImplicitShapes import *
from LevelsetFromGeometry import *



def LSParameters():
	parms = { "-clamp":10.0,  "-readLevelset":"" }
	return parms




#
#=========== obj construction ============================
#
parms = AllParameters(LSParameters())


parms["-fov"] = 20 
parms["-size"] = [ 960, 540 ]
parms["-near"] = 5 
parms["-maxpathlength"] = 10 
parms["-ambient"] = [ 0, 0, 0 ] 
parms["-scatter"] = 10 
parms["-ds"] = 0.005 
parms["-dsmNxNyNz"] = [ 200, 200, 200 ]
parms["-name"] = "levelset.exr"

parms = CmdLineFindValues( parms )
CmdLineHelp( "-h" )

print "Reading density from file " + parms["-readLevelset"]
ls = FloatVolumeGrid()
ReadFloatVolumeGrid( ls, parms["-readLevelset"] )
lssdf = GriddedVolume( ls )


llc = ls.llc()
urc = ls.urc()

center = (llc+urc)/2.0
dims = (urc-llc)


parms["-dsmllc"] = [ llc.X(), llc.Y(), llc.Z() ]
parms["-dsmlength"] = [ dims.X(), dims.Y(), dims.Z() ]


eye = center + Vector( 0, 0, dims.Z() )*3.0
view = Vector(0,0,-1)
parms["-eye"] = [ eye.X(), eye.Y(), eye.Z() ]
parms["-view"] = [ view.X(), view.Y(), view.Z() ]

tanfov = dims.X() / (2.0*dims.Z() )
fov = 2.0 * math.atan( tanfov ) * 180.0/3.14159265
parms["-fov"] = fov

parms["-near"] = 2.0*dims.Z()
parms["-maxpathlength"] = 3.0*dims.Z()
parms["-ds"] = dims.Z() / 30.0

parms = Turntable( parms )
PrintParameters( parms )

#
#========================= End of scene definition =========================
#


clamper = float( parms["-clamp"] )
multiplied = KeepVolume( MultiplyVolume( lssdf, clamper ) )
clamped = ClampVolume( multiplied, 0.0, 1.0 )
renderdensity = KeepVolume(  clamped )



StandardRender( renderdensity, parms )




LogIt("CLOSE")
