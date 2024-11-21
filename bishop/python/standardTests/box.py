#!/usr/bin/python



from bishoputils import *


beginJob()

csgparameters = {"-radius":2.5, "-center":[0.0, 0.0, 0.0], "-power":4.0,
		 "-clamp":0.05,
		 "-doturntable":0
		 }



parms = AllParameters( csgparameters )
CmdLineHelp("-h")
parms["-ds"] = 0.035
parms["-near"] = 4.0
parms["-maxpathlength"] = 10.0
parms["-fov"] = 55.0
parms["-scatter"] = 5.0
parms["-nbthreads"] = 2
parms["-dsmlenth"] = [ 10, 10, 10 ]
parms["-dsmcell"] = [ 0.05, 0.05, 0.05 ]
parms["-dsmllc"] = [ -5, -5, -5 ]
parms["-eye"] = [ 10, 3, 10 ]
parms = CmdLineFindValues( parms )
PrintParameters( parms )


##### Object 
radius = float(parms["-radius"])
center = Vector( float(parms["-center"][0]),  float(parms["-center"][1]),  float(parms["-center"][2]) )
pwr = float(parms["-power"])
v1 = CsgBox( center, radius, pwr ) 

clampv = float( parms["-clamp"])
v1= clamp( v1/constant(clampv), 0.0, 1.0 )

StandardRender( v1, parms )


endJob()

