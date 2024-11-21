#!/usr/bin/python



from bishoputils import *


beginJob()

csgparameters = {"-radius":2.5, "-center":[0.0, 0.0, 0.0],
		 "-clamp":0.1,
		 "-doturntable":0
		 }



parms = AllParameters( csgparameters )
CmdLineHelp("-h")
parms["-ds"] = 0.015
parms["-near"] = 4.0
parms["-maxpathlength"] = 10.0
parms["-fov"] = 65.0
parms["-scatter"] = 5.0
parms["-nbthreads"] = 2
parms["-dsmlenth"] = [ 10, 10, 10 ]
parms["-dsmcell"] = [ 0.02, 0.02, 0.02 ]
parms["-dsmllc"] = [ -5, -5, -5 ]
parms = CmdLineFindValues( parms )
PrintParameters( parms )


##### Sphere 
sphereRadius = float(parms["-radius"])
sphereCenter = Vector( float(parms["-center"][0]),  float(parms["-center"][1]),  float(parms["-center"][2]) )
v1 = Sphere( sphereCenter, sphereRadius ) 

clampv = float( parms["-clamp"])
v1= clamp( v1/constant(clampv), 0.0, 1.0 )

StandardRender( v1, parms )


endJob()

