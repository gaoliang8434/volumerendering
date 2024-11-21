#!/usr/bin/python



from bishoputils import *


beginJob()

csgparameters = {"-radius":2.5, "-center":[0.0, 0.0, 0.0], "-separation":2.0,
		 "-clamp":0.1,
		 "-csgoperation":"union",
		 "-blendstrength":4.0,
		 "-doturntable":0
		 }


csgOperationMenu = [ "union", "intersect", "blend", "cutout", "othercutout" ]

parms = AllParameters( csgparameters )
CmdLineHelp("-h")
parms["-ds"] = 0.015
parms["-near"] = 4.0
parms["-maxpathlength"] = 10.0
parms["-fov"] = 65.0
parms["-scatter"] = 5.0
parms["-dsmlenth"] = [ 10, 10, 10 ]
parms["-dsmcell"] = [ 0.02, 0.02, 0.02 ]
parms["-dsmllc"] = [ -5, -5, -5 ]
parms = CmdLineFindValues( parms )
PrintParameters( parms )


##### SPhere 1
sphereRadius = float(parms["-radius"])
sphereCenter = Vector( float(parms["-center"][0]) - float(parms["-separation"])/2.0,  float(parms["-center"][1]),  float(parms["-center"][2]) )
v1 = Sphere( sphereCenter, sphereRadius ) 

##### SPhere 2
sphereCenter = Vector( float(parms["-center"][0]) + float(parms["-separation"])/2.0,  float(parms["-center"][1]),  float(parms["-center"][2]) )
v2 = Sphere( sphereCenter, sphereRadius ) 


v3 = constant(0.0)
csgtype = str( parms["-csgoperation"] )




if csgtype == "intersect":
	v3 =  intersection( v1, v2 )
if csgtype == "blend":
	blendStrength = float( parms["-blendstrength"] )
	v1 = v1*constant(blendStrength)
	v2 = v2*constant(blendStrength )
	v3 = BlinnBlend( v1, v2 )
if csgtype == "cutout":
	v3 = cutout( v1, v2 )
if csgtype == "othercutout":
	v3 = cutout( v2, v1 )
if csgtype == "union":
	v3 = Union( v1, v2 )

clampv = float( parms["-clamp"])
v3= clamp( v3/constant(clampv), 0.0, 1.0 )

StandardRender( v3, parms )


endJob()

