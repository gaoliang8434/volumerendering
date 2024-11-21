#!/usr/bin/python



from bishopUtils import *


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
parms["dsmNxNyNz"] = [ 200, 200, 200 ]
parms = CmdLineFindValues( parms )
parms["-name"] = "sphereSphere_" + parms["-csgoperation"] + ".exr"
parms = Turntable(parms)
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
	v1 = multiply( v1, blendStrength )
	v2 = multiply( v2, blendStrength )
	v3 = BlinnBlendVolume( v1, v2 )
if csgtype == "cutout":
	v3 = cutout( v1, v2 )
if csgtype == "othercutout":
	v3 = cutout( v2, v1 )
if csgtype == "union":
	v3 = Union( v1, v2 )

clampv = float( parms["-clamp"])
v3= clamp( multiply(v3,1.0/clampv), 0.0, 1.0 )

if parms["-doturntable"] == 0:
	StandardRender( v3, parms )
else:
	TurntableRender( v3, parms )


endJob()

