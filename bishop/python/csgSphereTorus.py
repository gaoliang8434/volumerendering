#!/usr/bin/python



from bishoputils import *


beginJob()

csgparameters = {"-radius":2.5, "-center":[0.0, 0.0, 0.0], 
                 "-majorradius":3.0, "-minorradius":0.5, 
		 "-toruscenter":[1.5, 0.0, 0.0 ], "-torusaxis":[ 1.0, 1.0, 1.0 ],
		 "-clamp":0.1,
		 "-csgoperation":"union",
		 "-blendstrength":4.0
		 }


csgOperationMenu = [ "union", "intersect", "blend", "cutout", "othercutout" ]

parms = AllParameters( csgparameters )
CmdLineHelp("-h")
parms["-ds"] = 0.015
parms["-near"] = 4.0
parms["-maxpathlength"] = 10.0
parms["-fov"] = 65.0
parms["dsmNxNyNz"] = [ 200, 200, 200 ]
parms = CmdLineFindValues( parms )
parms["-name"] = "sphereTorus_" + parms["-csgoperation"] + ".exr"
PrintParameters( parms )


##### SPhere
sphereRadius = float(parms["-radius"])
sphereCenter = Vector( float(parms["-center"][0]),  float(parms["-center"][1]),  float(parms["-center"][2]) )
v1 = SphereVolume( sphereCenter, sphereRadius ) 

##### Torus
torusMajorRadius = float( parms["-majorradius"] )
torusMinorRadius = float( parms["-minorradius"] )
torusCenter = Vector( float(parms["-toruscenter"][0]),  float(parms["-toruscenter"][1]),  float(parms["-toruscenter"][2]) )
torusAxis = Vector( float(parms["-torusaxis"][0]),  float(parms["-torusaxis"][1]),  float(parms["-torusaxis"][2]) )
v2 = TorusVolume( torusCenter, torusAxis, torusMajorRadius, torusMinorRadius ) 

v3 = ConstantVolume(0.0)
csgtype = str( parms["-csgoperation"] )




if csgtype == "intersect":
	v3 =  IntersectionVolume( v1, v2 )
if csgtype == "blend":
	blendStrength = float( parms["-blendstrength"] )
	v1b = KeepVolume( MultiplyVolume( v1, blendStrength ) )
	v2b = KeepVolume( MultiplyVolume( v2, blendStrength ) )
	v3 = BlinnBlendVolume( v1b, v2b )
if csgtype == "cutout":
	v3 = CutoutVolume( v1, v2 )
if csgtype == "othercutout":
	v3 = CutoutVolume( v2, v1 )
if csgtype == "union":
	v3 = UnionVolume( v1, v2 )

clamp = float( parms["-clamp"])
v4 = MultiplyVolume( v3, 1.0/clamp )
v5 = ClampVolume( v4, 0.0, 1.0 )


StandardRender( v5, parms )


endJob()

