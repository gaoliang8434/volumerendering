#!/usr/bin/python



from bishoputils import *


beginJob()

csgparameters = {"-radius":2.5, "-center":[0.0, 0.0, 0.0], "-separation":2.0,
		 "-csgoperation":"union",
		 "-blendstrength":4.0
		 }


csgOperationMenu = [ "union", "intersect", "blend", "cutout", "othercutout" ]

parms = AllParameters( csgparameters )
CmdLineHelp("-h")
parms["-ds"] = 0.1
parms = CmdLineFindValues( parms )
parms["-name"] = "Spheres_" + parms["-csgoperation"] + ".vdb"
parms = Turntable(parms)
PrintParameters( parms )


##### SPhere 1
sphereRadius = float(parms["-radius"])
sphereCenter1 = Vector( float(parms["-center"][0]) - float(parms["-separation"])/2.0,  float(parms["-center"][1]),  float(parms["-center"][2]) )
v1 = Sphere( sphereCenter1, sphereRadius ) 

##### SPhere 2
sphereCenter2 = Vector( float(parms["-center"][0]) + float(parms["-separation"])/2.0,  float(parms["-center"][1]),  float(parms["-center"][2]) )
v2 = Sphere( sphereCenter2, sphereRadius ) 


v3 = constant(0.0)
csgtype = str( parms["-csgoperation"] )




if csgtype == "intersect":
	v3 =  intersection( v1, v2 )
if csgtype == "blend":
	blendStrength = float( parms["-blendstrength"] )
	v1 = multiply( v1, blendStrength )
	v2 = multiply( v2, blendStrength )
	v3 = BlinnBlend( v1, v2 )
if csgtype == "cutout":
	v3 = cutout( v1, v2 )
if csgtype == "othercutout":
	v3 = cutout( v2, v1 )
if csgtype == "union":
	v3 = Union( v1, v2 )



boxCenter = (sphereCenter1 + sphereCenter2)*0.5
boxLength = ( (sphereCenter1-sphereCenter2).magnitude() + 2.0*sphereRadius )*1.2
llc = boxCenter - Vector(1,1,1)*boxLength*0.5
urc = boxCenter + Vector(1,1,1)*boxLength*0.5
cellSize = float( parms["-ds"])
vdbgb = makeGridBox( llc, urc, Vector( cellSize, cellSize, cellSize) )
writeOpenVDB( parms["-name"], -v3, vdbgb )





endJob()

