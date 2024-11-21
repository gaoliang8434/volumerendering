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
parms["-name"] = "phasefunctionSphere.exr"
parms["-case"] = "uniform"
PrintParameters( parms )


##### SPhere
sphereRadius = float(parms["-radius"])
sphereCenter = Vector( float(parms["-center"][0]),  float(parms["-center"][1]),  float(parms["-center"][2]) )
v1 = Sphere( sphereCenter, sphereRadius ) 


v5 = clamp( v1/constant(0.1), 0.0, 1.0 )


renderdata = RenderSetup( v5, parms )
parms["-case"] = "uniform"
RenderFromData( renderdata, parms )
SetHenyeyGreensteinPhaseFunction( renderdata, 0.95 )
parms["-case"] = "henyeygreenstein"
RenderFromData( renderdata, parms )
SetDoubleHenyeyGreensteinPhaseFunction( renderdata, 0.95, -0.7, 0.5 )
parms["-case"] = "doublehenyeygreenstein"
RenderFromData( renderdata, parms )
SetFournierForandPhaseFunction( renderdata, 1.05, 4.5 )
parms["-case"] = "fournierforand"
RenderFromData( renderdata, parms )


endJob()

