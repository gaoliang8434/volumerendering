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
parms["-name"] = "blendVisual.vdb"
parms = Turntable(parms)
PrintParameters( parms )


##### Shape 1
v1 = CsgBox( Vector(0,0,0), 1.0, 3.0 ) 

##### Shape 2
v2 = Torus( Vector(-1.75,0,0), Vector(0,0,1), 2.0, 0.4 ) 

v3 =  Union( v1, v2 )
#blendStrength = float( parms["-blendstrength"] )
#v1 = multiply( v1, blendStrength )
#v2 = multiply( v2, blendStrength )
#v3 = BlinnBlend( v1, v2 )



llc = Vector(-5,-5,-5)
urc = Vector(5,5,5)
cellSize = float( parms["-ds"])
vdbgb = makeGridBox( llc, urc, Vector( cellSize, cellSize, cellSize) )
#writeOpenVDB( "unblended.vdb", -v3, vdbgb )



v3 = BlinnBlend( multiply(v1,1.0), multiply(v2, 1.0) )
writeOpenVDB( "blinnblended.vdb", -v3, vdbgb )



endJob()

