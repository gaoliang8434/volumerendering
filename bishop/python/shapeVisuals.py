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




llc = Vector(-1,-1,-1)
urc = Vector(1,1,1)
cellSize = float( parms["-ds"])
vdbgb = makeGridBox( llc, urc, Vector( cellSize, cellSize, cellSize) )

v3 = SteinerPatch()
#writeOpenVDB( "steinerpatch.vdb", -v3, vdbgb )

v3 = Cone( Vector(0,0,0), Vector(0,1,0), 0.75, 25.0 )
#writeOpenVDB( "cone.vdb", -v3, vdbgb )

v3 = CappedCylinder( Vector(0,0,0), Vector(0,1,0), 1.5, 0.5 )
#writeOpenVDB( "cappedcylinder.vdb", -v3, vdbgb )

v3 = Ellipse( Vector(0,0,0), Vector(1,0,0), 0.75, 0.35 )
writeOpenVDB( "ellipse.vdb", -v3, vdbgb )




endJob()

