#!/usr/bin/python

import os
import sys
from bishoputils import *


beginJob()



parms = AllParameters({})

parms = CmdLineFindValues(parms)


#curve = targetedPath( Vector(-1,-1,-1), Vector(1,1,1), Vector(0,-10,0), 0, 1.0 )

X0 = Vector(-1,-1,-1)
X1 = Vector(-1,1,-1)
X2 = Vector(1,0.5,1)

curve = splinePath( X0, X1, X2, 0, 0.5, 1.0 )

curveField = SpaceCurveField( curve, 0.25 )


X0 = Vector(-1,-1,-1)
X1 = Vector(1,-1,1)
X2 = Vector(1,0.5,1)

curve = splinePath( X0, X1, X2, 0, 0.5, 1.0 )

curveField = Union(curveField, SpaceCurveField( curve, 0.25 ))






density = clamp( curveField/constant(0.1), 0, 1.0 )

gb = makeGridBox( Vector(-2,-2,-2), Vector(4,4,4), Vector(0.025,0.025,0.025) )
dgrid = makeGrid( gb, 0.0 )
stamp( dgrid, density )
density = gridded(dgrid)

PrintParameters( parms )
StandardRender( density, parms )
#TurntableRender( density, parms )
 
endJob()
