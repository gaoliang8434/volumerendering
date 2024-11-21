#!/usr/bin/python

from bishoputils import *


from ImplicitShapes import *


parms = AllParameters( {} )

parms["-name"] = "particle.exr"
parms = CmdLineFindValues(parms)
CmdLineHelp( "-h" )
PrintParameters( parms )



p = Particle()

print "Particle worked"
