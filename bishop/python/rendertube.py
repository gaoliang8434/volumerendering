#!/usr/bin/python


import sys
import os
from bishoputils import *
from ImplicitShapes import *

print versionString()


#
#=========== tube construction ============================
#

allparms =  AllParameters(TubeParameters())
CmdLineHelp( "-h" )

allparms["-name"] = "tube.exr"

allparms = CmdLineFindValues( allparms )

print "Building Tube" 
volumeObject = Tube(allparms)


#,
#========================= End of scene definition =========================
#

StandardRender( volumeObject, allparms )
