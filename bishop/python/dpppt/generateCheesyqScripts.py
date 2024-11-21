#!/usr/bin/env python

import os
import sys
from bishoputils import *

beginJob()

parms = { "-script":"junk.py", "-start":1, "-end":10, "-increment":1, "-fname":"test" }


CmdLineHelp("-h")

parms = CmdLineFindValues( parms )


start = int( parms["-start"])
end   = int( parms["-end"])
increment = int( parms["-increment"])

frames = range( start, end+1, increment )

cwd = os.getcwd()

for f in frames:
    fname = parms["-fname"] + "." + formattedFrame(f) + ".sh"
    ff = open(fname, 'w')
    if ff:
        ff.write("#!/bin/bash\n")
        ff.write("source ~/.bashrc\n")
        ff.write("oldpipeup\n")
        ff.write("cd " + cwd + "\n")
        ff.write(cwd + "/" + parms["-script"] + " -turntableframe " + str(f) + "\n")
        ff.close()
        os.system("chmod 777 " + fname)


endJob()
