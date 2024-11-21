#!/usr/bin/python


import sys
import os
import math
from bishoputils import *
from ImplicitShapes import *

def TurntableParameters():
	parms = { "-first":1, "-last":1, "-nbframes":120, "-increment":1  }
	return parms

#
#=========== text volume construction ============================
#

textvolumeparms = { "-densityfile":"", "-colorfile":"",  "-textNxNyNz":[ 976, 529, 624 ], "-textsize":[ 9.76, 5.29, 6.24 ], "-textllc":[ -4.88, -2.645, -3.12 ],
                    "-densityscale":1.0, "-densityoffset":0.0 }


allparms =  AllParameters( textvolumeparms )
turntableparms = TurntableParameters()
allparms = dict( allparms.items() + turntableparms.items() )



#
# ======= Insert special setting HERE ==============
allparms["-name"] = "d3d11890.exr"
allparms["-fov"] = 90.0

# ==================================================
#





allparms = CmdLineFindValues( allparms )
CmdLineHelp( "-h" )

PrintParameters( allparms )

if allparms["-densityfile"] == "":
	print "No density file given"
	sys.exit()

tn = allparms["-textNxNyNz"]

mean = 0.0
stddev = 0.0
count = 0
maxvalue = -100000.0
minvalue = 1000000.0
if allparms["-densityfile"] != "":
	densityFile = open( allparms["-densityfile"], 'r')
	pm = ProgressMeter( long(int(tn[0])*int(tn[1])*int(tn[2])), "density volume" )
	for z in range(0,int(tn[2])):
		for y in range(0,int(tn[1])):
			for x in range(0,int(tn[0])):
				inputvalue = densityFile.readline()
				inputvalue = inputvalue.replace('D', 'E' )
				value = float( inputvalue ) 
				mean += value
				stddev += value*value
				count = count + 1
				if value > maxvalue:
					maxvalue = value
				if value < minvalue:
					minvalue = value
				pm.update()

mean = float(mean) / float(count)
stddev = math.sqrt( float(stddev)/float(count) -float( mean)*float(mean) )

print "\n\nMean value : " + str(mean)
print "Std dev    : " + str(stddev)
print "Max value  : " + str(maxvalue)
print "Min value  : " + str(minvalue)

