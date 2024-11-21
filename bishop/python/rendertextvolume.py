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


ambient = allparms["-ambient"]
ambientColor = ConstantColor( Color( float(ambient[0]), float( ambient[1] ), float( ambient[2] ), 1.0  ) ) 
colr = allparms["-color"]
litColor = ConstantColor( Color( float(colr[0]), float( colr[1]  ), float( colr[2] ), 1.0  ) )

print "Building Text Volume" 

llc = allparms["-textllc"]
tsize = allparms["-textsize"]
tn = allparms["-textNxNyNz"]
X = Vector( float(llc[0]), float(llc[1]), float(llc[2]) )
volumeObject = FloatVolumeGrid()
densityField = ConstantVolume( 0.0 )

scale = float( allparms["-densityscale"] )
offset = float( allparms["-densityoffset"] )
if allparms["-densityfile"] != "":
	volumeObject.init( int(tn[0]), int(tn[1]), int(tn[2]), float(tsize[0]), float(tsize[1]), float(tsize[2]), X )
	densityFile = open( allparms["-densityfile"], 'r')
	pm = ProgressMeter( long(int(tn[0])*int(tn[1])*int(tn[2])), "density volume" )
	for z in range(0,int(tn[2])):
		for y in range(0,int(tn[1])):
			for x in range(0,int(tn[0])):
				inputvalue = densityFile.readline()
				inputvalue = inputvalue.replace('D', 'E' )
				value = ( float( inputvalue ) + offset ) * scale
				volumeObject.set(x,y,z, value)
				pm.update()
	densityField = KeepVolume( GriddedVolume( volumeObject ) )


#
#========================= End of scene definition =========================
#

renderData = MakeRenderData( densityField, litColor, ambientColor, allparms )
dsms = MakeDSMs( renderData, allparms )
dsmVolumes = PackageDSMs( dsms, renderData )

first = int( allparms["-first"] )
last = int( allparms["-last"] )
increment = int( allparms["-increment"] )
nbframes = int(allparms["-nbframes"]) 
print "first frame: " + str(first) + "    last frame: " + str(last)
if nbframes <= 0 or last < first:
	print "Bad frame range"
	sys.exit()

cameraStart = Vector( float(allparms["-eye"][0]), float(allparms["-eye"][1]), float(allparms["-eye"][2] ) )
dtheta = 2.0*3.14159265/nbframes
theta = dtheta * (first-1)

frameNameComponents = os.path.splitext(str( allparms["-name"] ))
for frame in range ( first, last+1, increment ):
	padframe = str(frame)
	if frame < 1000:
		padframe = "0" + padframe
	if frame < 100:
		padframe = "0" + padframe
	if frame < 10:
		padframe = "0" + padframe
	print "FRAME " + padframe

	eyeX =  cameraStart.X() * math.cos( theta ) + cameraStart.Z() * math.sin( theta )
	eyeZ = -cameraStart.X() * math.sin( theta ) + cameraStart.Z() * math.cos( theta )
	allparms["-eye"][0] = eyeX
	allparms["-eye"][2] = eyeZ
	camera = MakeCamera( allparms )
	image = MakeImage( allparms )
	RenderLoop( camera, image, renderData, allparms )
	frameName = frameNameComponents[0] + "." + padframe + frameNameComponents[1]
	allparms["-name"] = frameName
	PrintParameters( allparms )
	OutputImage( image, allparms )
	theta += dtheta * increment
