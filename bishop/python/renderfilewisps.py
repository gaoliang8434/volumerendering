#!/usr/bin/python

import os
import sys
from bishoputils import *



user = os.getenv("USER")
homeDirectory = "/scratch/" + str(user)




#dataLocation = "/scratch/jeannettemyers/MWG_1p3M_S600_AvgVf2p5/"
#frameMin = 0
#frameMax = 232
dataLocationBase = "/scratch/jeannettemyers/"
dataSet = "Merged_2p0_1p5M_s600"
dataLocation = dataLocationBase + dataSet + "/"
frameMin = 0
frameMax = 474
dataFileBase = "snapshot_"
conversionRoutine = "/scratch/vlbyrd/jeannettemyers/filter_snapshot" 




baseCommandLine = " -volume filewisps -nodsm -clamp 1.0  -stampsize 360 560 160 -stamporigin -180 -280 -80 -fov 12 -near 10 -maxpathlength 800 -sparsesize 20 -stampNxNyNz 3240 5040 1440 -ds 0.075 -eye 0 0 1000 -view 0 0 0 -up 0 1 0 -scatter 1.0 -blurscale 0.0 -colorscale 0.0 1.0 1.0 -raysperpixel 20 -oiiobrightness 3 -oiiogamma 1.0 -magicknoowner  -magickbrightness 3 -magickgamma 1.0 -blurgrid 0 -nbpatches 10000 -useboundingbox -nbdensitysamples 1 " + res + " "






# remove the copied data file
cmd = "rm /scratch/jtessen/" + padded4Frame + "/" + rawDataFile
#print cmd
#os.system(cmd)

convertedFile1 =  "/scratch/jtessen/" + padded4Frame + "/milkywayBulge_particles.txt"
convertedFile2 =  "/scratch/jtessen/" + padded4Frame + "/milkywayDisk_particles.txt"
convertedFile3 =  "/scratch/jtessen/" + padded4Frame + "/sagittariusDisk_particles.txt"

version = "0010"

runTag = "0009"

filesets =             " -particlefile " + convertedFile1 + " -opacityscale  1000000.0 -particlecolor 0.5294 0.3882 0.3333 "
filesets = filesets +  " -particlefile " + convertedFile2 + " -opacityscale  1000000.0 -particlecolor 0.3412 0.4745 0.5843 "
filesets = filesets +  " -particlefile " + convertedFile3 + " -opacityscale  1000000.0 -particlecolor 0.5 0.5 0.0 "

cmd = volumeRenderer + " " + baseCommandLine + filesets + " -magickname " + outputLocation + "/" + outputFile + " -oiioname " + outputLocation + "/" + outputFile

cmd = "rm -rf /scratch/jtessen/" + padded4Frame
#print cmd
#os.system(cmd)











from ImplicitShapes import *

def TurntableParameters():
	parms = { "-first":1, "-last":1, "-nbframes":120, "-increment":1  }
	return parms


# Specific to render the temperature field for Miller's data

#
#=========== text volume construction ============================
#

textvolumeparms = { "-densityfile":"/Volumes/images/MillerData/d3d11890/d3d11890_HRadMass.dat", "-colorfile":"",  "-textNxNyNz":[ 976, 529, 624 ], "-textsize":[ 9.76, 5.29, 6.24 ], "-textllc":[ -4.88, -2.645, -3.12 ],
                    "-densityscale":39.0, "-densityoffset":0.0 }


allparms =  AllParameters( textvolumeparms )
turntableparms = TurntableParameters()
allparms = dict( allparms.items() + turntableparms.items() )



#
# ======= Insert special setting HERE ==============
allparms["-name"] = "hradmass/hradmass.exr"
allparms["-fov"] = 90.0
allparms["-color"] = [0,0,0]
allparms["-ambient"] = [1,1,1]
allparms["-near"] = 4.0
allparms["-maxpathlength"] = 12.0
allparms["-nodsm"] = 1
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

print "Building File Wisp Volume" 

llc = allparms["-textllc"]
tsize = allparms["-textsize"]
tn = allparms["-textNxNyNz"]
X = Vector( float(llc[0]), float(llc[1]), float(llc[2]) )
volumeObject = FloatVolumeGrid()
densityField = ConstantVolume( 0.0 )


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
