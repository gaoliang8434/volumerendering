#!/usr/bin/python


import sys
import os
from bishoputils import *
from ImplicitShapes import *
from LevelsetFromGeometry import *


tetmeter = ProgressMeter(1, "TOTAL EXECUTION TIME" )

def ObjParameters():
	parms = {  "-objfilename":"", "-objscale":1.0, "-objtranslate":[0,0,0], "-levelsetNxNyNz":[100, 100, 100], 
	           "-lsthresh":0.1, "-clamp":10.0, "-saveLevelset":"", "-readLevelset":"" }
	return parms




#
#=========== obj construction ============================
#
parms = CmdLineFindValues(Turntable(AllParameters(ObjParameters())))
CmdLineHelp( "-h" )


PrintParameters( parms )
print "Building obj levelset" 



p = ObjParser()
p.ParseFile( parms["-objfilename"] )
g = TriangleGeometry()
g.setScaling( parms["-objscale"] )
fillresult = p.Fill(g)
if fillresult == 0:
	print "Could not read geometry from file " + parms["-objfilename"]
	sys.exit()
print "Number of vertices = " + str( g.nbVertices() ) + "   number of faces = " + str( g.nbFaces() )
dims = g.URC() - g.LLC()
center = (g.URC() + g.LLC())*0.5;
dims *= 2.0;
llc = center - dims*0.5;
urc = llc + dims;
print "Obj BB:   ",
print str(llc.X()) + " " +  str(llc.Y()) + " " +  str(llc.Z()) + "  X  " +  str(urc.X()) + " " +  str(urc.Y()) + " " +  str(urc.Z())
print "Obj Size: " + str(dims.X()) + " X " + str(dims.Y()) + " X " + str(dims.Z()) 
levelsetNxNyNz = parms["-levelsetNxNyNz"]


lssdf = ConstantVolume(0.0)
if parms["-readLevelset"] != "":
	print "Reading density from file " + parms["-readLevelset"]
	ls = FloatVolumeGrid()
	ReadFloatVolumeGrid( ls, parms["-readLevelset"] )
	lssdf = GriddedVolume( ls )
else:
	lssdf = LevelsetFromGeometry( g, int(levelsetNxNyNz[0]),  int(levelsetNxNyNz[1]), int(levelsetNxNyNz[2]), dims.X(), dims.Y(), dims.Z() )
	if parms["-saveLevelset"] != "":
		print "Saving levelset grid to file " + parms["-saveLevelset"]
		WriteFloatVolumeGrid( lssdf.getGrid(), parms["-saveLevelset"] )





ambientColor = ConstantColor( Color( float(parms["-ambient"][0]), float( parms["-ambient"][1] ), float( parms["-ambient"][2] ), 1.0  ) )
objColor = ConstantColor( Color( float(parms["-color"][0]), float( parms["-color"][1]  ), float( parms["-color"][2] ), 1.0  ) )


#
#========================= End of scene definition =========================
#



print "Building render object" 
camera = MakeCamera( parms )
image = MakeImage( parms )

multiplied = KeepVolume( MultiplyVolume( lssdf, 10.0 ) )
clamped = ClampVolume( multiplied, 0.0, 1.0 )
renderdensity = KeepVolume(  clamped )
renderData = MakeRenderData( renderdensity, renderdensity, objColor, ambientColor, parms )


dsms = MakeDSMs( renderData, parms )
dsmVolumes = PackageDSMs( dsms, renderData )


RenderLoop( camera, image, renderData, parms )

OutputImage( image, parms )


tetmeter.update()

LogIt("CLOSE")
