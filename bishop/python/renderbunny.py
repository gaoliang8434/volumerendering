#!/usr/bin/python


from bishoputils import *
from ImplicitShapes import *
from LevelsetFromGeometry import *



def ObjParameters():
	parms = {  "-objfilename":"", "-objscale":1.0, "-objtranslate":[0,0,0], "-levelsetNxNyNz":[100, 100, 100], 
	           "-lsthresh":0.1, "-clamp":10.0, "-saveLevelset":"", "-readLevelset":"" }
	return parms




#
#=========== obj construction ============================
#
parms = AllParameters(ObjParameters())


parms["-objfilename"] = "../819/models/cleanbunny.obj"
parms["-levelsetNxNyNz"] = [ 400, 400, 400 ] 
parms["-fov"] = 20 
parms["-size"] = [ 960, 540 ]
parms["-near"] = 5 
parms["-maxpathlength"] = 10 
parms["-ambient"] = [ 0, 0, 0 ] 
parms["-scatter"] = 10 
parms["-readLevelset"] = "bunnyls400.grid" 
parms["-ds"] = 0.005 
parms["-dsmNxNyNz"] = [ 200, 200, 200 ]
parms["-dsmllc"] = [ -1.551592, -1.5354, -1.218761 ] 
parms["-dsmlength"] = [ 3.103184, 3.0708, 2.407876 ]
parms["-name"] = "bunnyTurntable.exr"

parms = CmdLineFindValues( parms )
parms = Turntable( parms )
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






#
#========================= End of scene definition =========================
#



multiplied = KeepVolume( MultiplyVolume( lssdf, 10.0 ) )
clamped = ClampVolume( multiplied, 0.0, 1.0 )
renderdensity = KeepVolume(  clamped )



StandardRender( renderdensity, parms )




LogIt("CLOSE")
