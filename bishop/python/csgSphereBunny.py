#!/usr/bin/python



from bishoputils import *
from LevelsetFromGeometry import *


beginJob()

csgparameters = {"-radius":0.42, "-center":[0.0, 0.15, 0.0], 
		 "-clamp":0.1,
		 "-csgoperation":"union",
		 "-blendstrength":4.0,
		 "-densityscale":1.0
		 }


csgOperationMenu = [ "union", "intersect", "blend", "cutout", "othercutout" ]



def ObjParameters():
	parms = {  "-objfilename":"", "-objscale":1.0, "-objtranslate":[0,0,0], "-levelsetNxNyNz":[100, 100, 100], 
	           "-lsthresh":0.1, "-clamp":10.0, "-saveLevelset":"", "-readLevelset":"" }
	return parms






parms = AllParameters( MergeParameters(csgparameters, ObjParameters() ) )
CmdLineHelp("-h")
parms["-ds"] = 0.005
parms["-near"] = 7.0
parms["-maxpathlength"] = 10.0
parms["-fov"] = 20.0
parms["-dsmNxNyNz"] = [ 200, 200, 200 ]
parms["-dsmllc"] = [ -1.551592, -1.5354, -1.218761 ] 
parms["-dsmlength"] = [ 3.103184, 3.0708, 2.407876 ]
parms["-ambient"] = [0.0, 0.0, 0.0]
parms["-color"] = [1,1,1]
parms["-scatter"] = 10.0
parms["-objfilename"] = "../819/models/cleanbunny.obj"
parms["-readLevelset"] = "bunnyls400.grid"
parms = CmdLineFindValues( parms )
parms["-name"] = "sphereBunny_" + parms["-csgoperation"] + ".exr"
PrintParameters( parms )



p = ObjParser()
p.ParseFile( parms["-objfilename"] )
g = TriangleGeometry()
g.setScaling( float(parms["-objscale"]) )
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


##### SPhere
sphereRadius = float(parms["-radius"])
sphereCenter = Vector( float(parms["-center"][0]),  float(parms["-center"][1]),  float(parms["-center"][2]) )
v1 = SphereVolume( sphereCenter, sphereRadius ) 

##### Bunny 
v2 = ConstantVolume(0.0)
if parms["-readLevelset"] != "":
	print "Reading density from file " + parms["-readLevelset"]
	ls = FloatVolumeGrid()
	ReadFloatVolumeGrid( ls, parms["-readLevelset"] )
	v2 = GriddedVolume( ls )
else:
	v2 = LevelsetFromGeometry( g, int(levelsetNxNyNz[0]),  int(levelsetNxNyNz[1]), int(levelsetNxNyNz[2]), dims.X(), dims.Y(), dims.Z() )
	if parms["-saveLevelset"] != "":
		print "Saving levelset grid to file " + parms["-saveLevelset"]
		WriteFloatVolumeGrid( lssdf.getGrid(), parms["-saveLevelset"] )






v3 = ConstantVolume(0.0)
csgtype = str( parms["-csgoperation"] )




if csgtype == "intersect":
	v3 =  IntersectionVolume( v1, v2 )
if csgtype == "blend":
	blendStrength = float( parms["-blendstrength"] )
	v1b = KeepVolume( MultiplyVolume( v1, blendStrength ) )
	v2b = KeepVolume( MultiplyVolume( v2, blendStrength ) )
	v3 = BlinnBlendVolume( v1b, v2b )
if csgtype == "cutout":
	v3 = CutoutVolume( v1, v2 )
if csgtype == "othercutout":
	v3 = CutoutVolume( v2, v1 )
if csgtype == "union":
	v3 = UnionVolume( v1, v2 )

clamp = float( parms["-clamp"])
v4 = MultiplyVolume( v3, 1.0/clamp )
v5 = ClampVolume( v4, 0.0, 1.0 )


StandardRender( v5, parms )


endJob()

