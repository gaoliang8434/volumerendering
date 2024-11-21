#!/usr/bin/python


from bishoputils import *
from ImplicitShapes import *
from LevelsetFromGeometry import *



def ObjParameters():
	parms = {  "-objfilename":"", "-objscale":1.0, "-objtranslate":[0,0,0], "-levelsetNxNyNz":[100, 100, 100], 
	           "-saveLevelset":"" }
	return parms




#
#=========== obj construction ============================
#
parms = AllParameters(ObjParameters())
parms = CmdLineFindValues( parms )
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


lssdf = constant(0.0)
lssdf = geom2ls( g, int(levelsetNxNyNz[0]),  int(levelsetNxNyNz[1]), int(levelsetNxNyNz[2]), dims.X(), dims.Y(), dims.Z() )
if parms["-saveLevelset"] != "":
	print "Saving levelset grid to file " + parms["-saveLevelset"]
	WriteFloatVolumeGrid( lssdf.getGrid(), parms["-saveLevelset"] )



