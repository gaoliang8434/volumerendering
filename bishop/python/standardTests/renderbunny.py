#!/usr/bin/python


from bishoputils import *
from ImplicitShapes import *

beginJob()

def ObjParameters():
	parms = {  "-objfilename":"", "-objscale":1.0, "-objtranslate":[0,0,0], "-levelsetNxNyNz":[100, 100, 100], 
	           "-lsthresh":0.1, "-clamp":10.0, "-saveLevelset":"", "-readLevelset":"" }
	return parms


#
#=========== obj construction ============================
#
parms = AllParameters(ObjParameters())


parms["-objfilename"] = "cleanbunny.obj"
parms["-fov"] = 20 
parms["-near"] = 5 
parms["-maxpathlength"] = 10 
parms["-ambient"] = [ 0, 0, 0 ] 
parms["-scatter"] = 10 
parms["-ds"] = 0.005 
parms["-dsmcell"] = [  0.03103184/3.0, 0.030708/3.0, 0.02407876/3.0 ]
parms["-dsmllc"] = [ -1.551592, -1.5354, -1.218761 ] 
parms["-dsmlength"] = [ 3.103184, 3.0708, 2.407876 ]
parms["-name"] = "bunnyTurntable.exr"
parms["-subsamples"] = 10

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
dims *= 1.2;
llc = center - dims*0.5;
urc = llc + dims;
print "Obj BB:   ",
print str(llc.X()) + " " +  str(llc.Y()) + " " +  str(llc.Z()) + "  X  " +  str(urc.X()) + " " +  str(urc.Y()) + " " +  str(urc.Z())
print "Obj Size: " + str(dims.X()) + " X " + str(dims.Y()) + " X " + str(dims.Z()) 

cell = Vector( float(parms["-dsmcell"][0]),  float(parms["-dsmcell"][1]),  float(parms["-dsmcell"][2])  )

gb = makeGridBox( llc, urc, cell )
sgrid = makeGrid( gb, 0.0 )

RayMarchLevelSet( g, sgrid )




lssdf = gridded(sgrid)
#
#========================= End of scene definition =========================
#

renderdensity = clamp( lssdf*constant(10.0), 0.0, 1.0 )

StandardRender( renderdensity, parms )

endJob()
