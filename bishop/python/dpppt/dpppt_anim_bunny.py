#!/usr/bin/python


import sys
import os
from bishoputils import *
from ImplicitShapes import *

def ObjParameters():
    parms = {  "-objfilename":"", "-objscale":1.0, "-objtranslate":[0,0,0], "-levelsetNxNyNz":[100, 100, 100], 
               "-lsthresh":0.1, "-clamp":10.0, "-saveLevelset":"", "-readLevelset":"" }
    return parms

beginJob()


allparms = AllParameters(MergeParameters(PyroSphereParameters(),ObjParameters()))

CmdLineHelp( "-h" )


allparms["-objfilename"] = "../standardTests/cleanbunny.obj"

allparms["-name"] = "dpppt_displ.exr"
allparms["-case"] = "bunnry"
allparms["-version"] = 1
allparms["-scatter"] = 2.5
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-size"] = [ 960, 540 ]
allparms["-dsmNxNyNz"] = [ 300, 300, 300 ]
allparms["-dsmllc"] = [ -2,-2,-2 ]
allparms["-dsmurc"] = [ 2,2,2 ]
allparms["-dsmlength"] = [ 4,4,4 ]
allparms["-dsmcell"] = [ 0.01,0.01,0.01 ]
allparms["-dsmsamples"] = 1
allparms["-dsmdxdydz"] = [ 0.01,0.01,0.01 ]
allparms["-fov"] = 30.0
allparms["-ds"] = 0.01
allparms["-nbthreads"] = 4
allparms["-subsamples"] = 1
allparms["-near"] = 0.1
allparms["-maxpathlength"] = 100.0
allparms["-freq"] = 1.124*1.5
allparms["-octaves"] = 3.0
allparms["-roughness"] = 0.75
allparms["-amp"] = 1.0
allparms["-backgroundcolor"] = [ 0.5,0.5,0.6,0 ]
allparms["-lightstyle"] = "whitekeyrimfill"

allparms = CmdLineFindValues( allparms )

noiseparms = Noise_t()
noiseparms.frequency = float(allparms["-freq"])
noiseparms.roughness = float(allparms["-roughness"] )
noiseparms.amplitude = float(allparms["-amp"])
noiseparms.octaves = float( allparms["-octaves"])

llc = Vector(-2,-2,-2)
urc = Vector(2,2,2)
cellSize = 0.05
vdbgb = makeGridBox( llc, urc, Vector( cellSize, cellSize, cellSize) )
vdbgbout = makeGridBox( llc, urc, Vector(0.005,0.005,0.005 ) )


allparms["-renderllc"] = [ llc.X(), llc.Y(), llc.Z()   ]
allparms["-renderurc"] = [ urc.X(), urc.Y(), urc.Z()   ]
PrintParameters( allparms )


radius = 1.5
#vol = Sphere( Vector(0,0,0), radius )



print "Building obj levelset" 

p = ObjParser()
p.ParseFile( allparms["-objfilename"] )
g = TriangleGeometry()
g.setScaling( allparms["-objscale"] )
fillresult = p.Fill(g)
if fillresult == 0:
    print "Could not read geometry from file " + allparms["-objfilename"]
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

cell = Vector( float(allparms["-dsmcell"][0]),  float(allparms["-dsmcell"][1]), float(allparms["-dsmcell"][2])  )

gb = makeGridBox( llc, urc, cell )
sgrid = makeGrid( gb, 0.0 )




#RayMarchLevelSet( g, sgrid )
#writeGrid( sgrid, "bunnylevelset_hi.sgrid" )
readGrid( sgrid, "bunnylevelset.sgrid" )
vol = gridded(sgrid)






boundingSphere = Sphere( Vector(0,0,0), 0.2 )

#for f in range(1,101):
#    allparms["-turntableframe"] = f
f = int(allparms["-turntableframe"])
noiseparms.translate = Vector(0, 0, -f*0.025 )
noise = perlin( noiseparms )
amp = 0.9
if f<10:
    amp = amp * float(f-1)/9.0
om1 = SFNoise( noise )*constant( amp )
om2 = translate( om1, Vector( 0.7563, -0.2029, 0.3 ) ) 
om3 = translate( om1, Vector( -0.48484, 0.1208, 0.75757 ) ) 
omega = component( om1, om2, om3 )
noiseM = rotation( omega )
XX = gradientDisplacement( noiseM, boundingSphere, cellSize )
Xgrid = makeGrid( vdbgb, Vector(0,0,0) )
stamp( Xgrid, XX, 1 )
XX = gridded( Xgrid ) + identity()
displ = warp( vol, XX )
vdbfilename =  "dpppt_displ_bunny_hi." + formattedFrame(f) + ".vdb"
writeOpenVDB( vdbfilename, -displ, vdbgbout )
#density = clamp( displ/constant(0.02), 0.0, 1.0 )
#dengrid = makeGrid( vdbgbout, 0.0 )
#stamp( dengrid, density, 4 )
#density = gridded( dengrid )
#StandardRender( density, allparms )







endJob()
