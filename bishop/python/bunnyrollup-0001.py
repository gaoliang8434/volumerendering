#!/usr/bin/python


import sys
import os
from bishoputils import *
from ImplicitShapes import *

def ObjParameters():
    parms = {  "-objfilename":"", "-objscale":1.0, "-objtranslate":[0,0,0], "-levelsetNxNyNz":[100, 100, 100], 
               "-lsthresh":0.1, "-clamp":10.0, "-saveLevelset":"", "-readLevelset":"" }
    return parms



def displaceGeometry( geo, vf ):
    nb = geo.nbVertices()
    pm = ProgressMeter( nb, "displace geometry")
    for i in range(0,nb):
        P = geo.getVertex(i)
        X = evaluate(vf, P )
        geo.setVertex(i,X)
        pm.update()



beginJob()


allparms = AllParameters(MergeParameters(PyroSphereParameters(),ObjParameters()))

CmdLineHelp( "-h" )


#allparms["-objfilename"] = "./tubearray_superhires.obj"
allparms["-objfilename"] = "./clean_hires_bunny.obj"
#allparms["-objfilename"] = "../../../ash/models/ajax/ajax.obj"

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
allparms["-freq"] = 1.124*1.5/2.0
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
cellSize = 0.025
vdbgb = makeGridBox( llc, urc, Vector( cellSize, cellSize, cellSize) )
vdbgbout = makeGridBox( llc, urc, Vector(0.005,0.005,0.005 ) )


allparms["-renderllc"] = [ llc.X(), llc.Y(), llc.Z()   ]
allparms["-renderurc"] = [ urc.X(), urc.Y(), urc.Z()   ]
PrintParameters( allparms )


nbfolds = 8
foldfactor = math.pow( 2.0, nbfolds )

boundingSphere = Sphere( Vector(0,0,0), 0.1 )

f = int(allparms["-turntableframe"])
noiseparms.translate = Vector(0, 0, -f*0.025 )
noise = perlin( noiseparms )
amp = (f-1)*0.1/foldfactor
omegadir = Vector( 0,0,1 ) 
xperp = identity() - constant(omegadir) * ( identity()*constant(omegadir) )
omega = constant( omegadir*amp ) * abs(xperp)
noiseM = rotation( omega ) 
XX = gradientDisplacement( noiseM, boundingSphere, 0.025 )
YY = ContinuedFractionDisplacement( XX, 10 )

Ygrid = makeGrid( vdbgb, Vector(0,0,0) )
stamp( Ygrid, YY - identity(), 1 )
YY = gridded( Ygrid ) + identity()
for fold in range(0,nbfolds):
    YY = warp( YY, YY )
    Ygrid = makeGrid( vdbgb, Vector(0,0,0) )
    stamp( Ygrid, YY - identity(), 1 )
    YY = gridded( Ygrid ) + identity()


p = ObjParser()
p.ParseFile( allparms["-objfilename"] )
g = TriangleGeometry()
g.setScaling( allparms["-objscale"] )
fillresult = p.Fill(g)
if fillresult == 0:
    print "Could not read geometry from file " + allparms["-objfilename"]
    sys.exit()
print "Number of vertices = " + str( g.nbVertices() ) + "   number of faces = " + str( g.nbFaces() )
llc = g.LLC()
urc = g.URC()
print "Obj BB:   ",
print str(llc.X()) + " " +  str(llc.Y()) + " " +  str(llc.Z()) + "  X  " +  str(urc.X()) + " " +  str(urc.Y()) + " " +  str(urc.Z())
         
        
displaceGeometry( g, YY )
outfilename = "./bunnyrollup-0001." + formattedFrame( f ) + ".obj"
writeObj( outfilename, g)









endJob()
