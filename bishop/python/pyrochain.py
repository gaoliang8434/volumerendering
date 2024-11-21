#!/usr/bin/python


import sys
import os
from bishoputils import *
from ImplicitShapes import *



beginJob()

#
#=========== pyro construction ============================
#

allparms =  AllParameters(PyroSphereParameters())
allparms = MergeParameters( allparms, {"-separationfraction":3.0, "-chainlength":15, "-chainroughness":0.98, "-pyrovariance":0.25, "-maxpyrovariance":2.0, "-usecolor":"no"} )

CmdLineHelp( "-h" )

allparms["-pyrogamma"] = 0.75
allparms["-amp"] = 1.5
allparms["-octaves"] = 6.0
allparms["-roughness"] = 0.45
allparms["-translate"] = [ 0, 0, 0.34 ]
allparms["-name"] = "pyrochain.exr"
allparms["-scatter"] = 5.0
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-fov"] = 70.0
allparms["-ds"] = 0.025
allparms["-near"] = 19.0
allparms["-radius"] = 2.5
allparms["-amp"] = 1.0
allparms["-maxpathlength"] = 10.0
allparms["-dsmcell"] = [ 0.05, 0.05, 0.05 ] 
allparms["-dsmlength"] = [ 28, 8, 8 ]
allparms["-dsmllc"] = [ -24, -4, -4 ]
allparms["-view"] = list(allparms["-center"])

allparms = CmdLineFindValues( allparms )


meanradius = allparms["-radius"]
lognormal = LognormalPRN()
noiseparms = Noise_t()
noiseparms.lognormalmean = float(1.0)
noiseparms.gaussianstandarddeviation = float( allparms["-pyrovariance"])
lognormal.setParameters(noiseparms)
center = list(allparms["-center"])
radii = []
centers = []
amp = allparms["-amp"]
amps = []
colors = []
for i in range(0,int(allparms["-chainlength"])):
    radius = lognormal.eval()
    if radius > float( allparms["-maxpyrovariance"]):
        radius = float( allparms["-maxpyrovariance"])
    print "Variance " + str(radius)
    radius = radius * meanradius
    centers.append( list(center) )
    radii.append(radius)
    amps.append(amp*lognormal.eval())
    separation = (radius + radius*allparms["-chainroughness"])/allparms["-separationfraction"]
    allparms["-view"][0] += center[0]
    allparms["-view"][1] += center[1]
    allparms["-view"][2] += center[2]
    meanradius = meanradius * allparms["-chainroughness"]
    center[0] -= separation
    amp = amp* allparms["-chainroughness"]
    color = Color( math.fabs(random.random()), math.fabs(random.random()), math.fabs(random.random()), 1 )
    colors.append(color)
allparms["-view"][0] /= int(allparms["-chainlength"])
allparms["-view"][1] /= int(allparms["-chainlength"])
allparms["-view"][2] /= int(allparms["-chainlength"])
allparms["-eye"] = list(allparms["-view"])
allparms["-eye"][2] += 25.0
allparms["-subsamples"] = 1

#allparms["-case"] = "Layout"

allparms = CmdLineFindValues( allparms )
PrintParameters( allparms )

dx = 0.1

llc = Vector( float(allparms["-dsmllc"][0]), float(allparms["-dsmllc"][1]), float(allparms["-dsmllc"][2]) )
urc = llc + Vector( float(allparms["-dsmlength"][0]), float(allparms["-dsmlength"][1]), float(allparms["-dsmlength"][2]) )
res = Vector( float(allparms["-dsmcell"][0]), float(allparms["-dsmcell"][1]), float(allparms["-dsmcell"][2]) )
gridb = makeGridBox( llc, urc, res )

setInterpOrder( gridb, 1 )


random.seed(484874)
for i in range(1,101):
    print "Frame " + str(i) + "  Building Pyro" 
    allparms["-turntableframe"] = i
    volumeObject = constant(-100000)
    volumeColor = constant( Color(0,0,0,0) )
    for c in range(0,len(centers)):
        allparms["-center"] = list( centers[c] )
        allparms["-radius"] = float(radii[c])
        allparms["-amp"] = amps[c]
        colorField = constant( Color(1,1,1,1) )
        if allparms["-usecolor"] == "yes":
            colorField = constant( colors[c] )
        if allparms["-case"] == "Layout":
            center = Vector( float(centers[c][0]), float(centers[c][1]), float(centers[c][2])    )
            radius = float( radii[c] )
            pyobj = Sphere( center, radius )
            volumeObject = Union(volumeObject, pyobj)
            volumeColor = which( volumeColor, colorField, pyobj )
        else:
            pyobj = PyroSphere(allparms)
            volumeObject = Union(volumeObject, pyobj)
            volumeColor = which( volumeColor, colorField, pyobj )
    densityField = clamp( volumeObject/constant(0.1), 0, 1.0 )
    densityGrid = makeGrid( gridb, 0.0 )
    stamp( densityGrid, densityField, 1 )
    densityField = gridded(densityGrid)
    colorGrid = makeGrid( gridb, Color(0,0,0,0) )
    stamp( colorGrid, volumeColor, 1 )
    volumeColor = gridded( colorGrid )
    if allparms["-usecolor"] == "yes":
        StandardTwoPartRender( constant(0.0), densityField, constant(Color(0,0,0,0)), volumeColor, allparms )
    else:
        StandardRender( densityField, allparms )
    allparms["-translate"][0] -= dx

#cmd = "dpaffmpeg -base " + GenerateStandardBasename(allparms) + " -usemask -useslate"
#os.system(cmd)
endJob()


