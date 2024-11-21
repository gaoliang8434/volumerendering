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
allparms = MergeParameters( allparms, { "-advection":1, "-smoke":1 }   )

CmdLineHelp( "-h" )


allparms["-pyrogamma"] = 0.5
allparms["-amp"] = 1.5
allparms["-octaves"] = 2.0
allparms["-roughness"] = 0.7
allparms["-name"] = "firepyrov1.exr"
allparms["-scatter"] = 0.01
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-size"] = [ 960, 540 ]
allparms["-fov"] = 70.0
allparms["-ds"] = 0.01
allparms["-ambient"] = [ 0.0, 0.0, 0.0 ]
allparms["-near"] = 4.0
allparms["-maxpathlength"] = 10.0
allparms["-size"] = [ 960, 540 ]
allparms["-nodsm"] = 0



allparms = CmdLineFindValues( allparms )


doAdvection = int(allparms["-advection"])
doSmoke = int(allparms["-smoke"])

setNbCores( int(allparms["-nbthreads"] ) )

fraction = ( float( allparms["-frame"] ) - 1.0 ) / float(  allparms["-nbframes"] )
allparms["-turntableframe"] = allparms["-frame"]


ambient = allparms["-ambient"]
ambientColor = constant( Color( float(ambient[0]), float( ambient[1] ), float( ambient[2] ), 1.0  ) ) 
colr = allparms["-color"]
litColor = constant( Color( float(colr[0]), float( colr[1]  ), float( colr[2] ), 1.0  ) )


Radius = (0.2*(1.0-fraction) + fraction)*float( allparms["-radius"] )
basesdf = constant(Radius) - abs( identity() )


dispfraction = math.tanh( fraction*3.0 )
dispgamma = 2.2 *( 1.0-dispfraction) + 0.24*dispfraction
print "Fraction " + str(fraction) + "    dispgamma " + str(dispgamma)

displacementNoise = PerlinFractalSum()
noiseparms = Noise_t()
displacementNoise.getParameters( noiseparms )
noiseparms.wavelength = 3.3 
noiseparms.octaves = 3.0
noiseparms.roughness = 0.8
displacementNoise.setParameters( noiseparms )
displacementNoise.getParameters( noiseparms )
dispN = pow(abs( SFNoise( displacementNoise ) ),dispgamma)

noiseamplitude = 1.1;


maxDisplacement = noiseamplitude * math.pow( (1.0 - math.pow(noiseparms.roughness,noiseparms.octaves+1) ) / (1.0 - noiseparms.roughness ) , 0.5 )/2.0


noiseradius = 1.0*(1.0-fraction) +  fraction*( 1.0*(1.0-0.78) +  0.78*0.02)

ix = identity()
surfaceOfSphere = unitvector(ix) * constant(noiseradius)
displacement = warp( dispN, surfaceOfSphere )*constant(noiseamplitude)



velocitynoise = PerlinFractalSum()
noiseparms.translate = Vector( 0.3, -0.1, 3.24 )
noiseparms.wavelength = 1.3 
noiseparms.octaves = 3.0
noiseparms.roughness = 0.4
velocitynoise.setParameters( noiseparms )
velocity = VFNoise( velocitynoise )



sdf = basesdf + displacement

if doAdvection:
	sdf = advect( sdf, velocity, 0.1 )

rg = makeGridBox( Vector(-6,-6,-6), Vector(6,6,6), Vector(0.02,0.02,0.02) )
densitygrid = makeGrid( rg, 0.0 )
stamp( densitygrid,  clamp(sdf/constant(0.1), 0, 1.0) )
emittedDensity = gridded( densitygrid )

tempgrid = makeGrid( rg, 0.0 )
stamp( tempgrid, mask(sdf)*constant(5000.0) )
temperature = gridded(tempgrid)
emittedCd = Blackbody( temperature )*constant(60.0)

scattered = constant(0.0)
scatteredCd = constant( Color(1,1,1,1) )
if doSmoke:
	smokenoise = displacement*constant(fraction)
	if doAdvection:
		smokenoise = advect(smokenoise, velocity, 0.1)
	smokesdf = cutout(sdf + smokenoise, sdf + constant(0.25))
	smokesdf = clamp( smokesdf/constant(0.002), 0.0, 1.0 ) * constant(1000.0)
	smokegrid = makeGrid(rg, 0.0)
	stamp( smokegrid, smokesdf )
	scattered = gridded(smokegrid)



#
#========================= End of scene definition =========================
#
PrintParameters( allparms )

StandardTwoPartRender( emittedDensity, scattered, emittedCd, scatteredCd, allparms )


