#!/usr/bin/python

from vrutils import *



# Blend:
#   Ramps on beginning at frame 6
#   Ramps off by frame 45
#   Spherical region ramps up to radius of 4
def GetBlendAlpha( f ):
	blendTimeAlpha = (f-6.0)/15.0
	if blendTimeAlpha > 1.0:
		blendTimeAlpha = 1.0
	blendTimeAlpha = blendTimeAlpha * (45.0-f)/15.0
	if blendTimeAlpha > 1.0:
		blendTimeAlpha = 1.0
	if blendTimeAlpha<0:
		blendTimeAlpha = 0.0
	blendIF = constant(0.0)
	if blendTimeAlpha > 0:
		blendIF = constant(9.0*blendTimeAlpha*blendTimeAlpha) - (identity()-constant(Vector(-0.25,3,0)))*(identity()-constant(Vector(-0.25,3,0)))
	blendSpaceAlpha = clamp( blendIF/constant(3.0),0.0, 1.0)
	blendAlpha = constant(blendTimeAlpha) * blendSpaceAlpha
	return blendAlpha


def GridCM( cm, gridbox ):
	dX = cm - identity()
	cmgrid = makeGrid( gridbox, Vector(0,0,0) )
	stamp( cmgrid, dX )
	cmX = gridded(cmgrid) + identity()
	return cmX

def GridSF( sf, gridbox ):
	sfgrid = makeGrid( gridbox, 0.0 )
	stamp( sfgrid, sf )
	return gridded(sfgrid)

def fftdivfree( vf, gridbox ):
	velgrid = makeGrid( gridbox, Vector(0,0,0) )
	FFTDivFree( velgrid, vf )
	return gridded(velgrid)


#
# Parameters that can be controlled from the command line
#
parms = AllParameters({ "-start":1, "-end":70 })
parms["-fov"] = 40.0
parms["-size"] = [480, 1080]
parms["-aspect"] = float( parms["-size"][0])/float( parms["-size"][1] )
parms["-nbthreads"] = int(20)
parms["-name"] = "blendimage.exr"
parms["-view"] = [0,4,0]
parms["-eye"][1] = float(parms["-eye"][1]) + float(parms["-view"][1])
parms["-ds"] = 0.001
parms["-dsmlength"] = [10, 15, 10]
parms["-dsmcell"] = [0.03, 0.03, 0.03]
parms["-scatter"] = 3.0
parms["-nbthreads"] = 20
parms = CmdLineFindValues(parms)
CmdLineHelp( "-h" )
PrintParameters( parms )
imageBaseName = parms["-name"]

#
# Construct base spherical implicit function and density
#
sphere = constant(1.0) - identity()*identity()
density0 = clamp( sphere/constant(0.1), 0.0, 1.0 )

#
# Grid definitions for CFD and CM selma
#
gb = makeGridBox( Vector(-3,-2,-3), Vector(3,10,3),Vector(0.015,0.015,0.015) )
cmgb = makeGridBox( Vector(-3,-2,-3), Vector(3,10,3),Vector(0.015,0.015,0.015) )

#
# Simulation parameters
#
dt = 0.1
gravity = constant(Vector(0,-10*dt,0))

#
#initialize simple flow
#
velocity  = constant(Vector(0,0,0))
cmX = identity()


#
# initialize noisy injection flow
#
nnvelocity = constant(0.0)
nncmX = identity()


#
# initialize velocity noise
#
noiseparm = Noise_t()
noiseparm.octaves = 3.0
vstrength = constant(0.1)
vnoise = SFNoise( perlin(noiseparm) )
velocityreference = grad( vnoise )
vnoise = vnoise*vstrength

#
# Loop over frames, run simualtions, blend, render
#
for f in range(1,int(parms["-end"])+1):
	print "Frame " + str(f)

	# CFD and CM for simple flow
	density = warp( density0, cmX )
	velocity += -density*gravity
	velocity = advect(velocity, velocity, dt )
	velocity = fftdivfree( velocity, gb )
	cmX = advect( cmX, velocity, dt )
	cmX = GridCM( cmX, cmgb ) 

	# CFD and CM for noisy injection flow
	noisevelocity = cross( grad(vnoise), velocityreference )
	nncmX = advect( nncmX, velocity+noisevelocity, dt )
	nncmX = GridCM( nncmX, cmgb )
	vnoise = GridSF(advect( vnoise, velocity + noisevelocity, dt ), cmgb )


	#Compute blend coefficient
	blendAlpha = GetBlendAlpha(f)

	# Blend operation
        #blendCMX= cmX*(constant(1.0)-blendAlpha) + nncmX*blendAlpha
        blendCMX= nncmX
        #blendCMX= cmX
	blenddensity = warp( density0, blendCMX )

	# Render
	if f>=int(parms["-start"]):
		parms["-name"] = imageBaseName
		parms["-turntableframe"] = int(f)
		StandardRender( blenddensity, parms )



	

