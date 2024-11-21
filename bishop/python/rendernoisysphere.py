#!/usr/bin/python


import sys
import os
from bishoputils import *
from ImplicitShapes import *
import random





def RandomDirection():
	D = Vector( random.random()-0.5, random.random()-0.5, random.random()-0.5 )
	D.normalize()
	return D


def RandomWalk( P, dt, corr ):
	Pnew = P*corr + RandomDirection()*dt*(1.0-corr)
	return Pnew


def ConstrainedIFRandomWalk( P, dt, corr, surfaceProjection ):
	return evaluate( surfaceProjection, RandomWalk( P, dt, corr ) )





beginJob()

#
#=========== pyro construction ============================
#

allparms =  AllParameters(NoisySphereParameters())

CmdLineHelp( "-h" )


allparms["-amp"] = 2.0
allparms["-octaves"] = 4.0
allparms["-roughness"] = 0.5
allparms["-name"] = "noisysphere.exr"
allparms["-case"] = "TranslateAnim"
allparms["-scatter"] = 2.5
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-size"] = [ 960, 540 ]
allparms["-dsmNxNyNz"] = [ 100, 100, 100 ]
allparms["-dsmllc"] = [ -4,-4,-4 ]
allparms["-dsmurc"] = [ 4,4,4 ]
allparms["-dsmdxdydz"] = [ 0.01,0.01,0.01 ]
allparms["-fov"] = 80.0
allparms["-ds"] = 0.02
allparms["-near"] = 0.1 
allparms["-maxpathlength"] = 100.0
allparms["-renderllc"] = [-4,-4,-4]
allparms["-renderurc"] = [4,4,4]

allparms = CmdLineFindValues( allparms )


nbframes = 200
#for i in range(0,nbframes):
#	print "Frame " + str(i+1)
#	allparms["-turntableframe"] = i+1
#	allparms["-translate"] = [ 0,0, i*0.1 ]
#	PrintParameters( allparms )
#	volumeObject = NoisySphere(allparms)
#	densityField =  clamp( volumeObject/constant(0.3), 0,1 )
#	StandardRender( densityField, allparms )

#for i in range(0,nbframes):
#	print "Frame " + str(i+1+nbframes)
#	allparms["-turntableframe"] = i+1+nbframes
#	angle = i*3.14159265*2.0/100.0
#	cs = math.cos(angle)
#	ss = math.sin(angle)
#	allparms["-translate"] = [ 0, ss*1.3, (nbframes-1)*0.1 + cs*1.3 ]
#	PrintParameters( allparms )
#	volumeObject = NoisySphere(allparms)
#	densityField =  clamp( volumeObject/constant(0.3), 0,1 )
#	StandardRender( densityField, allparms )

#for i in range(0,nbframes):
#	print "Frame " + str(i+1+2*nbframes)
#	allparms["-turntableframe"] = i+1+2*nbframes
#	allparms["-translate"] = [ i*0.1, i*0.1, (i + (nbframes-1))*0.1 ]
#	PrintParameters( allparms )
#	volumeObject = NoisySphere(allparms)
#	densityField =  clamp( volumeObject/constant(0.3), 0,1 )
#	StandardRender( densityField, allparms )



allparms["-case"] = "RandomWalk"
translateP = Vector(0,0,0)
translateV = Vector( random.random()-0.5, random.random()-0.5, random.random()-0.5 )
timestep = 0.5
correlation = 0.0
bouyancy = 0.1
gravity = Vector(0,1,0)
for i in range(0,nbframes):
	print "Frame " + str(i+1)
	allparms["-turntableframe"] = i+1
	allparms["-translate"] = [ translateP.X(), translateP.Y(), translateP.Z() ]
	PrintParameters( allparms )
	volumeObject = NoisySphere(allparms)
	densityField =  clamp( volumeObject/constant(0.3), 0,1 )
	StandardRender( densityField, allparms )
	translateP = translateP + translateV*timestep*0.5
	translateV = RandomWalk( translateV, timestep, correlation )
	translateP = translateP + translateV*timestep*0.5




endJob()
