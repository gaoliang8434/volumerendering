#!/usr/bin/python


import sys
import os
from vr import *
from vrutils import *

def PutInGrid( field, gb, defvalue ):
	grid = makeGrid( gb, defvalue )
	stamp( grid, field, 1 )
	return gridded( grid )



# From 
# Back and forth error compensation and correction mehtods for semi-lagrangian schemes with application to level set interface computations
# TODD F. DUPONT AND YINGJIE LIU
def advectBFECC( field, velocity, timestep ):
	step1 = advect( field, velocity, timestep )
	step2 = advect( step1, velocity, -timestep )
	step3 = field * constant(1.5) - step2 * constant(0.5)
	step4 = advect( step3, velocity, timestep )
	return step4


# From 
# An Unconditionally Stable MacCormack Method
# Andrew Selle	Ronald Fedkiw	ByungMoon Kim Yingjie Liu	Jarek Rossignac
def advectModifiedMacCormack( field, velocity, timestep ):
	step1 = advect( field, velocity, timestep )
	step2 = advect( step1, velocity, -timestep )
	step3 = step1 + (field - step2) * constant(0.5)
	return step3

#New type
def advectModifiedVelocity( field, velocity, timestep ):
	egrad = exp( grad(velocity) * constant(-timestep) )
	modifiedvelocity = velocity * egrad
	step1 = advect( field, modifiedvelocity, timestep )
	return step1

#New type with iteration
def advectModifiedVelocityIterated( field, velocity, timestep, nbiterations ):
	egrad = exp( grad(velocity) * constant(-timestep) )
	X = identity()
	for i in range(0,nbiterations):
		X = identity() - velocity*warp( egrad, X )*constant(timestep)
	step1 = warp( field, X )
	return step1

#New type with iteration, refined
def advectModifiedVelocityIteratedRefined( field, velocity, timestep, nbiterations ):
	halftimestep = timestep/2.0
	egrad = exp( grad(velocity) * constant(-halftimestep) )
	modifiedvelocity = velocity * egrad
	X = advect( identity(), velocity, timestep )
	for i in range(0,nbiterations):
		X = advect(identity(), velocity, halftimestep) - modifiedvelocity*warp( egrad, X )*constant(halftimestep)
	step1 = warp( field, X )
	return step1

#New type with evaluating approximately the integration
def advectModifiedVelocityIntegrated( field, velocity, timestep, nbiterations, gb ):
	X = identity()
	vv = velocity * constant(timestep/2.0)
	u = velocity * constant(timestep)
	A = grad(vv)
	for i in range(0,nbiterations):
		B = warp( A, X )
		X = identity() - u*orderedSinch(A, B)
		X = PutInGrid( X - identity(), gb, Vector(0,0,0) ) + identity()
	step1 = warp( field, X )
	return step1



#New type with iteratively handling the integration
def advectIntegrated( field, velocity, timestep, nbiterations, gb ):
	U = constant(unitMatrix())
	dt = timestep/nbiterations
	egrad = exp( grad(velocity) * constant(-dt) )
	X0 = advectBFECC( identity(), velocity, dt )
	X = advect( identity(), velocity, dt )
	for i in range(1,nbiterations):
		Q = warp( egrad, X0 )
		U = PutInGrid( U * Q, gb, Matrix() )
		X = X - velocity * U * constant(dt)
		X = PutInGrid( X - identity(), gb, Vector(0,0,0) ) + identity()
		X0 = X
	step1 = warp( field, X )
	return step1


#New type with iteratively handling the integration
def advectIntegratedIterated( field, velocity, timestep, nbiterations, gb ):
	U = constant(unitMatrix())
	dt = timestep/2.0
	egrad = exp( grad(velocity) * constant(-dt) )
	X = advectModifiedMacCormack( identity(), velocity, timestep )
	#X = advect( identity(), velocity, timestep )
	for i in range(1,nbiterations):
		Q = warp( egrad, X )
		X = identity() - velocity * constant(dt) * ( constant(unitMatrix()) + egrad*Q )
		X = PutInGrid( X - identity(), gb, Vector(0,0,0) ) + identity()
	step1 = warp( field, X )
	return step1



def TransectVolume( volume, dx ):
	x = -3.0
	nb = int(6.0/dx)
	for i in range( 0,nb ):
		P = Vector( x, 0, 0 )
		value = evaluate( volume, P )
		print str(x) + " " + str(value)
		x = x + dx





allparms = AllParameters( { "-nbadvectioniterations":1 } )

CmdLineHelp( "-h" )

allparms["-name"] = "advection.exr"
allparms["-scatter"] = 2.0
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-fov"] = 70.0
allparms["-ds"] = 0.01
allparms["-near"] = 3.0
allparms["-maxpathlength"] = 14.0
allparms["-dsmcell"] = [ 0.02, 0.02, 0.02 ] 
allparms["-dsmlength"] = [ 10, 10, 10 ]
allparms["-dsmllc"] = [ -5, -5, -5 ]
allparms["-case"] = "Comparison"
allparms["-version"] = 1
allparms["-subsamples"] = 10
allparms["-size"] = [1920, 540]
allparms["-aspect"] = 1920.0/540.0

allparms = CmdLineFindValues( allparms )


noiseParms1 = Noise_t()
noiseParms1.octaves = 2.0
noiseParms1.freq = 2.0
noiseParms1.translate = Vector(0,0,0)
perlin1 = perlin( noiseParms1 )
noisefield1 = SFNoise( perlin1 )


noiseParms2 = Noise_t()
noiseParms2.octaves = 2.0
noiseParms2.freq = 1.5
noiseParms2.translate = Vector(-0.3,1.5,0.6)
perlin2 = perlin( noiseParms2 )
noisefield2 = SFNoise( perlin2 )

T = 0.15
nbIterations = int( allparms["-nbadvectioniterations"] )
dt = T/float(nbIterations)

velocityField = cross( grad(noisefield1), grad(noisefield2) )

L = 3.0 
dL = 0.05
gridbox = makeGridBox( Vector(-L,-L,-L), Vector(L,L,L), Vector( dL, dL, dL ) )

baseVolume =  clamp( Sphere( Vector(0,0,0), 1.0 )/constant(0.1), 0.0, 1.0 )
baseMap = identity()
fields = [ PutInGrid(baseVolume,gridbox,0.0), 
           PutInGrid(baseVolume,gridbox,0.0), 
	   PutInGrid(baseVolume,gridbox,0.0), 
	   PutInGrid(baseVolume,gridbox,0.0) 
	 ]

maps = []
for f in fields:
	maps.append( identity() )

nbModifiedAdvectionIterations =  3
nbIntegrations = 6 


for i in range(0,nbIterations):
	print "Iteration " + str(i)
	fields[0] = PutInGrid( advect( fields[0], velocityField, dt ), gridbox, 0.0 )
	fields[1] = PutInGrid( advectBFECC( fields[1], velocityField, dt ), gridbox, 0.0 )
	fields[2] = PutInGrid( advectModifiedMacCormack( fields[2], velocityField, dt ), gridbox, 0.0 )
	fields[3] = PutInGrid( advectModifiedVelocityIntegrated( fields[3], velocityField, dt, nbIntegrations, gridbox ), gridbox, 0.0 )
	maps[0] = PutInGrid( advect( maps[0], velocityField, dt ) - identity(), gridbox, Vector(0.0,0.0,0.0) ) + identity()
	maps[1] = PutInGrid( advectBFECC( maps[1], velocityField, dt ) - identity(), gridbox, Vector(0.0,0.0,0.0) ) + identity()
	maps[2] = PutInGrid( advectModifiedMacCormack( maps[2], velocityField, dt ) - identity(), gridbox, Vector(0.0,0.0,0.0) ) + identity()
	maps[3] = PutInGrid( advectModifiedVelocityIntegrated( maps[3], velocityField, dt, nbIntegrations, gridbox ) - identity(), gridbox, Vector(0.0,0.0,0.0) ) + identity()


#Transects of density
dxTransect = dL/5.0
for f in fields:
	print "#Field:"
	TransectVolume( f, dxTransect )
	print "\n\n\n\n"

volumePreservation = []
for m in maps:
	volumePreservation.append( det(grad(m)) - constant(1.0) )


#Transects of determinant 
for f in volumePreservation:
	print "#Determinant:"
	TransectVolume( f, dxTransect )
	print "\n\n\n\n"

#Transect of exact determinant
print "#Exact Determinant:"
TransectVolume( det( exp( grad(velocityField) * constant(-T) )) - constant(1.0), dxTransect)
print "\n\n\n\n"



exit()

# Displace and Render
hdisplacement = Vector( 3.0, 0.0, 0.0 )

density = constant(0.0)
for i in range(0, len(fields) ):
	delx = 0.0
	if len(fields)%2 == 0:
		delx = 0.5
	density = density + translate( fields[i], hdisplacement * (float( i - int(len(fields)/2.0) ) + delx)  )
allparms["-case"] = "GriddedFields" + "Steps" + str(nbIterations)
#StandardRender( density, allparms )



allparms["-ds"] = 0.0015

density = constant(0.0)
for i in range(0, len(fields) ):
	delx = 0.0
	if len(fields)%2 == 0:
		delx = 0.5
	density = density + translate( warp( baseVolume, maps[i] ), hdisplacement * (float( i - int(len(fields)/2.0) ) + delx) ) 
allparms["-case"] = "GriddedMaps" + "Steps" + str(nbIterations)
#StandardRender( density, allparms )



gridbox = makeGridBox( Vector(-L,-L,-L), Vector(L,L,L), Vector( dL/3.0, dL/3.0, dL/3.0 ) )
vpfields = []
for i in range(0,len(maps)):
	vpf =  warp( baseVolume, maps[i] ) *  ( PutInGrid( abs(det(grad(maps[i])))-constant(1.0), gridbox, 0.0 ) + constant(1.0) )
	vpfields.append(vpf)

density = constant(0.0)
for i in range(0, len(fields) ):
	delx = 0.0
	if len(fields)%2 == 0:
		delx = 0.5
	density = density + translate( vpfields[i], hdisplacement * (float( i - int(len(fields)/2.0) ) + delx)  )
allparms["-case"] = "GriddedVPFields" + "Steps" + str(nbIterations)
#StandardRender( density, allparms )



allparms["-ds"] = 0.0015



vpShell = constant(1.0)
vpColors = []
for vp in volumePreservation:
	color = constant( Color(1,0,0,1) ) * mask( vp - vpShell ) + constant( Color(0,1,0,1) ) * mask( -vp + vpShell ) + constant( Color(0,0,1,1) ) * mask( -vp - vpShell ) * mask( vp + vpShell )
	vpColors.append( color )

density = constant(0.0)
color = constant(Color(0,0,0,0) )
for i in range(0, len(fields) ):
	delx = 0.0
	if len(fields)%2 == 0:
		delx = 0.5
	vpdensity = translate( warp( baseVolume, maps[i] ), hdisplacement * (float( i - int(len(fields)/2.0) ) + delx) )
	color = color + translate( vpColors[i], hdisplacement * (float( i - int(len(fields)/2.0) ) + delx) ) * mask(vpdensity)
	density = density + vpdensity
allparms["-case"] = "GriddedVolumePreservation" + "Steps" + str(nbIterations)
allparms["-nodsm"] = 1
#StandardTwoPartRender( density, constant(0.0), color, constant(Color(0,0,0,0)), allparms )
allparms["-nodsm"] = 0


