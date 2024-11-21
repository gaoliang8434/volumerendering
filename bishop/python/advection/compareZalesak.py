#!/usr/bin/python


import sys
import os
from bishoputils import *


#
#   Compare advection maps for rigid rotation (Zalesak)
#


beginJob()

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

def PutInGrid( field, gb, defvalue ):
	grid = makeGrid( gb, defvalue )
	stamp( grid, field, 1 )
	return gridded( grid )



def GetStats( vfield, gb ):
	maxvalue = -1.0
	minvalue = 1.0e6
	avgvalue = Vector(0,0,0)
	rmsvalue = 0.0
	ncount = 0
	NX = nx( gb )
	NY = ny( gb )
	NZ = nz( gb )
	pg = ProgressMeter( NZ, "stats" )
	for k in range(0,NZ):
		for j in range(0,NY):
			for i in range(0,NX):
				P = gb.evalP( i,j,k )
				vfvalue = evaluate( vfield, P )
				avgvalue += vfvalue
				rmsvalue += vfvalue*vfvalue
				absvalue = vfvalue.magnitude()
				if maxvalue < absvalue:
					maxvalue = absvalue
				if minvalue > absvalue:
					minvalue = absvalue
				ncount = ncount + 1
		pg.update()
	avgvalue = avgvalue / float(ncount)
	rmsvalue = rmsvalue / float(ncount)
	rmsvalue = rmsvalue - avgvalue*avgvalue
	rmsvalue = math.sqrt(rmsvalue)
	return [ avgvalue, rmsvalue, minvalue, maxvalue ]





allparms = AllParameters( { "-nbadvectioniterations":0 } )

CmdLineHelp( "-h" )

allparms["-name"] = "gs.exr"
allparms["-scatter"] = 5.0
allparms["-ambient"] = [ 0, 0, 0 ]
allparms["-fov"] = 50.0
allparms["-eye"] = [0, 0, 200]
allparms["-view"] = [0, 0, 0]
allparms["-up"] = [0, 1, 0]
allparms["-ds"] = 0.05
allparms["-near"] = 100.0
allparms["-maxpathlength"] = 200.0
allparms["-dsmcell"] = [ 0.05, 0.05, 0.05 ] 
allparms["-dsmlength"] = [ 100, 100, 100 ]
allparms["-dsmllc"] = [ -50,-50,-50 ]
allparms["-case"] = "ZalesakSphere"
allparms["-version"] = 1
allparms["-subsamples"] = 1
allparms["-size"] = [1920, 1080]
allparms["-aspect"] = 1920.0/1080.0
allparms["-nbthreads"] = 4

allparms = CmdLineFindValues( allparms )


nrange = 0
dilationScale = math.pow(2, nrange)

velocityField = component( constant(3.14159265/314.0)*( constant(00.0)-yIdentity() ) , constant(3.14159265/314.0)*( xIdentity() - constant(00.0) ) , constant(0.0) )



baseT = 6.28/10.0

L = 100 
dL = 10.0
gb = makeGridBox( Vector(-L,-L,-L), Vector(L,L,L), Vector(dL,dL,dL) )







results = []
for i in range(1,200):
	T = baseT*math.pow(10, i/50.0)
	gsmap = gradientStretchCM( velocityField, T, 256 )
	bfeccmap = advectBFECC( identity(), velocityField, T )
	mmmap = advectModifiedMacCormack( identity(), velocityField, T )
	slmap = advect( identity(), velocityField, T )
	bfeccstats = GetStats( bfeccmap-gsmap, gb ) 
	slstats = GetStats( slmap-gsmap, gb ) 
	mmstats = GetStats( mmmap-gsmap, gb ) 
	results.append( str(T) + "  " + str( slstats[1] ) + "  " + str(bfeccstats[1]) + "  " + str(mmstats[1])  )

for l in results:
	print l

endJob()
