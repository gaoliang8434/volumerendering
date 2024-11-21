#!/usr/bin/python

import os,sys
import NBWave as nbw
import bishoputils as bu

bu.beginJob()

#Set up input parameters
parms = bu.AllParameters( bu.MergeParameters( nbw.nbwave_parameters(), { "-source":0.1, "-sourceradius":0.2 } ) )
parms["-dt"] = 1.0/24.0
parms["-nbframes"] = 75 
parms["-name"] = "plane.vdb"
parms["-case"] = "8"
parms["-version"] = "0001"
parms = bu.CmdLineFindValues( parms )
bu.PrintParameters(parms)
bu.LogIt( bu.ShowParameters(parms) )


nbthreads = int( parms["-nbthreads"])
bu.setNbCores( nbthreads )




# Implicit function defining the surface S_B
base_surface =  bu.Plane( bu.Vector(0,0,0), bu.Vector(0,-1,0)  )


#Set the scale of the simulation domain and cell size
L = 1.5
dx = float(parms["-dx"])
gb = bu.makeGridBox( bu.Vector(-L,-0.5*L,-L), bu.Vector(L,0.5*L,L), bu.Vector(dx,dx,dx) )
print "Gridbox:"
print str(gb)

#Order of spatial gradient and time step
grad_width = int(parms["-grad_width"])
dt = float(parms["-dt"])


#induce a disturbance
radius = float(parms["-sourceradius"] )
source = float( parms["-source"] )
phi = bu.constant(0.0)
xe =  bu.constant(bu.Vector(0,-2.0*source,0)) * bu.exp( -(bu.identity() * bu.identity())/bu.constant(radius*radius)  )
se = xe *  bu.fdgrad( base_surface, grad_width, dx, dx, dx  )

# Gravity is along the normal of the levelset
#Initial the simultion state data
nbstate = nbw.build_state( phi, xe, se, base_surface, gb, grad_width, bu.constant(bu.Vector(0,-9.8,0)) )

erode = float( parms["-erode"] )

nbframes = int(parms["-nbframes"])
substeps = int(parms["-substeps"])
boundary_clamp = nbw.boundary_shell( nbstate["gridbox"], 0.25 )
nbstate["phi"] *= boundary_clamp
nbstate["SE"]  *= boundary_clamp
nbstate["XE"]  *= boundary_clamp

#Run simulation frames
for i in range(1,nbframes+1):
	pm = bu.ProgressMeter( 1, "EXEC TIME FRAME " + str(i) + " ---> " )
	parms["-turntableframe"] = i

	#Advance dynamics
	nbstate = nbw.solve( nbstate, dt, parms, substeps)

        #Clamp boundaries
	nbstate["phi"] *= boundary_clamp
	nbstate["SE"]  *= boundary_clamp
	nbstate["XE"]  *= boundary_clamp

	#Write openvdb file of the surface data
	fname = bu.GenerateStandardFilename( parms )
	surface = nbstate["SB"] + nbstate["SE"]
	if erode != 0.0:
		#Apply erosion to the output data.  This done not affect the simulation
		surface = nbstate["SB"] + nbw.cheesy_levelset_filter( nbstate["SE"], gb, erode )
	bu.writeOpenVDB( fname, -surface, gb )
	print "File " + fname + " written."
	pm.update()



bu.endJob()
