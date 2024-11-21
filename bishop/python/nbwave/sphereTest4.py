#!/usr/bin/python

import os,sys
import NBWave as nbw
import bishoputils as bu

bu.beginJob()

#Set up input parameters
parms = bu.AllParameters( bu.MergeParameters( nbw.nbwave_parameters(), { "-source":0.1, "-sourceradius":0.2 } ) )
parms["-dt"] = 1.0/24.0
parms["-nbframes"] = 100
parms["-name"] = "sphere.vdb"
parms["-case"] = "4"
parms["-version"] = "0001"
parms = bu.CmdLineFindValues( parms )
bu.PrintParameters(parms)
bu.LogIt( bu.ShowParameters(parms) )

nbthreads = int( parms["-nbthreads"])
bu.setNbCores( nbthreads )
 
# Implicit function defining the surface S_B
base_surface = bu.Sphere( bu.Vector(0,0.00132344,0), 1.0  )

#Set the scale of the simulation domain and cell size
dx = float(parms["-dx"])
L = 1.5
gb = bu.makeGridBox( bu.Vector(-L,-L,-L), bu.Vector(L,L,L), bu.Vector(dx,dx,dx) )
print "Gridbox:"
print str(gb)

#Order of spatial gradient and time step
grad_width = int(parms["-grad_width"])
print "Grad width " + str(grad_width)
dt = float(parms["-dt"])

#induce a disturbance
radius = float(parms["-sourceradius"] )
source = float( parms["-source"] )
phi = bu.constant(0.0)
xe =  bu.constant(bu.Vector(0,-source,0)) * bu.exp( -((bu.identity()-bu.constant(bu.Vector(0,1,0))) * (bu.identity()-bu.constant(bu.Vector(0,1,0))))/bu.constant(radius*radius)  )
se = xe *  bu.fdboundedgrad( base_surface, grad_width, gb  )



# Gravity is along the normal of the levelset
gravity = bu.constant(9.8)*bu.unitvector(bu.fdboundedgrad( base_surface, grad_width, gb ))

#Initial the simultion state data
nbstate = nbw.build_state( phi, xe, se, base_surface, gb, grad_width, gravity )

# Set up boundary clamping
boundary_clamp = nbw.boundary_shell( nbstate["gridbox"], 0.25 )
nbstate["phi"] *= boundary_clamp
nbstate["SE"]  *= boundary_clamp
nbstate["XE"]  *= boundary_clamp

erode = float( parms["-erode"] )

nbframes = int(parms["-nbframes"])
substeps = int(parms["-substeps"])

#Run simulation frames
for i in range(1,nbframes+1):
	pm = bu.ProgressMeter( 1, "EXEC TIME FRAME " + str(i) + " ---> " )
	parms["-turntableframe"] = i

	#Advance dynamics
	nbstate = nbw.solve( nbstate, dt, parms, substeps )

        #Clamp boundaries
	nbstate["phi"] *= boundary_clamp
	nbstate["SE"]  *= boundary_clamp
	nbstate["XE"]  *= boundary_clamp

	#Write openvdb file of the surface data
	fname = bu.GenerateStandardFilename( parms )
	surface = nbstate["SB"] + nbstate["SE"]
	if erode != 0.0:
		#Apply erosion to the output data.  This done not affect the simulation
		print "\n\tEROSION\n"
		surface = nbstate["SB"] + nbw.cheesy_levelset_filter( nbstate["SE"], gb, erode )
	bu.writeOpenVDB( fname, -surface, gb )
	print "File " + fname + " written."
	pm.update()



bu.endJob()
