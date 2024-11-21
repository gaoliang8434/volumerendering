#!/usr/bin/python

import os,sys
import NBWave as nbw
import bishoputils as bu

bu.beginJob()

parms = bu.AllParameters( nbw.nbwave_parameters() )
parms["-dt"] = 1.0/24.0
parms["-nbframes"] = 200
parms["-name"] = "plane"
parms["-case"] = "4"
parms["-version"] = "0001"
parms = bu.CmdLineFindValues( parms )
bu.PrintParameters(parms)
bu.LogIt( bu.ShowParameters(parms) )


nbthreads = int( parms["-nbthreads"])
bu.setNbCores( nbthreads )


base_surface = bu.Plane( bu.Vector(0,0,0), bu.Vector(0,-1,0)  )

dx = float(parms["-dx"])
L = 15.0
gb = bu.makeGridBox( bu.Vector(-L,-0.25*L,-L), bu.Vector(L,0.25*L,L), bu.Vector(dx,dx,dx) )
print "Gridbox:"
print str(gb)

grad_width = int(parms["-grad_width"])
dt = float(parms["-dt"])

#induce a disturbance in phi
radius = 5.0
disturbance =  bu.Sphere(bu.Vector(0,0,0), radius )
xe = bu.constant(bu.Vector(0,0,0))
se = bu.constant(0.0)
phi = bu.constant(5.0) * (bu.tanh( disturbance ) - bu.tanh( bu.constant( radius - L )))

nbstate = nbw.build_state( phi, xe, se, base_surface, gb, grad_width )
nbstate["dampen"] = float(parms["-dampen"])


nbframes = int(parms["-nbframes"])
substeps = int(parms["-substeps"])
print "Initial state: "
nbw.slice_state( nbstate )
for i in range(1,nbframes+1):
	pm = bu.ProgressMeter( 1, "EXEC TIME FRAME " + str(i) + " ---> " )
	parms["-turntableframe"] = i
	nbstate = nbw.solve( nbstate, dt, parms, substeps)
	fname = bu.formattedFrameName( parms, "vdb" )
	bu.writeOpenVDB( fname, -(nbstate["SB"] + bu.constant(float(parms["-amplify"]))*nbstate["SE"]), gb )
	print "File " + fname + " written."
	nbw.write_state_image( nbstate, parms, 0.0 )
	pm.update()



bu.endJob()
