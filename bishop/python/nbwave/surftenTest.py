#!/usr/bin/python

import os,sys
import NBWave as nbw
import bishoputils as bu

bu.beginJob()

parms = bu.AllParameters( bu.MergeParameters( nbw.nbwave_parameters(), { "-source":1.0, "-sourceradius":2.0 } ) )
parms["-dt"] = 1.0/24.0
parms["-nbframes"] = 75 
parms["-name"] = "surfaceTension.vdb"
parms["-case"] = "1"
parms["-version"] = "0001"
parms = bu.CmdLineFindValues( parms )
bu.PrintParameters(parms)
bu.LogIt( bu.ShowParameters(parms) )


nbthreads = int( parms["-nbthreads"])
bu.setNbCores( nbthreads )



L = 1.5

base_surface =  bu.Plane( bu.Vector(0,0,0), bu.Vector(0,-1,0)  )




dx = float(parms["-dx"])
gb = bu.makeGridBox( bu.Vector(-L,-0.5*L,-L), bu.Vector(L,0.5*L,L), bu.Vector(dx,dx,dx) )
print "Gridbox:"
print str(gb)

grad_width = int(parms["-grad_width"])
dt = float(parms["-dt"])


#induce a disturbance in phi
radius = float(parms["-sourceradius"] )
source = float( parms["-source"] )
phi = bu.constant(0.0)
xe =  bu.constant(bu.Vector(0,-2.0*source,0)) * bu.exp( -(bu.identity() * bu.identity())/bu.constant(radius*radius)  )
se = xe *  bu.fdgrad( base_surface, grad_width, dx, dx, dx  )

forceXe = bu.constant(0.01)*xe;
forceSe = bu.constant(0.01)*se;





nbstate = nbw.build_state( phi, xe, se, base_surface, gb, grad_width, bu.constant(bu.Vector(0,-9.8,0)) )
nbstate["dampen"] = float(parms["-dampen"])


nbframes = int(parms["-nbframes"])
substeps = int(parms["-substeps"])
stsubsteps = 1
boundary_clamp = nbw.boundary_shell( nbstate["gridbox"], 0.25 )
nbstate["phi"] *= boundary_clamp
nbstate["SE"]  *= boundary_clamp
nbstate["XE"]  *= boundary_clamp

erode = float( parms["-erode"] )

for i in range(1,nbframes+1):
        neold = bu.fdboundedgrad( nbstate["SE"] + nbstate["SB"], grad_width, gb )
	pm = bu.ProgressMeter( 1, "EXEC TIME FRAME " + str(i) + " ---> " )
	parms["-turntableframe"] = i
	nbstate = nbw.solve( nbstate, dt, parms, substeps)

	nbstate["phi"] *= boundary_clamp
	nbstate["SE"]  *= boundary_clamp
	nbstate["XE"]  *= boundary_clamp


	nbstate["phi"] *= boundary_clamp
	nbstate["SE"]  *= boundary_clamp
	nbstate["XE"]  *= boundary_clamp


	fname = bu.GenerateStandardFilename( parms )
	surface = nbstate["SB"] + nbstate["SE"]
	if erode != 0.0:
		surface = nbstate["SB"] + nbw.cheesy_levelset_filter( nbstate["SE"], gb, erode )
	bu.writeOpenVDB( fname, -surface, gb )
	print "File " + fname + " written."



	pm.update()



bu.endJob()
