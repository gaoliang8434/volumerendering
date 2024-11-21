#!/usr/bin/python

import os,sys
import NBWave as nbw
import bishoputils as bu

bu.beginJob()

parms = bu.AllParameters( nbw.nbwave_parameters() )
parms["-dt"] = 1.0/24.0
parms["-nbframes"] = 200
parms["-name"] = "plane"
parms["-case"] = "5"
parms["-version"] = "0001"
parms = bu.CmdLineFindValues( parms )
bu.PrintParameters(parms)
bu.LogIt( bu.ShowParameters(parms) )


nbthreads = int( parms["-nbthreads"])
bu.setNbCores( nbthreads )



L = 15.0

base_surface =  bu.Plane( bu.Vector(0,0,0), bu.Vector(0,-1,0)  )




dx = float(parms["-dx"])
gb = bu.makeGridBox( bu.Vector(-L,-0.25*L,-L), bu.Vector(L,0.25*L,L), bu.Vector(dx,dx,dx) )
print "Gridbox:"
print str(gb)

grad_width = int(parms["-grad_width"])
dt = float(parms["-dt"])


#induce a disturbance in phi
radius = 2.0
phi = bu.constant(0.0)
se =  bu.constant(3.0) * bu.exp( -(bu.identity() * bu.identity())/bu.constant(radius*radius)  )
xe = bu.fdboundedgrad( se, grad_width, gb  )
se = xe *  bu.fdgrad( base_surface, grad_width, dx, dx, dx  )

forceXe = bu.constant(0.01)*xe;
forceSe = bu.constant(0.01)*se;





nbstate = nbw.build_state( phi, xe, se, base_surface, gb, grad_width )
nbstate["dampen"] = float(parms["-dampen"])


nbframes = int(parms["-nbframes"])
substeps = int(parms["-substeps"])
for i in range(1,nbframes+1):
        neold = bu.fdboundedgrad( nbstate["SE"] + nbstate["SB"], grad_width, gb )
	pm = bu.ProgressMeter( 1, "EXEC TIME FRAME " + str(i) + " ---> " )
	parms["-turntableframe"] = i
	nbstate = nbw.solve( nbstate, dt, parms, substeps)
	fname = bu.formattedFrameName( parms, "vdb" )
	bu.writeOpenVDB( fname, -(nbstate["SB"] + bu.constant(float(parms["-amplify"]))*nbstate["SE"]), gb )
	print "File " + fname + " written."
        #nenew = bu.fdboundedgrad( nbstate["SE"] + nbstate["SB"], grad_width, gb )
        #sediff = neold * nenew
	#fname = "SeDiff" + bu.formattedFrameName( parms, "vdb" )
	#bu.writeOpenVDB( fname, sediff, gb )
	#print "File " + fname + " written."

        #forceXe = bu.translate(forceXe, bu.Vector(0,0,0.4) )
        #forceSe = bu.translate(forceSe, bu.Vector(0,0,0.4) )
        #nbstate["XE"] += forceXe;
        #nbstate["SE"] += forceSe;


	pm.update()



bu.endJob()
