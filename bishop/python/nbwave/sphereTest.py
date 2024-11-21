#!/usr/bin/python

import NBWave as nbw
import bishoputils as bu


def padFrame( f ):
	fstr = str(f)
	if f < 1000:
		fstr = "0" + fstr
	if f < 100:
		fstr = "0" + fstr
	if f < 10:
		fstr = "0" + fstr
	return fstr

nbw_parms = { "-grad_width":1, "-dx":0.1  }
parms = bu.AllParameters( nbw_parms )
parms["-dt"] = 0.04
parms["-nbframes"] = 30
parms = bu.CmdLineFindValues( parms )

sphere = bu.Sphere( bu.Vector(0,0,0), 1.0 )

dx = float(parms["-dx"])
gb = bu.makeGridBox( bu.Vector(-2,-2,-2), bu.Vector(2,2,2), bu.Vector(dx,dx,dx) )

grad_width = int(parms["-grad_width"])

nbstate = nbw.build_initial_state( sphere, gb, grad_width )

dt = float(parms["-dt"])

#induce a disturbance in phi

nbstate["phi"] += bu.clamp(  bu.intersection(bu.Sphere( bu.Vector(0,1,0), 0.2 ),  bu.Shell( nbstate["SB"], 0.1 )  ) , 0, 1 ) * bu.constant(0.1) 

bu.PrintParameters(parms)

nbframes = int(parms["-nbframes"])

outputCount = 1
count = 0
for i in range(1,nbframes+1):
	pm = bu.ProgressMeter( 1, "EXEC TIME FRAME " + str(i) + " ---> " )
	print "Frame " + str(i)
	nbstate = nbw.fourth_order_solver( nbstate, dt )
	#nbstate = nbw.blanes_moan_solver( nbstate, dt )
	#nbstate = nbw.leapfrog_solver( nbstate, dt )
	#nbstate = nbw.gridify_state( nbstate )
	count = count + 1
	if count == outputCount:
		fname = "surface-0006." + padFrame(i) + ".vdb"
		bu.writeOpenVDB( fname, -(nbstate["SB"] + nbstate["SE"]), gb )
		count = 0
	pm.update()




