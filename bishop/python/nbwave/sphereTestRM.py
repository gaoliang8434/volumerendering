#!/usr/bin/python

import NBWave as nbw
import bishoputils as bu

bu.beginJob()

def padFrame( f ):
	fstr = str(f)
	if f < 1000:
		fstr = "0" + fstr
	if f < 100:
		fstr = "0" + fstr
	if f < 10:
		fstr = "0" + fstr
	return fstr

parms = bu.AllParameters( nbw.nbwave_parameters() )
parms["-dt"] = 1.0/24.0
parms["-nbframes"] = 30
parms["-name"] = "sphere"
parms["-case"] = "2"
parms["-version"] = "0001"
parms = bu.CmdLineFindValues( parms )

base_surface = bu.Sphere( bu.Vector(0,0,0), 1.0 )

dx = float(parms["-dx"])
L = 1.5
gb = bu.makeGridBox( bu.Vector(-L,-L,-L), bu.Vector(L,L,L), bu.Vector(dx,dx,dx) )
print "GridBox:"
print str(gb)

grad_width = int(parms["-grad_width"])
nbstate = nbw.build_initial_state( base_surface, gb, grad_width )
nbstate["dampen"] = float(parms["-dampen"])

dt = float(parms["-dt"])

#induce a disturbance in phi
#nbstate["phi"] += bu.exp( -( (bu.identity() - bu.constant(bu.Vector(0,1,0)) )*( bu.identity() - bu.constant(bu.Vector(0,1,0) ) ) )/bu.constant( 0.01 )  ) * bu.clamp( bu.Shell( base_surface, 0.2 ), 0, 1.0 )   * bu.constant(3.0)
nbstate["XE"] += bu.exp( -( (bu.identity()-bu.constant(bu.Vector(0,1,0)))*(bu.identity()-bu.constant(bu.Vector(0,1,0))) )/bu.constant( 0.01 )  ) * bu.constant(bu.Vector(0,1,0))

bu.PrintParameters(parms)
bu.LogIt( bu.ShowParameters(parms) )

nbframes = int(parms["-nbframes"])
outputCount = int(parms["-output"])
count = 0
fname_base = parms["-name"] + parms["-case"] + "-" + parms["-version"] + "."
for i in range(1,nbframes+1):
	pm = bu.ProgressMeter( 1, "EXEC TIME FRAME " + str(i) + " ---> " )
	nbstate = nbw.solve( nbstate, dt, parms )
	count = count + 1
	if count == outputCount:
		fname = fname_base + padFrame(i) + ".vdb"
		bu.writeOpenVDB( fname, -(nbstate["SB"] + nbstate["SE"]), gb )
		print "File " + fname + " written."
		count = 0
	pm.update()



bu.endJob()
