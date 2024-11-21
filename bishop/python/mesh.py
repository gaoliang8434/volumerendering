#!/usr/bin/python

from bishoputils import *

beginJob()


parms = AllParameters({})
parms = CmdLineFindValues(parms)



#m = getMesh( "../819/models/cleanteapot.obj")
m = getMesh( "./standardTests/cleanbunny.obj")

lc = llc(m)
rc = urc(m)
print "LLC " + str(lc) + "\nURC " + str(rc)

origin = (lc+rc)/2.0
length = (rc-lc)*1.2
lc = origin-length/2.0
rc = origin+length/2.0
res = length.X()
if res > length.Y():
	res = length.Y()
if res > length.Z():
	res = length.Z()
totallength = res
res = res/50.0

print "Res " + str(res)


gb = makeGridBox( lc, rc, Vector(res,res,res) )

setInterpOrder( gb, 8 )

ls = mesh2ls( m, gb )

print ls


parms["-view"] = [ origin.X(), origin.Y(), origin.Z()  ]
parms["-ds"] = 0.014
parms["-eye"] = [ origin.X(), origin.Y(), rc.Z() + 3.0  ]
parms["-near"] = res 
parms["-maxpathlength"] = length.Z()*2.0
parms["-fov"] = 60.0
parms["-scatter"] = 10.0
parms["-nbthreads"] = 4
parms["-size"] = [1920,1080]
parms["-renderllc"] = [ lc.X(), lc.Y(), lc.Z()  ]
parms["-renderurc"] = [ rc.X(), rc.Y(), rc.Z()  ]
parms["-dsmllc"] = [ lc.X(), lc.Y(), lc.Z()  ]
parms["-dsmcell"] = [res,res,res]
parms["-dsmlength"] = [length.X(), length.Y(), length.Z()]

density = clamp( ls/constant(0.03), 0, 1 )


#dengrid = makeGrid( gb, 0.0 )
#stamp( dengrid, density, 1 )
#density = gridded(dengrid)

#TurntableRender( density, parms )
PrintParameters(parms)
StandardRender( density, parms )


endJob()
