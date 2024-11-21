#!/usr/bin/python


import sys
import os
from bishop import *



def CmdLineFindIndex( tag ):
	for i in range(len(sys.argv)):
		if sys.argv[i] == tag:
			return i
	return -1

def CmdLineFind( tag, defaultvalue ):
	i = CmdLineFindIndex(tag)
	if i > 0:
		if i < len(sys.argv)-1:
			return sys.argv[i+1]
	return defaultvalue





majorrad = float(CmdLineFind( "-majorradius", 3.0 ))
minorrad = float(CmdLineFind( "-minorradius", 1.0 ))
center = Vector(0,0,0)


axis0 = Vector( 1,1,1 )
axis1 = Vector( 1,-1,0 )
axis2 = axis0^axis1
axis0.normalize()
axis1.normalize()
axis2.normalize()

v1 = EllipseVolume( center, axis0, majorrad*1.3, minorrad )
v2 = EllipseVolume( center, axis1, majorrad, minorrad )
u1 = UnionVolume( v1, v2 )
v4 = EllipseVolume( center, axis2, majorrad, minorrad )
u2 = UnionVolume( u1, v4 )

s1 =  SphereVolume( center+axis1*majorrad, majorrad/3.0 )
u3 = UnionVolume( u2, s1 )
s2 =  SphereVolume( center-axis1*majorrad, majorrad/3.0 )
u4 = UnionVolume( u3, s2 )

s3 =  SphereVolume( center+axis2*majorrad, majorrad/3.0 )
u4 = UnionVolume( u3, s3 )
s4 =  SphereVolume( center-axis2*majorrad, majorrad/3.0 )
u5 = UnionVolume( u4, s4 )





y = 0
z = 0
x = -2.0*majorrad
dx = majorrad/50.0
while x < 2.0*majorrad:
	P = Vector( x, y, z )
	vValue = u5.eval(P)
	print str(x) + " " + str(vValue)
	x = x + dx
