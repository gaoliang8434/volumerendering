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


print "Starting" 


majorrad = float(CmdLineFind( "-majorradius", 3.0 ))
minorrad = float(CmdLineFind( "-minorradius", 1.0 ))
center = Vector(0,0,0)

print "Milestone 1"

axis0 = Vector( 1,1,1 )
axis1 = Vector( 1,-1,0 )
axis2 = axis0^axis1
axis0.normalize()
axis1.normalize()
axis2.normalize()

print "Milestone 2"
v1 =  Ellipse( center, axis0, majorrad*1.3, minorrad )
v2 =  Ellipse( center, axis1, majorrad, minorrad )
u1 =   Union(v1, v2 )
v4 =  Ellipse( center, axis2, majorrad, minorrad ) 
u2 =   Union(u1, v4 )

print "Milestone 3"
s1 =  Sphere( center+axis1*majorrad, majorrad/3.0 ) 
u3 =  Union(u2, s1 )
s2 =  Sphere( center-axis1*majorrad, majorrad/3.0 ) 
u4 =  Union(u3, s2 )

print "Milestone 4"
s3 =  Sphere( center+axis2*majorrad, majorrad/3.0 )
print "Milestone 4.1"
u4 = Union(u3, s3 )
print "Milestone 4.2"
s4 =  Sphere( center-axis2*majorrad, majorrad/3.0 ) 
print "Milestone 4.3"
u5 =  Union(u4, s4 )









print "Milestone 5"
y = 0
z = 0
x = -2.0*majorrad
dx = majorrad/50.0
while x < 2.0*majorrad:
	P = Vector( x, y, z )
	vValue = evaluate(u5,P)
	print str(x) + " " + str(vValue)
	x = x + dx
