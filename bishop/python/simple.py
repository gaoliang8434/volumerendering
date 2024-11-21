#!/usr/bin/python

from bishop import *


trans = Vector( 1.0, 0.334485, -32.84848 )
freq = 1.432234
amp = 1.5
rough = 0.5
octaves =  3.0
pyrogamma =  0.33333
pyrotime = 0.0
timescale = 1.0
pyrotime *= timescale
radius = 2.5
center = Vector(0,0,0)

v1 = PyroclasticVolume( center, radius, amp, octaves, freq, rough, trans, pyrotime, pyrogamma ) 
v2 = ClampVolume( v1, 0.0, 1.0 );

y = 0
z = 0
x = -2.0*radius
dx = radius/50.0
while x < 2.0*radius:
	P = Vector( x, y, z )
	vValue = v2.eval(P)
	print str(x) + " " + str(vValue)
	x = x + dx
