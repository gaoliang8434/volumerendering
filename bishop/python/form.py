#!/usr/bin/python


from bishoputils import *


f1 = Form( 0.0, Vector(1,0,0), Vector(0,0,0), 10.0 )
f2 = Form( 1.0, Vector(0,0,1), Vector(0,1,0), -10.0 )
f3 = wedge(f1,f2)
print f1
print f2
print f3
