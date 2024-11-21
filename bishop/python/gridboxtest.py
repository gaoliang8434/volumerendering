#!/usr/bin/python


from bishoputils import *


gb = makeGridBox( Vector(-1,-1,-1), Vector(1,1,1), Vector(0.1,0.1,0.1) )


dx = dx(gb)

print "dx = " + str(dx)
