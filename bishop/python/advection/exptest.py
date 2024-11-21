#!/usr/bin/python

import os,sys

from bishoputils import *

m = Matrix( 0,0,1,   0,0,0,   -1,0,0 )
#m = Matrix( 0,0,0,   0,1,0,   0,0,2 )

#m.setExpRange(4)

em = exp(m)







print str( em.Get(0,0) ) + "  " + str( em.Get(0,1) ) + "  " + str( em.Get(0,2) )
print str( em.Get(1,0) ) + "  " + str( em.Get(1,1) ) + "  " + str( em.Get(1,2) )
print str( em.Get(2,0) ) + "  " + str( em.Get(2,1) ) + "  " + str( em.Get(2,2) )

print "det " + str( det(em) )
