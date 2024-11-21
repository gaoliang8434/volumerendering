#!/usr/bin/python

import os
import sys
from bishoputils import *


beginJob()


pc = makePointCloud(10)

attrlist = show_all_attrs( pc )
print "Generic Point Cloud"
print str(attrlist)

create_attr( pc, "freq", 1.0 )
create_attr( pc, "newcolor", Color(0,1,0,0) )

attrlist = show_all_attrs( pc )
print "Generic Point Cloud with custom attrs"
print str(attrlist)

wc = makeWispCloud(2)
attrlist = show_all_attrs( wc )
print "Wisp Point Cloud"
print str(attrlist)




endJob()
