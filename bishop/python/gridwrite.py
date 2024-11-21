#!/usr/bin/python


from bishoputils import *

dom = RectangularGrid()
dom.init( Vector(-1,-1,-1), Vector(1,1,1), Vector(0.01, 0.01, 0.01) )
print "Dimensions " + str(dom.nx()) + " " + str(dom.ny()) + " " + str(dom.nz())
sgrid = makeScalarGrid( dom, float(0.0) )
print "Writing Scalar Grid ======================="
writeScalarGrid( sgrid, "sgrid.grid" )
print "Writing Scalar Grid ======================= DONE"
print "Dimensions " + str(dom.nx()) + " " + str(dom.ny()) + " " + str(dom.nz())
vgrid = makeVectorGrid( dom, Vector(0,0,0) )
print "Writing Vector Grid ======================="
writeVectorGrid( vgrid, "vgrid.grid" )
print "Writing Vector Grid ======================= DONE"
print "Dimensions " + str(dom.nx()) + " " + str(dom.ny()) + " " + str(dom.nz())
cgrid = makeColorGrid( dom, Color(0,0,0,1) )
print "Writing Color Grid ======================="
writeColorGrid( cgrid, "cgrid.grid" )
print "Writing Color Grid ======================= DONE"


