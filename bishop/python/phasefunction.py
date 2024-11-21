#!/usr/bin/python


from bishoputils import *


hgpf = HenyeyGreensteinPhaseFunction( 0.99 )
upf = UniformPhaseFunction( 1.0/(4.0*3.14159265) )
dhgpf = DoubleHenyeyGreensteinPhaseFunction( 0.99, -0.7, 0.5 )
ffpf = FournierForandPhaseFunction( 1.05, 4.5 )

nb = 180 
dtheta = 3.14159265/nb

for i in range(0,nb+1):
	theta = dtheta * float(i)
	uval = upf.eval(theta)
	hgval = hgpf.eval(theta)
	dhgval = dhgpf.eval(theta)
	ffval = ffpf.eval(theta)
	print str(theta*180.0/3.14159265) + "  " + str(uval) + " " + str(hgval) + " " + str(dhgval) + " " + str(ffval)

