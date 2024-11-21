#!/usr/bin/python


import sys
import os
from bishop import *
from bishopUtils import *


#
#=========== Jack construction ============================
#

def ImplicitJackParameters():
	parms = {  "-majorradius":3.0, "-minorradius":1 }
	return parms



def ImplicitJackSDF( parameters ):
	majorrad = float(parameters[ "-majorradius"])
	minorrad = float(parameters[ "-minorradius" ])
	center = Vector(0,0,0)
	axis0 = Vector( 1,1,1 )
	axis1 = Vector( 1,-1,0 )
	axis2 = axis0^axis1
	axis0.normalize()
	axis1.normalize()
	axis2.normalize()
	v1 = Ellipse( center, axis0, majorrad*1.3, minorrad )
	v2 = Ellipse( center, axis1, majorrad, minorrad )
	v3 = Ellipse( center, axis2, majorrad, minorrad )
	u1 = Union( Union( v1, v2 ), v3 )
	s1 = Sphere( center+axis1*majorrad, majorrad/3.0 )
	s2 = Sphere( center-axis1*majorrad, majorrad/3.0 )
	u1 = Union( Union( u1, s1 ), s2 )
	s3 = Sphere( center+axis2*majorrad, majorrad/3.0 )
	s4 = Sphere( center-axis2*majorrad, majorrad/3.0 )
	u1 = Union( Union( u1, s3 ), s4 )


	return u1



def ImplicitJack( parameters ):
	majorrad = float(parameters[ "-majorradius"])
	minorrad = float(parameters[ "-minorradius" ])
	center = Vector(0,0,0)
	axis0 = Vector( 1,1,1 )
	axis1 = Vector( 1,-1,0 )
	axis2 = axis0^axis1
	axis0.normalize()
	axis1.normalize()
	axis2.normalize()
	v1 = Ellipse( center, axis0, majorrad*1.3, minorrad )
	v2 = Ellipse( center, axis1, majorrad, minorrad )
	v3 = Ellipse( center, axis2, majorrad, minorrad )
	u1 = Union( Union( v1, v2 ), v3 )
	s1 = Sphere( center+axis1*majorrad, majorrad/3.0 )
	s2 = Sphere( center-axis1*majorrad, majorrad/3.0 )
	u1 = Union( Union( u1, s1 ), s2 )
	s3 = Sphere( center+axis2*majorrad, majorrad/3.0 )
	s4 = Sphere( center-axis2*majorrad, majorrad/3.0 )
	u1 = Union( Union( u1, s3 ), s4 )

	threshold = CmdLineFind("-thresh", 0.1 )
	u1 = multiply( u1, 1.0/threshold )
	jack = clamp( u1, 0.0, 1.0 )
	return jack




#
#=========== Pyroclastic sphere construction ============================
#

def PyroSphereParameters():
	parms = {  "-translate":[0,0,0], "-freq":1, "-roughness":0.5, "-octaves":2.0, "-pyrogamma":1.0, "-timescale":1.0, "-pyrotime":0.0, "-amp":1.0, "-radius":2.5, "-center":[0,0,0]  }
	return parms


def PyroSphere( parameters ):
	trans = Vector( 0,0,0 )
	tvalue = parameters["-translate"]
	if len(tvalue) >= 3:
		trans = Vector( tvalue[0], tvalue[1], tvalue[2] )
	center= Vector( 0,0,0 )
	cvalue = parameters["-center"]
	if len(cvalue) >= 3:
		center = Vector( cvalue[0], cvalue[1], cvalue[2] )
	freq = float(parameters["-freq"])
	amp = float(parameters["-amp"])
	rough = float(parameters["-roughness"] )
	octaves = float(parameters["-octaves"] )
	pyrogamma = float(parameters["-pyrogamma"] )
	pyrotime = float( parameters["-pyrotime" ] )
	timescale = float(parameters["-timescale" ] )
	pyrotime *= timescale
	radius = float(parameters["-radius"] )
	v1 = Pyroclast( center, radius, amp, octaves, freq, rough, trans, pyrotime, pyrogamma )
	return v1



def RadialPyroSphere( parameters ):
	trans = Vector( 0,0,0 )
	tvalue = parameters["-translate"]
	if len(tvalue) >= 3:
		trans = Vector( tvalue[0], tvalue[1], tvalue[2] )
	center= Vector( 0,0,0 )
	cvalue = parameters["-center"]
	if len(cvalue) >= 3:
		center = Vector( cvalue[0], cvalue[1], cvalue[2] )
	freq = float(parameters["-freq"])
	amp = float(parameters["-amp"])
	rough = float(parameters["-roughness"] )
	octaves = float(parameters["-octaves"] )
	pyrogamma = float(parameters["-pyrogamma"] )
	pyrotime = float( parameters["-pyrotime" ] )
	timescale = float(parameters["-timescale" ] )
	pyrotime *= timescale
	radius = float(parameters["-radius"] )
	v1 = RadialPyroclast( center, radius, amp, octaves, freq, rough, trans.magnitude(), pyrotime, pyrogamma )
	return v1


#
#=========== Tube construction ============================
#

def TubeParameters():
	parms = {  "-outerradius":1.0, "-innerradius":0.3, "-length":5.0, "-center":[0,0,0], "-axis":[0,1,0]  }
	return parms


def Tube( parameters ):
	rad = float(parameters[ "-outerradius"])
	innerrad = float(parameters[ "-innerradius"])
	length = float(parameters[ "-length"])
	ce = parameters["-center"]
	center = Vector( float(ce[0]), float(ce[1]), float(ce[2]) )
	ax = parameters["-axis"]
	axis = Vector( float(ax[0]), float(ax[1]), float(ax[2]) )
	v1 = CappedCylinder( center, axis, length, rad )
	v2 = CappedCylinder( center, axis, length, innerrad )
	v3 = cutout( v1, v2 )
	return v3


#
#=========== Ellipse construction ============================
#

def EllipseParameters():
	parms = {  "-majorradius":3.0, "-minorradius":1.0, "-gamma":0.0, "-center":[0,0,0], "-axis":[0,1,0]  }
	return parms


def ellipse( parameters ):
	majorrad = float(parameters[ "-majorradius"])
	minorrad = float(parameters[ "-minorradius"])
	gamma = float(parameters[ "-gamma"])
	tvalue =parameters[ "-axis"]
	axis = Vector( tvalue[0], tvalue[1], tvalue[2] )
	tvalue =parameters[ "-center"]
	center = Vector( tvalue[0], tvalue[1], tvalue[2] )
	v1 = Ellipse( center, axis, majorrad, minorrad )
	if gamma > 0.0:
            v1 = multiply( v1, pow( v1, gamma-1 ))
	return v1;

