#!/usr/bin/python


import sys
import os
from bishoputils import *


#
#=========== File wisp construction ============================
#

def FileWispParameters():
	parms = {  "-stampNxNyNz":[100, 100, 100], "-sparsesize":16, "-stampsize":[10,10,10], 
	           "-stampllc":[-5,-5,-5], "-blurscale":1.0, 
		   "-particlefiles":["","",""], "-particlecolor":[ [1,1,1], [1,1,1], [1,1,1] ],
                   "-opacityscale":[1, 1, 1], "-readdensity":""
	        }
	return parms


def FileWisps( parameters ):
	densityGrid = SparseGrid()
	colorGrid = SparseColorGrid()
	densityGrid.setPartitionSize( int(parameters["-sparsesize"]) )
	colorGrid.setPartitionSize( int(parameters["-sparsesize"]) )
	NNN = parameters["-stampNxNyNz"]
	lX = parameters["-stampsize"]
	X0 = Vector( float(parameters["-stampllc"][0]),  float(parameters["-stampllc"][1]),  float(parameters["-stampllc"][2])  )
	densityGrid.init( int(NNN[0]), int(NNN[1]), int(NNN[2]), float(lX[0]), float(lX[1]), float(lX[2]), X0 )
	colorGrid.init( int(NNN[0]), int(NNN[1]), int(NNN[2]), float(lX[0]), float(lX[1]), float(lX[2]), X0 )
	typeColors = ColorArray()
	for color in parameters["-particlecolor"]:
		col = Color( color[0], color[1], color[2], 1.0 )
		typeColors.push_back( col )
	opacities = FloatArray()
	for opac in parameters["-opacityscale"]:
		opacities.push_back( opac )
	files = StringArray(parameters["-particlefiles"])
	StampPointWisps( densityGrid, colorGrid, files, typeColors, opacities, float(parameters["-blurscale"]) )
	return [ densityGrid, colorGrid ]



