#!/usr/bin/python


import os
import sys

from bishoputils import *


bbe = BlackBodyEmission()
bbe.setNbLevels( 10 )

for temperature in range(700,30000):
	color = bbe.emission( float(temperature) )
	print str(temperature) + "      " + str(color[0]) + " "  + str(color[1]) + " "  + str(color[2])
