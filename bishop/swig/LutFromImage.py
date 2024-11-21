
from vr import *

def LutFromImage( img ):
	lut = ColorArray()
	nb = img.Width()
	y = img.Height()/2
	for x in range(0,nb):
		pix = img.pixel(x,y)
		red = float(pix[0])
		green = float(pix[1])
		blue = float(pix[2])
		col = Color( red, green, blue, 0.0 )
		lut.push_back( col )
	return lut
