#!/usr/bin/python

import os, sys


window = "0x3400006"

filebase = "surface-0006."

for i in range(1,101):
	raw_input("Ready for frame " + str(i) )
	padframe = str(i)
	if i < 1000:
		padframe = "0" + padframe
	if i < 100:
		padframe = "0" + padframe
	if i < 10:
		padframe = "0" + padframe
	cmd = "import -window " + window + " " + filebase + padframe + ".jpg"
	os.system(cmd)
 
