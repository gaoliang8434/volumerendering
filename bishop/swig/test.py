#!/usr/bin/python

import os
import sys
import bishop

eye = bishop.Vector(0,0,-10);
view = bishop.Vector(0,0,0);
up = bishop.Vector(0,1,0);


cam = bishop.Camera()
cam.setEyeViewUp( eye, view, up )
cam.setFov( 45.0 )
cam.setAspectRatio( 16.0/9.0 )

