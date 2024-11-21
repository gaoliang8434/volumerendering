#!/usr/bin/python

import os
import sys
from bishoputils import *

#INCLUDES = "-I../include -I/usr/include/python2.7 -I/usr/lib/python2.7/config/ "
INCLUDES = "-I../include -I/usr/include/python2.6 -I/usr/lib/python2.6/config/ "
SWIGLIB = "-shared"
#SWIGLIB = "-bundle -undefined suppress /usr/lib/libc.dylib  -L/usr/local/lib -L/opt/local/lib"

moduleIncludes = [
"#include \"AsciiParser.h\"",
"#include \"Ballistics.h\"",
"#include \"BlackBody.h\"",
"#include \"Camera.h\"",
"#include \"CmdLineFind.h\"",
"#include \"Color.h\"",
"#include \"FFTNoise.h\"",
"#include \"Fields.h\"",
"#include \"Forms.h\"",
"#include \"Frame.h\"",
"#include \"FrameViewer.h\"",
"#include \"GaussianPRN.h\"",
"#include \"GeometryVolumeShapes.h\"",
"#include \"GridVolumes.h\"",
"#include \"Grids.h\"",
"#include \"ImplicitColors.h\"",
"#include \"ImplicitFormShapes.h\"",
"#include \"ImplicitMatrixShapes.h\"",
"#include \"ImplicitVectorShapes.h\"",
"#include \"ImplicitVolumeShapes.h\"",
"#include \"LABLogo.h\"",
"#include \"LAFINGLogo.h\"",
"#include \"LinearAlgebra.h\"",
"#include \"Logger.h\"",
"#include \"LognormalPRN.h\"",
"#include \"Matrix.h\"",
"#include \"Meshes.h\"",
"#include \"MoreImplicitVolumes.h\"",
"#include \"Noise.h\"",
"#include \"NoiseMachine.h\"",
"#include \"ObjParser.h\"",
"#include \"Particle.h\"",
"#include \"PerlinNoise.h\"",
"#include \"PlyParser.h\"",
"#include \"PoissonSolvers.h\"",
"#include \"ProgressMeter.h\"",
"#include \"RayMarcher.h\"",
"#include \"RectangularGrid.h\"",
"#include \"RenderPatchAllocator.h\"",
"#include \"SparseGrid.h\"",
"#include \"Stamp.h\"",
"#include \"SurfaceFile.h\"",
"#include \"Tracking.h\"",
"#include \"UniformPRN.h\"",
"#include \"VShader.h\"",
"#include \"Vector.h\"",
"#include \"Volume.h\"",
"#include \"VolumeDots.h\"",
"#include \"VolumeGeometry.h\"",
"#include \"VolumeGrid.h\"",
"#include \"cfd.h\"",
"#include \"version.h\""
                 ]

swigIncludes = [
                  "%include \"std_vector.i\"",
                  "%include \"std_string.i\"",
                  "%include \"Vector.h\"",
                  "%include \"Matrix.h\"",
                  "%include \"Forms.h\"",
                  "%include \"Color.h\"",
		  "%include \"Volume.h\"",
		  "%include \"VolumeGrid.h\"",
                  "%include <boost_shared_ptr.i>",
                  "%shared_ptr(ScalarField);",
                  "%shared_ptr(VectorField);",
                  "%shared_ptr(ColorField);",
                  "%shared_ptr(MatrixField);",
                  "%shared_ptr(FormField);",
                  "%template(ScalarVolume) lux::Volume<float>;",
                  "%template(VectorVolume) lux::Volume<lux::Vector>;",
                  "%template(ColorVolume) lux::Volume<lux::Color>;",
                  "%template(ScalarVolumeGrid) lux::VolumeGrid<float>;",
                  "%template(VectorVolumeGrid) lux::VolumeGrid<lux::Vector>;",
                  "%template(ColorVolumeGrid) lux::VolumeGrid<lux::Color>;",
               ]


shaderList = CmdLineFindArray( "-s", "" )
IList = CmdLineFindArray("-I", "" )
LList = CmdLineFindArray("-L","" )

for inc in IList:
	INCLUDES = INCLUDES + " -I" + inc

for ll in LList:
	SWIGLIB = SWIGLIB + " -L" + ll



for vsfile in shaderList:
	shaderFilePath,shaderFileName = os.path.split(vsfile)
	shaderName,shaderExt = os.path.splitext( shaderFileName )
	print "\n\n\n\tWorking on shader " + shaderName + "\n\n\n"
	swigInterfaceFile =  shaderFilePath + shaderName+".i"
	shaderfile = open( swigInterfaceFile, 'w')
	shaderfile.write( "%module " + shaderName + "\n%{\n"  )
	for inc in moduleIncludes:
		shaderfile.write( inc + "\n" )
	shaderfile.write( "#include \"" + vsfile + "\"\n%}\n\n"  )
	for inc in swigIncludes:
		shaderfile.write( inc + "\n" )
	shaderfile.write( "%include \"" + vsfile + "\"" )
	shaderfile.close()
	swigWrapCFile =  shaderFilePath + shaderName + "_wrap.cxx"
	cmd = "swig  -c++ -python -shadow " + INCLUDES + " " + swigInterfaceFile
	#print cmd
	os.system(cmd)
	swigWrapOFile =  shaderFilePath + shaderName + "_wrap.o"
	cmd = "g++ -Wall -g -O2 -fPIC -c " + swigWrapCFile + " " + INCLUDES + " -o " + swigWrapOFile
	#print cmd
	os.system(cmd)
	swigPythonModule =  shaderFilePath + shaderName + ".py"
	swigSOModule =  "_" + shaderFilePath + shaderName + ".so"
	cmd = "g++ -flat_namespace " + swigWrapOFile + " " + SWIGLIB + " -o " + swigSOModule
	#print cmd
	os.system(cmd)
	# clean up
	os.system( "rm " + swigWrapOFile + " " + swigWrapCFile + " " + swigInterfaceFile )
	print "\n\n\n\tFinished shader " + shaderName + "\n\n\n"



