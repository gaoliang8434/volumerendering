
OFILES = base/Matrix.o \
	base/Forms.o \
	base/Volume.o \
	base/ImplicitVolumeShapes.o \
	base/ImplicitVectorShapes.o \
	base/ImplicitColors.o \
	base/ImplicitFormShapes.o \
	base/Fields.o \
	base/LinearAlgebra.o \
	base/PerlinNoise.o \
	base/Camera.o \
	base/RayMarcher.o \
	base/Color.o \
	base/RenderPatchAllocator.o \
	base/HighOrderInterpolator.o \
	base/RectangularGrid.o \
	base/VolumeDots.o \
	base/AsciiParser.o \
	base/GridVolumes.o \
	base/VolumeGeometry.o \
	base/ObjParser.o \
	base/GeometryVolumeShapes.o \
	base/VolumeGrid.o \
	base/PPMFiles.o  \
	base/FFTNoise.o \
	base/Stamp.o \
	base/Ballistics.o \
	base/Tracking.o \
	base/version.o \
	base/BlackBody.o \
	base/LABLogo.o \
	base/cfd.o \
	base/PoissonSolvers.o \
	base/Grids.o \
	base/ImplicitMatrixShapes.o \
	base/Image.o \
	base/NoiseMachine.o \
	base/Meshes.o \
	base/PhaseFunction.o \
	base/FrustumGrid.o \
	base/SpaceCurve.o \
	base/MoreImplicitVolumes.o \
	base/LognormalPRN.o \
	base/Curves.o \
	base/Noise.o \
	base/Incompressible.o \
	base/ParticleGroup.o  \
	base/PointCloud.o  \
	base/PointCloudVolumes.o  \
	base/LevelSetGeometry.o  \
	base/AARectangle.o  \
	base/VolumeLight.o  \
	base/IntervalTree.o  \
	base/SparseGrid.o \
	base/RayMarcher.o 
#	base/OIIOFiles.o \
	base/OpenVDBFiles.o

CUFILES = base/SparseGrid.o \
	base/cuRayMarcher.o

HFILES = include/AARectangle.h		include/GridVolumes.h		include/Noise.h			include/TextureViewer.h \
	include/AsciiParser.h		include/Grids.h			include/NoiseMachine.h		include/Tracking.h \
	include/Ballistics.h		include/ImageFiles.h		include/ObjParser.h		include/UniformPRN.h \
	include/BlackBody.h		include/ImplicitColors.h	include/Particle.h		include/VShader.h \
	include/Camera.h		include/ImplicitFormShapes.h	include/PerlinNoise.h		include/Vector.h \
	include/CmdLineFind.h		include/ImplicitMatrixShapes.h	include/PhaseFunction.h		include/Viewer.h \
	include/Color.h			include/ImplicitVectorShapes.h	include/PlyParser.h		include/Volume.h \
	include/Curves.h		include/ImplicitVolumeShapes.h	include/PoissonSolvers.h	include/VolumeDots.h \
	include/FFTNoise.h		include/LABLogo.h		include/ProgressMeter.h		include/VolumeGeometry.h \
	include/Fields.h		include/LAFINGLogo.h		include/RayMarcher.h		include/VolumeGrid.h \
	include/Forms.h			include/LinearAlgebra.h		include/RectangularGrid.h	include/cfd.h \
	include/Frame.h			include/Logger.h		include/RenderPatchAllocator.h	include/clSparseGrid.h \
	include/FrameViewer.h		include/LognormalPRN.h		include/SpaceCurve.h		include/version.h \
	include/FrustumGrid.h		include/Matrix.h		include/SparseGrid.h \
	include/GaussianPRN.h		include/Meshes.h	include/OIIOFIles.h		include/Stamp.h \
	include/GeometryVolumeShapes.h	include/MoreImplicitVolumes.h  include/SurfaceFile.h  include/HighOrderInterpolator.h

AFILES = $(OFILES)

CFDFILES =  base/cfd.o \
	    base/PoissonSolvers.o




ROOTDIR = .
LIB = -L$(ROOTDIR)/lib -lbishop  -lm
CFDLIB = $(ROOTDIR)/lib/libishopCFD.a

MAGICKLIBS =   -L/usr/local/lib  -L/opt/local/lib -lfftw3
#MAGICKLIBS =  `Magick++-config  --cppflags --cxxflags --ldflags --libs` -L/usr/local/lib  -L/opt/local/lib -lfftw3


ifeq ("$(shell uname)", "Darwin")
  GLLDFLAGS     = -framework Foundation -framework GLUT -framework OpenGL -lm
else
  ifeq ("$(shell uname)", "Linux")
    GLLDFLAGS     = -lglut -lGL -lm -lGLU
  endif
endif
OIIO_LIB_LOCATION = /home/liangg/Desktop/3rdparty/oiio-master/dist/linux64/lib
OIIOLIB = -L${OIIO_LIB_LOCATION} -lOpenImageIO

INSTALL_TOOL=


#CXX = /usr/local/cuda/bin/nvcc -g -O2 --compiler-options '-fPIC -fopenmp' $(DEFINES) -std=c++11 -D_FORCE_INLINES 
CXX = g++ -shared -g -O2 -fPIC -fopenmp $(DEFINES) -fopenmp -std=c++11
CFLAGS = -g -O2 --compiler-options '-fPIC -fopenmp' $(DEFINES) -std=c++11
DEFINES = -D__REVISION__="$(REVISION)"
#REVISION = $(shell svn info | grep Revision | sed "s/Revision: //g")
REVISION = "unknown"

SWIGCXX = g++ -shared -g -O2 -fPIC -fopenmp $(DEFINES) -fopenmp -std=c++11


ifeq ("$(shell uname)", "Darwin")
  CXX = clang-omp++ -Wall -g -O2 -fPIC $(DEFINES) -fopenmp -flat_namespace
  SWIGCXX = $(CXX)
  OIIOMACH = macosx
  OIIOLIB = -L/usr/local/lib -lOpenImageIO -L/opt/local/lib/libomp -lomp 
  #OIIOLIB = $(ROOTDIR)/../OpenImageIO/dist/$(OIIOMACH)/lib/libOpenImageIO.dylib  -L../imageTools/lib -lLaffImageTools
  SWIGLIB = -bundle -undefined suppress /usr/lib/libc.dylib
  CHROMELIB = -lfftw3
  PYTHONINCLUDE =  `python2.7-config --cflags --ldflags`
  #PYTHONINCLUDE = -I/usr/include/python3.3 -I/usr/lib/python3.3/config/
  SWIGLD = $(CXX)
  SWIGEXEC = swig
  #SWIGEXEC = ../swig-2.0.4/swig -I../swig/swig-2.0.4/Lib/
  OIIO_LIB_LOCATION = ../3rdbuild/lib
  OIIO_INCLUDE_LOCATION = ../3rdparty/3rdbuild/include -I /opt/local/include
  INSTALL_TOOL=sudo install_name_tool -change libOpenImageIO.1.7.dylib /usr/local/lib/libOpenImageIO.1.7.dylib swig/_bishop.so
  CUDALIB = 
  DATAINTERFACELIB = 
else
  ifeq ("$(shell uname)", "Linux")
  OIIOMACH = linux64
  OIIOLIB = -ldl -lm -L${OIIO_LIB_LOCATION} -lOpenImageIO -L/home/liangg/Desktop/3rdparty/openvdb -lopenvdb -lfftw3
  #OIIOLIB = -L/usr/local/lib/ -lOpenImageIO
  #OIIOLIB = -L/group/dpa/lib/ -lOpenImageIO
  #OIIOLIB = /group/dpa/local/oiio/lib/libOpenImageIO.so
  #OIIOLIB = /usr/local/lib/libOpenImageIO.so  -L../imageTools/lib -lLaffImageTools /usr/local/lib/libOpenImageIO.so
  #OIIOLIB = $(ROOTDIR)/../OpenImageIO/dist/$(OIIOMACH)/lib/libOpenImageIO.so
  SWIGLIB =
  CAPTUREALL = -Wl,-rpath,${ROOTDIR}/lib -L${ROOTDIR}/lib -Wl,--whole-archive -lbishop -Wl,--no-whole-archive
  CHROMELIB = -ldl
  DATAINTERFACELIB = -L${ROOTDIR}/../data_interface/lib -lDI
  CUDALIB = -L/usr/local/cuda/lib64 -lcudart
  PYTHONINCLUDE = -I/usr/include/python2.7 -I/usr/lib/python2.7/config
  #PYTHONINCLUDE = -I/usr/include/python2.4 -I/usr/lib/python2.4/config
  SWIGLD = $(CXX) -shared
  SWIGEXEC = swig
  #SWIGEXEC = ../swig/swig-2.0.4/swig -I../swig/swig-2.0.4/Lib/
  OIIO_LIB_LOCATION = /home/liangg/Desktop/3rdparty/oiio-master/dist/linux64/lib

  OIIO_INCLUDE_LOCATION =  /home/liangg/Desktop/3rdparty/oiio-master/dist/linux64/include
  endif
 endif



INCLUDES = -I ./include/ $(PYTHONINCLUDE) -I ${OIIO_INCLUDE_LOCATION} -I ../ -I ./
CU_INCLUDES = -I /usr/local/cuda/samples/common/inc
#INCLUDES = -I ./include/ -I ../Viewer/ -I ../OpenImageIO/dist/linux64/include/ `Magick++-config  --cppflags --cxxflags` -I /usr/local/include/ -I /opt/local/include/ $(PYTHONINCLUDE)
#INCLUDES = -I ./include/ -I ../imageTools/include -I ../Viewer/ `Magick++-config  --cppflags --cxxflags` -I /usr/local/include/ -I /opt/local/include/ $(PYTHONINCLUDE)




#all: $(AFILES) $(CUFILES)
all: $(AFILES) 
	ar rv ./lib/libbishop.a $?

.C.o: $<
	$(CXX) -x c++ -c $(INCLUDES) $< -o $@

#base/SparseGrid.o: base/SparseGrid.cu
#	$(CXX) -c $(INCLUDES) $(CU_INCLUDES) base/SparseGrid.cu -o $@


#base/cuRayMarcher.o: base/cuRayMarcher.cu cuda/ssRaymarchKernels.cu
#	$(CXX) -c $(INCLUDES) $(CU_INCLUDES) base/cuRayMarcher.cu -o $@


noiseimage:	exec/noiseimage.C $(AFILES)
	$(CXX) exec/noiseimage.C $(INCLUDES)  -L$(ROOTDIR)/lib -lVR $(OIIOLIB) $(MAGICKLIBS) -o bin/noiseimage

demo:	exec/demo.C MagickFiles.o
	$(CXX) exec/demo.C $(INCLUDES)  MagickFiles.o  $(MAGICKLIBS) -o bin/demo

combine:	exec/combine.C
	$(CXX) exec/combine.C $(INCLUDES)  $(MAGICKLIBS) -o bin/combine

doc:	$(HFILES)
	doxygen doc/Doxyfile
	rm doc/html/files.html
	rm doc/html/_*.html

cleandoc:
	rm -rf doc/html

clean:
	rm -rf *.o bin/* base/*.o cuda/*.o core ./lib/libbishop.a  *~ demo prntest noiseimage combine *.dSYM swig/*.cxx swig/*~ swig/*.so swig/*.o volumeShaders/*.i volumeShaders/*.so volumeShaders/*.py volumeShaders/*~ swig/bishop.py* doc/html doc/latex

#genswig:	swig/bishop.i	$(OFILES) $(CUFILES)
genswig:	swig/bishop.i	$(OFILES)
	$(SWIGEXEC) -c++ -python -shadow -I./include/ swig/bishop.i
	$(SWIGCXX) -c swig/bishop_wrap.cxx  $(INCLUDES) -o swig/bishop_wrap.o
	$(SWIGCXX) $(DEFINES)  swig/bishop_wrap.o $(SWIGLIB) $(LIB) $(DATAINTERFACELIB) $(CUDALIB) $(OIIOLIB) $(CHROMELIB) $(LIB) -L/opt/local/lib -ltbb -o swig/_bishop.so
	cat swig/bishopdoc.py >> swig/bishop.py
	$(INSTALL_TOOL)
