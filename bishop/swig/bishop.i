
%module bishop 



%include <std_string.i>
//%include <boost_shared_ptr.i>
%include <std_shared_ptr.i>

%include "std_vector.i"
%template(StringArray) std::vector<std::string>;
%include "Vector.i"
%include "Matrix.i"
%include "Camera.i"
%include "Image.i"
%include "Color.i"
%include "AARectangle.i"
%include "Volume.i"
%include "VolumeLight.i"
%include "ProgressMeter.i"
%include "Noise.i"
%include "NoiseMachine.i"
%include "LognormalPRN.i"
%include "RectangularGrid.i"
%include "FrustumGrid.i"
%include "VolumeGrid.i"
%include "SparseGrid.i"
%include "SpaceCurve.i"
%include "Curves.i"
%include "Ballistics.i"
%include "ImplicitVolumeShapes.i"
%include "ImplicitVectorShapes.i"
%include "ImplicitColors.i"
%include "Field.i"
%include "version.i"
%include "RayMarcher.i"
#%include "cuRayMarcher.i"
%include "Stamp.i"
%include "LABLogo.i"
%include "Logger.i"
%include "BlackBody.i"
%include "VolumeGeometry.i"
%include "LevelSetGeometry.i"
%include "ObjParser.i"
%include "GridVolumes.i"
%include "cfd.i"
%include "PoissonSolvers.i"
%include "Particles.i"
%include "Grids.i"
%include "Mesh.i"
%include "LinearAlgebra.i"
%include "Incompressible.i"
%include "ParticleGroup.i"
%include "PointCloud.i"

