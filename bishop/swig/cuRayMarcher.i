
%module bishop
%{
#include "AARectangle.h"
#include "PhaseFunction.h"
#include "cuRayMarcher.h"
%}

%include "std_vector.i"

%include "AARectangle.i"
namespace std
{
%template(FloatArray)       vector<float>;
%template(PixelArray)       vector< vector<float> >;
%template(ColorArray)       vector<lux::Color>;
%template(VectorArray)      vector<lux::Vector>;
%template(FloatVolumeArray) vector<lux::Volume<float>*>;
%template(AARectangleArray) vector<lux::AARectangle>;

}


%include "PhaseFunction.h"
%include "cuRayMarcher.h"
