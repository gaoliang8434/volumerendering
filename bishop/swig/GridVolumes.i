
%module bishop
%{
#include "GridVolumes.h"
%}

%include "GridVolumes.h"
%template(FloatSample) lux::Sample<float>;
%template(VectorSample) lux::Sample<lux::Vector>;
%template(ColorSample) lux::Sample<lux::Color>;
FloatSample( FloatVolumeGrid*, const FloatVolume* );
VectorSample( VectorVolumeGrid*, const VectorVolume* );
ColorSample( ColorVolumeGrid*, const ColorVolume* );


