
%module bishop
%{
#include "VolumeGrid.h"
%}

%include "VolumeGrid.h"

%shared_ptr(ScalarGridBase);
%shared_ptr(VectorGridBase);
%shared_ptr(ColorGridBase);


%template(FloatVolumeGrid) lux::VolumeGrid<float>;
%template(VectorVolumeGrid) lux::VolumeGrid<lux::Vector>;
%template(ColorVolumeGrid) lux::VolumeGrid<lux::Color>;


