
%module bishop
%{
#include "MoreImplicitVolumes.h"
%}

%include "MoreImplicitVolumes.h"
%template(FloatVolumeArray) vector<Volume<float>* >;

