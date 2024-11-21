
%module bishop
%{
#include "ParticleGroup.h"
%}

%include <std_shared_ptr.i>
//%include <boost_shared_ptr.i>
%shared_ptr(PointCloudBase);

%include "ParticleGroup.h"
%template(PointCloudBase)  std::shared_ptr<lux::ParticleGroup>;
