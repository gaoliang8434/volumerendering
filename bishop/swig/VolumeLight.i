
%module bishop
%{
#include "VolumeLight.h"
%}

//%include <boost_shared_ptr.i>
%include <std_shared_ptr.i>
%shared_ptr(VLBase);



%include "VolumeLight.h"

%template(VLBase)  std::shared_ptr<lux::VolumeLight >;
