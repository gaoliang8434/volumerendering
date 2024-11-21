
%module bishop
%{
#include "SpaceCurve.h"
%}

//%include <boost_shared_ptr.i>
%include <std_shared_ptr.i>
%shared_ptr(SpaceCurveBase);



%include "SpaceCurve.h"

%template(SpaceCurvedBase)  std::shared_ptr<lux::CurveFS>;
