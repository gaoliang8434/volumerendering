
%module bishop
%{
#include "AARectangle.h"
#include "IntervalTree.h"
%}
//%include <boost_shared_ptr.i>
%include <std_shared_ptr.i>
%shared_ptr(AABBBase);

%include "AARectangle.h"
%include "IntervalTree.h"

%template(AABBBase)  std::shared_ptr<lux::AARectangle >;
