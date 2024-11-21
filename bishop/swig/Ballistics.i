
%module bishop
%{
#include "Ballistics.h"
%}

namespace lux{
SpaceCurve ballisticPath( const Vector X0, const Vector V0, const Vector A0, const float startTime, const float targetTime );
SpaceCurve targetedPath( const Vector X0, const Vector X1, const Vector A0, const float startTime, const float targetTime );
SpaceCurve splinePath( const Vector X0, const Vector X1, const Vector A0, const float startTime, const float midTime, const float targetTime );

}
