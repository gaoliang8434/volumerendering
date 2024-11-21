
#include <cmath>
#include "PhaseFunction.h"

using namespace lux;





FournierForandPhaseFunction::FournierForandPhaseFunction( const float en, const float mu )
{
   nu = (3.0-mu)/2.0;
   delta180 = 4.0/( 3.0*(en-1.0)*(en-1.0) );
   double delta180nupower = std::pow( delta180, nu );
   p1factor =  (1.0-delta180nupower)/( 16.0*M_PI*delta180nupower*(delta180-1.0) );
}

const float FournierForandPhaseFunction::eval( const float theta ) const
{
   double s = std::sin( theta/2.0 );
   if( fabs(theta) < 0.00001 ){ s = sin(0.000001/2.0); }
   double s2 = s*s;
   double delta = delta180*s2;
   double deltanupower = std::pow( (double)delta, (double)nu );
   double c = std::cos(theta);
   if( fabs(theta) < 0.00001 ){ c = cos(0.000001); }
   double c2 = c*c;
   double p1 = p1factor * ( 3.0*c2 - 1.0 );
   double p0 = nu*(1.0-delta) - (1.0-deltanupower) + ( delta*(1.0-deltanupower) -nu*(1.0-delta) )/s2;
   p0 /= 4.0*3.14159265*(1.0-delta)*(1.0-delta)*deltanupower;
   float value = p0 + p1;
   return value;
}

