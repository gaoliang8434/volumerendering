
#include "Noise.h"

using namespace lux;



const Anchor lux::evaluateAnchorChain( const AnchorChain& anchorList, const double x )
{
   int anchorsize = (int)anchorList.size();
   int anchor = x * anchorsize;
   double w = x * anchorsize - (double)anchor;
   if( anchor < 0 )
   { 
      anchor = 0;
      w = 0;
   }
   if( anchor >= anchorsize-1 )
   {
      anchor = anchorsize-2;
      w = 1;
   }
   Anchor value = interpolateAnchors( anchorList[anchor], anchorList[anchor+1], w );
   return value;
}


const Anchor lux::interpolateAnchors( const Anchor& a1, const Anchor& a2, const double w )
{
double w0 = 1.0 - w;
Anchor value;
value.frequency                 = a1.frequency * w0                 + a2.frequency * w;
value.translate                 = a1.translate * w0                 + a2.translate * w; 
value.octaves                   = a1.octaves * w0                   + a2.octaves * w;
value.amplitude                 = a1.amplitude * w0                 + a2.amplitude * w;
value.offset                    = a1.offset * w0                    + a2.offset * w;
value.fjump                     = a1.fjump * w0                     + a2.fjump * w;
value.roughness                 = a1.roughness * w0                 + a2.roughness * w;
value.radius                    = a1.radius * w0                    + a2.radius * w;
value.capradius                 = a1.capradius * w0                 + a2.capradius * w;
value.pscale                    = a1.pscale * w0                    + a2.pscale * w;
value.amplitude                 = a1.amplitude * w0                 + a2.amplitude * w;
value.gamma                     = a1.gamma * w0                     + a2.gamma * w;
value.time                      = a1.time * w0                      + a2.time * w;
value.fftLowCutoff              = a1.fftLowCutoff * w0              + a2.fftLowCutoff * w;
value.fftHighCutoff             = a1.fftHighCutoff * w0             + a2.fftHighCutoff * w;
value.fftPower                  = a1.fftPower * w0                  + a2.fftPower * w;
value.fftNbGridPoints           = a1.fftNbGridPoints * w0           + a2.fftNbGridPoints * w;
value.fftLength                 = a1.fftLength * w0                 + a2.fftLength * w;
value.lognormalmean             = a1.lognormalmean * w0             + a2.lognormalmean * w;
value.gaussianstandarddeviation = a1.gaussianstandarddeviation * w0 + a2.gaussianstandarddeviation * w;
value.seed                      = a1.seed * w0                      + a2.seed * w;
value.tangent                   = (a1.tangent * w0                   + a2.tangent * w).unitvector();
value.normal                    = a1.normal * w0                    + a2.normal * w;
value.normal -= (value.normal*value.tangent) * value.tangent;
value.normal.normalize();
value.binormal                  = value.tangent ^ value.normal;
value.axis                      = (a1.axis * w0                      + a2.axis * w ).unitvector();
value.angle                     = a1.angle * w0                     + a2.angle * w;
value.P                         = a1.P * w0                         + a2.P * w;
value.v                         = a1.v * w0                         + a2.v * w;
value.A                         = a1.A * w0                         + a2.A * w;
   return value;
}



void lux::setAnchor( AnchorChain& a, int i, const Anchor& c )
{
   if( i >= 0 && i < (int)a.size() )
   {
      a[i] = c;
   }
}

Anchor lux::getAnchor( AnchorChain& a, int i )
{
   if( i >= 0 && i < (int)a.size() ){ return a[i]; }
}
