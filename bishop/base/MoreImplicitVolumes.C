
#include <iostream>
#include "MoreImplicitVolumes.h"

using namespace std;
using namespace lux;


SpaceCurveVolume::SpaceCurveVolume( const SpaceCurve& c, const float rad ) : curve(c), radius(rad) {}
SpaceCurveVolume::~SpaceCurveVolume(){}

const float SpaceCurveVolume::eval( const Vector& P ) const
{
     double qmin = curve->qMin();
     double qmax = curve->qMax();

     double tolerance = (qmax-qmin)/1000.0;
     double q = (qmin+qmax)*0.5;
     double qold = q;
     int count = 0;
     for( int n=0;n<1000;n++ )
     {
        ++count;
        Vector df = (P-curve->eval(q));
	double H = df*curve->fsT(q);
	double dHdq = -curve->speed(q) + df*curve->fsN(q)*curve->curvature(q);
	qold = q;
	q = qold - H/dHdq;
	if (fabs(q-qold) < tolerance){ break; }
     }
     //cout << "eval interations: " << count << endl;
     if( q < qmin ){ q = qmin; }
     if( q > qmax ){ q = qmax; }
     double distq = (P-curve->eval(q)).magnitude();
     double distqmin = (P-curve->eval(qmin)).magnitude();
     double distqmax = (P-curve->eval(qmax)).magnitude();

     if( distqmin < distq )
     {
        q = qmin;
	distq = distqmin;
     }
     if( distqmax < distq )
     {
        q = qmax;
	distq = distqmax;
     }
     return radius - distq;
}

const Vector SpaceCurveVolume::grad( const Vector& P ) const
{
   double dx = 0.01;
   float v0 = eval(P);
   float vx = eval(P+Vector(dx,0,0));
   float vy = eval(P+Vector(0,dx,0));
   float vz = eval(P+Vector(0,0,dx));
   return Vector( vx-v0, vy-v0, vz-v0  )/dx;
}





PiecewiseCurveVolume::PiecewiseCurveVolume( const AnchorChain& c ) : anchorList(c) {}
PiecewiseCurveVolume::~PiecewiseCurveVolume() {}

const float PiecewiseCurveVolume::eval( const Vector& P ) const
{
   int a = -1;
   double pos = 0;
   findClosestPoint( P, a, pos );
   if( a < 0 ){ return -1.0e10; }
   Anchor anch;
   interpolateAnchor( a, pos, anch );
   return anch.radius- distanceToSegment( P, a );
}

const Vector PiecewiseCurveVolume::grad( const Vector& P ) const
{
   double dx = 0.01;
   float v0 = eval(P);
   float vx = eval(P+Vector(dx,0,0));
   float vy = eval(P+Vector(0,dx,0));
   float vz = eval(P+Vector(0,0,dx));
   return Vector( vx-v0, vy-v0, vz-v0  )/dx;
}


void PiecewiseCurveVolume::findClosestPoint( const Vector& P, int& anchor, double& position ) const
{
   double distance = distanceToSegment( P, 0 );
   anchor = 0;
   position = positionOnSegment( P, 0 );
   for( size_t i=1;i<anchorList.size()-1;i++ )
   {
      double dist = distanceToSegment( P, i );
      if( dist < distance )
      { 
         distance = dist;
	 anchor = i;
	 position = positionOnSegment( P, i );
      }
   }
}


void PiecewiseCurveVolume::interpolateAnchor( const int anchor, const double w, Anchor& value ) const
{
double w0 = 1.0 - w;
value.frequency                 = anchorList[anchor].frequency * w0                 + anchorList[anchor+1].frequency * w;
value.translate                 = anchorList[anchor].translate * w0                 + anchorList[anchor+1].translate * w; 
value.octaves                   = anchorList[anchor].octaves * w0                   + anchorList[anchor+1].octaves * w;
value.amplitude                 = anchorList[anchor].amplitude * w0                 + anchorList[anchor+1].amplitude * w;
value.offset                    = anchorList[anchor].offset * w0                    + anchorList[anchor+1].offset * w;
value.fjump                     = anchorList[anchor].fjump * w0                     + anchorList[anchor+1].fjump * w;
value.roughness                 = anchorList[anchor].roughness * w0                 + anchorList[anchor+1].roughness * w;
value.radius                    = anchorList[anchor].radius * w0                    + anchorList[anchor+1].radius * w;
value.capradius                 = anchorList[anchor].capradius * w0                 + anchorList[anchor+1].capradius * w;
value.pscale                    = anchorList[anchor].pscale * w0                    + anchorList[anchor+1].pscale * w;
value.amplitude                 = anchorList[anchor].amplitude * w0                 + anchorList[anchor+1].amplitude * w;
value.gamma                     = anchorList[anchor].gamma * w0                     + anchorList[anchor+1].gamma * w;
value.time                      = anchorList[anchor].time * w0                      + anchorList[anchor+1].time * w;
value.fftLowCutoff              = anchorList[anchor].fftLowCutoff * w0              + anchorList[anchor+1].fftLowCutoff * w;
value.fftHighCutoff             = anchorList[anchor].fftHighCutoff * w0             + anchorList[anchor+1].fftHighCutoff * w;
value.fftPower                  = anchorList[anchor].fftPower * w0                  + anchorList[anchor+1].fftPower * w;
value.fftNbGridPoints           = anchorList[anchor].fftNbGridPoints * w0           + anchorList[anchor+1].fftNbGridPoints * w;
value.fftLength                 = anchorList[anchor].fftLength * w0                 + anchorList[anchor+1].fftLength * w;
value.lognormalmean             = anchorList[anchor].lognormalmean * w0             + anchorList[anchor+1].lognormalmean * w;
value.gaussianstandarddeviation = anchorList[anchor].gaussianstandarddeviation * w0 + anchorList[anchor+1].gaussianstandarddeviation * w;
value.seed                      = anchorList[anchor].seed * w0                      + anchorList[anchor+1].seed * w;
value.tangent                   = anchorList[anchor].tangent * w0                   + anchorList[anchor+1].tangent * w;
value.normal                    = anchorList[anchor].normal * w0                    + anchorList[anchor+1].normal * w;
value.binormal                  = anchorList[anchor].binormal * w0                  + anchorList[anchor+1].binormal * w;
value.axis                      = anchorList[anchor].axis * w0                      + anchorList[anchor+1].axis * w;
value.angle                     = anchorList[anchor].angle * w0                     + anchorList[anchor+1].angle * w;
value.P                         = anchorList[anchor].P * w0                         + anchorList[anchor+1].P * w;
value.v                         = anchorList[anchor].v * w0                         + anchorList[anchor+1].v * w;
value.A                         = anchorList[anchor].A * w0                         + anchorList[anchor+1].A * w;
}

double PiecewiseCurveVolume::distanceToSegment( const Vector& P, const size_t seg ) const
{
   if( seg >= anchorList.size()-1 ){ return -1.0e10; }

   // first, distant to endpoints
   double distance = (P-anchorList[seg  ].P).magnitude();
   double dist     = (P-anchorList[seg+1].P).magnitude();
   if( dist < distance ){ distance = dist; }

   // distance to line segment
   Vector D = (anchorList[seg+1].P - anchorList[seg].P);
   Vector del = P - anchorList[seg].P;
   double t = (del*D)/(D*D);
   if( t > 0 && t < 1 )
   {
      del -= D*t;
      dist = del.magnitude();
      if( dist < distance ){ distance = dist; }
   }
   return distance;
}

double PiecewiseCurveVolume::positionOnSegment( const Vector& P, const size_t seg ) const
{
   if( seg >= anchorList.size()-1 ){ return 1.0; }

   Vector D = (anchorList[seg+1].P - anchorList[seg].P);
   double t = ((P-anchorList[seg].P)*D)/(D*D);
   if( t < 0 ){ t = 0.0; }
   if( t > 1.0 ){ t = 1.0; }
   return t;
}





PiecewisePyroCurveVolume::PiecewisePyroCurveVolume( const AnchorChain& c ) : 
   PiecewiseCurveVolume(c)
   {}

PiecewisePyroCurveVolume::~PiecewisePyroCurveVolume(){}

const float PiecewisePyroCurveVolume::eval( const Vector& P ) const
{

   vector<float> segmentPositions;
   vector<float> segmentLSvalues;
   for( size_t i=0;i<anchorList.size()-1;i++ )
   {
      Vector D = (anchorList[i+1].P - anchorList[i].P);
      double t = ((P-anchorList[i].P)*D)/(D*D);
      Vector Pperp = P-anchorList[i].P - t*D;
      if( t < 0 )
      {
         float radius = anchorList[i].radius;
	 float value = radius - Pperp.magnitude(); 
	 float value2 = radius + t*D.magnitude();
	 if( value2 < value ){ value = value2; }
	 segmentPositions.push_back( t );
	 segmentLSvalues.push_back(value);
      }
      else if( t > 1 )
      {
         float radius = anchorList[i+1].radius;
	 float value = radius - Pperp.magnitude(); 
	 float value2 = radius + (1.0 - t)*D.magnitude();
	 if( value2 < value ){ value = value2; }
	 segmentPositions.push_back( t );
	 segmentLSvalues.push_back(value);
      }
      else
      {
      }
   }

   return 0;




/*
   int a = -1;
   double pos = 0;
   findClosestPoint( P, a, pos );
   if( a < 0 ){ return -1.0e10; }
   Anchor anch;
   interpolateAnchor( a, pos, anch );
   FractalSum<PerlinNoiseGustavson> fspn;
   fspn.setParameters( anch );
   Vector P0 = P - anch.P;
   Vector noiseP = Vector(  P0*anch.normal, 0, P0*anch.binormal ).unitvector() * anch.radius;

   Vector D = (anchorList[a+1].P - anchorList[a].P);
   double t = ((P-anchorList[a].P)*D)/(D*D);
   if( t < 0.0 )
   {
      return anch.capradius - distanceToSegment( P, a )  + std::pow(fabs( fspn.eval(noiseP) ), anch.gamma) * anch.amplitude;
      //return - distanceToSegment( P, a ) + std::pow(fabs( fspn.eval(noiseP) ), anch.gamma) * anch.amplitude;
   }
   if( t > 1.0 )
   {
      return anch.capradius - distanceToSegment( P, a )  + std::pow(fabs( fspn.eval(noiseP) ), anch.gamma) * anch.amplitude;
      //return - distanceToSegment( P, a ) + std::pow(fabs( fspn.eval(noiseP) ), anch.gamma) * anch.amplitude;
   }
   return anch.radius - distanceToSegment( P, a ) + std::pow(fabs( fspn.eval(noiseP) ), anch.gamma) * anch.amplitude;
*/
}




PiecewiseNoiseCurveVolume::PiecewiseNoiseCurveVolume( const AnchorChain& c ) :
   PiecewiseCurveVolume(c)
   {}

PiecewiseNoiseCurveVolume::~PiecewiseNoiseCurveVolume() {}

const float PiecewiseNoiseCurveVolume::eval( const Vector& P ) const
{
   int a = -1;
   double pos = 0;
   findClosestPoint( P, a, pos );
   if( a < 0 ){ return -1.0e10; }
   Anchor anch;
   interpolateAnchor( a, pos, anch );
   FractalSum<PerlinNoiseGustavson> fspn;
   fspn.setParameters( anch );
   Vector P0 = P - anch.P;
   Vector noiseP = Vector(  P0*anch.tangent, P0*anch.binormal, P0*anch.normal );
   float ff = 1.0 - P0.magnitude()*anch.falloff/anch.radius;
   ff = ( ff < 0 ) ? 0.0 : ff;
   ff = ( ff > 1 ) ? 1.0 : ff;
   ff = fade(ff);
   return fspn.eval(noiseP) * ff * anch.amplitude;
}


const float PiecewiseNoiseCurveVolume::fade( const float t ) const
{
   return t * t * t * ( t * ( t * 6.0 - 15.0 ) + 10.0 ); 
}





BoxedVolume::BoxedVolume( const ScalarField& f, const Vector& llc, const Vector& urc, const float def ) :
   AARectangle(llc,urc),
   elem(f),
   defValue(def)
{}

const float BoxedVolume::eval( const Vector& P ) const
{
   if( !isInside(P) ){ return defValue; }
   return elem->eval(P);
}




