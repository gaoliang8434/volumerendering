
#include "Curves.h"
using namespace lux;

FadeCurve::FadeCurve( const PointChain& Anchors, const float Qmin, const float Qmax  ) :
	CurveFS( Qmax, Qmin ),
	anchors (Anchors)
	{
	   qfactor = (double)(Anchors.size()-1)/(Qmax-Qmin);
	}

FadeCurve::~FadeCurve(){}

const Vector FadeCurve::eval( const float q ) const
{
   double t;
   int segment;
   getSegment( q, t, segment );
   double w = fade(t);
   return anchors[segment] * (1.0-w) + anchors[segment+1] * w;
}

const Vector FadeCurve::grad( const float q ) const
{
   double t;
   int segment;
   getSegment( q, t, segment );
   double w = fadedot(t)*qfactor;
   return -anchors[segment] * w + anchors[segment+1] * w;
}

const Vector FadeCurve::gradgrad( const float q ) const
{
   double t;
   int segment;
   getSegment( q, t, segment );
   double w = fadedotdot(t)*qfactor*qfactor;
   return -anchors[segment] * w + anchors[segment+1] * w;
}

const Vector FadeCurve::gradgradgrad( const float q ) const
{
   double t;
   int segment;
   getSegment( q, t, segment );
   double w = fadedotdotdot(t)*qfactor*qfactor*qfactor;
   return -anchors[segment] * w + anchors[segment+1] * w;
}

const Vector FadeCurve::fsTprime( const float q ) const
{
   Vector that = fsT(q);
   Vector xpp = gradgrad(q);
   Vector xp = grad(q);
   return ( xpp - that*(that*xpp) )/xp.magnitude();
}

const Vector FadeCurve::fsNprime( const float q ) const
{
   Vector that = fsT(q);
   Vector nhat = fsN(q);
   Vector thatprime = fsTprime(q);
   Vector xp = grad(q);
   Vector xpp = gradgrad(q);
   Vector xppp = gradgradgrad(q);
   Vector thatprimeprime = ( xppp - that*(that*xppp) - 2.0*thatprime*(that*xpp) - that*(thatprime*xpp) )/xp.magnitude();
   return ( thatprimeprime - nhat*(nhat*thatprimeprime) )/ thatprime.magnitude();
}

const Vector FadeCurve::fsN( const float q ) const
{
   return fsTprime(q).unitvector();
}

double FadeCurve::curvature( const float q ) const 
{
   return fsTprime(q).magnitude();
}

double FadeCurve::torsion  ( const float q ) const
{
   return fsNprime(q)*fsB(q);
}

void FadeCurve::getSegment( const double q, double& t, int& segment ) const
{
   t = (q-qmin)*qfactor;
   segment = (int)t;
   if (segment<0)
   {
      segment = 0;
      t = 0;
   }
   if (segment>=(int)anchors.size()-1 )
   {
      segment = (int)anchors.size() - 2;
      t = 1;
   }
}

double FadeCurve::fade( const double t ) const
{
    return t * t * t * ( t * ( t * 6.0 - 15.0 ) + 10.0 );
}

double FadeCurve::fadedot( const double t ) const
{
    return t * ( 20.0 + t * (24.0 * t - 45.0) );
}

double FadeCurve::fadedotdot( const double t ) const
{
    return 20.0 + t * (72.0 * t - 90.0);
}


double FadeCurve::fadedotdotdot( const double t ) const
{
    return 144.0 * t - 90.0;
}

SpaceCurve lux::fadePath( const PointChain& Anchors )
     { return SpaceCurve( new FadeCurve( Anchors, 0.0, 1.0 ) ); }




