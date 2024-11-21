

#include "Ballistics.h"

using namespace lux;

void lux::BallisticString( const Vector X0, const Vector V0, const Vector A0, const float startTime, const float currentTime, const float targetTime, const int maxNbParticles, ParticleGroupA& particles )
{
   float dt = (targetTime-startTime)/maxNbParticles;

   float t = startTime;
   int id = 0;
   while( t <=currentTime )
   {
      Particle p;
      p.P() = X0 + t*V0 + 0.5*t*t*A0;
      p.v() = V0 + t*A0;
      if( p.v().magnitude() != 0 )
      { 
         p.normal() = p.v().unitvector();
	 Vector up = A0 - (A0*p.normal())*p.normal();
	 if( up.magnitude() != 0 )
	 {
	    up.normalize();
	    p.up() = up;
	 }
	 else
	 {
	    p.up() = p.up() - (p.up() * p.normal()) * p.normal();
	    p.up().normalize();
	 }
	 p.right() = p.normal() ^ p.up();
      }

      p.accel() = A0;
      p.id() = id;
      p.lifetime() = targetTime - startTime;
      p.age() = currentTime - t;
      particles.push_back( p );
      t += dt;
      ++id;
   }
}


void lux::TargetedBallisticString( const Vector X0, const Vector X1, const Vector A0, const float startTime, const float currentTime, const float targetTime, const int maxNbParticles, ParticleGroupA& particles )
{
   float T = targetTime - startTime;
   Vector V0 = X1 - X0 - (0.5*T*T)*A0;
   V0 /= T;
   BallisticString( X0, V0, A0, startTime, currentTime, targetTime, maxNbParticles, particles );
}




void lux::ParticlesOnGeometry( TriangleGeometry& geom, ParticleGroupA& particles )
{
   geom.computeConnectivity();
   for( size_t i=0;i<geom.nbVertices();i++ )
   {
      Particle p;
      p.P() = geom.getVertex(i);
      p.pscale() = geom.averageNeighborDistance(i);
      p.id() = i;
   }
}


BallisticPath::BallisticPath( const Vector X0, const Vector V0, const Vector A0, const float startTime, const float targetTime  ) :
   CurveFS( targetTime, startTime ),
   Xstart      ( X0 ),
   Velocity    ( V0 ),
   Acceleration( A0 )
   {}

BallisticPath::~BallisticPath(){}

const Vector BallisticPath::eval( const float q ) const
{
   double dt = q - qmin;
   return Xstart + Velocity*dt + Acceleration*dt*dt/2.0;
}

const Vector BallisticPath::grad( const float q ) const
{
   double dt = q - qmin;
   return Velocity + Acceleration*dt;
}

const Vector BallisticPath::fsN( const float q ) const 
{
   Vector V = fsT(q);
   Vector dTdt = (Acceleration - V*(V*Acceleration))/speed(q);
   return dTdt.unitvector();
}


double BallisticPath::curvature( const float q ) const 
{
   Vector V = fsT(q);
   Vector dTdt = (Acceleration - V*(V*Acceleration))/speed(q);
   return dTdt.magnitude();
}

double BallisticPath::torsion  ( const float q ) const 
{
  return 0.0; // ballistic trajectories dont twist
}




TargetedBallisticPath::TargetedBallisticPath( const Vector X0, const Vector X1, const Vector A0, const float startTime, const float targetTime  ) :
         BallisticPath( X0,  (X1 - X0 - (0.5*(targetTime-startTime)*(targetTime-startTime))*A0)/(targetTime-startTime)  , A0, startTime, targetTime )
      {}

TargetedBallisticPath::~TargetedBallisticPath(){}



SpaceCurve lux::ballisticPath( const Vector X0, const Vector V0, const Vector A0, const float startTime, const float targetTime )
   {  return SpaceCurve( new BallisticPath(X0, V0, A0, startTime, targetTime ) ); }

SpaceCurve lux::targetedPath( const Vector X0, const Vector X1, const Vector A0, const float startTime, const float targetTime )
   {  return SpaceCurve( new TargetedBallisticPath( X0, X1, A0, startTime, targetTime ) ); }

SpaceCurve lux::splinePath( const Vector X0, const Vector X1, const Vector X2, const float startTime, const float midTime, const float endTime )
{
   float tbar = ( startTime + midTime + endTime )/3.0;
   float t2bar = ( startTime*startTime + midTime*midTime + endTime*endTime )/3.0;
   float dt0 = startTime - tbar;
   float dt1 = midTime - tbar;
   float dt20 = startTime*startTime - t2bar;
   float dt21 = midTime*midTime - t2bar;
   Vector Xbar = ( X0 + X1 + X2 )/ 3.0;
   Vector C = ( X1-Xbar - (X0-Xbar)*dt1/dt0 );
   C = C/( dt21 - dt20*dt1/dt0 );
   Vector BB = ( X0 - Xbar - C*dt20 )/dt0;
   Vector A = Xbar - BB*tbar - C*t2bar;
   C = C*2.0;
   return lux::ballisticPath( A, BB, C, startTime, endTime );
}
