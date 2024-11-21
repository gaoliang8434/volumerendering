

#ifndef __BALLISTICS_H__
#define __BALLISTICS_H__

#include "Particle.h"
#include "VolumeGeometry.h"
#include "SpaceCurve.h"

namespace lux{

// generates a string of particles along a trajectory, filling in lifetime, age, id, P, v, accel, age
void BallisticString( const Vector X0, const Vector V0, const Vector A0, const float startTime, const float currentTime, const float targetTime, const int maxNbParticles, ParticleGroupA& particles );

void TargetedBallisticString( const Vector X0, const Vector X1, const Vector A0, const float startTime, const float currentTime, const float targetTime, const int maxNbParticles, ParticleGroupA& particles );

void ParticlesOnGeometry( TriangleGeometry& geom, ParticleGroupA& particles );


class BallisticPath : public CurveFS
{
    public:

      BallisticPath( const Vector X0, const Vector V0, const Vector A0, const float startTime, const float targetTime  );
     ~BallisticPath();

     const Vector eval( const float t ) const;
     const Vector grad( const float q ) const;
     const Vector fsN( const float q ) const ;
     double curvature( const float q ) const ;
     double torsion  ( const float q ) const ;

  private:

    Vector Xstart, Velocity, Acceleration;
};



class TargetedBallisticPath : public BallisticPath
{
    public:

      TargetedBallisticPath( const Vector X0, const Vector X1, const Vector A0, const float startTime, const float targetTime  );
     ~TargetedBallisticPath();

};


SpaceCurve ballisticPath( const Vector X0, const Vector V0, const Vector A0, const float startTime, const float targetTime );
SpaceCurve targetedPath ( const Vector X0, const Vector X1, const Vector A0, const float startTime, const float targetTime );
SpaceCurve splinePath   ( const Vector X0, const Vector X1, const Vector X2, const float startTime, const float midTime, const float endTime );



}
#endif
