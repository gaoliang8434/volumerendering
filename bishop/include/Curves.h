
#ifndef __CURVES_H__
#define __CURVES_H__

#include <vector>
#include "SpaceCurve.h"
#include "Noise.h"


namespace lux
{

typedef std::vector<Vector> PointChain;

class FadeCurve : public CurveFS
{
    public:

      FadeCurve( const PointChain& Anchors, const float Qmin, const float Qmax  );
     ~FadeCurve();

     const Vector eval( const float t ) const;
     const Vector grad( const float q ) const;
     const Vector fsN( const float q ) const ;
     double curvature( const float q ) const ;
     double torsion  ( const float q ) const ;

  private:

    PointChain anchors;
    double qfactor;
    int nbSegments;

    void getSegment( const double q, double& t, int& segment ) const;
    double fade( const double t ) const;
    double fadedot( const double t ) const;
    double fadedotdot( const double t ) const;
    double fadedotdotdot( const double t ) const;

    const Vector fsTprime( const float q ) const;
    const Vector fsNprime( const float q ) const;
    const Vector gradgrad( const float q ) const;
    const Vector gradgradgrad( const float q ) const;
};

SpaceCurve fadePath( const PointChain& Anchors );


}
#endif
