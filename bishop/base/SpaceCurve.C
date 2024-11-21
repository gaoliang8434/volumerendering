
#include "SpaceCurve.h"


using namespace lux;

  
CurveFS::CurveFS(const double Qmax, const double Qmin) : qmax(Qmax), qmin(Qmin) {}
CurveFS::~CurveFS(){}
    
const Vector CurveFS::eval( const float q ) const { return Vector(1,0,0)*q; }
const Vector CurveFS::grad( const float q ) const { return Vector(1,0,0); }

const Vector CurveFS::fsT( const float q ) const { return grad(q).unitvector(); }
const Vector CurveFS::fsN( const float q ) const { return Vector(0,1,0); }
const Vector CurveFS::fsB( const float q ) const { return fsT(q)^fsN(q); }

double CurveFS::speed( const float q ) const { return grad(q).magnitude(); }
double CurveFS::curvature( const float q ) const { return 0.0; }
double CurveFS::torsion  ( const float q ) const { return 0.0; }

double CurveFS::qMax() const { return qmax; }
double CurveFS::qMin() const { return qmin; }



SpaceCurve::SpaceCurve() :  std::shared_ptr<CurveFS>() {}
SpaceCurve::SpaceCurve( CurveFS* f ) :  std::shared_ptr<CurveFS>( f ) {}
SpaceCurve::~SpaceCurve(){}

