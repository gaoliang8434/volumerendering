
#ifndef __FIELDS_H__
#define __FIELDS_H__

#include "Volume.h"
#include "VolumeGrid.h"
#include "SparseGrid.h"
#include "SpaceCurve.h"
//#include "OVDBGrid.h"
#include "NoiseMachine.h"


namespace lux
{


const float  evaluate( const ScalarField& v, const Vector& P );
const Vector evaluate( const VectorField& v, const Vector& P );
const Color  evaluate( const ColorField&  v, const Vector& P );
const Matrix evaluate( const MatrixField&  v, const Vector& P );
const Form   evaluate( const FormField&  v, const Vector& P );

const VectorField grad( const ScalarField& e ); 
const VectorField fdgrad( const ScalarField& e, const int N, const float dx, const float dy, const float dz ); 
const VectorField fdboundedgrad( const ScalarField& e, const int N, const GridBox& gb ); 
const VectorField fdinteriorgrad( const ScalarField& e, const int N, const float dx, const float dy, const float dz, const ScalarField& gb ); 
const MatrixField grad( const VectorField& e ); 
const FormField grad( const FormField& e ); 

ScalarField constant( const float v );
VectorField constant( const Vector& v );
ColorField  constant( const Color& v );
MatrixField constant( const Matrix& v );
FormField   constant( const Form& v );

ScalarField scale( const ScalarField& , const Vector& s );
VectorField scale( const VectorField& , const Vector& s );
ColorField  scale( const ColorField& , const Vector& s );
//MatrixField scale( const MatrixField& , const Vector& s );
FormField   scale( const FormField& , const Vector& s );

ScalarField translate( const ScalarField& , const Vector& s );
VectorField translate( const VectorField& , const Vector& s );
ColorField  translate( const ColorField& , const Vector& s );
//MatrixField translate( const MatrixField& , const Vector& s );
FormField   translate( const FormField& , const Vector& s );

ScalarField rotate( const ScalarField& , const Vector& s );
VectorField rotate( const VectorField& , const Vector& s );
ColorField  rotate( const ColorField& , const Vector& s );
//MatrixField rotate( const MatrixField& , const Vector& s );
FormField   rotate( const FormField& f, const Vector& s );
MatrixField rotation( const VectorField& v );

ScalarField exp( const ScalarField& v ); 
MatrixField exp( const MatrixField& v ); 
ScalarField det( const MatrixField& v ); 

MatrixField sinch( const MatrixField& v ); 
MatrixField orderedSinch( const MatrixField& v1, const MatrixField& v2 ); 

ScalarField report( const ScalarField& v, const string& tag ); 
VectorField report( const VectorField& v, const string& tag ); 
ColorField  report( const ColorField& v, const string& tag ); 
//MatrixField report( const MatrixField& v, const string& tag ); 
FormField   report( const FormField& v, const string& tag ); 

ScalarField negate( const ScalarField& v );
VectorField negate( const VectorField& v );
ColorField  negate( const ColorField& v );
MatrixField negate( const MatrixField& v );
FormField   negate( const FormField& v );

ScalarField abs( const ScalarField& v ); 
ScalarField abs( const VectorField& v ); 

ScalarField which( const ScalarField& v1, const ScalarField& v2, const ScalarField& swtch );
VectorField which( const VectorField& v1, const VectorField& v2, const ScalarField& swtch );
ColorField  which( const ColorField& v1, const ColorField& v2, const ScalarField& swtch );
//MatrixField which( const MatrixField& v1, const MatrixField& v2, const ScalarField& swtch );
FormField   which( const FormField& v1, const FormField& v2, const ScalarField& swtch );

ScalarField multiply( const ScalarField& v, const float a ); 
ScalarField multiply( const ScalarField& v, const ScalarField& u ); 
VectorField multiply( const VectorField& v, const float a ); 
VectorField multiply( const VectorField& v, const ScalarField& u ); 
ColorField  multiply( const ColorField& v, const float a ); 
ColorField  multiply( const ColorField& v, const ScalarField& u ); 
ColorField  multiply( const ColorField& v, const ColorField& u ); 
MatrixField multiply( const MatrixField& v, const float a ); 
MatrixField multiply( const MatrixField& v, const ScalarField& u ); 
FormField   multiply( const FormField& v, const float a ); 
FormField   multiply( const FormField& v, const ScalarField& u ); 

ScalarField divide( const ScalarField& v, const float a ); 
ScalarField divide( const ScalarField& v, const ScalarField& u ); 
VectorField divide( const VectorField& v, const float a ); 
VectorField divide( const VectorField& v, const ScalarField& u ); 
ColorField  divide( const ColorField& v, const float a ); 
ColorField  divide( const ColorField& v, const ScalarField& u ); 
MatrixField divide( const MatrixField& v, const float a ); 
MatrixField divide( const MatrixField& v, const ScalarField& u ); 
FormField   divide( const FormField& v, const float a ); 
FormField   divide( const FormField& v, const ScalarField& u ); 

ScalarField add( const ScalarField&  v1, const ScalarField& v2 );
VectorField add( const VectorField&  v1, const VectorField& v2 );
ColorField  add( const ColorField&  v1, const ColorField& v2 );
MatrixField add( const MatrixField&  v1, const MatrixField& v2 );
FormField   add( const FormField&  v1, const FormField& v2 );

ScalarField subtract( const ScalarField&  v1, const ScalarField& v2 );
VectorField subtract( const VectorField&  v1, const VectorField& v2 );
ColorField  subtract( const ColorField&  v1, const ColorField& v2 );
MatrixField subtract( const MatrixField&  v1, const MatrixField& v2 );
FormField   subtract( const FormField&  v1, const FormField& v2 );


ScalarField Sphere( const Vector& cen, const float rad );
ScalarField Ellipse( const Vector& cen, const Vector& axs, const float majorrad, const float minorrad ); 
ScalarField CsgBox( const Vector& cen, const float rad, const float pwr );
ScalarField CsgRectangleBox( const Vector& cen, const float rad, const Vector& asp, const float pwr );
ScalarField HardBox( const Vector _llc, const Vector _urc );
ScalarField Cone( const Vector& cen, const Vector& ax, const float h, const float theta );
ScalarField Plane( const Vector cen, const Vector norm );
ScalarField Torus( const Vector& cen, const Vector& axis, const float majorRad, const float minorRad );
ScalarField SteinerPatch(); 
ScalarField Icosahedron();
ScalarField Cylinder( const Vector axis, const float rad );
ScalarField CappedCylinder( const Vector cen, const Vector axis, const float length, const float radius );
ScalarField Shell( const ScalarField& v, const float thickness );

ScalarField SpaceCurveField( const SpaceCurve& curve, float radius );
ScalarField PiecewiseCurveField( const AnchorChain& list );
ScalarField PiecewisePyroCurveField( const AnchorChain& list );
ScalarField PiecewiseNoiseCurveField( const AnchorChain& list );

ScalarField mask( const ScalarField& v );
ScalarField clamp( const ScalarField& v, float minv, float maxv );
ScalarField pow( const ScalarField& v, float gam );
ScalarField pow( const ScalarField& v, const ScalarField& gam );
ColorField  pow( const ColorField& v, float gam );
ColorField  pow( const ColorField& v, const ScalarField& gam );
ScalarField BlinnBlend( const ScalarField& v1, const ScalarField& v2, const float _alpha = 1.0 );
ScalarField MultiBlend( std::vector<ScalarField>& vs, const float a = 1.0  );
ScalarField Union( const ScalarField& v1, const ScalarField& v2 );
ScalarField intersection( const ScalarField& v1, const ScalarField& v2 );
ScalarField cutout( const ScalarField& v1, const ScalarField& v2 );
ScalarField boxed( const ScalarField& v1, const Vector& llc, const Vector& urc, const float defvalue );

ScalarField sin( const ScalarField& v1 );
ScalarField cos( const ScalarField& v1 );
ScalarField tan( const ScalarField& v1 );
ScalarField acos( const ScalarField& v1 );
ScalarField asin( const ScalarField& v1 );
ScalarField atan( const ScalarField& v1 );
ScalarField sinh( const ScalarField& v1 );
ScalarField cosh( const ScalarField& v1 );
ScalarField tanh( const ScalarField& v1 );

MatrixField outer( const VectorField& v1, const VectorField& v2 );
MatrixField inverse( const MatrixField& m );

ScalarField Pyroclast( const Vector& Center, const float Radius, const float Amp, 
                         const float octaves, const float freq, const float rough, 
                         const Vector trans, const float time, const float Gamma = 1.0/3.0 );

ScalarField RadialPyroclast( const Vector& Center, const float Radius, const float Amp, 
                             const float octaves, const float freq, const float rough, 
                             const float trans, const float time, const float Gamma = 1.0/3.0 );

ScalarField SFFFTNoise( const float power, const float low, const float high, const float length, const int sz );
ScalarField SFNoise( NoiseMachine n, const float d = 0.01 ); 
VectorField VFNoise( NoiseMachine n, const float d = 0.01 ); 
ScalarField JitterSample( const ScalarField& v, const float rad, const int nb, const int seed = 2847573 ); 
//VectorField JitterSample( const VectorField& v, const float rad, const int nb, const int seed = 2847573 ); 
//ColorField JitterSample( const ColorField& v, const float rad, const int nb, const int seed = 2847573 ); 
//MatrixField JitterSample( const MatrixField& v, const float rad, const int nb, const int seed = 2847573 ); 

ScalarField gridded( const VolumeGrid<float>* g );
ScalarField gridded( const SparseGrid* g );
ScalarField gridded( const ScalarGrid& g );
ScalarField gridded( const ScalarFrustumGrid& g );
//ScalarField gridded( const ScalarOGrid& g );
VectorField gridded( const VolumeGrid<Vector>* g );
VectorField gridded( const SparseVectorGrid* g );
VectorField gridded( const VectorGrid& g );
VectorField gridded( const VectorFrustumGrid& g );
//VectorField gridded( const VectorOGrid& g );
ColorField  gridded( const VolumeGrid<Color>* g );
ColorField  gridded( const SparseColorGrid* g );
ColorField  gridded( const ColorGrid& g );
ColorField  gridded( const ColorFrustumGrid& g );
//ColorField  gridded( const ColorOGrid& g );
MatrixField gridded( const MatrixGrid& g );


ScalarField advect( const ScalarField& v, const VectorField& u, const float delt ); 
VectorField advect( const VectorField& v, const VectorField& u, const float delt ); 
ColorField  advect( const ColorField& v, const VectorField& u, const float delt ); 
MatrixField advect( const MatrixField& v, VectorField& u, const float delt ); 
FormField   advect( const FormField& v, VectorField& u, const float delt ); 

VectorField gradientStretchCM( const VectorField& v, const float T, const int nb );

VectorField gradientDisplacement( const MatrixField& M, const ScalarField& LS, const double step );

ScalarField warp( const ScalarField& v, VectorField& map );
VectorField warp( const VectorField& v, VectorField& map );
ColorField  warp( const ColorField& v, VectorField& map );
MatrixField warp( const MatrixField& v, VectorField& map );
FormField   warp( const FormField& v, VectorField& map );

ScalarField Periodic( const ScalarField v, const Vector& o, const Vector& L );
VectorField Periodic( const VectorField v, const Vector& o, const Vector& L );
ColorField  Periodic( const ColorField v, const Vector& o, const Vector& L );
//MatrixField Periodic( const MatrixField v, const Vector& o, const Vector& L );
FormField   Periodic( const FormField v, const Vector& o, const Vector& L );


FormField wedge( const FormField& e1, const FormField& e2 );
FormField star( const FormField& e1 );
FormField contraction( const VectorField& v, const FormField f );

ScalarField dot( const VectorField& e1, const VectorField& e2 );
VectorField unitvector( const VectorField& e );
VectorField identity();
ScalarField xIdentity();
ScalarField yIdentity();
ScalarField zIdentity();
VectorField ImplicitSurfacePoint( const ScalarField& e, const float step, const int nbIterations );
VectorField cross( const VectorField& e1, const VectorField& e2 );
ScalarField div( const VectorField& e );
ScalarField fddiv( const VectorField& e, const int N, const float dx, const float dy, const float dz );
ScalarField fdboundeddiv( const VectorField& e, const int N, const GridBox& gb );
ScalarField fdinteriordiv( const VectorField& e, const int N, const float dx, const float dy, const float dz, const ScalarField& gb );
VectorField curl( const VectorField& e );
VectorField ContinuedFractionDisplacement( const VectorField& dX, const int nbIterations );
VectorField XYZ( const ScalarField& x, const ScalarField& y, const ScalarField& z );

ColorField Chroma( const ColorField& e );
ColorField Blackbody( const ScalarField& e );
ColorField RGB( const ScalarField& red, const ScalarField& green, const ScalarField& blue );
ColorField LUTColor( const ScalarField& field, const vector<Color> lut, const float maxValue, const float minValue );


ScalarField DetGrad( const VectorField& e );

ScalarField rComponent( const ColorField& C  );
ScalarField gComponent( const ColorField& C  );
ScalarField bComponent( const ColorField& C  );

ColorField rgbComponent( const ScalarField& red, const ScalarField& green, const ScalarField& blue );

ScalarField xComponent( const VectorField& X );
ScalarField yComponent( const VectorField& X );
ScalarField zComponent( const VectorField& X );

VectorField component( const ScalarField& X, const ScalarField& Y, const ScalarField& Z );

ScalarField zeroComponent( const FormField& f );
VectorField oneComponent( const FormField& f );
VectorField twoComponent( const FormField& f );
ScalarField threeComponent( const FormField& f );

FormField component( const ScalarField& X, const VectorField& Y, const VectorField& Z, const ScalarField& W );

FormField lie( VectorField& X, FormField& a );


//VectorField bracket( VectorField& X, VectorField& Y );



//
//  Rendering tools
//

float transmissivity( const ScalarField& density, const Vector& start, const Vector& end, const double step, const double scatter );
Color lineIntegral( const ScalarField& density, const ColorField& color, const Vector& start, const Vector& end, const double step, const double scatter );




void fieldStatistics( const ScalarField& f, const GridBox& g );
void fieldStatistics( const VectorField& f, const GridBox& g );


void setFDSize( ScalarField& f, int n );
void setFDSize( VectorField& f, int n );
void setFDSize( MatrixField& f, int n );
void setFDSize( ColorField& f, int n );
void setFDSize( FormField& f, int n );

void setFDStep( ScalarField& f, double dx );
void setFDStep( VectorField& f, double dx );
void setFDStep( MatrixField& f, double dx );
void setFDStep( ColorField& f, double dx );
void setFDStep( FormField& f, double dx );

void setFDStep( ScalarField& f, double dx, double dy, double dz );
void setFDStep( VectorField& f, double dx, double dy, double dz );
void setFDStep( MatrixField& f, double dx, double dy, double dz );
void setFDStep( ColorField& f, double dx, double dy, double dz );
void setFDStep( FormField& f, double dx, double dy, double dz );



}





#endif



