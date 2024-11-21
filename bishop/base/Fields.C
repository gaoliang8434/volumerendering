
#include "MoreImplicitVolumes.h"
#include "ImplicitVolumeShapes.h"
#include "ImplicitVectorShapes.h"
#include "ImplicitMatrixShapes.h"
#include "ImplicitFormShapes.h"
#include "ImplicitColors.h"
#include "BlackBody.h"
#include "Fields.h"


using namespace std;
namespace lux
{


const float  evaluate( const ScalarField& v, const Vector& P ) { return v->eval(P); }
const Vector evaluate( const VectorField& v, const Vector& P ) { return v->eval(P); }
const Color  evaluate( const ColorField&  v, const Vector& P ) { return v->eval(P); }
const Matrix evaluate( const MatrixField& v, const Vector& P ) { return v->eval(P); }
const Form   evaluate( const FormField& v,   const Vector& P ) { return v->eval(P); }

const ScalarField ScalarField::operator+( const ScalarField& e2 ) { return add(*this,e2); }
const ScalarField ScalarField::operator-( const ScalarField& e2 ) { return subtract(*this,e2); }
const ScalarField ScalarField::operator-() { return negate(*this); }
const ScalarField ScalarField::operator*( const ScalarField& e2 ) { return multiply( *this, e2 ); }
const ScalarField ScalarField::operator/( const ScalarField& e2 ) { return divide(*this,e2); }
const ScalarField ScalarField::operator^( const ScalarField& e2 ) { return cutout(*this,e2); }
const ScalarField ScalarField::operator&&( const ScalarField& e2 ) { return Union(*this,e2); }
const ScalarField ScalarField::operator||( const ScalarField& e2 ) { return intersection(*this,e2); }
const VectorField ScalarField::operator*( const VectorField& e1 ) { return multiply(e1,*this); }
const ColorField  ScalarField::operator*( const ColorField& e1 ) { return multiply(e1,*this); }

const ColorField ColorField::operator+( const ColorField& e2 ) { return add(*this,e2); }
const ColorField ColorField::operator-( const ColorField& e2 ) { return subtract(*this,e2); }
const ColorField ColorField::operator-( ) { return negate(*this); }
const ColorField ColorField::operator*( const ColorField& e2 ) { return multiply(*this,e2); }
const ColorField ColorField::operator*( const ScalarField& e2 ) { return multiply(*this,e2); }
const ColorField ColorField::operator/( const ScalarField& e2 ) { return divide(*this,e2); }

const VectorField VectorField::operator+( const VectorField& e2 ) { return add(*this,e2); }
const VectorField VectorField::operator-( const VectorField& e2 ) { return subtract(*this,e2); }
const VectorField VectorField::operator-( ) { return negate(*this); }
const VectorField VectorField::operator*( const ScalarField& e2 ) { return multiply(*this,e2); }
const ScalarField VectorField::operator*( const VectorField& e2 ) { return dot(*this,e2); }
const VectorField VectorField::operator/( const ScalarField& e2 ) { return divide(*this,e2); }
const VectorField VectorField::operator^( const VectorField& e2 ) { return cross(*this,e2); }
const VectorField grad( const ScalarField& e ) { return VectorField(new GradientVectorVolume(e)); }
const VectorField fdgrad( const ScalarField& e, const int N, const float dx, const float dy, const float dz ) 
        { return VectorField(new FiniteDifferenceGradientVectorVolume(e,N,dx,dy,dz)); }
const VectorField fdboundedgrad( const ScalarField& e, const int N, const GridBox& gb )
        { return VectorField(new FiniteDifferenceBoundedGradientVectorVolume(e,N,gb)); }
const VectorField fdinteriorgrad( const ScalarField& e, const int N, const float dx, const float dy, const float dz, const ScalarField& gb )
        { return VectorField(new FiniteDifferenceInteriorGradientVectorVolume(e,N,dx,dy,dz,gb)); }

const MatrixField grad( const VectorField& e ) { return MatrixField(new GradientMatrixVolume(e)); }
const FormField   grad( const FormField& e )   
{ 
   if( e.isGrad )
   {
      FormField f = constant(Form(0,Vector(0,0,0), Vector(0,0,0), 0) );
      f.isGrad = true;
      return f;
   }
   FormField f = FormField(new GradientFormVolume(e)); 
   f.isGrad = true; 
   return f; 
}


const FormField FormField::operator+( const FormField& e2 ) { return add(*this,e2); }
const FormField FormField::operator-( const FormField& e2 ) { return subtract(*this,e2); }
const FormField FormField::operator-() { return negate(*this); }
const FormField FormField::operator/( const ScalarField& e2 ) { return divide(*this,e2); }
const FormField FormField::operator*( const ScalarField& e1 ) { return multiply(*this,e1); }
const FormField FormField::operator^(  const FormField& e2 )  { return wedge(*this,e2); }
const FormField ScalarField::operator*( const FormField& e1 ) { return multiply(e1,*this); }
const FormField VectorField::operator*( const FormField& e1 ) { return contraction(*this,e1); }


const MatrixField MatrixField::operator+( const MatrixField& e2 ) { return MatrixField( new AddMatrixVolume(*this,e2) ); }
const MatrixField MatrixField::operator-( const MatrixField& e2 ) { return MatrixField( new SubtractMatrixVolume(*this,e2) ); }
const MatrixField MatrixField::operator-() { return MatrixField( new NegateMatrixVolume(*this) ); }
const MatrixField MatrixField::operator*( const ScalarField& e2 ) { return MatrixField( new ScalarMatrixProductVolume( *this, e2 ) ); }
const VectorField MatrixField::operator*( const VectorField& e2 ) { return VectorField( new MatrixVectorProductVolume( e2, *this ) ); }
const MatrixField MatrixField::operator*( const MatrixField& e2 ) { return MatrixField( new MatrixMatrixProductVolume( e2, *this ) ); }
const MatrixField MatrixField::operator/( const ScalarField& e2 ) { return MatrixField( ); }
const VectorField VectorField::operator*( const MatrixField& e2 ) { return VectorField( new VectorMatrixProductVolume( e2, *this ) ); }
const MatrixField ScalarField::operator*( const MatrixField& e1 ) { return MatrixField( new ScalarMatrixProductVolume( e1, *this )  ); }
ScalarField det( const MatrixField& m ) { return ScalarField( new DetMatrixVolume( m ) ); }
MatrixField outer( const VectorField& v1, const VectorField& v2 ) { return MatrixField( new OuterProductMatrixVolume(v1,v2)); }



ScalarField constant( const float v ) { return ScalarField( new ConstantVolume(v) ); }
VectorField constant( const Vector& v ) { return VectorField(new ConstantVectorVolume(v)); }
ColorField  constant( const Color& v ) { return ColorField(new ConstantColor(v)); }
MatrixField constant( const Matrix& v ) { return MatrixField( new ConstantMatrixVolume(v) ); }
FormField   constant( const Form& v ) { return FormField( new ConstantFormVolume(v) ); }


ScalarField which( const ScalarField& v1, const ScalarField& v2, const ScalarField& swtch )
   { return ScalarField( new SwitchVolume(v1,v2,swtch) ); }
VectorField which( const VectorField& v1, const VectorField& v2, const ScalarField& swtch )
   { return VectorField( new SwitchVectorVolume(v1,v2,swtch) ); }
ColorField  which( const ColorField& v1, const ColorField& v2, const ScalarField& swtch )
   { return ColorField( new SwitchColorVolume(v1,v2,swtch) ); }
FormField  which( const FormField& v1, const FormField& v2, const ScalarField& swtch )
   { return FormField( new SwitchFormVolume(v1,v2,swtch) ); }


ScalarField scale( const ScalarField& v, const Vector& s ) { return ScalarField(new ScaleVolume(v,s)); }
VectorField scale( const VectorField& v, const Vector& s ) { return VectorField(new ScaleVectorVolume(v,s)); }
ColorField scale( const ColorField& v, const Vector& s ) { return ColorField(new ScaleColor(v,s)); }
//MatrixField scale( const MatrixField& v, const Vector& s );
FormField scale( const FormField& v, const Vector& s ) { return FormField(new ScaleFormVolume(v,s)); }

ScalarField translate( const ScalarField& v, const Vector& s ) { return ScalarField(new TranslateVolume(v,s)); }
VectorField translate( const VectorField& v, const Vector& s ) { return VectorField(new TranslateVectorVolume(v,s)); }
ColorField translate( const ColorField& v, const Vector& s ) { return ColorField(new TranslateColor(v,s)); }
//MatrixField translate( const MatrixField& v, const Vector& s );
FormField translate( const FormField& v, const Vector& s ) { return FormField(new TranslateFormVolume(v,s)); }

ScalarField rotate( const ScalarField& v, const Vector& s ) { return ScalarField(new RotateVolume(v,s)); }
VectorField rotate( const VectorField& v, const Vector& s ) { return VectorField(new RotateVectorVolume(v,s)); }
ColorField rotate( const ColorField& v, const Vector& s ) { return ColorField(new RotateColor(v,s)); }
//MatrixField rotate( const MatrixField& v, const Vector& s );
FormField rotate( const FormField& v, const Vector& s ) { return FormField(new RotateFormVolume(v,s)); }

MatrixField rotation( const VectorField& v ) { return MatrixField( new RotatorMatrixVolume(v) ); }




ScalarField exp( const ScalarField& v ) { return ScalarField(new ExpVolume(v)); }
MatrixField exp( const MatrixField& m ) { return MatrixField( new ExpMatrixVolume( m ) ); }
MatrixField sinch( const MatrixField& m ) { return MatrixField( new SinchMatrixVolume( m ) ); }
MatrixField orderedSinch( const MatrixField& v1, const MatrixField& v2 ) { return MatrixField( new OrderedSinchMatrixVolume( v1, v2 ) ); }

ScalarField report( const ScalarField& v, const string& tag ) { return ScalarField(new ReportVolume(v,tag)); } 
VectorField report( const VectorField& v, const string& tag ) { return VectorField(new ReportVectorVolume(v,tag)); }
ColorField report( const ColorField& v, const string& tag )   { return ColorField(new ReportColorVolume(v,tag)); }
//MatrixField Report( const MatrixField& v );
FormField report( const FormField& v, const string& tag )   { return FormField(new ReportFormVolume(v,tag)); }

ScalarField negate( const ScalarField& v ) { return ScalarField(new NegateVolume(v)); }
VectorField negate( const VectorField& v ) { return VectorField(new NegateVectorVolume(v)); }
ColorField  negate( const ColorField& v )  { return ColorField(new NegateColor(v)); }
MatrixField negate( const MatrixField& v ) { return MatrixField( new NegateMatrixVolume(v) ); }
FormField   negate( const FormField& v )   { return FormField( new NegateFormVolume(v) ); }

ScalarField abs( const ScalarField& v ) { return ScalarField(new AbsoluteVolume(v)); }
ScalarField abs( const VectorField& v ) { return ScalarField(new MagnitudeVectorVolume(v)); }

ScalarField multiply( const ScalarField& v, const float a ) { return ScalarField(new MultiplyVolume(v,a)); }
ScalarField multiply( const ScalarField& v, const ScalarField& u ) { return ScalarField(new MultiplyVolume(v,u)); }
VectorField multiply( const VectorField& v, const float a ) { return VectorField(new MultiplyVectorVolume(v,a)); }
VectorField multiply( const VectorField& v, const ScalarField& u ) { return VectorField(new MultiplyVectorVolume(v,u)); }
ColorField multiply( const ColorField& v, const float a ) { return ColorField(new FloatMultiplyColor(v,a)); }
ColorField multiply( const ColorField& v, const ScalarField& u ) { return ColorField(new VolumeMultiplyColor(v,u)); }
ColorField multiply( const ColorField& v, const ColorField& u ) { return ColorField(new MultiplyColor(v,u)); }
MatrixField multiply( const MatrixField& v, const float a )  { return MatrixField( new FloatMatrixProductVolume(v,a)); }
MatrixField multiply( const MatrixField& v, const ScalarField& u ) { return MatrixField( new ScalarMatrixProductVolume(v,u) ); }
FormField multiply( const FormField& v, const float a )  { return FormField( new MultiplyFormVolume(v,a)); }
FormField multiply( const FormField& v, const ScalarField& u ) { return FormField( new MultiplyFormVolume(v,u) ); }

ScalarField divide( const ScalarField& v, const float a ) { return ScalarField(new DivideVolume(v,a)); }
ScalarField divide( const ScalarField& v, const ScalarField& u ) { return ScalarField(new DivideVolume(v,u)); }
VectorField divide( const VectorField& v, const float a ) { return VectorField(new DivideVectorVolume(v,a)); }
VectorField divide( const VectorField& v, const ScalarField& u ) { return VectorField(new DivideVectorVolume(v,u)); }
ColorField divide( const ColorField& v, const float a ) { return ColorField(new FloatDivideColor(v,a)); }
ColorField divide( const ColorField& v, const ScalarField& u ) { return ColorField(new VolumeDivideColor(v,u)); }
MatrixField divide( const MatrixField& v, const float a ) { return MatrixField(new FloatMatrixDivideVolume(v,a)); }
MatrixField divide( const MatrixField& v, const ScalarField& u ) { return MatrixField( new ScalarMatrixDivideVolume( v, u ) ); }
FormField divide( const FormField& v, const float a ) { return FormField(new DivideFormVolume(v,a)); }
FormField divide( const FormField& v, const ScalarField& u ) { return FormField(new DivideFormVolume(v,u)); }

ScalarField add( const ScalarField&  v1, const ScalarField& v2 ) { return ScalarField(new AddVolume(v1,v2));}
VectorField add( const VectorField&  v1, const VectorField& v2 ) { return VectorField(new AddVectorVolume(v1,v2)); }
ColorField  add( const ColorField&  v1, const ColorField& v2 )   { return ColorField(new AddColor(v1,v2)); }
MatrixField add( const MatrixField&  v1, const MatrixField& v2 ) { return MatrixField( new AddMatrixVolume( v1, v2 )); }
FormField   add( const FormField&  v1, const FormField& v2 )     { return FormField( new AddFormVolume( v1, v2 )); }

ScalarField subtract( const ScalarField&  v1, const ScalarField& v2 ) { return ScalarField(new SubtractVolume(v1,v2)); }
VectorField subtract( const VectorField&  v1, const VectorField& v2 ) { return VectorField(new SubtractVectorVolume(v1,v2)); }
ColorField  subtract( const ColorField&  v1, const ColorField& v2 )   { return ColorField(new SubtractColor(v1,v2)); }
MatrixField subtract( const MatrixField&  v1, const MatrixField& v2 ) { return MatrixField( new SubtractMatrixVolume(v1,v2)); }
FormField   subtract( const FormField&  v1, const FormField& v2 )     { return FormField( new SubtractFormVolume(v1,v2)); }

ScalarField Sphere( const Vector& cen, const float rad ) { return ScalarField( new SphereVolume(cen, rad) ); }
ScalarField Ellipse( const Vector& cen, const Vector& axs, const float majorrad, const float minorrad )
             { return ScalarField(new EllipseVolume(cen, axs, majorrad, minorrad)); }
ScalarField CsgBox( const Vector& cen, const float rad, const float pwr )
             { return ScalarField(new CsgBoxVolume(cen, rad, pwr)); }
ScalarField CsgRectangleBox( const Vector& cen, const float rad, const Vector& asp, const float pwr )
             {  return ScalarField(new CsgRectangleVolume(cen,rad,asp,pwr));    }
ScalarField HardBox( const Vector _llc, const Vector _urc )
             { return ScalarField(new  HardBoxVolume(_llc,_urc)); }
ScalarField Cone( const Vector& cen, const Vector& ax, const float h, const float theta )
             {  return ScalarField(new ConeVolume(cen,ax,h,theta));  }
ScalarField Plane( const Vector cen, const Vector norm )
             {  return ScalarField(new PlaneVolume(cen,norm)); }
ScalarField Torus( const Vector& cen, const Vector& axis, const float majorRad, const float minorRad )
             { return ScalarField(new TorusVolume(cen,axis,majorRad,minorRad)); }
ScalarField SteinerPatch() { return ScalarField(new SteinerPatchVolume()); }
ScalarField Icosahedron() { return ScalarField(new IcosahedronVolume()); }
ScalarField Cylinder( const Vector axis, const float rad ) { return ScalarField(new InfiniteCylinder(axis,rad)); }
ScalarField CappedCylinder( const Vector cen, const Vector axis, const float length, const float radius )
             { return ScalarField(new ImplicitCylinder(cen,axis,length,radius)); }
ScalarField Shell( const ScalarField& v, const float thickness )
             { return ScalarField(new ImplicitShell(v,thickness)); }

ScalarField SpaceCurveField( const SpaceCurve& curve, float radius )
             { return ScalarField( new SpaceCurveVolume( curve, radius ) ); } 
ScalarField PiecewiseCurveField( const AnchorChain& list )
             { return ScalarField( new PiecewiseCurveVolume( list ) ); }
ScalarField PiecewisePyroCurveField( const AnchorChain& list )
             { return ScalarField( new PiecewisePyroCurveVolume( list ) ); }
ScalarField PiecewiseNoiseCurveField( const AnchorChain& list )
             { return ScalarField( new PiecewiseNoiseCurveVolume( list ) ); }



ScalarField mask( const ScalarField& v ) { return ScalarField(new MaskVolume(v)); }
ScalarField clamp( const ScalarField& v, float minv, float maxv ) { return ScalarField(new ClampVolume(v, minv, maxv));}
ScalarField pow( const ScalarField& v, float gam ) { return ScalarField(new GammaVolume(v,gam)); }
ScalarField pow( const ScalarField& v, const ScalarField& gam ) { return ScalarField(new VolumeGammaVolume(v,gam)); }
ColorField pow( const ColorField& v, float gam ) { return ColorField(new FloatGammaColor(v,gam)); }
ColorField pow( const ColorField& v, const ScalarField& gam ) { return ColorField(new VolumeGammaColor(v,gam)); }
ScalarField BlinnBlend( const ScalarField& v1, const ScalarField& v2, const float _alpha )
             { return ScalarField(new BlinnBlendVolume(v1,v2,_alpha)); }
ScalarField MultiBlend( std::vector<ScalarField>& vs, const float a )
             { return ScalarField(new MultiBlendVolume(vs,a) ); }
ScalarField Union( const ScalarField& v1, const ScalarField& v2 ) { return ScalarField(new UnionVolume(v1,v2)); }
ScalarField intersection( const ScalarField& v1, const ScalarField& v2 ) { return ScalarField(new IntersectionVolume(v1,v2)); }
ScalarField cutout( const ScalarField& v1, const ScalarField& v2 ) { return ScalarField(new CutoutVolume(v1,v2)); }
ScalarField boxed( const ScalarField& v1, const Vector& llc, const Vector& urc, const float defvalue ){ return  ScalarField(new BoxedVolume(v1,llc,urc,defvalue));  }


ScalarField Pyroclast( const Vector& Center, const float Radius, const float Amp, 
                         const float octaves, const float freq, const float rough, 
                         const Vector trans, const float time, const float Gamma )
             { return ScalarField(new PyroclasticVolume(Center,Radius,Amp,octaves,freq,rough,trans,time,Gamma)); }

ScalarField RadialPyroclast( const Vector& Center, const float Radius, const float Amp, 
                             const float octaves, const float freq, const float rough, 
                             const float trans, const float time, const float Gamma )
             { return ScalarField(new RadialPyroclasticVolume(Center,Radius,Amp,octaves,freq,rough,trans,time,Gamma)); }

ScalarField SFFFTNoise( const float power, const float low, const float high, const float length, const int sz )
             { return ScalarField(new FFTNoiseVolume(power,low,high,length,sz)); }
ScalarField SFNoise( NoiseMachine n, const float d ) { return ScalarField(new NoiseVolume(n,d)); }
VectorField VFNoise( NoiseMachine n, const float d ) { return VectorField(new NoiseVectorVolume(n,d) ); }
ScalarField JitterSample( const ScalarField& v, const float rad, const int nb, const int seed )
             { return ScalarField(new JitterSampleVolume(v,rad,nb,seed) ); } 
//VectorField JitterSample( const VectorField& v, const float rad, const int nb, const int seed )
//             { return VectorField(new JitterSampleVectorVolume(v,rad,nb,seed) ); } 
//ColorField JitterSample( const ColorField& v, const float rad, const int nb, const int seed ); 
//             { return ColorField(new JitterSampleColor(v,rad,nb,seed) ); } 
//MatrixField JitterSample( const MatrixField& v, const float rad, const int nb, const int seed ); 

ScalarField gridded( const VolumeGrid<float>* g )
              { return ScalarField(new GriddedVolume(g)); }
ScalarField gridded( const SparseGrid* g )
              { return ScalarField(new GriddedVolume(g)); }
ScalarField gridded( const ScalarGrid& g )
              { return ScalarField(new GriddedSGridVolume(g)); }
ScalarField gridded( const ScalarFrustumGrid& g )
              { return ScalarField(new GriddedFrustumVolume(g)); }
VectorField gridded( const VolumeGrid<Vector>* g )
              { return VectorField(new GriddedVectorVolume(g)); }
VectorField gridded( const SparseVectorGrid* g )
              { return VectorField(new GriddedVectorVolume(g)); }
VectorField gridded( const VectorGrid& g )
              { return VectorField(new GriddedSGridVectorVolume(g)); }
VectorField gridded( const VectorFrustumGrid& g )
              { return VectorField(new GriddedFrustumVectorVolume(g)); }
ColorField gridded( const VolumeGrid<Color>* g )
              { return ColorField(new GriddedColor(g)); }
ColorField gridded( const SparseColorGrid* g )
              { return ColorField(new GriddedColor(g)); }
ColorField gridded( const ColorGrid& g )
              { return ColorField(new GriddedSGridColor(g)); }
ColorField gridded( const ColorFrustumGrid& g )
              { return ColorField(new GriddedFrustumColor(g)); }


MatrixField gridded( const MatrixGrid& g ) { return MatrixField( new GriddedMatrixVolume(g) ); }

ScalarField advect( const ScalarField& v, const VectorField& u, const float delt )
             { return ScalarField(new AdvectVolume(v,u,delt)); }
VectorField advect( const VectorField& v, const VectorField& u, const float delt ) 
             { return VectorField(new AdvectVectorVolume(v,u,delt)); }
ColorField advect( const ColorField& v, const VectorField& u, const float delt )
             { return ColorField(new AdvectColorVolume(v,u,delt)); }
MatrixField advect( const MatrixField& v, VectorField& u, const float delt ) 
             { return MatrixField(new AdvectMatrixVolume(v,u,delt)); }
FormField advect( const FormField& v, VectorField& u, const float delt ) 
             { return FormField(new AdvectFormVolume(v,u,delt)); }

VectorField gradientStretchCM( const VectorField& v, const float T, const int nb )
             { return VectorField( new GradientStretchCMVolume( v, T, nb ) ); }


VectorField gradientDisplacement( const MatrixField& M, const ScalarField& LS, const double step )
             { return VectorField( new GradDisplaceVectorVolume( M, LS, step ) ); }


ScalarField warp( const ScalarField& v, VectorField& map ) { return ScalarField(new WarpVolume(v,map)); }
VectorField warp( const VectorField& v, VectorField& map ) { return VectorField(new WarpVectorVolume(v,map)); }
ColorField  warp( const ColorField& v, VectorField& map )  { return ColorField(new WarpColorVolume(v,map)); }
MatrixField warp( const MatrixField& v, VectorField& map ) { return MatrixField(new WarpMatrixVolume(v,map)); }
FormField   warp( const FormField& v, VectorField& map )   { return FormField(new WarpFormVolume(v,map)); }

ScalarField Periodic( const ScalarField v, const Vector& o, const Vector& L ) { return ScalarField(new PeriodicVolume(v,o,L)); }
VectorField Periodic( const VectorField v, const Vector& o, const Vector& L ) { return VectorField(new PeriodicVectorVolume(v,o,L)); }
ColorField Periodic( const ColorField v, const Vector& o, const Vector& L ) { return ColorField(new PeriodicColorVolume(v,o,L)); }
//MatrixField Periodic( const MatrixField v, const Vector& o, const Vector& L );
FormField Periodic( const FormField v, const Vector& o, const Vector& L ) { return FormField(new PeriodicFormVolume(v,o,L)); }


FormField wedge( const FormField& f1, const FormField& f2 ) { return FormField(new WedgeFormVolume(f1,f2)); }


FormField contraction( const VectorField& v, const FormField f ) { return FormField( new ContractionFormVolume( v, f ) ); }

ScalarField dot( const VectorField& e1, const VectorField& e2 )
              { return ScalarField(new DotProductVectorVolume(e1,e2)); }
VectorField unitvector( const VectorField& e ) { return VectorField(new UnitVectorVolume(e)); }
VectorField identity() { return VectorField(new IdentityVectorVolume()); }
ScalarField xIdentity() { return ScalarField(new XIdentityVolume()); }
ScalarField yIdentity() { return ScalarField(new YIdentityVolume()); }
ScalarField zIdentity() { return ScalarField(new ZIdentityVolume()); }
VectorField ImplicitSurfacePoint( const ScalarField& e, const float step, const int nbIterations )
              { return VectorField(new ImplicitPointVectorVolume(e,step,nbIterations)); }
VectorField cross( const VectorField& e1, const VectorField& e2 )
              { return VectorField(new CrossProductVectorVolume(e1,e2)); }
ScalarField div( const VectorField& e ) { return ScalarField(new DivergenceVectorVolume(e)); }

ScalarField fddiv( const VectorField& e, const int N, const float dx, const float dy, const float dz )
     { return ScalarField(new FiniteDifferenceDivergenceVectorVolume(e,N,dx,dy,dz)); }

ScalarField fdboundeddiv( const VectorField& e, const int N, const GridBox& gb )
     { return ScalarField(new FiniteDifferenceBoundedDivergenceVectorVolume(e,N,gb)); }

ScalarField fdinteriordiv( const VectorField& e, const int N, const float dx, const float dy, const float dz, const ScalarField& gb )
     { return ScalarField(new FiniteDifferenceInteriorDivergenceVectorVolume(e,N,dx,dy,dz,gb)); }

VectorField curl( const VectorField& e ) { return VectorField(new CurlVectorVolume(e)); }
VectorField ContinuedFractionDisplacement( const VectorField& dX, const int nbIterations )
               { return VectorField(new ContinuedFractionDisplacementVectorVolume(dX,nbIterations)); }
VectorField XYZ( const ScalarField& x, const ScalarField& y, const ScalarField& z )
            { return VectorField(new ComponentVectorVolume(x,y,z)); }

ColorField Chroma( const ColorField& e ) { return ColorField(new ChromaColor(e)); }
ColorField Blackbody( const ScalarField& e ) { return ColorField(new BlackBodyVolume(e)); }
ColorField RGB( const ScalarField& red, const ScalarField& green, const ScalarField& blue )
            { return ColorField(new ComponentColor(red,green,blue)); }
ColorField LUTColor( const ScalarField& field, const vector<Color> lut, const float maxValue, const float minValue )
            { return ColorField(new VolumeLUTColor(field,lut,maxValue,minValue)); }




ScalarField rComponent( const ColorField& C  ) { return ScalarField( new RColorVolume(C)  );  }
ScalarField gComponent( const ColorField& C  ) { return ScalarField( new GColorVolume(C)  );  }
ScalarField bComponent( const ColorField& C  ) { return ScalarField( new BColorVolume(C)  );  }

ColorField rgbComponent( const ScalarField& red, const ScalarField& green, const ScalarField& blue )
   {  return ColorField( new ComponentColor( red, green, blue )  ); }

ScalarField xComponent( const VectorField& X ) { return ScalarField( new XVectorVolume(X) );  }
ScalarField yComponent( const VectorField& X ) { return ScalarField( new YVectorVolume(X) );  }
ScalarField zComponent( const VectorField& X ) { return ScalarField( new ZVectorVolume(X) );  }


VectorField component( const ScalarField& X, const ScalarField& Y, const ScalarField& Z )
   {  return VectorField( new ComponentVectorVolume( X, Y, Z )  ); }





ScalarField DetGrad( const VectorField& e )
{
   return ScalarField( new DetGradMapVolume(e) );
}

MatrixField inverse( const MatrixField& m )
{
   return MatrixField( new InverseMatrixVolume(m) ); 
}



ScalarField zeroComponent( const FormField& f ) { return ScalarField(new ZeroFormVolume(f)); } 
VectorField oneComponent( const FormField& f )  { return VectorField(new OneFormVolume(f)); }
VectorField twoComponent( const FormField& f )  { return VectorField(new TwoFormVolume(f)); }
ScalarField threeComponent( const FormField& f ){ return ScalarField(new ThreeFormVolume(f)); }

FormField component( const ScalarField& X, const VectorField& Y, const VectorField& Z, const ScalarField& W )
     { return FormField(new ComponentFormVolume( X, Y, Z, W )); } 



FormField star( const FormField& e1 ) { return FormField(new StarFormVolume( e1 ) ); }


FormField lie( VectorField& X, FormField& a )
{
   return ( contraction(X,grad(a)) + grad( contraction(X,a) ) );
}

//VectorField bracket( VectorField& X, VectorField& Y )
//{
//   return ( (X*grad(Y))-(Y*grad(X))  );
//}







ScalarField sin( const ScalarField& v1 ) { return ScalarField( new SineVolume(v1)); }
ScalarField cos( const ScalarField& v1 ) { return ScalarField( new CosineVolume(v1)); }
ScalarField tan( const ScalarField& v1 ) { return ScalarField( new TangentVolume(v1)); }
ScalarField acos( const ScalarField& v1 ) { return ScalarField( new ArcsineVolume(v1)); }
ScalarField asin( const ScalarField& v1 ) { return ScalarField( new ArccosineVolume(v1)); }
ScalarField atan( const ScalarField& v1 ) { return ScalarField( new ArctangentVolume(v1)); }
ScalarField sinh( const ScalarField& v1 ) { return ScalarField( new HyperbolicSineVolume(v1)); }
ScalarField cosh( const ScalarField& v1 ) { return ScalarField( new HyperbolicCosineVolume(v1)); }
ScalarField tanh( const ScalarField& v1 ) { return ScalarField( new HyperbolicTangentVolume(v1)); }




float transmissivity( const ScalarField& density, const Vector& start, const Vector& end, const double step, const double scatter )
{
   double T = 0;
   Vector D = end-start;
   double distance = D.magnitude();
   D *= step/distance;

   int nsteps = (int)(distance/step);
   Vector X = start;
   for(size_t i=0;i<nsteps;i++)
   {
      T += density->eval(X);
      X += D;
   }
   double fraction = (distance - nsteps*step);
   if(fraction > 0)
   {
      X -= D*(1.0-fraction);
      T += density->eval(X)*fraction;
   }
   T = std::exp( -T*step*scatter );
   return (float)T;
}

Color lineIntegral( const ScalarField& density, const ColorField& color, const Vector& start, const Vector& end, const double step, const double scatter )
{
   double T = 1;
   Color Ci;
   Vector D = end-start;
   double distance = D.magnitude();
   D *= step/distance;
   int nsteps = (int)(distance/step);
   Vector X = start;
   for(size_t i=0;i<nsteps;i++)
   {
      double dT = std::exp( -step*scatter*density->eval(X) );
      Ci += color->eval(X) * (1.0-dT) * T;
      T *= dT;
      X += D;
   }
   double fraction = (distance - nsteps*step);
   if(fraction > 0)
   {
      X -= D*(1.0-fraction);
      double dT = std::exp( -step*fraction*scatter*density->eval(X) );
      Ci += color->eval(X) * (1.0-dT) * T;
      T *= dT;
   }
   Ci /= scatter;
   Ci[3] = 1.0 - T;
   return Ci;
}









void fieldStatistics( const ScalarField& f, const GridBox& g )
{
   double mean = 0.0, stddev = 0.0;
   double max = f->eval( g->evalP(0,0,0) );
   double min = max;
   for( int k=0;k<g->nz();k++ )
   {
      for( int j=0;j<g->ny();j++ )
      {
         for( int i=0;i<g->nx();i++ )
         {
            float data = f->eval( g->evalP(i,j,k) );
            mean += data;
            stddev += data*data;
            if (max < data) { max = data; }
            if (min > data) { min = data; }
         }
      }
   }
   mean /= (g->nx())*(g->ny())*(g->nz());
   stddev /= (g->nx())*(g->ny())*(g->nz());
   stddev -= mean*mean;
   cout << "\tMean:   " << mean << endl;
   cout << "\tStddev: " << stddev << endl;
   cout << "\tMax:    " << max << endl;
   cout << "\tMin:    " << min << endl;
}

void fieldStatistics( const VectorField& f, const GridBox& g )
{
   Vector mean;
   double stddev = 0.0;
   Vector max = f->eval( g->evalP(0,0,0) );
   Vector min = max;
   for( int k=0;k<g->nz();k++ )
   {
      for( int j=0;j<g->ny();j++ )
      {
         for( int i=0;i<g->nx();i++ )
         {
            Vector data = f->eval( g->evalP(i,j,k) );
            mean += data;
            stddev += data*data;
            if (max < data) { max = data; }
            if (min > data) { min = data; }
         }
      }
   }
   mean /= (g->nx())*(g->ny())*(g->nz());
   stddev /= (g->nx())*(g->ny())*(g->nz());
   stddev -= mean*mean;
   cout << "\tMean:   " << mean.X() << " " << mean.Y() << " " << mean.Z() << endl;
   cout << "\tStddev: " << stddev << endl;
   cout << "\tMax:    " << max.X() << " " << max.Y() << " " << max.Z() << endl;
   cout << "\tMin:    " << min.X() << " " << min.Y() << " " << min.Z() << endl;
}






void setFDSize( ScalarField& f, int n ){ f->setFDSize(n); }
void setFDSize( VectorField& f, int n ){ f->setFDSize(n); }
void setFDSize( MatrixField& f, int n ){ f->setFDSize(n); }
void setFDSize( ColorField& f, int n ){ f->setFDSize(n); }
void setFDSize( FormField& f, int n ){ f->setFDSize(n); }

void setFDStep( ScalarField& f, double dx ){ f->setFDStep(dx); }
void setFDStep( VectorField& f, double dx ){ f->setFDStep(dx); }
void setFDStep( MatrixField& f, double dx ){ f->setFDStep(dx); }
void setFDStep( ColorField& f, double dx ){ f->setFDStep(dx); }
void setFDStep( FormField& f, double dx ){ f->setFDStep(dx); }

void setFDStep( ScalarField& f, double dx, double dy, double dz ){ f->setFDStep(dx,dy,dz); }
void setFDStep( VectorField& f, double dx, double dy, double dz ){ f->setFDStep(dx,dy,dz); }
void setFDStep( MatrixField& f, double dx, double dy, double dz ){ f->setFDStep(dx,dy,dz); }
void setFDStep( ColorField& f, double dx, double dy, double dz ){ f->setFDStep(dx,dy,dz); }
void setFDStep( FormField& f, double dx, double dy, double dz ){ f->setFDStep(dx,dy,dz); }









}
