
#include "ImplicitVolumeShapes.h"
#include <iostream>
using namespace std;

using namespace lux;




ReportVolume::ReportVolume( Volume<float>* v, const string tag ) :
  elem (v),
  label (tag)
{}

ReportVolume::ReportVolume( const ScalarField& v, const string tag ) :
  elem (v),
  label (tag)
{}



const float ReportVolume::eval( const Vector& P ) const 
{
   float value =  elem->eval(P);
   cout << label << " eval: " << value << " at position " << P[0] << " " << P[1] << " " << P[2] << endl; 
   return value; 
}

const Vector ReportVolume::grad( const Vector& P ) const 
{ 
   Vector value =  elem->grad(P);
   cout << label << " grad: " << value[0] << " " << value[1] << " " << value[2] << endl; 
   return value; 
}








PyroclasticVolume::PyroclasticVolume( const Vector& Center, const float Radius, const float Amp, 
                                      const float octaves, const float freq, const float rough, const Vector trans, const float time, const float Gamma  ) :
   amplitude (Amp),
   center    (Center),
   gamma     (Gamma),
   elem      ( new SphereVolume( Center, Radius ) )
{
   Noise_t parms;
   parms.octaves = octaves;
   parms.frequency = freq;
   parms.roughness = rough;
   parms.translate = trans;
   parms.time = time;
   noise.setParameters( parms );
   dx = 0.1/( freq * pow( (double)parms.fjump, (double)octaves) );
}


const float PyroclasticVolume::eval( const Vector& P ) const
{
   Vector X = (P-center).unitvector();
   float value = elem->eval(P) + pow(fabs( noise.eval(X) ), gamma )*amplitude;
   return value;
}

/*
const Vector PyroclasticVolume::grad( const Vector& P ) const
{
   float g1 = (eval(P+Vector(dx,0,0))-eval(P-Vector(dx,0,0)))/(2.0*dx);
   float g2 = (eval(P+Vector(0,dx,0))-eval(P-Vector(0,dx,0)))/(2.0*dx);
   float g3 = (eval(P+Vector(0,0,dx))-eval(P-Vector(0,0,dx)))/(2.0*dx);
   return Vector(g1,g2,g3);
}
*/




FFTNoiseVolume::FFTNoiseVolume( const float power, const float low, const float high, const float length, const int sz )
{
   dx = length/sz;
   Noise_t parms;
   parms.fftPower = power;
   parms.fftLowCutoff = low;
   parms.fftHighCutoff = high;
   parms.fftLength = length;
   noise = new FFTNoise();
   noise->setParameters( parms );
}


const float FFTNoiseVolume::eval( const Vector& P ) const { return noise->eval(P); }

/*
const Vector FFTNoiseVolume::grad( const Vector& P ) const
{
   Vector g( noise->eval( P+Vector(dx,0,0) )-noise->eval( P-Vector(dx,0,0)), noise->eval(P+Vector(0,dx,0))-noise->eval(P-Vector(0,dx,0)), noise->eval(P+Vector(0,0,dx))-noise->eval(P-Vector(0,0,dx)) );
   g /= 2.0*dx;
   return g;
}
*/




RadialPyroclasticVolume::RadialPyroclasticVolume( const Vector& Center, const float Radius, const float Amp, 
                                      const float octaves, const float freq, const float rough, const float trans, const float time, const float Gamma  ) :
   amplitude (Amp),
   center    (Center),
   gamma     (Gamma),
   radius    (Radius),
   translate (trans)
{
   parms.octaves = octaves;
   parms.frequency = freq;
   parms.roughness = rough;
   parms.translate = Vector(trans,trans,trans);
   parms.time = time;
   noise.setParameters( parms );
   dx = 0.1/( freq * pow( (double)parms.fjump, (double)octaves) );
   gradParams.setStep( dx );
}


const float RadialPyroclasticVolume::eval( const Vector& P ) const
{
   Vector X = (P-center);
   float rad = X.magnitude();
   X.normalize();
   parms.translate = X * translate;
   noise.setParameters( parms );
   float value = radius - rad + pow(fabs( noise.eval(X) ), gamma )*amplitude;
   return value;
}

/*
const Vector RadialPyroclasticVolume::grad( const Vector& P ) const
{
   float e0 = eval(P);
   float g1 = (eval(P+Vector(dx,0,0))-e0)/dx;
   float g2 = (eval(P+Vector(0,dx,0))-e0)/dx;
   float g3 = (eval(P+Vector(0,0,dx))-e0)/dx;
   return Vector(g1,g2,g3);
}
*/

ConstantVolume::ConstantVolume( const float v ) :
       value (v),
       gradvalue (Vector(0,0,0))
    {}


const float ConstantVolume::eval( const Vector& P ) const 
{
   return value; 
}

const Vector ConstantVolume::grad( const Vector& P ) const { return gradvalue; }


ExpVolume::ExpVolume( Volume<float>* v ) :
      elem (v)
    {}

ExpVolume::ExpVolume( const ScalarField& v ) :
      elem (v)
    {}

const float ExpVolume::eval( const Vector& P ) const 
{
   return std::exp( elem->eval(P) ); 
}

const Vector ExpVolume::grad( const Vector& P ) const { return eval(P) * elem->grad(P); }




SphereVolume::SphereVolume( const Vector& cen, const float rad ) :
       center (cen),
       radius (rad)
    {}


const float SphereVolume::eval( const Vector& P ) const 
{
   Vector x = P-center;
   float rad = x.magnitude();
   return radius-rad;
}

const Vector SphereVolume::grad( const Vector& P ) const 
{
   Vector x = P-center;
   if( x.magnitude() != 0 ){ return -x.unitvector(); }
   return Vector(0,1,0);
}





CsgBoxVolume::CsgBoxVolume( const Vector& cen, const float rad, const float pwr ) :
       center (cen),
       radius (pow((double)rad,(double)(2.0*pwr))),
       power  (2.0*pwr)
    {}



const float CsgBoxVolume::eval( const Vector& P ) const 
{
   Vector x = P-center;
   float rad = pow( fabs(x[0]), power ) + pow( fabs(x[1]), power) + pow( fabs(x[2]), power );
   return radius-rad;
}

const Vector CsgBoxVolume::grad( const Vector& P ) const 
{
   Vector x = P-center;
   return -Vector( pow( fabs(x[0]), power-1.0), pow( fabs(x[1]), power-1.0), pow( fabs(x[2]), power-1.0 )  )*power;
}



CsgRectangleVolume::CsgRectangleVolume( const Vector& cen, const float rad, const Vector& asp, const float pwr ) :
       center (cen),
       radius (pow((double)rad,(double)pwr)),
       aspect (asp),
       power  (2.0*pwr)
    {}

const float CsgRectangleVolume::eval( const Vector& P ) const 
{
   Vector x = P-center;
   float rad = pow( fabs(x[0]/aspect[0]), power ) + pow( fabs(x[1]/aspect[1]), power) + pow( fabs(x[2]/aspect[2]), power );
   return radius-rad;
}

const Vector CsgRectangleVolume::grad( const Vector& P ) const 
{
   Vector x = P-center;
   return -Vector( pow( fabs(x[0]), power-1.0), pow( fabs(x[1]), power-1.0), pow( fabs(x[2]), power-1.0 )  )*power;
}







ConeVolume::ConeVolume( const Vector& cen, const Vector& ax, const float h, const float theta ) :
       center (cen),
       axis   (ax.unitvector()),
       height (h),
       angle  (theta*M_PI/180.0)
    {}

const float ConeVolume::eval( const Vector& P ) const 
{
   if( P == center ){ return 0.0; }
   Vector x = P-center;
   float Y = x*axis;
   if( Y > height ){ return height-Y; }
   if( Y < 0 ){ return Y; }
   Vector test = x.unitvector();
   return angle - acos( x.unitvector() * axis );
}

const Vector ConeVolume::grad( const Vector& P ) const 
{
   Vector x = P-center;
   Vector xhat = x.unitvector();
   float Y = x*axis;
   Vector value;
   if( Y > height ){ value = Vector( 0, -1, 0); }
   else if( Y < 0 ){ value = Vector(0,1,0); }
   else { value = xhat/sqrt( 1.0 - pow( xhat * axis, 2.0 ) ); }
   return value;
}







PlaneVolume::PlaneVolume( const Vector cen, const Vector norm ) :
       center (cen),
       normal (norm.unitvector())
    {}



const float PlaneVolume::eval( const Vector& P ) const 
{
   Vector x = P-center;
   return x*normal;
}

const Vector PlaneVolume::grad( const Vector& P ) const 
{
   return normal;
}



HardBoxVolume::HardBoxVolume( const Vector _llc, const Vector _urc ) :
    llc(_llc),
    urc(_urc)
    {}




const float HardBoxVolume::eval( const Vector& P ) const 
{
   // check that it is inside
   Vector X = P-llc;
   Vector XX = urc-P;
   float val = X[0];
   val = ( val > X[1] ) ? X[1] : val;
   val = ( val > X[2] ) ? X[2] : val;
   val = ( val > XX[0] ) ? XX[0] : val;
   val = ( val > XX[1] ) ? XX[1] : val;
   val = ( val > XX[2] ) ? XX[2] : val;
   return val;
}




TorusVolume::TorusVolume( const Vector& cen, const Vector& axis, const float majorRad, const float minorRad ) :
       tcenter (cen),
       taxis (axis),
       majorRadiusSquare (majorRad*majorRad),
       minorRadiusSquare (minorRad*minorRad)
    {
       taxis.normalize();
    }

const float TorusVolume::eval( const Vector& P ) const 
{
   Vector x = P-tcenter;
   float y = x*taxis;
   Vector xperp = x - y*taxis;
   float rperpsquare = xperp*xperp;
   float xsquare = x*x;
   return 4.0*majorRadiusSquare*rperpsquare - pow( ( xsquare + majorRadiusSquare - minorRadiusSquare ), 2 );
}

const Vector TorusVolume::grad( const Vector& P ) const 
{
   Vector x = P-tcenter;
   float y = x*taxis;
   Vector xperp = x - y*taxis;
   float xsquare = x*x;
   Vector g = 8.0*majorRadiusSquare*xperp - 4.0*x*( xsquare + majorRadiusSquare - minorRadiusSquare );
   return g;
}




MobiusStripVolume::MobiusStripVolume( const Vector& cen, const Vector& axis, const float rad, const float thick ) :
       tcenter (cen),
       taxis (axis),
       radius (rad),
       thickness (thick)
    {
       taxis.normalize();
    }



const float MobiusStripVolume::eval( const Vector& P ) const 
{
   Vector x = P-tcenter;
   float z = x*taxis - thickness;
   Vector xperp = x - z*taxis;
   float r = xperp.magnitude()/radius;
   float theta = asin((xperp^taxis).magnitude())/2.0;
   return  -log(r) * sin( theta ) + z*cos(theta);
}







SteinerPatchVolume::SteinerPatchVolume(){}

const float SteinerPatchVolume::eval( const Vector& P ) const 
{
   float xx = P[0]*P[0];
   float yy = P[1]*P[1];
   float zz = P[2]*P[2];
   return -( xx*yy + xx*zz + yy*zz - P[0]*P[1]*P[2] );
}

const Vector SteinerPatchVolume::grad( const Vector& P ) const 
{
   float xx = P[0]*P[0];
   float yy = P[1]*P[1];
   float zz = P[2]*P[2];
   Vector g( -2.0*P[0]*(yy+zz) + P[1]*P[2], -2.0*P[1]*(xx+zz) + P[0]*P[2], -2.0*P[2]*(xx+yy) + P[0]*P[1]  );
   return g;
}


IcosahedronVolume::IcosahedronVolume(){}

const float IcosahedronVolume::eval( const Vector& P ) const 
{
   float threshold = 1.8 * M_PI;
   float x = P[0];
   float y = P[1];
   float z = P[2];
   float T = 1.61803399; // golden ratio
   float rad = P.magnitude();
   if( rad > threshold ){ return -threshold; }
   return (cos(x + T*y) + cos(x - T*y) + cos(y + T*z) + cos(y - T*z) + cos(z - T*x) + cos(z + T*x)) - 2.0;
}




ScaleVolume::ScaleVolume( Volume<float> * v, const Vector& s ) :
       elem (v),
       scale (s)
    {}

ScaleVolume::ScaleVolume( Volume<float> * v, const float& s ) :
       elem (v),
       scale (Vector(s,s,s))
    {}


ScaleVolume::ScaleVolume( const ScalarField& v, const Vector& s ) :
       elem (v),
       scale (s)
    {}

ScaleVolume::ScaleVolume( const ScalarField& v, const float& s ) :
       elem (v),
       scale (Vector(s,s,s))
    {}


const float ScaleVolume::eval( const Vector& P ) const
{
   Vector X = P;
   X[0] /= scale[0];
   X[1] /= scale[1];
   X[2] /= scale[2];
   return  (elem->eval(X));
}

const Vector ScaleVolume::grad( const Vector& P ) const
{
   Vector X = P;
   X[0] /= scale[0];
   X[1] /= scale[1];
   X[2] /= scale[2];
   Vector d =  (elem->grad(X));
   d[0] *= scale[0];
   d[1] *= scale[1];
   d[2] *= scale[2];
   return d;
}





TranslateVolume::TranslateVolume( Volume<float> * v, const Vector& s ) :
       elem (v),
       translate (s)
    {}

TranslateVolume::TranslateVolume( const ScalarField& v, const Vector& s ) :
       elem (v),
       translate (s)
    {}

const float TranslateVolume::eval( const Vector& P ) const
{
   Vector X = P - translate;
   return  (elem->eval(X));
}

const Vector TranslateVolume::grad( const Vector& P ) const
{
   Vector X = P - translate;
   return  (elem->grad(X));
}



RotateVolume::RotateVolume( Volume<float> * v, const Vector& s ) :
        elem (v),
	axis( s.unitvector()),
	sina (cos(s.magnitude()*M_PI/180.0)),
	cosa (sin(s.magnitude()*M_PI/180.0))
    {
       R = inverse(rotation( s.unitvector() , s.magnitude() ));
    }

RotateVolume::RotateVolume( const ScalarField& v, const Vector& s ) :
        elem (v),
	axis( s.unitvector()),
	sina (cos(s.magnitude()*M_PI/180.0)),
	cosa (sin(s.magnitude()*M_PI/180.0))
    {
       R = inverse(rotation( s.unitvector() , s.magnitude() ));
    }

const float RotateVolume::eval( const Vector& P ) const
{
   Vector X = P*cosa + (axis*P)*(1.0-cosa)*axis - (axis^P)*sina;
   return  (elem->eval(X));
}

const Vector RotateVolume::grad( const Vector& P ) const
{
   Vector X = P*cosa + (axis*P)*(1.0-cosa)*axis - (axis^P)*sina;
   Vector v = elem->grad(X);
   v = v*cosa + (axis*v)*(1.0-cosa)*axis + (axis^v)*sina;
   return v;
}




NegateVolume::NegateVolume( Volume<float> * v ) :
      elem(v)
    {}

NegateVolume::NegateVolume( const ScalarField& v ) :
      elem(v)
    {}


const float NegateVolume::eval( const Vector& P ) const
{
       return  -(elem->eval(P));
}

const Vector NegateVolume::grad( const Vector& P ) const
{
   return  -(elem->grad(P));
}



AbsoluteVolume::AbsoluteVolume( Volume<float> * v )  : elem (v) {}

AbsoluteVolume::AbsoluteVolume( const ScalarField& v )  : elem (v) {}

const float AbsoluteVolume::eval( const Vector& P ) const
{
   return  fabs(elem->eval(P));
}

const Vector AbsoluteVolume::grad( const Vector& P ) const
{
   float val = elem->eval(P);
   Vector g =  elem->grad(P);
   return  (val > 0 ) ? g : -g;
}







MultiplyVolume::MultiplyVolume( Volume<float> * v, const float a ) : 
       elem (v),
       factor( new ConstantVolume(1.0) ),
       constant (a) 
    {}

MultiplyVolume::MultiplyVolume( Volume<float> * v, Volume<float>* u ) : 
       elem (v),
       factor(u),
       constant (1.0) 
    {}

MultiplyVolume::MultiplyVolume( const ScalarField& v, const float a ) : 
       elem (v),
       factor( new ConstantVolume(1.0) ),
       constant (a) 
    {}

MultiplyVolume::MultiplyVolume( const ScalarField& v, const ScalarField& u ) : 
       elem (v),
       factor(u),
       constant (1.0) 
    {}

const float MultiplyVolume::eval( const Vector& P ) const
{
   return  elem->eval(P) * factor->eval(P) * constant;
}

const Vector MultiplyVolume::grad( const Vector& P ) const
{
   return elem->grad(P) * (factor->eval(P) * constant ) + factor->grad(P) * ( elem->eval(P) * constant);
}




DivideVolume::DivideVolume( Volume<float> * v, const float a ) : 
       elem (v),
       factor( new ConstantVolume(1.0) ),
       constant (a) 
    {}

DivideVolume::DivideVolume( Volume<float> * v, Volume<float>* u ) : 
       elem (v),
       factor(u),
       constant (1.0) 
    {}

DivideVolume::DivideVolume( const ScalarField& v, const float a ) : 
       elem (v),
       factor( new ConstantVolume(1.0) ),
       constant (a) 
    {}

DivideVolume::DivideVolume( const ScalarField& v, const ScalarField& u ) : 
       elem (v),
       factor(u),
       constant (1.0) 
    {}

const float DivideVolume::eval( const Vector& P ) const
{
   return  elem->eval(P) / (factor->eval(P) * constant);
}

const Vector DivideVolume::grad( const Vector& P ) const
{
   float elemv = elem->eval(P);
   float factv = factor->eval(P);
   return ( elem->grad(P) /factv - factor->grad(P)*elemv/(factv*factv) ) / constant;
}





AddVolume::AddVolume( Volume<float> * v1, Volume<float> * v2 ) :
      e1 (v1),
      e2 (v2)
    {}

AddVolume::AddVolume( const ScalarField&  v1, const ScalarField& v2 ) :
      e1 (v1),
      e2 (v2)
    {}

const float AddVolume::eval( const Vector& P ) const
{
   return  e1->eval(P) + e2->eval(P);
}

const Vector AddVolume::grad( const Vector& P ) const
{
   return  e1->grad(P) + e2->grad(P);
}



SubtractVolume::SubtractVolume( Volume<float> * v1, Volume<float> * v2 ) :
      e1 (v1),
      e2 (v2)
    {}

SubtractVolume::SubtractVolume( const ScalarField& v1, const ScalarField& v2 ) :
      e1 (v1),
      e2 (v2)
    {}

const float SubtractVolume::eval( const Vector& P ) const
{
   return  e1->eval(P) - e2->eval(P);
}

const Vector SubtractVolume::grad( const Vector& P ) const
{
   return  e1->grad(P) - e2->grad(P);
}




MaskVolume::MaskVolume( Volume<float> * v ) :
      elem (v)
    {}


MaskVolume::MaskVolume( const ScalarField& v ) :
      elem (v)
    {}

const float MaskVolume::eval( const Vector& P ) const
{
   const float e = elem->eval(P);
   if( e > 0 ){ return 1.0; }
   if( e <= 0 ){ return 0.0; }
   return 0.0;
}

const Vector MaskVolume::grad( const Vector& P ) const
{
   const float e = elem->eval(P);
   if( e == 0 ){ return elem->grad(P); }
   return Vector(0,0,0);
}



ClampVolume::ClampVolume( Volume<float> * v, float minv, float maxv ) :
       elem (v),
       vmin (minv),
       vmax (maxv)
    {}

ClampVolume::ClampVolume( const ScalarField& v, float minv, float maxv ) :
       elem (v),
       vmin (minv),
       vmax (maxv)
    {}

const float ClampVolume::eval( const Vector& P ) const
{
   const float e = elem->eval(P);
   if( e > vmax ){ return vmax; }
   if( e < vmin ){ return vmin; }
   return e;
}

const Vector ClampVolume::grad( const Vector& P ) const
{
   const float e = elem->eval(P);
   if( e > vmax ){ return Vector(0,0,0); }
   if( e < vmin ){ return Vector(0,0,0); }
   return elem->grad(P);
}



GammaVolume::GammaVolume( Volume<float> * v, float gam ) :
       elem (v),
       gamma (gam)
    {}

GammaVolume::GammaVolume( const ScalarField& v, float gam ) :
       elem (v),
       gamma (gam)
    {}

const float GammaVolume::eval( const Vector& P ) const
{
   const float e = elem->eval(P);
   return pow( (double)(e), (double)gamma );
}

const Vector GammaVolume::grad( const Vector& P ) const
{
   const float e = elem->eval(P);
   Vector g = elem->grad(P);
   g = g * ( gamma * pow( (double)e, (double)gamma - 1.0 ) );
   return g;
}




VolumeGammaVolume::VolumeGammaVolume( Volume<float> * v, Volume<float>* gam ) :
       elem (v),
       gamma (gam)
    {}

VolumeGammaVolume::VolumeGammaVolume( const ScalarField& v, const ScalarField& gam ) :
       elem (v),
       gamma (gam)
    {}

const float VolumeGammaVolume::eval( const Vector& P ) const
{
   const float e = elem->eval(P);
   const float g = gamma->eval(P);
   return pow( (double)(e), (double)g );
}

const Vector VolumeGammaVolume::grad( const Vector& P ) const
{
   const float e = elem->eval(P);
   Vector eg = elem->grad(P);
   const float g = gamma->eval(P);
   Vector gg = gamma->grad(P);
   eg = g * pow( (double)(e), (double)g - 1.0 ) * eg;
   eg += gg * pow( (double)(e), (double)g ) * log(e);
   return eg;
}







BlinnBlendVolume::BlinnBlendVolume( Volume<float> * v1, Volume<float> * v2, const float _alpha ) :
     elem1 (v1),
     elem2 (v2),
     alpha (_alpha)
    {}

BlinnBlendVolume::BlinnBlendVolume( const ScalarField& v1, const ScalarField& v2, const float _alpha ) :
     elem1 (v1),
     elem2 (v2),
     alpha (_alpha)
    {}

const float BlinnBlendVolume::eval( const Vector& P ) const
{
   const float e0 = elem1->eval(P);
   const float e1 = elem2->eval(P);
   return (std::exp(e0) + std::exp(e1) - 2.0*alpha );
}

const Vector BlinnBlendVolume::grad( const Vector& P ) const
{
   const float e0 = elem1->eval(P);
   const float e1 = elem2->eval(P);
   const Vector g0 = elem1->grad(P);
   const Vector g1 = elem2->grad(P);
   return (std::exp(e0)*g0 + std::exp(e1)*g1 );
}




MultiBlendVolume::MultiBlendVolume( std::vector<Volume<float>*>& vs, const float a  ) :
     alpha (a)
    {
       for( size_t i=0;i<vs.size();i++ )
       {
          const ScalarField e( vs[i] );
          elems.push_back(e);
       }
    }

MultiBlendVolume::MultiBlendVolume( std::vector<ScalarField>& vs, const float a  ) :
      alpha(a),
      elems (vs)
    {}

const float MultiBlendVolume::eval( const Vector& P ) const
{
   float sum = 0;
   for( size_t i=0;i<elems.size();i++ )
   {
          sum += std::exp( elems[i]->eval(P) );
   }
   sum -= (float)elems.size() * alpha;
   return sum;
}

const Vector MultiBlendVolume::grad( const Vector& P ) const
{
   Vector sum(0,0,0);
   for( size_t i=0;i<elems.size();i++ )
   {
      sum += ( elems[i]->grad(P) ) * std::exp( elems[i]->eval(P) );
   }
   return sum;
}

  

UnionVolume::UnionVolume( Volume<float> * v1, Volume<float> * v2 ) :
      elem1(v1),
      elem2(v2)
    {}

UnionVolume::UnionVolume( const ScalarField& v1, const ScalarField& v2 ) :
      elem1(v1),
      elem2(v2)
    {}

const float UnionVolume::eval( const Vector& P ) const
{
   const float e0 = elem1->eval(P);
   const float e1 = elem2->eval(P);
   return ( (e0 > e1) ? e0 : e1  );
}

const Vector UnionVolume::grad( const Vector& P ) const
{
   const float e0 = elem1->eval(P);
   const float e1 = elem2->eval(P);
   const Vector g0 = elem1->grad(P);
   const Vector g1 = elem2->grad(P);
   return ( (e0 > e1) ? g0 : g1  );
}

  




IntersectionVolume::IntersectionVolume( Volume<float> * v1, Volume<float> * v2 ) :
       elem1(v1),
       elem2(v2)
    {}

IntersectionVolume::IntersectionVolume( const ScalarField& v1, const ScalarField& v2 ) :
       elem1(v1),
       elem2(v2)
    {}

const float IntersectionVolume::eval( const Vector& P ) const
{
       const float e0 = elem1->eval(P);
       const float e1 = elem2->eval(P);
       return ( (e0 < e1) ? e0 : e1  );
}

const Vector IntersectionVolume::grad( const Vector& P ) const
{
   const float e0 = elem1->eval(P);
   const float e1 = elem2->eval(P);
   const Vector g0 = elem1->grad(P);
   const Vector g1 = elem2->grad(P);
   return ( (e0 < e1) ? g0 : g1  );
}

  




CutoutVolume::CutoutVolume( Volume<float> * v1, Volume<float> * v2 ):
      elem1(v1),
      elem2(v2)
    {}

CutoutVolume::CutoutVolume( const ScalarField& v1, const ScalarField& v2 ):
      elem1(v1),
      elem2(v2)
    {}

const float CutoutVolume::eval( const Vector& P ) const
{
   const float e0 = elem1->eval(P);
   const float e1 = -elem2->eval(P);
   return ( (e0 < e1) ? e0 : e1  );
}

const Vector CutoutVolume::grad( const Vector& P ) const
{
   const float e0 = elem1->eval(P);
   const float e1 = -elem2->eval(P);
   const Vector g0 = elem1->grad(P);
   const Vector g1 = -elem2->grad(P);
   return ( (e0 < e1) ? g0 : g1  );
}


  



NoiseVolume::NoiseVolume( Noise* n, const float d ) : noise (n), dx (d) { gradParams.setStep(dx); }
NoiseVolume::NoiseVolume( NoiseMachine n, const float d ) : noise (n), dx (d) { gradParams.setStep(dx); }

const float NoiseVolume::eval( const Vector& P ) const { return noise->eval(P); }

/*
const Vector NoiseVolume::grad( const Vector& P ) const
{
 //  float n = noise->eval(P);
   Vector g( noise->eval( P+Vector(dx,0,0) )-noise->eval( P-Vector(dx,0,0)), noise->eval(P+Vector(0,dx,0))-noise->eval(P-Vector(0,dx,0)), noise->eval(P+Vector(0,0,dx))-noise->eval(P-Vector(0,0,dx)) );
   g /= 2.0*dx;
   return g;
}
*/

GriddedVolume::GriddedVolume( const VolumeGrid<float>* g ) :
       grid (g),
       sparsegrid(0)
    {
       dx = g->dx();
       dy = g->dy();
       dz = g->dz();
       gradParams.setStep(dx,dy,dz);
    }

GriddedVolume::GriddedVolume( const SparseGrid* g ) :
       grid (0),
       sparsegrid(g)
    {
       dx = g->dx();
       dy = g->dy();
       dz = g->dz();
       gradParams.setStep(dx,dy,dz);
    }

/*
const Vector GriddedVolume::grad( const Vector& P ) const
{
   Vector value( (eval(P+Vector(dx,0,0))-eval(P-Vector(dx,0,0)))/dx,
                 (eval(P+Vector(0,dy,0))-eval(P-Vector(0,dy,0)))/dy,
                 (eval(P+Vector(0,0,dz))-eval(P-Vector(0,0,dz)))/dz    );
   return value*0.5;
}
*/

const float GriddedVolume::eval( const Vector& P ) const 
{
   if( grid ){ return grid->eval(P); }
   if( sparsegrid ){ return sparsegrid->eval(P); }
   return 0.0;
}


GriddedSGridVolume::GriddedSGridVolume( const ScalarGrid& g ) :
       scgrid(g)
    {
       dx = g->dx();
       dy = g->dy();
       dz = g->dz();
       gradParams.setStep(dx,dy,dz);
    }


const float GriddedSGridVolume::eval( const Vector& P ) const 
{
   return scgrid->eval(P);
}

/*
const Vector GriddedSGridVolume::grad( const Vector& P ) const
{
   Vector value( (eval(P+Vector(dx,0,0))-eval(P-Vector(dx,0,0)))/dx,
                 (eval(P+Vector(0,dy,0))-eval(P-Vector(0,dy,0)))/dy,
                 (eval(P+Vector(0,0,dz))-eval(P-Vector(0,0,dz)))/dz    );
   return value*0.5;
}
*/



GriddedFrustumVolume::GriddedFrustumVolume( const ScalarFrustumGrid& g ) :
       scgrid(g)
    {
       dx = g->dx();
       dy = g->dy();
       dz = g->dz();
       gradParams.setStep(dx,dy,dz);
    }

const float GriddedFrustumVolume::eval( const Vector& P ) const 
{
   return scgrid->eval(P);
}

/*
const Vector GriddedFrustumVolume::grad( const Vector& P ) const
{
   Vector value( (eval(P+Vector(dx,0,0))-eval(P-Vector(dx,0,0)))/dx,
                 (eval(P+Vector(0,dy,0))-eval(P-Vector(0,dy,0)))/dy,
                 (eval(P+Vector(0,0,dz))-eval(P-Vector(0,0,dz)))/dz    );
   return value*0.5;
}
*/
  
/*
GriddedOGridVolume::GriddedOGridVolume( const ScalarOGrid& g ) :
       scgrid(g)
    {
       dx = g->dx();
       dy = g->dy();
       dz = g->dz();
       gradParams.setStep(dx,dy,dz);
    }


const float GriddedOGridVolume::eval( const Vector& P ) const 
{
   return scgrid->eval(P);
}
*/



EllipseVolume::EllipseVolume( const Vector& cen, const Vector& axs, const float majorrad, const float minorrad ) : 
       center (cen), 
       axis (axs.unitvector()), 
       majorRadius (majorrad*majorrad),
       minorRadius (minorrad*minorrad)
       {}

const float EllipseVolume::eval( const Vector& P ) const 
{
   Vector X = P-center;
   float z = X*axis;
   X -= z*axis;
   return 1.0 - z*z/majorRadius - (X*X)/minorRadius;
}

const Vector EllipseVolume::grad( const Vector& P ) const
{
   Vector X = P-center;
   float z = X*axis;
   X -= z*axis;
   Vector g = -axis * (2.0*z/majorRadius) - X * ( 2.0/minorRadius);
   return g;
}

  




JitterSampleVolume::JitterSampleVolume( Volume<float>* v, const float rad, const int nb, const int seed  ) : 
       radius (rad), 
       nbSamples (nb),
       prnSeed (seed),
       elem(v)
       {}

JitterSampleVolume::JitterSampleVolume( const ScalarField& v, const float rad, const int nb, const int seed ) : 
       radius (rad), 
       nbSamples (nb),
       prnSeed (seed),
       elem(v)
       {}

const float JitterSampleVolume::eval( const Vector& P ) const 
{
   float nvalue = (P.magnitude()) * prnSeed;
   parms.seed = (int)nvalue;
   prn.setParameters(parms);
   float result = 0;
   for( int i=0;i<nbSamples;i++ )
   {
      Vector X = P + Vector( prn.eval()-0.5, prn.eval()-0.5, prn.eval()-0.5 )*radius;
      result += elem->eval(X);
   }
   result /= nbSamples;
   return result;
}

const Vector JitterSampleVolume::grad( const Vector& P ) const
{
   float nvalue = fabs(perlin.eval( P )) * prnSeed;
   parms.seed = (int)nvalue;
   prn.setParameters(parms);
   Vector result(0,0,0);
   for( int i=0;i<nbSamples;i++ )
   {
      Vector X = P + Vector( prn.eval()-0.5, prn.eval()-0.5, prn.eval()-0.5 )*radius;
      result += elem->grad(X);
   }
   result /= nbSamples;
   return result;
}

  




AdvectVolume::AdvectVolume( Volume<float>* v, Volume<Vector>* u, const float delt ) : 
       elem(v),
       velocity(u),
       dt (delt)
    {}

AdvectVolume::AdvectVolume( const ScalarField& v, const VectorField& u, const float delt ) : 
       elem(v),
       velocity(u),
       dt (delt)
    {}


const float AdvectVolume::eval( const Vector& P ) const 
{ 
   Vector X = P - ( velocity->eval(P) )*dt;
   return elem->eval(X); 
}

const Vector AdvectVolume::grad( const Vector& P ) const
{
   Vector X = P - ( velocity->eval(P) )*dt;
   Matrix M = unitMatrix() - ( velocity->grad(P) )*dt;
   return M * (elem->grad(X)); 
}

/*
GradStretchAdvectVolume::GradStretchAdvectVolume( Volume<float>* v, Volume<Vector>* u, Volume<Matrix>* gu,  const float delt, const int nb ) :
    elem(v),
    velocity(u),
    gradvelocity(gu),
    dt (delt),
    nbIterations (nb)
{}

GradStretchAdvectVolume::GradStretchAdvectVolume( const ScalarField& v, const VectorField& u, const MatrixField& gu, const float delt, const int nb ) :
    elem(v),
    velocity(u),
    gradvelocity(gu),
    dt (delt),
    nbIterations (nb)
{}

const float GradStretchAdvectVolume::eval( const Vector& P ) const 
{ 
   Vector X = P - ( velocity->eval(P) )*dt;
   return elem->eval(X); 
}
*/

/*
class BFECCAdvectVolume : public Volume<float> 
{
  public:

    BFECCAdvectVolume( const Volume<float>* v, const Volume<Vector>* u, const float dt, const int nb = 1 ) :
      base(v),
      velocity(u)
    {
       if( nb == 0 )
       {
	  const ScalarField e( new AdvectVolume( v, u, dt ) );
	  elem = e;
       }
       else
       {
          Volume<float>* advected = new BFECCAdvectVolume( v, u, dt, nb-1 );
	  Volume<Vector>* negu = new NegateVectorVolume( u );
	  Volume<float>* backadvected = new BFECCAdvectVolume( advected, negu, dt, nb-1 );
	  Volume<float>* error = new MultiplyVolume( new SubtractVolume( v, backadvected ) , 0.5);
	  advected = new AddVolume( advected, error );
	  const ScalarField e( advected );
          elem = e;
       }
    }

   ~BFECCAdvectVolume(){}

    const float eval( const Vector& P ) const 
    {
       return elem->eval(P);
    }

    const Vector grad( const Vector& P ) const
    {
       return elem->grad(P); 
    }

  

    const ScalarField base;
    VectorField velocity;
    const ScalarField elem;

*/





WarpVolume::WarpVolume( Volume<float>* v, Volume<Vector>* map ):
      elem(v),
      mapX(map)
    {}

WarpVolume::WarpVolume( const ScalarField& v, VectorField& map ):
      elem(v),
      mapX(map)
    {}

const float WarpVolume::eval( const Vector& P ) const 
{ 
   Vector X = mapX->eval(P);
   return elem->eval(X); 
}
const Vector WarpVolume::grad( const Vector& P ) const
{
   Vector X = mapX->eval(P);
   Matrix M = mapX->grad(P);
   return M * (elem->grad(X)); 
}

  


PeriodicVolume::PeriodicVolume( Volume<float>* v, const Vector& o, const Vector& L ) :
      origin (o),
      periods (L),
      elem(v)
   {}

PeriodicVolume::PeriodicVolume( const ScalarField& v, const Vector& o, const Vector& L ) :
      origin (o),
      periods (L),
      elem(v)
   {}

const float PeriodicVolume::eval( const Vector& P ) const
{
   return elem->eval( periodP(P) );
}
   
const Vector PeriodicVolume::grad( const Vector& P ) const
{
   return elem->grad( periodP(P) );
}
   
const Vector PeriodicVolume::periodP( const Vector& P ) const
{
      Vector X = P-origin;
      int  x = (int)(X[0]/periods[0]);
      int  y = (int)(X[1]/periods[1]);
      int  z = (int)(X[2]/periods[2]);
      X[0] -= x*periods[0];
      X[1] -= y*periods[1];
      X[2] -= z*periods[2];
      if( X[0] < 0 ){ X[0] += periods[0]; }
      if( X[1] < 0 ){ X[1] += periods[1]; }
      if( X[2] < 0 ){ X[2] += periods[2]; }
      return X+origin;
}




SwitchVolume::SwitchVolume( Volume<float>* v1, Volume<float>* v2, Volume<float>* swtch ) :
   elem1(v1),
   elem2(v2),
   swtchelem(swtch)
{}

SwitchVolume::SwitchVolume( const ScalarField& v1, const ScalarField& v2, const ScalarField& swtch ) :
   elem1(v1),
   elem2(v2),
   swtchelem(swtch)
{}

const float SwitchVolume::eval( const Vector& P ) const
{
   float swvalue = swtchelem->eval(P);
   if( swvalue <= 0 )
   {
      return elem1->eval(P);
   }
   return elem2->eval(P);
}
   
const Vector SwitchVolume::grad( const Vector& P ) const
{
   float swvalue = swtchelem->eval(P);
   if( swvalue <= 0 )
   {
      return elem1->grad(P);
   }
   return elem2->grad(P);
}

   



SineVolume::SineVolume( Volume<float>* v1 ) : elem1(v1) {}

SineVolume::SineVolume( const ScalarField& v1 ) : elem1(v1) {}


const float SineVolume::eval( const Vector& P ) const
{
   return sin( elem1->eval(P) );
}
   
const Vector SineVolume::grad( const Vector& P ) const
{
   return cos( elem1->eval(P) ) * elem1->grad(P);
}
   





CosineVolume::CosineVolume( Volume<float>* v1 ) : elem1(v1) {}

CosineVolume::CosineVolume( const ScalarField& v1 ) : elem1(v1) {}


const float CosineVolume::eval( const Vector& P ) const
{
   return cos( elem1->eval(P) );
}
   
const Vector CosineVolume::grad( const Vector& P ) const
{
   return -sin(elem1->eval(P) ) * elem1->grad(P);
}
   





TangentVolume::TangentVolume( Volume<float>* v1 ) : elem1(v1) {}

TangentVolume::TangentVolume( const ScalarField& v1 ) : elem1(v1) {}

const float TangentVolume::eval( const Vector& P ) const
{
   return tan( elem1->eval(P) );
}
   
const Vector TangentVolume::grad( const Vector& P ) const
{
   float cp = cos( elem1->eval(P) );
   return elem1->grad(P) /( cp*cp );
}
   






ArcsineVolume::ArcsineVolume( Volume<float>* v1 ) : elem1(v1) {}

ArcsineVolume::ArcsineVolume( const ScalarField& v1 ) : elem1(v1) {}

const float ArcsineVolume::eval( const Vector& P ) const
{
   return asin( elem1->eval(P) );
}
   
const Vector ArcsineVolume::grad( const Vector& P ) const
{
   float p = elem1->eval(P);
   return elem1->grad(P) / sqrt( 1.0 - p*p );
}
   



ArccosineVolume::ArccosineVolume( Volume<float>* v1 ) : elem1(v1) {}

ArccosineVolume::ArccosineVolume( const ScalarField& v1 ) : elem1(v1) {}

const float ArccosineVolume::eval( const Vector& P ) const
{
   return acos( elem1->eval(P) );
}
   
const Vector ArccosineVolume::grad( const Vector& P ) const
{
   float p = elem1->eval(P);
   return -(elem1->grad(P))/sqrt( 1.0 - p*p ); 
}
   

ArctangentVolume::ArctangentVolume( Volume<float>* v1 ) : elem1(v1) {}

ArctangentVolume::ArctangentVolume( const ScalarField& v1 ) : elem1(v1) {}

const float ArctangentVolume::eval( const Vector& P ) const
{
   return atan( elem1->eval(P) );
}
   
const Vector ArctangentVolume::grad( const Vector& P ) const
{
   float p = elem1->eval(P);
   return elem1->grad(P) / (1.0 + p*p);
}

   



HyperbolicSineVolume::HyperbolicSineVolume( Volume<float>* v1 ) : elem1(v1) {}

HyperbolicSineVolume::HyperbolicSineVolume( const ScalarField& v1 ) : elem1(v1) {}

const float HyperbolicSineVolume::eval( const Vector& P ) const
{
   return sinh( elem1->eval(P) );
}
   
const Vector HyperbolicSineVolume::grad( const Vector& P ) const
{
   return elem1->grad(P) * cosh( elem1->eval(P) );
}
   



HyperbolicCosineVolume::HyperbolicCosineVolume( Volume<float>* v1 ) : elem1(v1) {}

HyperbolicCosineVolume::HyperbolicCosineVolume( const ScalarField& v1 ) : elem1(v1) {}

const float HyperbolicCosineVolume::eval( const Vector& P ) const
{
   return cosh( elem1->eval(P) );
}
   
const Vector HyperbolicCosineVolume::grad( const Vector& P ) const
{
   return elem1->grad(P) * sinh( elem1->eval(P) );
}
   



HyperbolicTangentVolume::HyperbolicTangentVolume( Volume<float>* v1 ) : elem1(v1) {}

HyperbolicTangentVolume::HyperbolicTangentVolume( const ScalarField& v1 ) : elem1(v1) {}

const float HyperbolicTangentVolume::eval( const Vector& P ) const
{
   return tanh( elem1->eval(P) );
}
   
const Vector HyperbolicTangentVolume::grad( const Vector& P ) const
{
   float cp = cosh(elem1->eval(P));
   return elem1->grad(P) / ( cp*cp );
}
   


XIdentityVolume::XIdentityVolume() : gradvalue(Vector(1,0,0)) {}
const float XIdentityVolume::eval( const Vector& P ) const { return P[0]; }
const Vector XIdentityVolume::grad( const Vector& P ) const { return gradvalue; }
   
YIdentityVolume::YIdentityVolume() : gradvalue(Vector(0,1,0)) {}
const float YIdentityVolume::eval( const Vector& P ) const { return P[1]; }
const Vector YIdentityVolume::grad( const Vector& P ) const { return gradvalue; }
   
ZIdentityVolume::ZIdentityVolume() : gradvalue(Vector(0,0,1)) {}
const float ZIdentityVolume::eval( const Vector& P ) const { return P[2]; }
const Vector ZIdentityVolume::grad( const Vector& P ) const { return gradvalue; }
   


