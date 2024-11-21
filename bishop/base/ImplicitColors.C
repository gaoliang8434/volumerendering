
#include "ImplicitColors.h"

using namespace lux;


ScaleColor::ScaleColor( Volume<Color> * v, const Vector& s ) :
       elem(v),
       scale (s)
    {}

ScaleColor::ScaleColor( const ColorField& v, const Vector& s ) :
       elem(v),
       scale (s)
    {}

ScaleColor::ScaleColor( Volume<Color> * v, const float& s ) :
       elem(v),
       scale (Vector(s,s,s))
    {}

ScaleColor::ScaleColor( const ColorField& v, const float& s ) :
       elem(v),
       scale (Vector(s,s,s))
    {}

const Color ScaleColor::eval( const Vector& P ) const
{
   Vector X = P;
   X[0] /= scale[0];
   X[1] /= scale[1];
   X[2] /= scale[2];
   return elem->eval(X);
}





TranslateColor::TranslateColor( Volume<Color> * v, const Vector& s ) :
       elem(v),
       translate (s)
    {}

TranslateColor::TranslateColor( const ColorField& v, const Vector& s ) :
       elem(v),
       translate (s)
    {}



const Color TranslateColor::eval( const Vector& P ) const
{
   Vector X = P - translate;
   return  elem->eval(X);
}

NegateColor::NegateColor( Volume<Color>* v ) : elem(v) {}

NegateColor::NegateColor( const ColorField& v ) : elem(v) {}

const Color NegateColor::eval( const Vector& P ) { return -(elem->eval(P)); }





VolumeLUTColor::VolumeLUTColor( Volume<float>* field, const vector<Color> lut, const float minValue, const float maxValue, const float bright , const float gam  ) :
      elem(field),
      LUT (lut),
      minV (minValue),
      maxV (maxValue),
      brightness (bright),
      gamma (gam)
    {}

VolumeLUTColor::VolumeLUTColor( const ScalarField& field, const vector<Color> lut, const float minValue, const float maxValue, const float bright , const float gam ) :
      elem(field),
      LUT (lut),
      minV (minValue),
      maxV (maxValue),
      brightness (bright),
      gamma (gam)
    {}


const Color VolumeLUTColor::eval( const Vector& P ) const
{
       float fvalue = elem->eval(P);
       fvalue -= minV;
       fvalue /= (maxV-minV);
       if (fvalue < 0 ){ fvalue = 0; }
       if (fvalue > 1.0 ){ fvalue = 1.0; }
       fvalue *= LUT.size();
       int i0 = (int)fvalue;
       int i1 = i0+1;
       float w1 = fvalue-i0;
       float w0 = 1.0 - w1;
       if( i1 >= LUT.size() ){ i1 = i0; }
       Color value = LUT[i0]*w0 + LUT[i1]*w1;
       value[0] = std::pow((double)value[0], (double)gamma) * brightness;
       value[1] = std::pow((double)value[1], (double)gamma) * brightness;
       value[2] = std::pow((double)value[2], (double)gamma) * brightness;
       return value;
}





AdvectColorVolume::AdvectColorVolume( Volume<Color>* v, Volume<Vector>* u, const float delt ) :
       elem(v),
       velocity(u),
       dt (delt)
    {}

AdvectColorVolume::AdvectColorVolume( const ColorField& v, const VectorField& u, const float delt ) :
       elem(v),
       velocity(u),
       dt (delt)
    {}


const Color AdvectColorVolume::eval( const Vector& P ) const 
{ 
   Vector X = P - ( velocity->eval(P) )*dt;
   return elem->eval(X); 
}




SwitchColorVolume::SwitchColorVolume( Volume<Color>* v1, Volume<Color>* v2, Volume<float>* swtch ) :
   elem1(v1),
   elem2(v2),
   swtchelem(swtch)
{}

SwitchColorVolume::SwitchColorVolume( const ColorField& v1, const ColorField& v2, const ScalarField& swtch ) :
   elem1(v1),
   elem2(v2),
   swtchelem(swtch)
{}

const Color SwitchColorVolume::eval( const Vector& P ) const
{
   float swvalue = swtchelem->eval(P);
   if( swvalue <= 0 )
   {
      return elem1->eval(P);
   }
   return elem2->eval(P);
}

