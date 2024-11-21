

#ifndef __IMPLICITCOLORS_H__
#define __IMPLICITCOLORS_H__


#include "Volume.h"
#include "Noise.h"
#include "UniformPRN.h"
#include "LinearAlgebra.h"
#include "SparseGrid.h"
//#include "OVDBGrid.h"
#include "PerlinNoise.h"

namespace lux
{

class ConstantColor : public Volume<Color> 
{
  public:

    ConstantColor( const Color& v ) :
       value (v)
    {}

   ~ConstantColor(){}


    const Color eval( const Vector& P ) const 
    {
       return value; 
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Constant";
       return lbl;
    }


  private:

    Color value;
};


class FloatMultiplyColor : public Volume<Color> 
{
  public:

    FloatMultiplyColor( Volume<Color> * e, const float v ) :
       elem(e),
       value (v)
    {}

    FloatMultiplyColor( const ColorField& e, const float v ) :
       elem(e),
       value (v)
    {}

   ~FloatMultiplyColor(){}


    const Color eval( const Vector& P ) const 
    {
       return (elem->eval(P))*value; 
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Multiply";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ColorField elem;
    float value;
};


class VolumeMultiplyColor : public Volume<Color> 
{
  public:

    VolumeMultiplyColor( Volume<Color> * e, Volume<float>* v ) :
       elem(e),
       factor(v)
    {}

    VolumeMultiplyColor( const ColorField& e, const ScalarField& v ) :
       elem(e),
       factor(v)
    {}

   ~VolumeMultiplyColor(){}


    const Color eval( const Vector& P ) const 
    {
       return (elem->eval(P)) * (factor->eval(P)); 
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Multiply";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ",";
       lbl = lbl + factor->typelabel();
       lbl = lbl + ")";
       return lbl;
    }


  private:

    ColorField elem;
    ScalarField factor;
};







class FloatDivideColor : public Volume<Color> 
{
  public:

    FloatDivideColor( Volume<Color> * e, const float v ) :
       elem(e),
       value (v)
    {}

    FloatDivideColor( const ColorField& e, const float v ) :
       elem(e),
       value (v)
    {}

   ~FloatDivideColor(){}


    const Color eval( const Vector& P ) const 
    {
       return (elem->eval(P))/value; 
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Divide";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ColorField elem;
    float value;
};


class VolumeDivideColor : public Volume<Color> 
{
  public:

    VolumeDivideColor( Volume<Color> * e, Volume<float>* v ) :
       elem(e),
       factor(v)
    {}

    VolumeDivideColor( const ColorField& e, const ScalarField& v ) :
       elem(e),
       factor(v)
    {}

   ~VolumeDivideColor(){}


    const Color eval( const Vector& P ) const 
    {
       return (elem->eval(P)) / (factor->eval(P)); 
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Divide";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ",";
       lbl = lbl + factor->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ColorField elem;
    ScalarField factor;
};














class MultiplyColor : public Volume<Color> 
{
  public:

    MultiplyColor( Volume<Color> * e, Volume<Color>* v ) :
       elem(e),
       factor(v)
    {}

    MultiplyColor( const ColorField& e, const ColorField& v ) :
       elem(e),
       factor(v)
    {}

   ~MultiplyColor(){}


    const Color eval( const Vector& P ) const 
    {
       return (elem->eval(P)) * (factor->eval(P)); 
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Multiply";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ",";
       lbl = lbl + factor->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ColorField elem;
    ColorField factor;
};


class ChromaColor : public Volume<Color> 
{
  public:

    ChromaColor( Volume<Color> * e ) :
      elem(e)
    {}

    ChromaColor( const ColorField& e ) :
      elem(e)
    {}

   ~ChromaColor(){}


    const Color eval( const Vector& P ) const 
    {
       Color cd = elem->eval(P);
       float denom = cd[0] + cd[1] + cd[2];
       if( denom > 0.0 )
       {
          cd /= denom;
       }
       return cd;
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Chroma";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ColorField elem;
};




class FloatGammaColor : public Volume<Color> 
{
  public:

    FloatGammaColor( Volume<Color> * e, const float v ) :
       elem (e),
       value (v)
    {}

    FloatGammaColor( const ColorField& e, const float v ) :
       elem (e),
       value (v)
    {}

   ~FloatGammaColor(){}


    const Color eval( const Vector& P ) const 
    {
       Color cd = elem->eval(P);
       cd[0] = std::pow( (double)cd[0], (double)value );
       cd[1] = std::pow( (double)cd[1], (double)value );
       cd[2] = std::pow( (double)cd[2], (double)value );
       return cd;
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Gamma";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ColorField elem;
    float value;
};


class VolumeGammaColor : public Volume<Color> 
{
  public:

    VolumeGammaColor( Volume<Color> * e, Volume<float>* v ) :
      elem (e),
      gam  (v)
    {}

    VolumeGammaColor( const ColorField& e, const ScalarField& v ) :
      elem (e),
      gam  (v)
    {}

   ~VolumeGammaColor(){}


    const Color eval( const Vector& P ) const 
    {
       Color cd = elem->eval(P);
       float value = gam->eval(P); 
       cd[0] = std::pow( (double)cd[0], (double)value );
       cd[1] = std::pow( (double)cd[1], (double)value );
       cd[2] = std::pow( (double)cd[2], (double)value );
       return cd;
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Gamma";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ",";
       lbl = lbl + gam->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ColorField elem;
    ScalarField gam;
};







class AddColor : public Volume<Color> 
{
  public:

    AddColor( Volume<Color> * e, Volume<Color>* v ) :
      elem1(e),
      elem2(v)
    {}

    AddColor( const ColorField& e, const ColorField& v ) :
      elem1(e),
      elem2(v)
    {}

   ~AddColor(){}


    const Color eval( const Vector& P ) const 
    {
       return ( elem1->eval(P) + elem2->eval(P) ); 
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Add";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ColorField elem1, elem2;
};


class SubtractColor : public Volume<Color> 
{
  public:

    SubtractColor( Volume<Color> * e, Volume<Color>* v ) :
      elem1(e),
      elem2(v)
    {}

    SubtractColor( const ColorField& e, const ColorField& v ) :
      elem1(e),
      elem2(v)
    {}

   ~SubtractColor(){}


    const Color eval( const Vector& P ) const 
    {
       return ( elem1->eval(P) - elem2->eval(P) ); 
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Subtract";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ColorField elem1, elem2;
};





class GriddedColor : public Volume<Color> 
{
  public:

    GriddedColor( const VolumeGrid<Color>* g ) :
       grid (g),
       sgrid(0)
    {}

    GriddedColor( const SparseColorGrid* g ) :
       grid (0),
       sgrid(g)
    {}

   ~GriddedColor(){}

    const Color eval( const Vector& P ) const 
    {
       if( grid ){ return grid->eval(P); }
       if( sgrid ){ return sgrid->eval(P); }
       return Color(0,0,0,0);
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Gridded";
       return lbl;
    }

  private:

    const VolumeGrid<Color>* grid;
    const SparseColorGrid* sgrid;
};





class GriddedSGridColor : public Volume<Color> 
{
  public:

    GriddedSGridColor( const ColorGrid& g ) :
       elem(g)
    {}


   ~GriddedSGridColor(){}

    const Color eval( const Vector& P ) const 
    {
       return elem->eval(P);
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Gridded";
       return lbl;
    }

  private:

    ColorGrid elem;
};



class GriddedFrustumColor : public Volume<Color> 
{
  public:

    GriddedFrustumColor( const ColorFrustumGrid& g ) :
       felem ( g )
    {}


   ~GriddedFrustumColor(){}

    const Color eval( const Vector& P ) const 
    {
       return felem->eval(P);
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Gridded";
       return lbl;
    }

  private:

    ColorFrustumGrid felem;
};




/*
class GriddedOGridColor : public Volume<Color> 
{
  public:

    GriddedOGridColor( const ColorOGrid& g ) :
       elem(g)
    {}


   ~GriddedOGridColor(){}

    const Color eval( const Vector& P ) const 
    {
       return elem->eval(P);
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Gridded";
       return lbl;
    }

  private:

    ColorOGrid elem;
};
*/













class JitterSampleColor : public Volume<Color>
{
  public:

    JitterSampleColor( Volume<Color>* v, const float rad, const int nb, const int seed = 2847573 ) : 
       elem(v),
       radius (rad), 
       nbSamples (nb),
       prnSeed (seed)
       {}

   ~JitterSampleColor(){}

    const Color eval( const Vector& P ) const 
    {
       float nvalue = fabs(perlin.eval( P )) * prnSeed;
       parms.seed = (int)nvalue;
       prn.setParameters(parms);
       Color result;
       for( int i=0;i<nbSamples;i++ )
       {
          Vector X = P + Vector( prn.eval()-0.5, prn.eval()-0.5, prn.eval()-0.5 )*radius;
	  result += elem->eval(X);
       }
       result /= nbSamples;
       return result;
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Jitter";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }


  private:

    ColorField elem;
    float radius;
    int nbSamples;
    int prnSeed;

    mutable Noise_t parms;
    mutable UniformPRN prn;
    PerlinNoise perlin;

};

class RotateColor: public Volume<Color> 
{
  public:

    RotateColor( Volume<Color> * v, const Vector& s ) :
      elem (v),
      R ( inverse(rotation( s.unitvector() , s.magnitude() )) )
    {}

    RotateColor( const ColorField& v, const Vector& s ) :
      elem (v),
      R ( inverse(rotation( s.unitvector() , s.magnitude() )) )
    {}

    ~RotateColor(){}


    const Color eval( const Vector& P ) const
    {
       Vector X = R*P;
       return  elem->eval(X);
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Rotate";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }


  private:

    ColorField elem;
    Matrix R;
};



class ScaleColor : public Volume<Color> 
{
  public:

    ScaleColor( Volume<Color> * v, const Vector& s ) ;

    ScaleColor( const ColorField& v, const Vector& s ) ;

    ScaleColor( Volume<Color> * v, const float& s ) ;

    ScaleColor( const ColorField& v, const float& s ) ;



    ~ScaleColor(){}


    const Color eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Scale";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }


  private:

    ColorField elem;
    Vector scale;
};





class TranslateColor: public Volume<Color> 
{
  public:

    TranslateColor( Volume<Color> * v, const Vector& s ) ;

    TranslateColor( const ColorField& v, const Vector& s ) ;

    ~TranslateColor(){}


    const Color eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Translate";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }



  private:

    ColorField elem;
    Vector translate;
};


class NegateColor : public Volume<Color>
{
  public:

    NegateColor( Volume<Color>* v );
    NegateColor( const ColorField& v );

   ~NegateColor(){}


    const Color eval( const Vector& P );


    virtual std::string typelabel() 
    { 
       std::string lbl = "Negate";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ColorField elem;

};


class ComponentColor: public Volume<Color> 
{
  public:

    ComponentColor( Volume<float>* r, Volume<float>* g, Volume<float>* b ) :
       red(r),
       green(g),
       blue(b)
    {}

    ComponentColor( const ScalarField& r, const ScalarField& g, const ScalarField& b ) :
       red(r),
       green(g),
       blue(b)
    {}

    ~ComponentColor(){}


    const Color eval( const Vector& P ) const
    {
       float r = red->eval(P);
       float g = green->eval(P);
       float b = blue->eval(P);

       if( r < 0 ){ r = 0; }
       if( g < 0 ){ g = 0; }
       if( b < 0 ){ b = 0; }
       return  Color(  r, g, b, 1.0 );
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Component";
       lbl = lbl + "(";
       lbl = lbl + red->typelabel();
       lbl = lbl + ",";
       lbl = lbl + green->typelabel();
       lbl = lbl + ",";
       lbl = lbl + blue->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ScalarField red, green, blue;
};

class RColorVolume : public Volume<float> 
{
  public:

    RColorVolume( Volume<Color>* v ) : elem(v) {}
    RColorVolume( const ColorField& v ) : elem(v) {}

   ~RColorVolume(){}

    const float eval( const Vector& P ) const
    {
      const Color value = elem->eval(P);
      return value.red();
    }

    virtual std::string typelabel() 
    { 
       std::string lbl = "RComponent";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }


  private:

    ColorField elem;
};



class GColorVolume : public Volume<float> 
{
  public:

    GColorVolume( Volume<Color>* v ) : elem(v) {}
    GColorVolume( const ColorField& v ) : elem(v) {}

   ~GColorVolume(){}


    const float eval( const Vector& P ) const
    {
      const Color value = elem->eval(P);
      return value.green();
    }

    virtual std::string typelabel() 
    { 
       std::string lbl = "GComponent";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }



  private:

    ColorField elem;
};



class BColorVolume : public Volume<float> 
{
  public:

    BColorVolume( Volume<Color>* v ) : elem(v) {}
    BColorVolume( const ColorField& v ) : elem(v) {}

   ~BColorVolume(){}


    const float eval( const Vector& P ) const
    {
      const Color value = elem->eval(P);
      return value.blue();
    }

    virtual std::string typelabel() 
    { 
       std::string lbl = "BComponent";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }



  private:

    ColorField elem;
};



class VolumeLUTColor: public Volume<Color> 
{
  public:

    VolumeLUTColor( Volume<float>* field, const vector<Color> lut, const float minValue, const float maxValue, const float bright = 1, const float gam = 1 ) ;

    VolumeLUTColor( const ScalarField& field, const vector<Color> lut, const float minValue, const float maxValue, const float bright = 1, const float gam = 1 ) ;

    ~VolumeLUTColor(){}


    const Color eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "LUT";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

   ScalarField elem;
   vector<Color> LUT;
   float minV, maxV;
   float brightness, gamma;
};



class ReportColorVolume: public Volume<Color> 
{
  public:

    ReportColorVolume( Volume<Color>* v, const string tag ) :
       elem(v),
       label (tag)
    {}

    ReportColorVolume( const ColorField& v, const string tag ) :
       elem(v),
       label (tag)
    {}

   ~ReportColorVolume(){}


    const Color eval( const Vector& P ) const 
    {
       Color value =  elem->eval(P);
       cout << label << " eval: " << value[0] << " " << value[1] << " " << value[2] << " " << value[3] << " at position " << P[0] << " " << P[1] << " " << P[2] << endl;
       return value; 
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Report";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }


  private:

    ColorField elem;
    string label;

};


class WarpColorVolume : public Volume<Color> 
{
  public:

    WarpColorVolume( Volume<Color>* v, Volume<Vector>* u ) :
      elem(v),
      warp(u)
    {}

    WarpColorVolume( const ColorField& v, const VectorField& u ) :
      elem(v),
      warp(u)
    {}

   ~WarpColorVolume(){}

    const Color eval( const Vector& P ) const 
    { 
       Vector X = warp->eval(P);
       return elem->eval(X); 
    }


    virtual std::string typelabel() 
    { 
       std::string lbl = "Warp";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ",";
       lbl = lbl + warp->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ColorField elem;
    VectorField warp;
};


class PeriodicColorVolume : public Volume<Color>
{
  public:
   PeriodicColorVolume( Volume<Color>* v, const Vector& o, const Vector& L ) :
      elem(v),
      origin (o),
      periods (L)
   {}

   PeriodicColorVolume( const ColorField& v, const Vector& o, const Vector& L ) :
      elem(v),
      origin (o),
      periods (L)
   {}

   ~PeriodicColorVolume(){}

   const Color eval( const Vector& P ) const
   {
      return elem->eval( periodP(P) );
   }
   
  

    virtual std::string typelabel() 
    { 
       std::string lbl = "Periodic";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }


  private:

   ColorField elem;
   Vector origin;
   Vector periods;

   const Vector periodP( const Vector& P ) const
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

};

class AdvectColorVolume : public Volume<Color> 
{
  public:

    AdvectColorVolume( Volume<Color>* v, Volume<Vector>* u, const float delt ) ;

    AdvectColorVolume( const ColorField& v, const VectorField& u, const float delt ) ;

   ~AdvectColorVolume(){}

    const Color eval( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Advect";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ",";
       lbl = lbl + velocity->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ColorField elem;
    VectorField velocity;
    float dt;
};



class SwitchColorVolume : public Volume<Color>
{
  public:

   SwitchColorVolume( Volume<Color>* v1, Volume<Color>* v2, Volume<float>* swtch );

   SwitchColorVolume( const ColorField& v1, const ColorField& v2, const ScalarField& swtch );

   ~SwitchColorVolume(){}

   const Color eval( const Vector& P ) const;
   

    virtual std::string typelabel() 
    { 
       std::string lbl = "Switch";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ",";
       lbl = lbl + swtchelem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }
 

  private:

   const ColorField elem1;
   const ColorField elem2;
   const ScalarField swtchelem;


};





}
#endif
