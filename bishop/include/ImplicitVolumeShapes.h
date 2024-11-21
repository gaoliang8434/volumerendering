
#ifndef __IMPLICITVOLUMESHAPES_H__
#define __IMPLICITVOLUMESHAPES_H__

#include "Volume.h"
#include "LinearAlgebra.h"
#include "Noise.h"
#include "FFTNoise.h"
#include "UniformPRN.h"
#include "VolumeGrid.h"
#include "SparseGrid.h"
//#include "OVDBGrid.h"
#include "ImplicitVectorShapes.h"
#include <cmath>
#include <vector>
#include <iostream>
#include "PerlinNoise.h"



namespace lux
{


class ConstantVolume : public Volume<float> 
{
  public:

    ConstantVolume( const float v ) ;

   ~ConstantVolume(){}


    const float eval( const Vector& P ) const; 

    const Vector grad( const Vector& P ) const; 

   virtual std::string typelabel() { return "Constant"; }


  private:

    float value;
    Vector gradvalue;
};


class ExpVolume: public Volume<float> 
{
  public:

    ExpVolume( Volume<float>* v );

    ExpVolume( const ScalarField& v );

   ~ExpVolume(){}


    const float eval( const Vector& P ) const; 

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Exp";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField elem;

};


class ReportVolume: public Volume<float> 
{
  public:

    ReportVolume( Volume<float>* v, const string tag );

    ReportVolume( const ScalarField& v, const string tag );

   ~ReportVolume(){}


    const float eval( const Vector& P ) const;
    const Vector grad( const Vector& P ) const;


    virtual std::string typelabel() 
    {
       std::string lbl = "Report";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    const ScalarField elem;
    string label;

};



class SphereVolume : public Volume<float> 
{
  public:

    SphereVolume( const Vector& cen, const float rad );

   ~SphereVolume(){}


    const float eval( const Vector& P ) const; 

    const Vector grad(  const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Sphere"); }


  private:

    Vector center;
    float radius;
};



class CsgBoxVolume : public Volume<float> 
{
  public:

    CsgBoxVolume( const Vector& cen, const float rad, const float pwr );

   ~CsgBoxVolume(){}


    const float eval( const Vector& P ) const; 

    const Vector grad(  const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("CsgBox"); }

  private:

    Vector center;
    float radius;
    double power;
};





class CsgRectangleVolume : public Volume<float> 
{
  public:

    CsgRectangleVolume( const Vector& cen, const float rad, const Vector& asp, const float pwr );

   ~CsgRectangleVolume(){}


    const float eval( const Vector& P ) const; 

    const Vector grad(  const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("CsgRectangle"); }

  private:

    Vector center;
    float radius;
    Vector aspect;
    double power;
};







class ConeVolume : public Volume<float> 
{
  public:

    ConeVolume( const Vector& cen, const Vector& ax, const float h, const float theta );

   ~ConeVolume(){}


    const float eval( const Vector& P ) const; 

    const Vector grad(  const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Cone"); }

  private:

    Vector center;
    Vector axis;
    float height;
    float angle;
};






class PlaneVolume : public Volume<float> 
{
  public:

    PlaneVolume( const Vector cen, const Vector norm );

   ~PlaneVolume(){}


    const float eval( const Vector& P ) const; 

    const Vector grad(  const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Plane"); }

  private:

    Vector center, normal;
};



class HardBoxVolume : public Volume<float>
{
  public:

    HardBoxVolume( const Vector _llc, const Vector _urc );

   ~HardBoxVolume(){}


    const float eval( const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("HardBox"); }

  private:

  Vector llc, urc;

};


class TorusVolume : public Volume<float> 
{
  public:

    TorusVolume( const Vector& cen, const Vector& axis, const float majorRad, const float minorRad );

   ~TorusVolume(){}


    const float eval( const Vector& P ) const; 

    const Vector grad(  const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Torus"); }

  private:

    Vector tcenter, taxis;
    float majorRadiusSquare, minorRadiusSquare;
};






class MobiusStripVolume : public Volume<float> 
{
  public:

    MobiusStripVolume( const Vector& cen, const Vector& axis, const float rad, const float thick );

   ~MobiusStripVolume(){}


    const float eval( const Vector& P ) const; 

  private:

    Vector tcenter, taxis;
    float radius, thickness;
};






class SteinerPatchVolume : public Volume<float> 
{
  public:

    SteinerPatchVolume();

   ~SteinerPatchVolume(){}


    const float eval( const Vector& P ) const; 

    const Vector grad(  const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("SteinerPatch"); }
};


class IcosahedronVolume : public Volume<float> 
{
  public:

    IcosahedronVolume();

   ~IcosahedronVolume(){}


    const float eval( const Vector& P ) const; 

    virtual std::string typelabel() { return std::string("Icosahedron"); }
};




class ScaleVolume : public Volume<float> 
{
  public:

    ScaleVolume( Volume<float> * v, const Vector& s );

    ScaleVolume( Volume<float> * v, const float& s );


    ScaleVolume( const ScalarField& v, const Vector& s );

    ScaleVolume( const ScalarField& v, const float& s );


    ~ScaleVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Scale";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField elem;
    Vector scale;
};



class TranslateVolume : public Volume<float> 
{
  public:

    TranslateVolume( Volume<float> * v, const Vector& s );

    TranslateVolume( const ScalarField& v, const Vector& s );

    ~TranslateVolume(){}


    const float eval( const Vector& P ) const;


    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Translate";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }



  private:

    const ScalarField elem;
    Vector translate;
};


class RotateVolume: public Volume<float> 
{
  public:

    RotateVolume( Volume<float> * v, const Vector& s );

    RotateVolume( const ScalarField& v, const Vector& s );

    ~RotateVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Rotate";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

    const ScalarField elem;
    Vector axis;
    double sina, cosa;
    Matrix R;
};






class NegateVolume : public Volume<float> 
{
  public:

    NegateVolume( Volume<float> * v );

    NegateVolume( const ScalarField& v );

    ~NegateVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Negate";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:
  
    const ScalarField elem;
};



class AbsoluteVolume : public Volume<float> 
{
  public:

    AbsoluteVolume( Volume<float> * v )  ; 

    AbsoluteVolume( const ScalarField& v )  ; 

    ~AbsoluteVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Absolute";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField elem;

};







class MultiplyVolume : public Volume<float> 
{
  public:

    MultiplyVolume( Volume<float> * v, const float a ); 

    MultiplyVolume( Volume<float> * v, Volume<float>* u ); 

    MultiplyVolume( const ScalarField& v, const float a ); 

    MultiplyVolume( const ScalarField& v, const ScalarField& u ); 


    ~MultiplyVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;

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

    const ScalarField elem;
    const ScalarField factor;
    float constant;
};





class DivideVolume : public Volume<float> 
{
  public:

    DivideVolume( Volume<float> * v, const float a ); 

    DivideVolume( Volume<float> * v, Volume<float>* u ); 

    DivideVolume( const ScalarField& v, const float a ); 

    DivideVolume( const ScalarField& v, const ScalarField& u ); 


    ~DivideVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;

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

    const ScalarField elem;
    const ScalarField factor;
    float constant;
};




class AddVolume : public Volume<float> 
{
  public:

    AddVolume( Volume<float> * v1, Volume<float> * v2 );

    AddVolume( const ScalarField&  v1, const ScalarField& v2 );

    ~AddVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Add";
      lbl = lbl + "(";
      lbl = lbl + e1->typelabel();
      lbl = lbl + ",";
      lbl = lbl + e2->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField e1, e2;
};


class SubtractVolume : public Volume<float> 
{
  public:

    SubtractVolume( Volume<float> * v1, Volume<float> * v2 );

    SubtractVolume( const ScalarField& v1, const ScalarField& v2 );

    ~SubtractVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Subtract";
      lbl = lbl + "(";
      lbl = lbl + e1->typelabel();
      lbl = lbl + ",";
      lbl = lbl + e2->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField e1, e2;
};



class MaskVolume : public Volume<float> 
{
  public:

    MaskVolume( Volume<float> * v );

    MaskVolume( const ScalarField& v );

    ~MaskVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Mask";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField elem;

};

class ClampVolume : public Volume<float> 
{
  public:

    ClampVolume( Volume<float> * v, float minv, float maxv );

    ClampVolume( const ScalarField& v, float minv, float maxv );


    ~ClampVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Clamp";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

   const ScalarField elem;
    float vmin, vmax;
};


class GammaVolume : public Volume<float> 
{
  public:

    GammaVolume( Volume<float> * v, float gam );

    GammaVolume( const ScalarField& v, float gam );



    ~GammaVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Gamma";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField elem;
    float gamma;
};



class VolumeGammaVolume : public Volume<float> 
{
  public:

    VolumeGammaVolume( Volume<float> * v, Volume<float>* gam );

    VolumeGammaVolume( const ScalarField& v, const ScalarField& gam );



    ~VolumeGammaVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Gamma";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ",";
      lbl = lbl + gamma->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField elem;
    const ScalarField gamma;
};




class BlinnBlendVolume : public Volume<float> 
{
  public:

    BlinnBlendVolume( Volume<float> * v1, Volume<float> * v2, const float _alpha = 1.0 );

    BlinnBlendVolume( const ScalarField& v1, const ScalarField& v2, const float _alpha = 1.0 );


    ~BlinnBlendVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "BlinnBlend";
      lbl = lbl + "(";
      lbl = lbl + elem1->typelabel();
      lbl = lbl + ",";
      lbl = lbl + elem2->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  protected:

    const ScalarField elem1, elem2;
    float alpha;
};


class MultiBlendVolume : public Volume<float> 
{
  public:

    MultiBlendVolume( std::vector<Volume<float>*>& vs, const float a = 1.0  );

    MultiBlendVolume( std::vector<ScalarField>& vs, const float a = 1.0  );

    ~MultiBlendVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "MultiBlend";
      lbl = lbl + "(";
      for(size_t i=0;i<elems.size();i++)
      {
         lbl = lbl + elems[i]->typelabel();
         if(i<elems.size()-1) {lbl = lbl + ",";}
      }
      lbl = lbl + ")";
      return lbl;
   }

  private:

    float alpha;
    std::vector<ScalarField> elems;
};




class UnionVolume : public Volume<float> 
{
  public:

    UnionVolume( Volume<float> * v1, Volume<float> * v2 );

    UnionVolume( const ScalarField& v1, const ScalarField& v2 );

    ~UnionVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Union";
      lbl = lbl + "(";
      lbl = lbl + elem1->typelabel();
      lbl = lbl + ",";
      lbl = lbl + elem2->typelabel();
      lbl = lbl + ")";
      return lbl;
   }
  private:

  const ScalarField elem1, elem2;
};



class IntersectionVolume : public Volume<float> 
{
  public:

    IntersectionVolume( Volume<float> * v1, Volume<float> * v2 );

    IntersectionVolume( const ScalarField& v1, const ScalarField& v2 );

    ~IntersectionVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Intersection";
      lbl = lbl + "(";
      lbl = lbl + elem1->typelabel();
      lbl = lbl + ",";
      lbl = lbl + elem2->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const ScalarField elem1, elem2;
};


class CutoutVolume : public Volume<float> 
{
  public:

    CutoutVolume( Volume<float> * v1, Volume<float> * v2 );

    CutoutVolume( const ScalarField& v1, const ScalarField& v2 );

    ~CutoutVolume(){}


    const float eval( const Vector& P ) const;

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Cutout";
      lbl = lbl + "(";
      lbl = lbl + elem1->typelabel();
      lbl = lbl + ",";
      lbl = lbl + elem2->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

    const ScalarField elem1, elem2;

};




class PyroclasticVolume : public Volume<float> 
{
  public:

    PyroclasticVolume( const Vector& Center, const float Radius, const float Amp, 
                       const float octaves, const float freq, const float rough, const Vector trans, const float time, const float Gamma = 1.0/3.0 );

   ~PyroclasticVolume(){}

    const float eval( const Vector& P ) const;
    //const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Pyroclastic";
      return lbl;
   }

  private:

    FractalSum<PerlinNoiseGustavson> noise;
    float amplitude;
    Vector center;
    float gamma;
    const ScalarField elem;

    float dx;
};




class RadialPyroclasticVolume : public Volume<float> 
{
  public:

    RadialPyroclasticVolume( const Vector& Center, const float Radius, const float Amp, 
                       const float octaves, const float freq, const float rough, const float trans, const float time, const float Gamma = 1.0/3.0 );

   ~RadialPyroclasticVolume(){}

    const float eval( const Vector& P ) const;
    //const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "RadialPyroclastic";
      return lbl;
   }

  private:

    mutable FractalSum<PerlinNoiseGustavson> noise;
    float amplitude;
    Vector center;
    float gamma;
    float translate;
    float radius;
    mutable Noise_t parms;

    float dx;
};



class FFTNoise;

class FFTNoiseVolume : public Volume<float> 
{
  public:

    FFTNoiseVolume( const float power, const float low, const float high, const float length, const int sz );

   ~FFTNoiseVolume(){ delete noise;}

    const float eval( const Vector& P ) const; 
    //const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "FFTNoise";
      return lbl;
   }

  private:

    FFTNoise * noise;
    float dx;
};

class NoiseVolume : public Volume<float> 
{
  public:

    NoiseVolume( Noise* n, const float d = 0.01 ); 
    NoiseVolume( NoiseMachine n, const float d = 0.01 ); 

   ~NoiseVolume(){}

    const float eval( const Vector& P ) const; 
    //const Vector grad(  const Vector& P ) const;  

   virtual std::string typelabel() 
   { 
      std::string lbl = "Noise";
      return lbl;
   }

  private:

    NoiseMachine noise;
    float dx;
};


class GriddedVolume : public Volume<float> 
{
  public:

    GriddedVolume( const VolumeGrid<float>* g );

    GriddedVolume( const SparseGrid* g );

   ~GriddedVolume(){}

    const float eval( const Vector& P ) const;
    //const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Gridded";
      return lbl;
   }

  private:

    const VolumeGrid<float>* grid;
    const SparseGrid* sparsegrid;
    float dx, dy, dz;
};


class GriddedSGridVolume : public Volume<float> 
{
  public:

    GriddedSGridVolume( const ScalarGrid& g );

   ~GriddedSGridVolume(){}

    const float eval( const Vector& P ) const;
    //const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Gridded";
      return lbl;
   }

  private:

    ScalarGrid scgrid;
    float dx, dy, dz;
};


class GriddedFrustumVolume : public Volume<float> 
{
  public:

    GriddedFrustumVolume( const ScalarFrustumGrid& g );

   ~GriddedFrustumVolume(){}

    const float eval( const Vector& P ) const;
    //const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Gridded";
      return lbl;
   }

  private:

    ScalarFrustumGrid scgrid;
    float dx, dy, dz;
};

/*
class GriddedOGridVolume : public Volume<float> 
{
  public:

    GriddedOGridVolume( const ScalarOGrid& g );

   ~GriddedOGridVolume(){}

    const float eval( const Vector& P ) const;
    //const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Gridded";
      return lbl;
   }

  private:

    ScalarOGrid scgrid;
    float dx, dy, dz;
};
*/


// This is NOT a signed distance function, although it is a signed implicit function
class EllipseVolume : public Volume<float>
{
  public:

    EllipseVolume( const Vector& cen, const Vector& axs, const float majorrad, const float minorrad ); 

   ~EllipseVolume(){}

    const float eval( const Vector& P ) const; 

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Ellipse";
      return lbl;
   }

  private:

    Vector center, axis;
    float majorRadius, minorRadius;

};

class JitterSampleVolume : public Volume<float>
{
  public:

    JitterSampleVolume( Volume<float>* v, const float rad, const int nb, const int seed = 2847573 ); 

    JitterSampleVolume( const ScalarField& v, const float rad, const int nb, const int seed = 2847573 ); 

   ~JitterSampleVolume(){}

    const float eval( const Vector& P ) const; 

    const Vector grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Jitter";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    float radius;
    int nbSamples;
    int prnSeed;

    const ScalarField elem;

    mutable Noise_t parms;
    mutable UniformPRN prn;
    PerlinNoise perlin;

};



class AdvectVolume : public Volume<float> 
{
  public:

    AdvectVolume( Volume<float>* v, Volume<Vector>* u, const float delt ); 

    AdvectVolume( const ScalarField& v, const VectorField& u, const float delt ); 

   ~AdvectVolume(){}

    const float eval( const Vector& P ) const; 
    const Vector grad(  const Vector& P ) const;

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

    const ScalarField elem;
    VectorField velocity;
    float dt;
};







/*
class GradStretchAdvectVolume : public Volume<float> 
{
  public:

    GradStretchAdvectVolume( Volume<float>* v, Volume<Vector>* u, Volume<Matrix>* gu,  const float delt, const int nb ); 

    GradStretchAdvectVolume( const ScalarField& v, const VectorField& u, const MatrixField& gu, const float delt, const int nb ); 

   ~GradStretchAdvectVolume(){}

    const float eval( const Vector& P ) const; 
    const Vector grad(  const Vector& P ) const;

    virtual std::string typelabel() 
    { 
       std::string lbl = "GradStretchAdvect";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ",";
       lbl = lbl + velocity->typelabel();
       lbl = lbl + ",";
       lbl = lbl + gradvelocity->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    const ScalarField elem;
    VectorField velocity;
    MatrixField gradvelocity;
    float dt;
    int nbIterations;
};
*/





/*
class BFECCAdvectVolume : public Volume<float> 
{
  public:

    BFECCAdvectVolume( const Volume<float>* v, const Volume<Vector>* u, const float dt, const int nb = 1 );
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

    const float eval( const Vector& P ) const; 
    {
       return elem->eval(P);
    }

    const Vector grad(  const Vector& P ) const;
    {
       return elem->grad(P); 
    }

  private:

    const ScalarField base;
    VectorField velocity;
    const ScalarField elem;
};
*/





class WarpVolume : public Volume<float> 
{
  public:

    WarpVolume( Volume<float>* v, Volume<Vector>* map );

    WarpVolume( const ScalarField& v, VectorField& map );

   ~WarpVolume(){}

    const float eval( const Vector& P ) const; 
    const Vector grad(  const Vector& P ) const;

    virtual std::string typelabel() 
    { 
       std::string lbl = "Warp";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ",";
       lbl = lbl + mapX->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    const ScalarField elem;
    VectorField mapX;
};

class PeriodicVolume : public Volume<float>
{
  public:

   PeriodicVolume( Volume<float>* v, const Vector& o, const Vector& L );

   PeriodicVolume( const ScalarField& v, const Vector& o, const Vector& L );

   ~PeriodicVolume(){}

   const float eval( const Vector& P ) const;
   
   const Vector grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "Periodic";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

   Vector origin;
   Vector periods;
   const ScalarField elem;

   const Vector periodP( const Vector& P ) const;

};


class DetGradMapVolume : public Volume<float> 
{
  public:

    DetGradMapVolume( Volume<Vector>* map, float d = 0.01 ) :
       elem(map)
    {
       dx = d;
    }

    DetGradMapVolume( const VectorField& map, float d = 0.01 ) :
       elem(map)
    {
       dx = d;
    }

   ~DetGradMapVolume(){}

    const float eval( const Vector& P ) const
    { 
       Matrix M = elem->grad(P);
       return det(M);
    }
    const Vector grad(  const Vector& P ) const
    {
      Vector value( (eval(P+Vector(dx,0,0))-eval(P-Vector(dx,0,0)))/dx,
                    (eval(P+Vector(0,dx,0))-eval(P-Vector(0,dx,0)))/dx,
                    (eval(P+Vector(0,0,dx))-eval(P-Vector(0,0,dx)))/dx    );
       return value*0.5;

    }

  private:

     float dx;
     VectorField elem;
};





class SwitchVolume : public Volume<float>
{
  public:

   SwitchVolume( Volume<float>* v1, Volume<float>* v2, Volume<float>* swtch );

   SwitchVolume( const ScalarField& v1, const ScalarField& v2, const ScalarField& swtch );

   ~SwitchVolume(){}

   const float eval( const Vector& P ) const;
   
   const Vector grad( const Vector& P ) const;
   
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

   const ScalarField elem1;
   const ScalarField elem2;
   const ScalarField swtchelem;


};




class SineVolume : public Volume<float>
{
  public:

   SineVolume( Volume<float>* v1 );

   SineVolume( const ScalarField& v1 );

   ~SineVolume(){}

   const float eval( const Vector& P ) const;
   
   const Vector grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "Sine";
      lbl = lbl + "(";
      lbl = lbl + elem1->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

   const ScalarField elem1;


};




class CosineVolume : public Volume<float>
{
  public:

   CosineVolume( Volume<float>* v1 );

   CosineVolume( const ScalarField& v1 );

   ~CosineVolume(){}

   const float eval( const Vector& P ) const;
   
   const Vector grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "Cosine";
      lbl = lbl + "(";
      lbl = lbl + elem1->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

   const ScalarField elem1;


};




class TangentVolume : public Volume<float>
{
  public:

   TangentVolume( Volume<float>* v1 );

   TangentVolume( const ScalarField& v1 );

   ~TangentVolume(){}

   const float eval( const Vector& P ) const;
   
   const Vector grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "Tangent";
      lbl = lbl + "(";
      lbl = lbl + elem1->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

   const ScalarField elem1;


};






class ArcsineVolume : public Volume<float>
{
  public:

   ArcsineVolume( Volume<float>* v1 );

   ArcsineVolume( const ScalarField& v1 );

   ~ArcsineVolume(){}

   const float eval( const Vector& P ) const;
   
   const Vector grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "Arcsine";
      lbl = lbl + "(";
      lbl = lbl + elem1->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

   const ScalarField elem1;


};



class ArccosineVolume : public Volume<float>
{
  public:

   ArccosineVolume( Volume<float>* v1 );

   ArccosineVolume( const ScalarField& v1 );

   ~ArccosineVolume(){}

   const float eval( const Vector& P ) const;
   
   const Vector grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "Arccosine";
      lbl = lbl + "(";
      lbl = lbl + elem1->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

   const ScalarField elem1;


};


class ArctangentVolume : public Volume<float>
{
  public:

   ArctangentVolume( Volume<float>* v1 );

   ArctangentVolume( const ScalarField& v1 );

   ~ArctangentVolume(){}

   const float eval( const Vector& P ) const;
   
   const Vector grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "Arctangent";
      lbl = lbl + "(";
      lbl = lbl + elem1->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

   const ScalarField elem1;


};



class HyperbolicSineVolume : public Volume<float>
{
  public:

   HyperbolicSineVolume( Volume<float>* v1 );

   HyperbolicSineVolume( const ScalarField& v1 );

   ~HyperbolicSineVolume(){}

   const float eval( const Vector& P ) const;
   
   const Vector grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "HyperbolicSine";
      lbl = lbl + "(";
      lbl = lbl + elem1->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

   const ScalarField elem1;


};



class HyperbolicCosineVolume : public Volume<float>
{
  public:

   HyperbolicCosineVolume( Volume<float>* v1 );

   HyperbolicCosineVolume( const ScalarField& v1 );

   ~HyperbolicCosineVolume(){}

   const float eval( const Vector& P ) const;
   
   const Vector grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "HyperbolicCosine";
      lbl = lbl + "(";
      lbl = lbl + elem1->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

   const ScalarField elem1;


};



class HyperbolicTangentVolume : public Volume<float>
{
  public:

   HyperbolicTangentVolume( Volume<float>* v1 );

   HyperbolicTangentVolume( const ScalarField& v1 );

   ~HyperbolicTangentVolume(){}

   const float eval( const Vector& P ) const;
   
   const Vector grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "HyperbolicTangent";
      lbl = lbl + "(";
      lbl = lbl + elem1->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

   const ScalarField elem1;


};




class XIdentityVolume : public Volume<float>
{
  public:

   XIdentityVolume();

   ~XIdentityVolume(){}

   const float eval( const Vector& P ) const;
   
   const Vector grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "XIdentity";
      lbl = lbl + "(";
      lbl = lbl + ")";
      return lbl;
   }


  private:

   const Vector gradvalue;

};



class YIdentityVolume : public Volume<float>
{
  public:

   YIdentityVolume();

   ~YIdentityVolume(){}

   const float eval( const Vector& P ) const;
   
   const Vector grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "YIdentity";
      lbl = lbl + "(";
      lbl = lbl + ")";
      return lbl;
   }


  private:

   const Vector gradvalue;

};



class ZIdentityVolume : public Volume<float>
{
  public:

   ZIdentityVolume();

   ~ZIdentityVolume(){}

   const float eval( const Vector& P ) const;
   
   const Vector grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "ZIdentity";
      lbl = lbl + "(";
      lbl = lbl + ")";
      return lbl;
   }


  private:

   const Vector gradvalue;

};





}





#endif



