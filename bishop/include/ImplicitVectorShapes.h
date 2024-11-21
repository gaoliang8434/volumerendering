
#ifndef __IMPLICITVECTORSHAPES_H__
#define __IMPLICITVECTORSHAPES_H__

#include "Volume.h"
#include "LinearAlgebra.h"
#include "Noise.h"
#include "NoiseMachine.h"
#include "FFTNoise.h"
#include "UniformPRN.h"
#include <cmath>
#include <vector>
#include <iostream>
#include "VolumeGrid.h"
#include "RectangularGrid.h"
//#include "OVDBGrid.h"
#include "PerlinNoise.h"


namespace lux
{



class ConstantVectorVolume : public Volume<Vector> 
{
  public:

    ConstantVectorVolume( const Vector v ) ;

   ~ConstantVectorVolume(){}


    const Vector eval( const Vector& P ) const; 

    const Matrix grad( const Vector& P ) const;

    virtual std::string typelabel() 
    { 
       std::string lbl = "Constant";
       return lbl;
    }


  private:

    Vector value;
    Matrix gradvalue;
};


class ComponentVectorVolume : public Volume<Vector> 
{
  public:

    ComponentVectorVolume( Volume<float>* x, Volume<float>* y, Volume<float>* z ) ;
    ComponentVectorVolume( const ScalarField& x, const ScalarField& y, const ScalarField& z ) ;

   ~ComponentVectorVolume(){}


    const Vector eval( const Vector& P ) const; 

    const Matrix grad( const Vector& P ) const;

    virtual std::string typelabel() 
    { 
       std::string lbl = "Component";
       lbl = lbl + "(";
       lbl = lbl + X->typelabel();
       lbl = lbl + ",";
       lbl = lbl + Y->typelabel();
       lbl = lbl + ",";
       lbl = lbl + Z->typelabel();
       lbl = lbl + ")";
       return lbl;
    }


  private:

    ScalarField X;
    ScalarField Y;
    ScalarField Z;
};



class XVectorVolume : public Volume<float> 
{
  public:

    XVectorVolume( Volume<Vector>* v ) ;
    XVectorVolume( const VectorField& v ) ;

   ~XVectorVolume(){}


    const float eval( const Vector& P ) const; 

    const Vector grad( const Vector& P ) const;

    virtual std::string typelabel() 
    { 
       std::string lbl = "XComponent";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }


  private:

    VectorField elem;
};



class YVectorVolume : public Volume<float> 
{
  public:

    YVectorVolume( Volume<Vector>* v ) ;
    YVectorVolume( const VectorField& v ) ;

   ~YVectorVolume(){}


    const float eval( const Vector& P ) const; 

    const Vector grad( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "YComponent";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }



  private:

    VectorField elem;
};



class ZVectorVolume : public Volume<float> 
{
  public:

    ZVectorVolume( Volume<Vector>* v ) ;
    ZVectorVolume( const VectorField& v ) ;

   ~ZVectorVolume(){}


    const float eval( const Vector& P ) const; 

    const Vector grad( const Vector& P ) const;

    virtual std::string typelabel() 
    { 
       std::string lbl = "ZComponent";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }



  private:

    VectorField elem;
};








class TranslateVectorVolume : public Volume<Vector> 
{
  public:

    TranslateVectorVolume( Volume<Vector> * v, const Vector& s ) ;

    TranslateVectorVolume( const VectorField& v, const Vector& s ) ;

    ~TranslateVectorVolume(){}


    const Vector eval( const Vector& P ) const;

    const Matrix grad( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Translate";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }


  private:

    VectorField elem;
    Vector translate;
};



class ScaleVectorVolume : public Volume<Vector> 
{
  public:

    ScaleVectorVolume( Volume<Vector> * v, const Vector& s ) ;

    ScaleVectorVolume( const VectorField& v, const Vector& s ) ;

    ScaleVectorVolume( Volume<Vector> * v, const float& s ) ;

    ScaleVectorVolume( const VectorField& v, const float& s ) ;



    ~ScaleVectorVolume(){}


    const Vector eval( const Vector& P ) const;

    const Matrix grad( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Scale";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }


  private:

    VectorField elem;
    Vector scale;
};







class RotateVectorVolume: public Volume<Vector> 
{
  public:

    RotateVectorVolume( Volume<Vector> * v, const Vector& s ) ;

    RotateVectorVolume( const VectorField& v, const Vector& s ) ;

    ~RotateVectorVolume(){}


    const Vector eval( const Vector& P ) const;

    const Matrix grad( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Rotate";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }


  private:

    VectorField elem;
    Matrix R, Rinv;
};






class NegateVectorVolume : public Volume<Vector> 
{
  public:

    NegateVectorVolume( Volume<Vector> * v ) ;

    NegateVectorVolume( const VectorField& v ) ;

    ~NegateVectorVolume(){}


    const Vector eval( const Vector& P ) const;

    const Matrix grad( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Negate";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    VectorField elem;
};

class DotProductVectorVolume : public Volume<float>
{
  public:

    DotProductVectorVolume( Volume<Vector>* v1, Volume<Vector>* v2 ) ;

    DotProductVectorVolume( const VectorField& v1, const VectorField& v2 ) ;

    ~DotProductVectorVolume(){}

    const float eval( const Vector& P ) const;

    const Vector grad( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "DotProduct";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }


  private:

    VectorField elem1, elem2;
};

class MagnitudeVectorVolume : public Volume<float>
{
  public:

    MagnitudeVectorVolume( Volume<Vector>* v1 ) ;

    MagnitudeVectorVolume( const VectorField& v1 ) ;

    ~MagnitudeVectorVolume(){}

    const float eval( const Vector& P ) const;

    const Vector grad( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Magnitude";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    VectorField elem;
};


class UnitVectorVolume : public Volume<Vector>
{
  public:

    UnitVectorVolume( Volume<Vector>* v1 ) ;

    UnitVectorVolume( const VectorField& v1 ) ;

    ~UnitVectorVolume(){}

    const Vector eval( const Vector& P ) const;

    const Matrix grad( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "UnitVector";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    VectorField elem;
};




class MultiplyVectorVolume : public Volume<Vector> 
{
  public:

    MultiplyVectorVolume( Volume<Vector> * v, const float a );

    MultiplyVectorVolume( Volume<Vector> * v, Volume<float>* u );

    MultiplyVectorVolume( const VectorField& v, const float a );

    MultiplyVectorVolume( const VectorField& v, const ScalarField& u );


    ~MultiplyVectorVolume(){}


    const Vector eval( const Vector& P ) const;

    const Matrix grad( const Vector& P ) const;


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

    VectorField elem;
    ScalarField factor;
    float constant;
};



class DivideVectorVolume : public Volume<Vector> 
{
  public:

    DivideVectorVolume( Volume<Vector> * v, const float a );

    DivideVectorVolume( Volume<Vector> * v, Volume<float>* u );

    DivideVectorVolume( const VectorField& v, const float a );

    DivideVectorVolume( const VectorField& v, const ScalarField& u );


    ~DivideVectorVolume(){}


    const Vector eval( const Vector& P ) const;

    const Matrix grad( const Vector& P ) const;


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

    VectorField elem;
    ScalarField factor;
    float constant;
};




class AddVectorVolume : public Volume<Vector> 
{
  public:

    AddVectorVolume( Volume<Vector> * v1, Volume<Vector> * v2 ) ;

    AddVectorVolume( const VectorField& v1, const VectorField& v2 ) ;

    ~AddVectorVolume(){}


    const Vector eval( const Vector& P ) const;

    const Matrix grad( const Vector& P ) const;


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
  
    VectorField elem1, elem2;
};


class SubtractVectorVolume : public Volume<Vector> 
{
  public:

    SubtractVectorVolume( Volume<Vector> * v1, Volume<Vector> * v2 ) ;

    SubtractVectorVolume( const VectorField& v1, const VectorField& v2 ) ;

    ~SubtractVectorVolume(){}


    const Vector eval( const Vector& P ) const;

    const Matrix grad( const Vector& P ) const;


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
  
    VectorField elem1, elem2;
};

class IdentityVectorVolume : public Volume<Vector> 
{
  public:

    IdentityVectorVolume();

    ~IdentityVectorVolume(){}


    const Vector eval( const Vector& P ) const ;

    const Matrix grad( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Identity";
       return lbl;
    }


  private:

    Matrix gradvalue;
};



class ImplicitPointVectorVolume : public Volume<Vector> 
{
  public:

    ImplicitPointVectorVolume( Volume<float>* v, const float dx, const int nb = 3, const float thresh = 1.0e6 ) ;

    ImplicitPointVectorVolume( const ScalarField& v, const float dx, const int nb = 3, const float thresh = 1.0e6 ) ;

    ~ImplicitPointVectorVolume(){}


    const Vector eval( const Vector& P ) const;

    // evaluate gradient by finite difference
    // size of difference is controlled via the "dx" attribute
    const Matrix grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "ImplicitPoint";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ScalarField elem;
    float step;
    int nbIterations;
    float threshold;
};


class GradientVectorVolume : public Volume<Vector> 
{
  public:

    GradientVectorVolume( Volume<float>* v, const float dx = 0.001 ) ;

    GradientVectorVolume( const ScalarField& v, const float dx = 0.001 ) ;

    ~GradientVectorVolume(){}

    const Vector eval( const Vector& P ) const;

    // evaluate gradient by finite difference
    // size of difference is controlled via the "dx" attribute
    const Matrix grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Grad";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ScalarField elem;
    float step;
};

class FiniteDifferenceGradientVectorVolume : public Volume<Vector> 
{
  public:

    // N is the length of the difference kernel
    FiniteDifferenceGradientVectorVolume( Volume<float>* v, const int N, const float dx, const float dy, const float dz ) ;

    FiniteDifferenceGradientVectorVolume( const ScalarField& v, const int N, const float dx, const float dy, const float dz ) ;

    ~FiniteDifferenceGradientVectorVolume(){}

    const Vector eval( const Vector& P ) const;

    const Matrix grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "FDGrad";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ScalarField elem;
    float dX, dY, dZ;
    vector<double> coefficients;
};


class FiniteDifferenceBoundedGradientVectorVolume : public Volume<Vector> 
{
  public:

    // N is the length of the difference kernel
    FiniteDifferenceBoundedGradientVectorVolume( Volume<float>* v, const int N,  const GridBox& gb ) ;

    FiniteDifferenceBoundedGradientVectorVolume( const ScalarField& v, const int N, const GridBox& gb ) ;

    ~FiniteDifferenceBoundedGradientVectorVolume(){}

    const Vector eval( const Vector& P ) const;

    const Matrix grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "FDBoundedGrad";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ScalarField elem;
    int maxN;
    vector< vector<double> > coefficients;
    GridBox gB;
};





class FiniteDifferenceInteriorGradientVectorVolume : public Volume<Vector> 
{
  public:

    // N is the length of the difference kernel
    FiniteDifferenceInteriorGradientVectorVolume( Volume<float>* v, const int N, const double dx, const double dy, const double dz,  Volume<float>* gb ) ;

    FiniteDifferenceInteriorGradientVectorVolume( const ScalarField& v, const int N, const double dx, const double dy, const double dz, const ScalarField& gb ) ;

    ~FiniteDifferenceInteriorGradientVectorVolume(){}

    const Vector eval( const Vector& P ) const;

    const Matrix grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "FDInteriorGrad";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    ScalarField elem;
    int maxN;
    double dX, dY, dZ;
    ScalarField gB;
    vector< vector<double> > coefficients;
};






class CrossProductVectorVolume : public Volume<Vector> 
{
  public:

    CrossProductVectorVolume( Volume<Vector>* v1, Volume<Vector>* v2 ) ;

    CrossProductVectorVolume( const VectorField& v1, const VectorField& v2 ) ;

    ~CrossProductVectorVolume(){}

    const Vector eval( const Vector& P ) const ;

    // evaluate gradient by finite difference
    // size of difference is controlled via the "dx" attribute
    const Matrix grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "CrossProduct";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    VectorField elem1, elem2;
    float step;
};







class DivergenceVectorVolume : public Volume<float> 
{
  public:

    DivergenceVectorVolume( Volume<Vector>* v, const float dx = 0.001 ) ;

    DivergenceVectorVolume( const VectorField& v, const float dx = 0.001 ) ;

    ~DivergenceVectorVolume(){}

    const float eval( const Vector& P ) const ;

    // evaluate gradient by finite difference
    // size of difference is controlled via the "dx" attribute
    const Vector grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Div";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    VectorField elem;
    float step;
};



class FiniteDifferenceDivergenceVectorVolume : public Volume<float> 
{
  public:

    FiniteDifferenceDivergenceVectorVolume( Volume<Vector>* v, const int N, const float dx, const float dy, const float dz ) ;

    FiniteDifferenceDivergenceVectorVolume( const VectorField& v, const int N, const float dx, const float dy, const float dz ) ;

    ~FiniteDifferenceDivergenceVectorVolume(){}

    const float eval( const Vector& P ) const ;

    // evaluate gradient by finite difference
    // size of difference is controlled via the "dx" attribute
    const Vector grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "FDDiv";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    VectorField elem;
    float dX, dY, dZ;
    vector<double> coefficients;

};



class FiniteDifferenceBoundedDivergenceVectorVolume : public Volume<float> 
{
  public:

    FiniteDifferenceBoundedDivergenceVectorVolume( Volume<Vector>* v, const int N, const GridBox& gb ) ;

    FiniteDifferenceBoundedDivergenceVectorVolume( const VectorField& v, const int N, const GridBox& gb ) ;

    ~FiniteDifferenceBoundedDivergenceVectorVolume(){}

    const float eval( const Vector& P ) const ;

    // evaluate gradient by finite difference
    // size of difference is controlled via the "dx" attribute
    const Vector grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "FDBoundedDiv";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    VectorField elem;
    int maxN;
    GridBox gB;
    vector< vector<double> > coefficients;

};





class FiniteDifferenceInteriorDivergenceVectorVolume : public Volume<float> 
{
  public:

    FiniteDifferenceInteriorDivergenceVectorVolume( Volume<Vector>* v, const int N, const double dx, const double dy, const double dz, Volume<float>* gb ) ;

    FiniteDifferenceInteriorDivergenceVectorVolume( const VectorField& v, const int N, const double dx, const double dy, const double dz, const ScalarField& gb ) ;

    ~FiniteDifferenceInteriorDivergenceVectorVolume(){}

    const float eval( const Vector& P ) const ;

    // evaluate gradient by finite difference
    // size of difference is controlled via the "dx" attribute
    const Vector grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "FDInteriorDiv";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    VectorField elem;
    int maxN;
    double dX, dY, dZ;
    ScalarField gB;
    vector< vector<double> > coefficients;

};










class CurlVectorVolume : public Volume<Vector> 
{
  public:

    CurlVectorVolume( Volume<Vector>* v1, const float dx = 0.001 ) ;

    CurlVectorVolume( const VectorField& v1, const float dx = 0.001 ) ;

    ~CurlVectorVolume(){}

    const Vector eval( const Vector& P ) const ;

    // evaluate gradient by finite difference
    // size of difference is controlled via the "dx" attribute
    const Matrix grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Curl";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

     VectorField elem;
     float step;
};






class ReportVectorVolume : public Volume<Vector> 
{
  public:

    ReportVectorVolume( Volume<Vector>* v, const string tag ) ;

    ReportVectorVolume( const VectorField& v, const string tag ) ;

    ~ReportVectorVolume(){}

    const Vector eval( const Vector& P ) const ;

    const Matrix grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Report";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    VectorField elem;
    string label;
};








class AdvectVectorVolume : public Volume<Vector> 
{
  public:

    AdvectVectorVolume( Volume<Vector>* v, Volume<Vector>* u, const float delt ) ;

    AdvectVectorVolume( const VectorField& v, const VectorField& u, const float delt ) ;

   ~AdvectVectorVolume(){}

    const Vector eval( const Vector& P ) const ;
    const Matrix grad( const Vector& P ) const ;


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

    VectorField elem, velocity;
    float dt;
};

/*
class BFECCAdvectVectorVolume : public Volume<Vector> 
{
  public:

    BFECCAdvectVectorVolume( Volume<Vector>* v, Volume<Vector>* u, const float dt, const int nb = 1 )  
    {
       if( nb == 0 )
       {
	  elem = VectorField( new AdvectVectorVolume( v, u, dt )  );
       }
       else
       {
          Volume<Vector>* advected = new BFECCAdvectVectorVolume( v, u, dt, nb-1 );
	  Volume<Vector>* negu = new NegateVectorVolume( u );
	  Volume<Vector>* backadvected = new BFECCAdvectVectorVolume( advected, negu, dt, nb-1 );
	  Volume<Vector>* error = new MultiplyVectorVolume( new SubtractVectorVolume( v, backadvected ) , 0.5);
	  advected = new AddVectorVolume( advected, error );
	  elem = VectorField(advected);
       }
    }

   ~BFECCAdvectVectorVolume(){}

    const Vector eval( const Vector& P ) const 
    {
       return elem->eval(P);
    }

    const Matrix grad( const Vector& P ) const
    {
       return elem->grad(P); 
    }

  private:

    VectorField elem;
};
*/

class DisplaceVectorVolume : public Volume<Vector> 
{
  public:

    DisplaceVectorVolume( Volume<Vector>* v, Volume<Vector>* u ) ;

    DisplaceVectorVolume( const VectorField& v, const VectorField& u ) ;

   ~DisplaceVectorVolume(){}

    const Vector eval( const Vector& P ) const ;
    const Matrix grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Displace";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    VectorField elem, disp;
};



class ContinuedFractionDisplacementVectorVolume : public Volume<Vector> 
{
  public:

    ContinuedFractionDisplacementVectorVolume( Volume<Vector>* u, int _iterations ) ;

    ContinuedFractionDisplacementVectorVolume( const VectorField& u, int _iterations ) ;

   ~ContinuedFractionDisplacementVectorVolume(){}

    const Vector eval( const Vector& P ) const ;

    const Matrix grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "ContinutedFractionDisplacement";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    VectorField elem;
    int iterations;
};








class WarpVectorVolume : public Volume<Vector> 
{
  public:

    WarpVectorVolume( Volume<Vector>* v, Volume<Vector>* u ) ;

    WarpVectorVolume( const VectorField& v, const VectorField& u ) ;

   ~WarpVectorVolume(){}

    const Vector eval( const Vector& P ) const ;
    const Matrix grad( const Vector& P ) const ;


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

    VectorField elem, warp;
};













class NoiseVectorVolume : public Volume<Vector> 
{
  public:

    NoiseVectorVolume( Noise* n, const float d = 0.01 ) ; 
    NoiseVectorVolume( NoiseMachine n, const float d = 0.01 ) ; 

   ~NoiseVectorVolume(){}

    const Vector eval( const Vector& P ) const;
    const Matrix grad( const Vector& P ) const;  


    virtual std::string typelabel() 
    { 
       std::string lbl = "Noise";
       return lbl;
    }

  private:

    NoiseMachine noise;
    float dx;
};


class NoiseSampleVectorVolume : public Volume<Vector> 
{
  public:

    NoiseSampleVectorVolume( Noise* n, const float d = 0.01 ) ;
    NoiseSampleVectorVolume( NoiseMachine n, const float d = 0.01 ) ;

   ~NoiseSampleVectorVolume(){}

    const Vector eval( const Vector& P ) const;
    const Matrix grad( const Vector& P ) const;  


    virtual std::string typelabel() 
    { 
       std::string lbl = "NoiseSample";
       return lbl;
    }

  private:

    NoiseMachine noise;
    float dx;
};




class GriddedVectorVolume : public Volume<Vector> 
{
  public:

    GriddedVectorVolume( const VolumeGrid<Vector>* g ) ;

    GriddedVectorVolume( const SparseVectorGrid* g ) ;

   ~GriddedVectorVolume(){}

    const Vector eval( const Vector& P ) const ;
    const Matrix grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Gridded";
       return lbl;
    }

  private:

    const VolumeGrid<Vector>* grid;
    const SparseVectorGrid * sparsegrid;
    float dx, dy, dz;
};




class GriddedSGridVectorVolume : public Volume<Vector> 
{
  public:

    GriddedSGridVectorVolume( const VectorGrid& g ) ;

   ~GriddedSGridVectorVolume(){}

    const Vector eval( const Vector& P ) const ;
    const Matrix grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Gridded";
       return lbl;
    }

  private:

    VectorGrid elem;
    float dx, dy, dz;
};




class GriddedFrustumVectorVolume : public Volume<Vector> 
{
  public:

    GriddedFrustumVectorVolume( const VectorFrustumGrid& g ) ;

   ~GriddedFrustumVectorVolume(){}

    const Vector eval( const Vector& P ) const ;
    const Matrix grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Gridded";
       return lbl;
    }

  private:

    VectorFrustumGrid felem;
    float dx, dy, dz;
};


/*
class GriddedOGridVectorVolume : public Volume<Vector> 
{
  public:

    GriddedOGridVectorVolume( const VectorOGrid& g ) ;

   ~GriddedOGridVectorVolume(){}

    const Vector eval( const Vector& P ) const ;
    const Matrix grad( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Gridded";
       return lbl;
    }

  private:

    VectorOGrid elem;
    float dx, dy, dz;
};
*/


class PeriodicVectorVolume : public Volume<Vector>
{
  public:

   PeriodicVectorVolume( Volume<Vector>* v, const Vector& o, const Vector& L ) ;
   PeriodicVectorVolume( const VectorField& v, const Vector& o, const Vector& L ) ;

   ~PeriodicVectorVolume(){}

   const Vector eval( const Vector& P ) const;
   
   const Matrix grad( const Vector& P ) const;
  
  
    virtual std::string typelabel() 
    { 
       std::string lbl = "Periodic";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }


  private:

   VectorField elem;
   Vector origin;
   Vector periods;

   const Vector periodP( const Vector& P ) const;

};


class SwitchVectorVolume : public Volume<Vector>
{
  public:

   SwitchVectorVolume( Volume<Vector>* v1, Volume<Vector>* v2, Volume<float>* swtch );

   SwitchVectorVolume( const VectorField& v1, const VectorField& v2, const ScalarField& swtch );

   ~SwitchVectorVolume(){}

   const Vector eval( const Vector& P ) const;
   
   const Matrix grad( const Vector& P ) const;
  
  
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

   const VectorField elem1;
   const VectorField elem2;
   const ScalarField swtchelem;


};



class GradientStretchCMVolume: public Volume<Vector>
{
   public:
      GradientStretchCMVolume( const VectorField& vel, const float T, const int nbiter ) :
          elem         (vel),
	  totalTime    (T),
	  nbiterations (nbiter)
	  {}

      ~GradientStretchCMVolume(){}

       const Vector eval( const Vector& P ) const ;
       const Matrix grad( const Vector& P ) const ;

       virtual std::string typelabel()
       {
           std::string lbl = "GradientStretchCM";
           lbl = lbl + "(";
           lbl = lbl + elem->typelabel();
           lbl = lbl + ",";
           lbl = lbl + "float";
           lbl = lbl + ",";
           lbl = lbl + "int";
           lbl = lbl + ")";
           return lbl;
       }

     private:
        VectorField elem;
	float totalTime;
	int nbiterations;
};





class FFTDivFreeNoiseVectorVolume: public Volume<Vector>
{
   public:
      FFTDivFreeNoiseVectorVolume( const float powerLaw, const float largeScale, const float smallScale, const float length, const int dim, const int seed ); 

      ~FFTDivFreeNoiseVectorVolume(){}

       const Vector eval( const Vector& P ) const ;
       const Matrix grad( const Vector& P ) const ;

       virtual std::string typelabel()
       {
           std::string lbl = "FFTDivFreeNoiseVector";
           lbl = lbl + "(";
           lbl = lbl + "float";
           lbl = lbl + ",";
           lbl = lbl + "float";
           lbl = lbl + ",";
           lbl = lbl + "float";
           lbl = lbl + ",";
           lbl = lbl + "float";
           lbl = lbl + ",";
           lbl = lbl + "int";
           lbl = lbl + ",";
           lbl = lbl + "int";
           lbl = lbl + ")";
           return lbl;
       }

     private:

	VolumeGrid<Vector> data;
	float dx;
       
};




class GradDisplaceVectorVolume: public Volume<Vector> 
{
  public:

    GradDisplaceVectorVolume( Volume<Matrix> * m, Volume<float>* ls, const double step ) ;

    GradDisplaceVectorVolume( const MatrixField& m, const ScalarField& ls, const double step ) ;

    ~GradDisplaceVectorVolume(){}


    const Vector eval( const Vector& P ) const;

    const Matrix grad( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "GradDisplace";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ", ";
       lbl = lbl + sdf->typelabel();
       lbl = lbl + ", float";
       lbl = lbl + ")";
       return lbl;
    }


  private:

    MatrixField elem;
    ScalarField sdf;
    double ds;
};




/*
class RayMarchDivFreeZeroNormalVectorVolume : public Volume<Vector>
{
  public:
    RayMarchDivFreeZeroNormalVectorVolume( const VectorField& W, const ScalarField& LS, const float resolution, const VectorField& BC );
    RayMarchDivFreeZeroNormalVectorVolume( Volume<Vector>* W, Volume<float>* LS, const float resolution, Volume<Vector>* BC );
   ~RayMarchDivFreeZeroNormalVectorVolume(){}


    const Vector eval( const Vector& P ) const;



    virtual std::string typelabel() 
    { 
       std::string lbl = "RayMarchDivFreeZeroNormal";
       lbl = lbl + "(";
       lbl = lbl + W->typelabel();
       lbl = lbl + ", ";
       lbl = lbl + sdf->typelabel();
       lbl = lbl + ", float";
       lbl = lbl + ")";
       return lbl;
    }


  private:

    VectorField W;
    ScalarField sdf;
    double ds;
    VectorField boundary_vel;
    ScalarField divW;
    VectorField normal;
    ScalarField curvature;
    VectorField CPT;

};


class RayMarchDivFreePlanarVectorVolume : public Volume<Vector>
{
  public:
    RayMarchDivFreePlanarVectorVolume( const VectorField& W, const Vector& planeP, const Vector& norm, const float resolution, const VectorField& BC );
    RayMarchDivFreePlanarVectorVolume( Volume<Vector>* W, const Vector& planeP, const Vector& norm, const float resolution, Volume<Vector>* BC );
   ~RayMarchDivFreePlanarVectorVolume(){}


    const Vector eval( const Vector& P ) const;



    virtual std::string typelabel() 
    { 
       std::string lbl = "RayMarchDivFreePlanar";
       lbl = lbl + "(";
       lbl = lbl + W->typelabel();
       lbl = lbl + ", Vector, Vector, float, ";
       lbl = lbl + boundary_vel->typelabel();
       lbl = lbl + ")";
       return lbl;
    }


  private:

    VectorField W;
    Vector Pplane;
    Vector normal;
    double ds;
    VectorField boundary_vel;

    Vector xhat, yhat, zhat;
    Vector xperphat, yperphat, zperphat;

};


*/








}
#endif
