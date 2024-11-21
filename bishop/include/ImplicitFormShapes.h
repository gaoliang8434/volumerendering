
#ifndef __IMPLICITFORMSHAPES_H__
#define __IMPLICITFORMSHAPES_H__

#include "Volume.h"
#include "LinearAlgebra.h"
#include "Noise.h"
#include "FFTNoise.h"
#include "UniformPRN.h"
#include "VolumeGrid.h"
#include "SparseGrid.h"
#include "ImplicitVolumeShapes.h"
#include "ImplicitVectorShapes.h"
#include <cmath>
#include <vector>
#include <iostream>
#include "PerlinNoise.h"



namespace lux
{


class ConstantFormVolume : public Volume<Form> 
{
  public:

    ConstantFormVolume( const Form v ) ;

   ~ConstantFormVolume(){}


    const Form eval( const Vector& P ) const; 

    const Form grad( const Vector& P ) const; 

   virtual std::string typelabel() { return "Constant"; }


  private:

    Form value;
    Form gradvalue;
};



class ReportFormVolume: public Volume<Form> 
{
  public:

    ReportFormVolume( Volume<Form>* v, const string tag );

    ReportFormVolume( const FormField& v, const string tag );

   ~ReportFormVolume(){}


    const Form eval( const Vector& P ) const;
    const Form grad( const Vector& P ) const;


    virtual std::string typelabel() 
    {
       std::string lbl = "Report";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    const FormField elem;
    string label;

};




class ScaleFormVolume : public Volume<Form> 
{
  public:

    ScaleFormVolume( Volume<Form> * v, const Vector& s );

    ScaleFormVolume( Volume<Form> * v, const float& s );


    ScaleFormVolume( const FormField& v, const Vector& s );

    ScaleFormVolume( const FormField& v, const float& s );


    ~ScaleFormVolume(){}


    const Form eval( const Vector& P ) const;

    const Form grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Scale";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const FormField elem;
    Vector scale;
};



class TranslateFormVolume : public Volume<Form> 
{
  public:

    TranslateFormVolume( Volume<Form> * v, const Vector& s );

    TranslateFormVolume( const FormField& v, const Vector& s );

    ~TranslateFormVolume(){}


    const Form eval( const Vector& P ) const;


    const Form grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Translate";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }



  private:

    const FormField elem;
    Vector translate;
};


class RotateFormVolume: public Volume<Form> 
{
  public:

    RotateFormVolume( Volume<Form> * v, const Vector& s );

    RotateFormVolume( const FormField& v, const Vector& s );

    ~RotateFormVolume(){}


    const Form eval( const Vector& P ) const;

    const Form grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Rotate";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

    const FormField elem;
    Vector axis;
    double sina, cosa;
    Matrix R;
};






class NegateFormVolume : public Volume<Form> 
{
  public:

    NegateFormVolume( Volume<Form> * v );

    NegateFormVolume( const FormField& v );

    ~NegateFormVolume(){}


    const Form eval( const Vector& P ) const;

    const Form grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Negate";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:
  
    const FormField elem;
};




class GradientFormVolume : public Volume<Form> 
{
  public:

    GradientFormVolume( Volume<Form> * v );

    GradientFormVolume( const FormField& v );

    ~GradientFormVolume(){}


    const Form eval( const Vector& P ) const;

    const Form grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Grad";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:
  
    const FormField elem;
};





class MultiplyFormVolume : public Volume<Form> 
{
  public:

    MultiplyFormVolume( Volume<Form> * v, const float a ); 

    MultiplyFormVolume( Volume<Form> * v, Volume<float>* u ); 

    MultiplyFormVolume( const FormField& v, const float a ); 

    MultiplyFormVolume( const FormField& v, const ScalarField& u ); 


    ~MultiplyFormVolume(){}


    const Form eval( const Vector& P ) const;

    const Form grad(  const Vector& P ) const;

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

    const FormField elem;
    const ScalarField factor;
    float constant;
};





class DivideFormVolume : public Volume<Form> 
{
  public:

    DivideFormVolume( Volume<Form> * v, const float a ); 

    DivideFormVolume( Volume<Form> * v, Volume<float>* u ); 

    DivideFormVolume( const FormField& v, const float a ); 

    DivideFormVolume( const FormField& v, const ScalarField& u ); 


    ~DivideFormVolume(){}


    const Form eval( const Vector& P ) const;

    const Form grad(  const Vector& P ) const;

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

    const FormField elem;
    const ScalarField factor;
    float constant;
};




class AddFormVolume : public Volume<Form> 
{
  public:

    AddFormVolume( Volume<Form> * v1, Volume<Form> * v2 );

    AddFormVolume( const FormField&  v1, const FormField& v2 );

    ~AddFormVolume(){}


    const Form eval( const Vector& P ) const;

    const Form grad(  const Vector& P ) const;

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

    const FormField e1, e2;
};


class SubtractFormVolume : public Volume<Form> 
{
  public:

    SubtractFormVolume( Volume<Form> * v1, Volume<Form> * v2 );

    SubtractFormVolume( const FormField& v1, const FormField& v2 );

    ~SubtractFormVolume(){}


    const Form eval( const Vector& P ) const;

    const Form grad(  const Vector& P ) const;

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

    const FormField e1, e2;
};



class WedgeFormVolume : public Volume<Form> 
{
  public:

    WedgeFormVolume( Volume<Form> * v1, Volume<Form> * v2 );

    WedgeFormVolume( const FormField& v1, const FormField& v2 );

    ~WedgeFormVolume(){}


    const Form eval( const Vector& P ) const;

    const Form grad(  const Vector& P ) const;

   virtual std::string typelabel() 
   { 
      std::string lbl = "Wedge";
      lbl = lbl + "(";
      lbl = lbl + elem1->typelabel();
      lbl = lbl + ",";
      lbl = lbl + elem2->typelabel();
      lbl = lbl + ")";
      return lbl;
   }

  private:

    const FormField elem1, elem2;
};


class AdvectFormVolume : public Volume<Form> 
{
  public:

    AdvectFormVolume( Volume<Form>* v, Volume<Vector>* u, const float delt ); 

    AdvectFormVolume( const FormField& v, const VectorField& u, const float delt ); 

   ~AdvectFormVolume(){}

    const Form eval( const Vector& P ) const; 
    const Form grad(  const Vector& P ) const;

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

    const FormField elem;
    VectorField velocity;
    float dt;
};





class WarpFormVolume : public Volume<Form> 
{
  public:

    WarpFormVolume( Volume<Form>* v, Volume<Vector>* map );

    WarpFormVolume( const FormField& v, VectorField& map );

   ~WarpFormVolume(){}

    const Form eval( const Vector& P ) const; 
    const Form grad(  const Vector& P ) const;

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

    const FormField elem;
    VectorField mapX;
};


class PeriodicFormVolume : public Volume<Form>
{
  public:

   PeriodicFormVolume( Volume<Form>* v, const Vector& o, const Vector& L );

   PeriodicFormVolume( const FormField& v, const Vector& o, const Vector& L );

   ~PeriodicFormVolume(){}

   const Form eval( const Vector& P ) const;
   
   const Form grad( const Vector& P ) const;
   
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
   const FormField elem;

   const Vector periodP( const Vector& P ) const;

};




class SwitchFormVolume : public Volume<Form>
{
  public:

   SwitchFormVolume( Volume<Form>* v1, Volume<Form>* v2, Volume<float>* swtch );

   SwitchFormVolume( const FormField& v1, const FormField& v2, const ScalarField& swtch );

   ~SwitchFormVolume(){}

   const Form eval( const Vector& P ) const;
   
   const Form grad( const Vector& P ) const;
   
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

   const FormField elem1;
   const FormField elem2;
   const ScalarField swtchelem;
};




class ZeroFormVolume : public Volume<float>
{
  public:

   ZeroFormVolume( Volume<Form>* v1 );

   ZeroFormVolume( const FormField& v1 );

   ~ZeroFormVolume(){}

   const float eval( const Vector& P ) const;
   
   const Vector grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "ZeroForm";
      lbl = lbl + "(";
      lbl = lbl + elem1->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

   const FormField elem1;
};



class ThreeFormVolume : public Volume<float>
{
  public:

   ThreeFormVolume( Volume<Form>* v1 );

   ThreeFormVolume( const FormField& v1 );

   ~ThreeFormVolume(){}

   const float eval( const Vector& P ) const;
   
   const Vector grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "ThreeForm";
      lbl = lbl + "(";
      lbl = lbl + elem1->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

   const FormField elem1;
};



class OneFormVolume : public Volume<Vector>
{
  public:

   OneFormVolume( Volume<Form>* v1 );

   OneFormVolume( const FormField& v1 );

   ~OneFormVolume(){}

   const Vector eval( const Vector& P ) const;
   
   const Matrix grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "OneForm";
      lbl = lbl + "(";
      lbl = lbl + elem1->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

   const FormField elem1;
};



class TwoFormVolume : public Volume<Vector>
{
  public:

   TwoFormVolume( Volume<Form>* v1 );

   TwoFormVolume( const FormField& v1 );

   ~TwoFormVolume(){}

   const Vector eval( const Vector& P ) const;
   
   const Matrix grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "TwoForm";
      lbl = lbl + "(";
      lbl = lbl + elem1->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

   const FormField elem1;
};




class ComponentFormVolume : public Volume<Form>
{
  public:

   ComponentFormVolume( Volume<float>* v0, Volume<Vector>* v1, Volume<Vector>* v2, Volume<float>* v3 );

   ComponentFormVolume( const ScalarField& v0, const VectorField& v1, const VectorField& v2, const ScalarField& v3 );

   ~ComponentFormVolume(){}

   const Form eval( const Vector& P ) const;
   
   const Form grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "ComponentForm";
      lbl = lbl + "(";
      lbl = lbl + elem0->typelabel();
      lbl = lbl + ",";
      lbl = lbl + elem1->typelabel();
      lbl = lbl + ",";
      lbl = lbl + elem2->typelabel();
      lbl = lbl + ",";
      lbl = lbl + elem3->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

   const ScalarField elem0;
   const VectorField elem1, elem2;
   const ScalarField elem3;
   const VectorField cur;
   const ScalarField dive;
};




class StarFormVolume : public Volume<Form>
{
  public:

   StarFormVolume( Volume<Form>* v1 );

   StarFormVolume( const FormField& v1 );

   ~StarFormVolume(){}

   const Form eval( const Vector& P ) const;
   
   const Form grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "HodgeStar";
      lbl = lbl + "(";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

   const FormField elem;
};


class ContractionFormVolume : public Volume<Form>
{
  public:

    ContractionFormVolume( Volume<Vector>* v, Volume<Form>* f );

    ContractionFormVolume( const VectorField& v, const FormField& f );

   ~ContractionFormVolume(){}

   const Form eval( const Vector& P ) const;
   
   const Form grad( const Vector& P ) const;
   
   virtual std::string typelabel() 
   { 
      std::string lbl = "Contraction";
      lbl = lbl + "(";
      lbl = lbl + X->typelabel();
      lbl = lbl + ",";
      lbl = lbl + elem->typelabel();
      lbl = lbl + ")";
      return lbl;
   }


  private:

   VectorField X;
   FormField elem;
   FormField contraxion;

};

}





#endif



