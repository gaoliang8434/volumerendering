
#ifndef __IMPLICITMATRIXSHAPES_H__
#define __IMPLICITMATRIXSHAPES_H__

#include "Volume.h"
#include <cmath>
#include <vector>
#include <iostream>



namespace lux
{

class ConstantMatrixVolume: public Volume<Matrix>
{
  public:

    ConstantMatrixVolume( const Matrix& m );
   ~ConstantMatrixVolume(){}

    const Matrix eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Constant";
       return lbl;
    }

  private:

    Matrix elem;

};

class GriddedMatrixVolume: public Volume<Matrix>
{
  public:

    GriddedMatrixVolume( const MatrixGrid& m );
   ~GriddedMatrixVolume(){}

    const Matrix eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Gridded";
       return lbl;
    }


  private:

    MatrixGrid elem;

};




class VectorMatrixProductVolume: public Volume<Vector>
{
  public:

    VectorMatrixProductVolume( Volume<Matrix>* m, Volume<Vector>* v );
    VectorMatrixProductVolume( const MatrixField& m, const VectorField& v );
   ~VectorMatrixProductVolume(){}

    const Vector eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Multiply";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    MatrixField elem1;
    VectorField elem2;

};



class MatrixVectorProductVolume: public Volume<Vector>
{
  public:

    MatrixVectorProductVolume( Volume<Vector>* v, Volume<Matrix>* m );
    MatrixVectorProductVolume( const VectorField& v, const MatrixField& m );
   ~MatrixVectorProductVolume(){}

    const Vector eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Multiply";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    VectorField elem1;
    MatrixField elem2;

};




class ScalarMatrixProductVolume: public Volume<Matrix>
{
  public:

    ScalarMatrixProductVolume( Volume<Matrix>* m, Volume<float>* v );
    ScalarMatrixProductVolume( const MatrixField& m, const ScalarField& v );
   ~ScalarMatrixProductVolume(){}

    const Matrix eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Multiply";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    MatrixField elem1;
    ScalarField elem2;

};


class MatrixMatrixProductVolume: public Volume<Matrix>
{
  public:

    MatrixMatrixProductVolume( Volume<Matrix>* m, Volume<Matrix>* v );
    MatrixMatrixProductVolume( const MatrixField& m, const MatrixField& v );
   ~MatrixMatrixProductVolume(){}

    const Matrix eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Multiply";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    MatrixField elem1;
    MatrixField elem2;

};





class FloatMatrixProductVolume: public Volume<Matrix>
{
  public:

    FloatMatrixProductVolume( Volume<Matrix>* m, const float v );
    FloatMatrixProductVolume( const MatrixField& m, const float v );
   ~FloatMatrixProductVolume(){}

    const Matrix eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Multiply";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    MatrixField elem1;
    float elem2;

};




class ScalarMatrixDivideVolume: public Volume<Matrix>
{
  public:

    ScalarMatrixDivideVolume( Volume<Matrix>* m, Volume<float>* v );
    ScalarMatrixDivideVolume( const MatrixField& m, const ScalarField& v );
   ~ScalarMatrixDivideVolume(){}

    const Matrix eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Divide";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    MatrixField elem1;
    ScalarField elem2;

};


class FloatMatrixDivideVolume: public Volume<Matrix>
{
  public:

    FloatMatrixDivideVolume( Volume<Matrix>* m, const float v );
    FloatMatrixDivideVolume( const MatrixField& m, const float v );
   ~FloatMatrixDivideVolume(){}

    const Matrix eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Divide";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    MatrixField elem1;
    float elem2;

};





class AddMatrixVolume: public Volume<Matrix>
{
  public:

    AddMatrixVolume( Volume<Matrix>* m1, Volume<Matrix>* m2 );
    AddMatrixVolume( const MatrixField& m1, const MatrixField& m2 );
   ~AddMatrixVolume(){}

    const Matrix eval( const Vector& P ) const;


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

    MatrixField elem1, elem2;

};


class SubtractMatrixVolume: public Volume<Matrix>
{
  public:

    SubtractMatrixVolume( Volume<Matrix>* m1, Volume<Matrix>* m2 );
    SubtractMatrixVolume( const MatrixField& m1, const MatrixField& m2 );
   ~SubtractMatrixVolume(){}

    const Matrix eval( const Vector& P ) const;


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

    MatrixField elem1, elem2;

};


class ProductMatrixVolume: public Volume<Matrix>
{
  public:

    ProductMatrixVolume( Volume<Matrix>* m1, Volume<Matrix>* m2 );
    ProductMatrixVolume( const MatrixField& m1, const MatrixField& m2 );
   ~ProductMatrixVolume(){}

    const Matrix eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Multiply";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ",";
       lbl = lbl + elem2->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    MatrixField elem1, elem2;

};



class NegateMatrixVolume: public Volume<Matrix>
{
  public:

    NegateMatrixVolume( Volume<Matrix>* m1 );
    NegateMatrixVolume( const MatrixField& m1 );
   ~NegateMatrixVolume(){}

    const Matrix eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Negate";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    MatrixField elem1;

};


class OuterProductMatrixVolume: public Volume<Matrix>
{
  public:

    OuterProductMatrixVolume( Volume<Vector>* v1, Volume<Vector>* v2 );
    OuterProductMatrixVolume( const VectorField& v1, const VectorField& v2 );
   ~OuterProductMatrixVolume(){}

    const Matrix eval( const Vector& P ) const;

  private:

    VectorField elem1, elem2;

};



class ExpMatrixVolume: public Volume<Matrix>
{
  public:

    ExpMatrixVolume( Volume<Matrix>* v1 );
    ExpMatrixVolume( const MatrixField& v1 );
   ~ExpMatrixVolume(){}

    const Matrix eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Exp";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    MatrixField elem1;

};


class SinchMatrixVolume: public Volume<Matrix>
{
  public:

    SinchMatrixVolume( Volume<Matrix>* v1 );
    SinchMatrixVolume( const MatrixField& v1 );
   ~SinchMatrixVolume(){}

    const Matrix eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Sinch";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    MatrixField elem1;

};

class OrderedSinchMatrixVolume: public Volume<Matrix>
{
  public:

    OrderedSinchMatrixVolume( Volume<Matrix>* v1, Volume<Matrix>* v2 );
    OrderedSinchMatrixVolume( const MatrixField& v1, const MatrixField& v2 );
   ~OrderedSinchMatrixVolume(){}

    const Matrix eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "OrderedSinch";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    MatrixField elem1;
    MatrixField elem2;

};





class InverseMatrixVolume: public Volume<Matrix>
{
  public:

    InverseMatrixVolume( Volume<Matrix>* v1 );
    InverseMatrixVolume( const MatrixField& v1 );
   ~InverseMatrixVolume(){}

    const Matrix eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Inverse";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    MatrixField elem1;

};



class DetMatrixVolume: public Volume<float>
{
  public:

    DetMatrixVolume( Volume<Matrix>* v1 );
    DetMatrixVolume( const MatrixField& v1 );
   ~DetMatrixVolume(){}

    const float eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Det";
       lbl = lbl + "(";
       lbl = lbl + elem1->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    MatrixField elem1;

};



class GradientMatrixVolume : public Volume<Matrix> 
{
  public:

    GradientMatrixVolume( Volume<Vector>* v ) ;

    GradientMatrixVolume( const VectorField& v ) ;

    ~GradientMatrixVolume(){}

    const Matrix eval( const Vector& P ) const;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Grad";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }


  private:

    VectorField elem;
};





class WarpMatrixVolume : public Volume<Matrix> 
{
  public:

    WarpMatrixVolume( Volume<Matrix>* v, Volume<Vector>* u ) ;

    WarpMatrixVolume( const MatrixField& v, const VectorField& u ) ;

   ~WarpMatrixVolume(){}

    const Matrix eval( const Vector& P ) const ;


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

    MatrixField elem;
    VectorField warp;
};




class AdvectMatrixVolume : public Volume<Matrix> 
{
  public:

    AdvectMatrixVolume( Volume<Matrix>* v, Volume<Vector>* u, float _dt ) ;

    AdvectMatrixVolume( const MatrixField& v, const VectorField& u, float _dt ) ;

   ~AdvectMatrixVolume(){}

    const Matrix eval( const Vector& P ) const ;


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

    MatrixField elem;
    VectorField velocity;
    float dt;
};





class RotatorMatrixVolume : public Volume<Matrix> 
{
  public:

    RotatorMatrixVolume( Volume<Vector>* v ) ;

    RotatorMatrixVolume( const VectorField& v ) ;

   ~RotatorMatrixVolume(){}

    const Matrix eval( const Vector& P ) const ;


    virtual std::string typelabel() 
    { 
       std::string lbl = "Rotator";
       lbl = lbl + "(";
       lbl = lbl + elem->typelabel();
       lbl = lbl + ")";
       return lbl;
    }

  private:

    VectorField elem;
};




}


#endif
