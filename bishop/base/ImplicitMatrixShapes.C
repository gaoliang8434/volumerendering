
#include "ImplicitMatrixShapes.h"
#include "LinearAlgebra.h"


using namespace lux;
using namespace std;

ConstantMatrixVolume::ConstantMatrixVolume( const Matrix& m ) :
   elem(m)
   {}

const Matrix ConstantMatrixVolume::eval( const Vector& P ) const
{
   return elem;
}


GriddedMatrixVolume::GriddedMatrixVolume( const MatrixGrid& m ) :
  elem(m)
  {}

const Matrix GriddedMatrixVolume::eval( const Vector& P ) const
{
   return elem->eval(P);
}





VectorMatrixProductVolume::VectorMatrixProductVolume( Volume<Matrix>* m, Volume<Vector>* v ) :
   elem1(m),
   elem2(v)
   {}

VectorMatrixProductVolume::VectorMatrixProductVolume( const MatrixField& m, const VectorField& v ) :
   elem1(m),
   elem2(v)
   {}


const Vector VectorMatrixProductVolume::eval( const Vector& P ) const
{
   return (elem2->eval(P))*(elem1->eval(P));
}





MatrixVectorProductVolume::MatrixVectorProductVolume( Volume<Vector>* v, Volume<Matrix>* m )  :
   elem1(v),
   elem2(m)
   {}

MatrixVectorProductVolume::MatrixVectorProductVolume( const VectorField& v, const MatrixField& m ) :
   elem1(v),
   elem2(m)
   {}

const Vector MatrixVectorProductVolume::eval( const Vector& P ) const
{
   return (elem2->eval(P))*(elem1->eval(P));
}






ScalarMatrixProductVolume::ScalarMatrixProductVolume( Volume<Matrix>* m, Volume<float>* v )  :
   elem1(m),
   elem2(v)
   {}

ScalarMatrixProductVolume::ScalarMatrixProductVolume( const MatrixField& m, const ScalarField& v ) :
   elem1(m),
   elem2(v)
   {}

const Matrix ScalarMatrixProductVolume::eval( const Vector& P ) const
{
   return (elem1->eval(P))*(elem2->eval(P));
}




MatrixMatrixProductVolume::MatrixMatrixProductVolume( Volume<Matrix>* m, Volume<Matrix>* v )  :
   elem1(m),
   elem2(v)
   {}

MatrixMatrixProductVolume::MatrixMatrixProductVolume( const MatrixField& m, const MatrixField& v ) :
   elem1(m),
   elem2(v)
   {}

const Matrix MatrixMatrixProductVolume::eval( const Vector& P ) const
{
   return (elem1->eval(P))*(elem2->eval(P));
}









ScalarMatrixDivideVolume::ScalarMatrixDivideVolume( Volume<Matrix>* m, Volume<float>* v )  :
   elem1(m),
   elem2(v)
   {}

ScalarMatrixDivideVolume::ScalarMatrixDivideVolume( const MatrixField& m, const ScalarField& v ) :
   elem1(m),
   elem2(v)
   {}

const Matrix ScalarMatrixDivideVolume::eval( const Vector& P ) const
{
   return (elem1->eval(P))/(elem2->eval(P));
}




AddMatrixVolume::AddMatrixVolume( Volume<Matrix>* m1, Volume<Matrix>* m2 ) :
  elem1(m1),
  elem2(m2)
  {}

AddMatrixVolume::AddMatrixVolume( const MatrixField& m1, const MatrixField& m2 ) :
  elem1(m1),
  elem2(m2)
  {}

const Matrix AddMatrixVolume::eval( const Vector& P ) const
{
   return elem1->eval(P) + elem2->eval(P);
}



SubtractMatrixVolume::SubtractMatrixVolume( Volume<Matrix>* m1, Volume<Matrix>* m2 ) :
   elem1(m1),
   elem2(m2)
   {}

SubtractMatrixVolume::SubtractMatrixVolume( const MatrixField& m1, const MatrixField& m2 ) :
   elem1(m1),
   elem2(m2)
   {}

const Matrix SubtractMatrixVolume::eval( const Vector& P ) const
{
   return elem1->eval(P) - elem2->eval(P);
}



ProductMatrixVolume::ProductMatrixVolume( Volume<Matrix>* m1, Volume<Matrix>* m2 ) :
   elem1(m1),
   elem2(m2)
   {}

ProductMatrixVolume::ProductMatrixVolume( const MatrixField& m1, const MatrixField& m2 ) :
   elem1(m1),
   elem2(m2)
   {}


const Matrix ProductMatrixVolume::eval( const Vector& P ) const
{
    return (elem1->eval(P))*(elem2->eval(P));
}




NegateMatrixVolume::NegateMatrixVolume( Volume<Matrix>* m1 ) :
   elem1(m1)
   {}

NegateMatrixVolume::NegateMatrixVolume( const MatrixField& m1 ) :
   elem1(m1)
   {}


const Matrix NegateMatrixVolume::eval( const Vector& P ) const
{
   return (-1.0)*(elem1->eval(P));
}




OuterProductMatrixVolume::OuterProductMatrixVolume( Volume<Vector>* v1, Volume<Vector>* v2 ) :
   elem1(v1),
   elem2(v2)
   {}

OuterProductMatrixVolume::OuterProductMatrixVolume( const VectorField& v1, const VectorField& v2 ) :
   elem1(v1),
   elem2(v2)
   {}


const Matrix OuterProductMatrixVolume::eval( const Vector& P ) const
{
   Vector v1 = elem1->eval(P);
   Vector v2 = elem2->eval(P);
   Matrix result;
   outer_product( v1, v2, result );
   return result;
}




ExpMatrixVolume::ExpMatrixVolume( Volume<Matrix>* v1 ) : 
   elem1(v1)
   {}

ExpMatrixVolume::ExpMatrixVolume( const MatrixField& v1 ) :
   elem1(v1)
   {}

const Matrix ExpMatrixVolume::eval( const Vector& P ) const
{
   return lux::exp( elem1->eval(P) );
}




SinchMatrixVolume::SinchMatrixVolume( Volume<Matrix>* v1 ) : 
   elem1(v1)
   {}

SinchMatrixVolume::SinchMatrixVolume( const MatrixField& v1 ) :
   elem1(v1)
   {}

const Matrix SinchMatrixVolume::eval( const Vector& P ) const
{
   return lux::sinch( elem1->eval(P) );
}




OrderedSinchMatrixVolume::OrderedSinchMatrixVolume( Volume<Matrix>* v1, Volume<Matrix>* v2 ) :
    elem1(v1),
    elem2(v2)
    {}

OrderedSinchMatrixVolume::OrderedSinchMatrixVolume( const MatrixField& v1, const MatrixField& v2 ) :
    elem1(v1),
    elem2(v2)
    {}

const Matrix OrderedSinchMatrixVolume::eval( const Vector& P ) const
{
   Matrix a = elem1->eval(P);
   Matrix b = elem2->eval(P);
   return orderedSinch( a, b );
}












InverseMatrixVolume::InverseMatrixVolume( Volume<Matrix>* v1 ) : 
   elem1(v1)
   {}

InverseMatrixVolume::InverseMatrixVolume( const MatrixField& v1 ) :
   elem1(v1)
   {}

const Matrix InverseMatrixVolume::eval( const Vector& P ) const
{
   return lux::inverse( elem1->eval(P) );
}






DetMatrixVolume::DetMatrixVolume( Volume<Matrix>* v1 ) :
   elem1(v1)
   {}

DetMatrixVolume::DetMatrixVolume( const MatrixField& v1 ) :
   elem1(v1)
   {}

const float DetMatrixVolume::eval( const Vector& P ) const
{
   return lux::det( elem1->eval(P) );
}






GradientMatrixVolume::GradientMatrixVolume( Volume<Vector>* v ) :
  elem(v)
  {}

GradientMatrixVolume::GradientMatrixVolume( const VectorField& v ) :
  elem(v)
  {}


const Matrix GradientMatrixVolume::eval( const Vector& P ) const { return elem->grad(P); }





WarpMatrixVolume::WarpMatrixVolume( Volume<Matrix>* v, Volume<Vector>* u ) :
   elem(v),
   warp(u)
   {}

WarpMatrixVolume::WarpMatrixVolume( const MatrixField& v, const VectorField& u ) :
   elem(v),
   warp(u)
   {}


const Matrix WarpMatrixVolume::eval( const Vector& P ) const 
{
   Vector X = warp->eval(P);
   return elem->eval(X);
}




FloatMatrixDivideVolume::FloatMatrixDivideVolume( Volume<Matrix>* m, float v ) :
   elem1(m),
   elem2(v)
   {}

FloatMatrixDivideVolume::FloatMatrixDivideVolume( const MatrixField& m, float v ) :
   elem1(m),
   elem2(v)
   {}

const Matrix FloatMatrixDivideVolume::eval( const Vector& P ) const
{
   return (elem1->eval(P))/elem2;
}



FloatMatrixProductVolume::FloatMatrixProductVolume( Volume<Matrix>* m, const float v ) :
   elem1(m),
   elem2(v)
   {}

FloatMatrixProductVolume::FloatMatrixProductVolume( const MatrixField& m, const float v ) :
   elem1(m),
   elem2(v)
   {}

const Matrix FloatMatrixProductVolume::eval( const Vector& P ) const
{
   return  (elem1->eval(P))*elem2;
}




AdvectMatrixVolume::AdvectMatrixVolume( Volume<Matrix>* v, Volume<Vector>* u, const float delt ) :
       elem(v),
       velocity(u),
       dt (delt)
    {}

AdvectMatrixVolume::AdvectMatrixVolume( const MatrixField& v, const VectorField& u, const float delt ) :
       elem(v),
       velocity(u),
       dt (delt)
    {}


const Matrix AdvectMatrixVolume::eval( const Vector& P ) const 
{ 
   Vector X = P - ( velocity->eval(P) )*dt;
   return elem->eval(X); 
}



RotatorMatrixVolume::RotatorMatrixVolume( Volume<Vector>* v ) :
   elem(v)
   {}

RotatorMatrixVolume::RotatorMatrixVolume( const VectorField& v ) :
   elem(v)
   {}


const Matrix RotatorMatrixVolume::eval( const Vector& P ) const 
{
   Vector e = elem->eval(P);
   return rotation( e.unitvector(), e.magnitude() );
}



