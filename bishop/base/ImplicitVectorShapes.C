
#include "ImplicitVolumeShapes.h"
#include "ImplicitVectorShapes.h"
#include "RayMarcher.h"
#include "Fields.h"

#include <iostream>
using namespace std;
using namespace lux;

MultiplyVectorVolume::MultiplyVectorVolume( Volume<Vector> * v, const float a ) : 
       elem(v),
       factor( new ConstantVolume(1.0) ),
       constant (a) 
    {}

MultiplyVectorVolume::MultiplyVectorVolume( Volume<Vector> * v, Volume<float>* u ) : 
      elem(v),
      factor(u),
      constant (1.0) 
    {}

MultiplyVectorVolume::MultiplyVectorVolume( const VectorField& v, const float a ) : 
       elem(v),
       factor( new ConstantVolume(1.0) ),
       constant (a) 
    {}

MultiplyVectorVolume::MultiplyVectorVolume( const VectorField& v, const ScalarField& u ) : 
      elem(v),
      factor(u),
      constant (1.0) 
    {}







const Vector NoiseVectorVolume::eval( const Vector& P ) const
{
   float nx = noise->eval(P+Vector(dx,0,0));
   float ny = noise->eval(P+Vector(0,dx,0));
   float nz = noise->eval(P+Vector(0,0,dx));
   Vector g( ny-nz, nz-nx, nx-ny );
   g /= dx;
   return g;
}


const Matrix NoiseVectorVolume::grad( const Vector& P ) const
{
   Vector n = eval(P);
   Vector nx = (eval(P+Vector(dx,0,0))-n)/dx;
   Vector ny = (eval(P+Vector(0,dx,0))-n)/dx;
   Vector nz = (eval(P+Vector(0,0,dx))-n)/dx;
   Matrix result(nx,ny,nz);
   return result.transpose();
}




const Vector NoiseSampleVectorVolume::eval( const Vector& P ) const
{
   float nx = noise->eval(P+Vector(dx,0,0));
   float ny = noise->eval(P+Vector(0,dx,0));
   float nz = noise->eval(P+Vector(0,0,dx));
   Vector g( nx, ny, nz );
   return g;
}


const Matrix NoiseSampleVectorVolume::grad( const Vector& P ) const
{
   Vector n = eval(P);
   Vector nx = (eval(P+Vector(dx,0,0))-n)/dx;
   Vector ny = (eval(P+Vector(0,dx,0))-n)/dx;
   Vector nz = (eval(P+Vector(0,0,dx))-n)/dx;
   Matrix result(nx,ny,nz);
   return result.transpose();
}









ConstantVectorVolume::ConstantVectorVolume( const Vector v ) :
       value (v),
       gradvalue (Vector(0,0,0), Vector(0,0,0), Vector(0,0,0) )
    {}



const Vector ConstantVectorVolume::eval( const Vector& P ) const 
{
   return value; 
}

const Matrix ConstantVectorVolume::grad( const Vector& P ) const { return gradvalue; }







TranslateVectorVolume::TranslateVectorVolume( Volume<Vector> * v, const Vector& s ) :
       elem(v),
       translate (s)
    {}

TranslateVectorVolume::TranslateVectorVolume( const VectorField& v, const Vector& s ) :
       elem(v),
       translate (s)
    {}



const Vector TranslateVectorVolume::eval( const Vector& P ) const
{
   Vector X = P - translate;
   return  (elem->eval(X));
}


const Matrix TranslateVectorVolume::grad( const Vector& P ) const
{
   Vector X = P - translate;
   return  (elem->grad(X));
}





ScaleVectorVolume::ScaleVectorVolume( Volume<Vector> * v, const Vector& s ) :
       elem(v),
       scale (s)
    {}

ScaleVectorVolume::ScaleVectorVolume( const VectorField& v, const Vector& s ) :
       elem(v),
       scale (s)
    {}

ScaleVectorVolume::ScaleVectorVolume( Volume<Vector> * v, const float& s ) :
       elem(v),
       scale (Vector(s,s,s))
    {}

ScaleVectorVolume::ScaleVectorVolume( const VectorField& v, const float& s ) :
       elem(v),
       scale (Vector(s,s,s))
    {}





const Vector ScaleVectorVolume::eval( const Vector& P ) const
{
   Vector X = P;
   X[0] /= scale[0];
   X[1] /= scale[1];
   X[2] /= scale[2];
   Vector result = (elem->eval(X));
   result[0] *= scale[0];
   result[1] *= scale[1];
   result[2] *= scale[2];
   return result;
}

const Matrix ScaleVectorVolume::grad( const Vector& P ) const
{
   Vector X = P;
   X[0] /= scale[0];
   X[1] /= scale[1];
   X[2] /= scale[2];
   Matrix d =  (elem->grad(X));
   for( int a=0;a<3;a++)
   {
      for( int b=0;b<3;b++)
      {
         d.Set( a, b, d.Get(a,b)*scale[a]*scale[a] );
      }
   }
   return d;
}








RotateVectorVolume::RotateVectorVolume( Volume<Vector> * v, const Vector& s ) :
       elem(v)
    {
       R = inverse(rotation( s.unitvector() , s.magnitude() ));
       Rinv = rotation( s.unitvector(), s.magnitude() );
    }

RotateVectorVolume::RotateVectorVolume( const VectorField& v, const Vector& s ) :
       elem(v)
    {
       R = inverse(rotation( s.unitvector() , s.magnitude() ));
       Rinv = rotation( s.unitvector(), s.magnitude() );
    }



const Vector RotateVectorVolume::eval( const Vector& P ) const
{
   Vector X = R*P;
   return Rinv*(elem->eval(X));
}

const Matrix RotateVectorVolume::grad( const Vector& P ) const
{
   Vector X = R*P;
   return  Rinv*(elem->grad(X))*R;
}








NegateVectorVolume::NegateVectorVolume( Volume<Vector> * v ) :
       elem(v)
    {}

NegateVectorVolume::NegateVectorVolume( const VectorField& v ) :
       elem(v)
    {}



const Vector NegateVectorVolume::eval( const Vector& P ) const
{
   return  -(elem->eval(P));
}

const Matrix NegateVectorVolume::grad( const Vector& P ) const
{
   return  -(elem->grad(P));
}





DotProductVectorVolume::DotProductVectorVolume( Volume<Vector>* v1, Volume<Vector>* v2 ) :
      elem1(v1),
      elem2(v2)
    {}

DotProductVectorVolume::DotProductVectorVolume( const VectorField& v1, const VectorField& v2 ) :
      elem1(v1),
      elem2(v2)
    {}


const float DotProductVectorVolume::eval( const Vector& P ) const
{
   Vector v1 = elem1->eval(P);
   Vector v2 = elem2->eval(P);
   return v1*v2;
}

const Vector DotProductVectorVolume::grad( const Vector& P ) const
{
   Vector v1 = elem1->eval(P);
   Vector v2 = elem2->eval(P);
   return v1*(elem2->grad(P)) + (elem1->grad(P))*v2; 
}





MagnitudeVectorVolume::MagnitudeVectorVolume( Volume<Vector>* v1 ) :
       elem(v1)
    {}

MagnitudeVectorVolume::MagnitudeVectorVolume( const VectorField& v1 ) :
       elem(v1)
    {}


const float MagnitudeVectorVolume::eval( const Vector& P ) const
    {
       Vector v1 = elem->eval(P);
       return v1.magnitude();
    }

const Vector MagnitudeVectorVolume::grad( const Vector& P ) const
    {
       Vector v1 = elem->eval(P);
       return v1.unitvector() * (elem->grad(P));
    }




UnitVectorVolume::UnitVectorVolume( Volume<Vector>* v1 ) :
      elem(v1)
    {}

UnitVectorVolume::UnitVectorVolume( const VectorField& v1 ) :
      elem(v1)
    {}


const Vector UnitVectorVolume::eval( const Vector& P ) const
    {
       Vector v1 = elem->eval(P);
       return v1.unitvector();
    }

const Matrix UnitVectorVolume::grad( const Vector& P ) const
    {
       Vector v1 = elem->eval(P);
       float vmag = v1.magnitude();
       if( vmag > 0 )
       {
          Vector vu = v1/vmag;
          Matrix g = elem->grad(P);
          Matrix gg;
          outer_product( vu, vu, gg);
          gg = g*gg;
          return (g - gg)/vmag;
       }
       return unitMatrix();
    }






const Vector MultiplyVectorVolume::eval( const Vector& P ) const
    {
       Vector result = elem->eval(P);
       result *= constant;
       result *= factor->eval(P);
       return  result;
    }

    const Matrix MultiplyVectorVolume::grad( const Vector& P ) const
    {
       Matrix result = elem->grad(P);
       result *= factor->eval(P);
       Vector g = factor->grad(P);
       Vector f = elem->eval(P);
       Matrix fg;
       outer_product( f, g, fg );
       result += fg;
       return result*constant;
    }



AddVectorVolume::AddVectorVolume( Volume<Vector> * v1, Volume<Vector> * v2 ) :
      elem1(v1),
      elem2(v2)
    {}

AddVectorVolume::AddVectorVolume( const VectorField& v1, const VectorField& v2 ) :
      elem1(v1),
      elem2(v2)
    {}



    const Vector AddVectorVolume::eval( const Vector& P ) const
    {
       return  elem1->eval(P) + elem2->eval(P);
    }

    const Matrix AddVectorVolume::grad( const Vector& P ) const
    {
       return  elem1->grad(P) + elem2->grad(P);
    }


SubtractVectorVolume::SubtractVectorVolume( Volume<Vector> * v1, Volume<Vector> * v2 ) :
      elem1(v1),
      elem2(v2)
    {}

SubtractVectorVolume::SubtractVectorVolume( const VectorField& v1, const VectorField& v2 ) :
      elem1(v1),
      elem2(v2)
    {}



    const Vector SubtractVectorVolume::eval( const Vector& P ) const
    {
       return  elem1->eval(P) - elem2->eval(P);
    }

    const Matrix SubtractVectorVolume::grad( const Vector& P ) const
    {
       return  elem1->grad(P) - elem2->grad(P);
    }




IdentityVectorVolume::IdentityVectorVolume()
    {
       gradvalue = unitMatrix();
    }


    const Vector IdentityVectorVolume::eval( const Vector& P ) const { return P; }

    const Matrix IdentityVectorVolume::grad( const Vector& P ) const { return gradvalue;  }




ImplicitPointVectorVolume::ImplicitPointVectorVolume( Volume<float>* v, const float dx, const int nb, const float thresh ) :
       elem(v),
       step (dx),
       nbIterations (nb),
       threshold (thresh)
    {}

ImplicitPointVectorVolume::ImplicitPointVectorVolume( const ScalarField& v, const float dx, const int nb, const float thresh ) :
       elem(v),
       step (dx),
       nbIterations (nb),
       threshold (thresh)
    {}



    const Vector ImplicitPointVectorVolume::eval( const Vector& P ) const
    {
       float dx = step;
       Vector X = P;
       float previous = elem->eval(P);
       if( fabs(previous) > threshold ){ return X; } // If you are far out, dont worry about it
       Vector D = elem->grad(P);
       if( D.magnitude() > 0 )
       { 
          D.normalize(); 
       }
       else
       {
          return X;
       }
       int flipcount = 0;
       if( previous == 0 ){ return X; }

       int totalcount = 0;
       while( flipcount < nbIterations  && totalcount < 1000 )
       {
          X += D * dx;
	  float value = elem->eval(X);
	  if( fabs(value) < 1.0e-10 ){ return X; }
	  if( value * previous  < 0 )
	  {
	     dx /= 2.0;
	     ++flipcount;
	     totalcount = 0;
	  } 
	  else if( fabs(value) > fabs(previous) )
	  {
	     dx *= -1;
	  }
	  previous = value;
	  D =  elem->grad(X);
          if( D.magnitude() > 0 )
          { 
             D.normalize(); 
          }
          else
          {
             return X;
          }
	  ++totalcount;
       }
       return X;
    }

    // evaluate gradient by finite difference
    // size of difference is controlled via the "dx" attribute
    const Matrix ImplicitPointVectorVolume::grad( const Vector& P ) const 
    {
       float dx = step;
       Vector e = eval(P);
       Vector ex = eval(P+Vector(dx,0,0));
       Vector ey = eval(P+Vector(0,dx,0));
       Vector ez = eval(P+Vector(0,0,dx));
       Matrix result( ex-e, ey-e, ez-e );
       result /= dx;
       return result.transpose();
    }




GradientVectorVolume::GradientVectorVolume( Volume<float>* v, const float dx ) :
      elem(v),
      step(dx)
    {}

GradientVectorVolume::GradientVectorVolume( const ScalarField& v, const float dx ) :
      elem(v),
      step(dx)
    {}


    const Vector GradientVectorVolume::eval( const Vector& P ) const { return elem->grad(P); }

    // evaluate gradient by finite difference
    // size of difference is controlled via the "dx" attribute
    const Matrix GradientVectorVolume::grad( const Vector& P ) const 
    {
       float dx = step;
       Vector exp = eval(P+Vector(dx,0,0));
       Vector eyp = eval(P+Vector(0,dx,0));
       Vector ezp = eval(P+Vector(0,0,dx));
       Vector exm = eval(P-Vector(dx,0,0));
       Vector eym = eval(P-Vector(0,dx,0));
       Vector ezm = eval(P-Vector(0,0,dx));
       Matrix result( exp-exm, eyp-eym, ezp-ezm );
       result /= 2.0*dx;
       return result.transpose();
    }




// N is the length of the difference kernel
FiniteDifferenceGradientVectorVolume::FiniteDifferenceGradientVectorVolume( Volume<float>* v, const int N, const float dx, const float dy, const float dz ) :
elem(v),
dX(dx), dY(dy), dZ(dz)
{
   if( N <= 1 )
   {
      coefficients.push_back(0.5);
   }
   else if( N == 2 )
   {
      coefficients.push_back( 0.57142857142857151); 
      coefficients.push_back(-0.071428571428571438);
   }
   else if( N==3 )
   {
      coefficients.push_back( 0.60810810810810767); 
      coefficients.push_back(-0.1188063063063063);
      coefficients.push_back( 0.010698198198198195);
   }
   else if( N==4 )
   {
      coefficients.push_back( 0.63039399624765435);
      coefficients.push_back(-0.15064623723160311);
      coefficients.push_back( 0.02101573900354391);
      coefficients.push_back(-0.00076349801959558095);
   }
   else if( N==5 )
   {
      coefficients.push_back( 0.64535955746772966);
      coefficients.push_back(-0.17321459058935865);
      coefficients.push_back( 0.029559879293859895); 
      coefficients.push_back(-0.0017378184115276803); 
      coefficients.push_back( 3.2972239295226059e-05);
   }
   else if( N==6 )
   {
      coefficients.push_back( 0.65609985088633171);
      coefficients.push_back(-0.18996590374149147);
      coefficients.push_back( 0.036481469793162628);
      coefficients.push_back(-0.0026999821278442119);
      coefficients.push_back( 8.5528288191602389e-05);
      coefficients.push_back(-9.6309840124388086e-07);
   }
   else if( N==7 )
   {
      coefficients.push_back( 0.66418180822321082);
      coefficients.push_back(-0.20286437997462345);
      coefficients.push_back( 0.042119156095069443);
      coefficients.push_back(-0.0035797086604024858);
      coefficients.push_back( 0.00014591197349313843);
      coefficients.push_back(-2.8080248121598451e-06);
      coefficients.push_back( 2.0367364328490342e-08);
   }
   else if( N==8 )
   {
      coefficients.push_back( 0.67048320765246205);
      coefficients.push_back(-0.213090576806653);
      coefficients.push_back( 0.046767457952830598);
      coefficients.push_back(-0.004362453117772256);
      coefficients.push_back( 0.00020752664200811052);
      coefficients.push_back(-5.2279046772863686e-06);
      coefficients.push_back( 6.598736368669243e-08);
      coefficients.push_back(-3.26811256996574e-10);
   }
   else if( N==9 )
   {
      coefficients.push_back( 0.67553384734048427);
      coefficients.push_back(-0.22139161101119745);
      coefficients.push_back( 0.050651292408117839);
      coefficients.push_back(-0.0050528139760881441);
      coefficients.push_back( 0.00026710733705086263);
      coefficients.push_back(-7.9566998459801905e-06);
      coefficients.push_back( 1.3328639394368043e-07);
      coefficients.push_back(-1.165378723448066e-09);
      coefficients.push_back( 4.1187471834005603e-12);
   }
   else if( N>=10 )
   {
      coefficients.push_back( 0.67967266644170954);
      coefficients.push_back(-0.22826167257103971);
      coefficients.push_back( 0.053937736604997032);
      coefficients.push_back(-0.005661092015920628);
      coefficients.push_back( 0.00032320263464067625);
      coefficients.push_back(-1.0809865109373177e-05);
      coefficients.push_back( 2.1645208248117699e-07);
      coefficients.push_back(-2.5401252316582755e-09);
      coefficients.push_back( 1.6035414536489318e-11);
      coefficients.push_back(-4.1871266927217782e-14);
   }
}

FiniteDifferenceGradientVectorVolume::FiniteDifferenceGradientVectorVolume( const ScalarField& v, const int N, const float dx, const float dy, const float dz ) :
elem(v),
dX(dx), dY(dy), dZ(dz)
{
   if( N <= 1 )
   {
      coefficients.push_back(0.5);
   }
   else if( N == 2 )
   {
      coefficients.push_back( 0.57142857142857151); 
      coefficients.push_back(-0.071428571428571438);
   }
   else if( N==3 )
   {
      coefficients.push_back( 0.60810810810810767); 
      coefficients.push_back(-0.1188063063063063);
      coefficients.push_back( 0.010698198198198195);
   }
   else if( N==4 )
   {
      coefficients.push_back( 0.63039399624765435);
      coefficients.push_back(-0.15064623723160311);
      coefficients.push_back( 0.02101573900354391);
      coefficients.push_back(-0.00076349801959558095);
   }
   else if( N==5 )
   {
      coefficients.push_back( 0.64535955746772966);
      coefficients.push_back(-0.17321459058935865);
      coefficients.push_back( 0.029559879293859895); 
      coefficients.push_back(-0.0017378184115276803); 
      coefficients.push_back( 3.2972239295226059e-05);
   }
   else if( N==6 )
   {
      coefficients.push_back( 0.65609985088633171);
      coefficients.push_back(-0.18996590374149147);
      coefficients.push_back( 0.036481469793162628);
      coefficients.push_back(-0.0026999821278442119);
      coefficients.push_back( 8.5528288191602389e-05);
      coefficients.push_back(-9.6309840124388086e-07);
   }
   else if( N==7 )
   {
      coefficients.push_back( 0.66418180822321082);
      coefficients.push_back(-0.20286437997462345);
      coefficients.push_back( 0.042119156095069443);
      coefficients.push_back(-0.0035797086604024858);
      coefficients.push_back( 0.00014591197349313843);
      coefficients.push_back(-2.8080248121598451e-06);
      coefficients.push_back( 2.0367364328490342e-08);
   }
   else if( N==8 )
   {
      coefficients.push_back( 0.67048320765246205);
      coefficients.push_back(-0.213090576806653);
      coefficients.push_back( 0.046767457952830598);
      coefficients.push_back(-0.004362453117772256);
      coefficients.push_back( 0.00020752664200811052);
      coefficients.push_back(-5.2279046772863686e-06);
      coefficients.push_back( 6.598736368669243e-08);
      coefficients.push_back(-3.26811256996574e-10);
   }
   else if( N==9 )
   {
      coefficients.push_back( 0.67553384734048427);
      coefficients.push_back(-0.22139161101119745);
      coefficients.push_back( 0.050651292408117839);
      coefficients.push_back(-0.0050528139760881441);
      coefficients.push_back( 0.00026710733705086263);
      coefficients.push_back(-7.9566998459801905e-06);
      coefficients.push_back( 1.3328639394368043e-07);
      coefficients.push_back(-1.165378723448066e-09);
      coefficients.push_back( 4.1187471834005603e-12);
   }
   else if( N>=10 )
   {
      coefficients.push_back( 0.67967266644170954);
      coefficients.push_back(-0.22826167257103971);
      coefficients.push_back( 0.053937736604997032);
      coefficients.push_back(-0.005661092015920628);
      coefficients.push_back( 0.00032320263464067625);
      coefficients.push_back(-1.0809865109373177e-05);
      coefficients.push_back( 2.1645208248117699e-07);
      coefficients.push_back(-2.5401252316582755e-09);
      coefficients.push_back( 1.6035414536489318e-11);
      coefficients.push_back(-4.1871266927217782e-14);
   }
}


const Vector FiniteDifferenceGradientVectorVolume::eval( const Vector& P ) const
{
   double fx = 0;
   double fy = 0;
   double fz = 0;
   int ii = 0;
   for( size_t i=0;i<coefficients.size();i++)
   {
      ii++;
      Vector DX(ii*dX,0,0);
      fx += elem->eval(P+DX)*coefficients[i];
      fx -= elem->eval(P-DX)*coefficients[i];
      Vector DY(0,ii*dY,0);
      fy += elem->eval(P+DY)*coefficients[i];
      fy -= elem->eval(P-DY)*coefficients[i];
      Vector DZ(0,0,ii*dZ);
      fz += elem->eval(P+DZ)*coefficients[i];
      fz -= elem->eval(P-DZ)*coefficients[i];
   }
   return Vector( fx/dX, fy/dY, fz/dZ );
}

const Matrix FiniteDifferenceGradientVectorVolume::grad( const Vector& P ) const 
{
   Vector fx;
   Vector fy;
   Vector fz;
   int ii = 0;
   for( size_t i=0;i<coefficients.size();i++)
   {
      ii++;
      Vector DX(ii*dX,0,0);
      fx += eval(P+DX)*coefficients[i];
      fx -= eval(P-DX)*coefficients[i];
      Vector DY(0,ii*dY,0);
      fy += eval(P+DY)*coefficients[i];
      fy -= eval(P-DY)*coefficients[i];
      Vector DZ(0,0,ii*dZ);
      fz += eval(P+DZ)*coefficients[i];
      fz -= eval(P-DZ)*coefficients[i];
   }
   Matrix result( fx/dX, fy/dY, fz/dZ );
   return result;
}





















// N is the length of the difference kernel
FiniteDifferenceBoundedGradientVectorVolume::FiniteDifferenceBoundedGradientVectorVolume( Volume<float>* v, const int N, const GridBox& gb ) :
elem(v),
maxN(N),
gB(gb)
{
   for(int nn=1;nn<11;nn++)
   {
      vector<double> coeffs;
      if( nn==1 )
      {
         coeffs.push_back(0.5);
      }
      else if( nn==2 )
      {
         coeffs.push_back( 0.57142857142857151); 
         coeffs.push_back(-0.071428571428571438);
      }
      else if( nn==3 )
      {
         coeffs.push_back( 0.60810810810810767); 
         coeffs.push_back(-0.1188063063063063);
         coeffs.push_back( 0.010698198198198195);
      }
      else if( nn==4 )
      {
         coeffs.push_back( 0.63039399624765435);
         coeffs.push_back(-0.15064623723160311);
         coeffs.push_back( 0.02101573900354391);
         coeffs.push_back(-0.00076349801959558095);
      }
      else if( nn==5 )
      {
         coeffs.push_back( 0.64535955746772966);
         coeffs.push_back(-0.17321459058935865);
         coeffs.push_back( 0.029559879293859895); 
         coeffs.push_back(-0.0017378184115276803); 
         coeffs.push_back( 3.2972239295226059e-05);
      }
      else if( nn==6 )
      {
         coeffs.push_back( 0.65609985088633171);
         coeffs.push_back(-0.18996590374149147);
         coeffs.push_back( 0.036481469793162628);
         coeffs.push_back(-0.0026999821278442119);
         coeffs.push_back( 8.5528288191602389e-05);
         coeffs.push_back(-9.6309840124388086e-07);
      }
      else if( nn==7 )
      {
         coeffs.push_back( 0.66418180822321082);
         coeffs.push_back(-0.20286437997462345);
         coeffs.push_back( 0.042119156095069443);
         coeffs.push_back(-0.0035797086604024858);
         coeffs.push_back( 0.00014591197349313843);
         coeffs.push_back(-2.8080248121598451e-06);
         coeffs.push_back( 2.0367364328490342e-08);
      }
      else if( nn==8 )
      { 
         coeffs.push_back( 0.67048320765246205);
         coeffs.push_back(-0.213090576806653);
         coeffs.push_back( 0.046767457952830598);
         coeffs.push_back(-0.004362453117772256);
         coeffs.push_back( 0.00020752664200811052);
         coeffs.push_back(-5.2279046772863686e-06);
         coeffs.push_back( 6.598736368669243e-08);
         coeffs.push_back(-3.26811256996574e-10);
      }
      else if( nn==9 )
      {
         coeffs.push_back( 0.67553384734048427);
         coeffs.push_back(-0.22139161101119745);
         coeffs.push_back( 0.050651292408117839);
         coeffs.push_back(-0.0050528139760881441);
         coeffs.push_back( 0.00026710733705086263);
         coeffs.push_back(-7.9566998459801905e-06);
         coeffs.push_back( 1.3328639394368043e-07);
         coeffs.push_back(-1.165378723448066e-09);
         coeffs.push_back( 4.1187471834005603e-12);
      }
      else if( nn==10 )
      {
         coeffs.push_back( 0.67967266644170954);
         coeffs.push_back(-0.22826167257103971);
         coeffs.push_back( 0.053937736604997032);
         coeffs.push_back(-0.005661092015920628);
         coeffs.push_back( 0.00032320263464067625);
         coeffs.push_back(-1.0809865109373177e-05);
         coeffs.push_back( 2.1645208248117699e-07);
         coeffs.push_back(-2.5401252316582755e-09);
         coeffs.push_back( 1.6035414536489318e-11);
         coeffs.push_back(-4.1871266927217782e-14);
      }
      coefficients.push_back( coeffs );
   }
}

FiniteDifferenceBoundedGradientVectorVolume::FiniteDifferenceBoundedGradientVectorVolume( const ScalarField& v, const int N, const GridBox& gb ) :
elem(v),
maxN(N),
gB(gb)
{
   for(int nn=1;nn<11;nn++)
   {
      vector<double> coeffs;
      if( nn==1 )
      {
         coeffs.push_back(0.5);
      }
      else if( nn==2 )
      {
         coeffs.push_back( 0.57142857142857151); 
         coeffs.push_back(-0.071428571428571438);
      }
      else if( nn==3 )
      {
         coeffs.push_back( 0.60810810810810767); 
         coeffs.push_back(-0.1188063063063063);
         coeffs.push_back( 0.010698198198198195);
      }
      else if( nn==4 )
      {
         coeffs.push_back( 0.63039399624765435);
         coeffs.push_back(-0.15064623723160311);
         coeffs.push_back( 0.02101573900354391);
         coeffs.push_back(-0.00076349801959558095);
      }
      else if( nn==5 )
      {
         coeffs.push_back( 0.64535955746772966);
         coeffs.push_back(-0.17321459058935865);
         coeffs.push_back( 0.029559879293859895); 
         coeffs.push_back(-0.0017378184115276803); 
         coeffs.push_back( 3.2972239295226059e-05);
      }
      else if( nn==6 )
      {
         coeffs.push_back( 0.65609985088633171);
         coeffs.push_back(-0.18996590374149147);
         coeffs.push_back( 0.036481469793162628);
         coeffs.push_back(-0.0026999821278442119);
         coeffs.push_back( 8.5528288191602389e-05);
         coeffs.push_back(-9.6309840124388086e-07);
      }
      else if( nn==7 )
      {
         coeffs.push_back( 0.66418180822321082);
         coeffs.push_back(-0.20286437997462345);
         coeffs.push_back( 0.042119156095069443);
         coeffs.push_back(-0.0035797086604024858);
         coeffs.push_back( 0.00014591197349313843);
         coeffs.push_back(-2.8080248121598451e-06);
         coeffs.push_back( 2.0367364328490342e-08);
      }
      else if( nn==8 )
      { 
         coeffs.push_back( 0.67048320765246205);
         coeffs.push_back(-0.213090576806653);
         coeffs.push_back( 0.046767457952830598);
         coeffs.push_back(-0.004362453117772256);
         coeffs.push_back( 0.00020752664200811052);
         coeffs.push_back(-5.2279046772863686e-06);
         coeffs.push_back( 6.598736368669243e-08);
         coeffs.push_back(-3.26811256996574e-10);
      }
      else if( nn==9 )
      {
         coeffs.push_back( 0.67553384734048427);
         coeffs.push_back(-0.22139161101119745);
         coeffs.push_back( 0.050651292408117839);
         coeffs.push_back(-0.0050528139760881441);
         coeffs.push_back( 0.00026710733705086263);
         coeffs.push_back(-7.9566998459801905e-06);
         coeffs.push_back( 1.3328639394368043e-07);
         coeffs.push_back(-1.165378723448066e-09);
         coeffs.push_back( 4.1187471834005603e-12);
      }
      else if( nn==10 )
      {
         coeffs.push_back( 0.67967266644170954);
         coeffs.push_back(-0.22826167257103971);
         coeffs.push_back( 0.053937736604997032);
         coeffs.push_back(-0.005661092015920628);
         coeffs.push_back( 0.00032320263464067625);
         coeffs.push_back(-1.0809865109373177e-05);
         coeffs.push_back( 2.1645208248117699e-07);
         coeffs.push_back(-2.5401252316582755e-09);
         coeffs.push_back( 1.6035414536489318e-11);
         coeffs.push_back(-4.1871266927217782e-14);
      }
      coefficients.push_back( coeffs );
   }
}


const Vector FiniteDifferenceBoundedGradientVectorVolume::eval( const Vector& P ) const
{
   int i,j,k, Nx=maxN, Ny=maxN, Nz=maxN;
   if( gB->isInside(P) )
   {
      gB->getGridIndex( P,i,j,k );
      if( i<Nx ) { Nx = i; }
      if( i+Nx>=(gB->nx()) ){ Nx = gB->nx() - i - 1; }
      if( j<Ny ) { Ny = j; }
      if( j+Ny>=(gB->ny()) ){ Ny = gB->ny() - j - 1; }
      if( k<Nz ) { Nz = k; }
      if( k+Nz>=(gB->nz()) ){ Nz = gB->nz() - k - 1; }
   }
   double fx = 0;
   double fy = 0;
   double fz = 0;
   if( Nx > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Nx-1].size();i++)
      {
         ii++;
         Vector DX(ii*gB->dx(),0,0);
         fx += elem->eval(P+DX)*coefficients[Nx-1][i];
         fx -= elem->eval(P-DX)*coefficients[Nx-1][i];
      }
   }
   else
   {
      if( i == 0 )
      {
         fx = elem->eval(P + Vector(gB->dx(),0,0)) - elem->eval(P);
      }
      else
      {
         fx = elem->eval(P) - elem->eval(P - Vector(gB->dx(),0,0));
      }
   }
   if( Ny > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Ny-1].size();i++)
      {
         ii++;
         Vector DX(0,ii*gB->dy(),0);
         fy += elem->eval(P+DX)*coefficients[Ny-1][i];
         fy -= elem->eval(P-DX)*coefficients[Ny-1][i];
      }
   }
   else
   {
      if( j == 0 )
      {
         fy = elem->eval(P + Vector(0,gB->dy(),0)) - elem->eval(P);
      }
      else
      {
         fy = elem->eval(P) - elem->eval(P - Vector(0,gB->dy(),0));
      }
   }
   if( Nz > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Nz-1].size();i++)
      {
         ii++;
         Vector DX(0,0,ii*gB->dz());
         fz += elem->eval(P+DX)*coefficients[Nz-1][i];
         fz -= elem->eval(P-DX)*coefficients[Nz-1][i];
      }
   }
   else
   {
      if( k == 0 )
      {
         fz = elem->eval(P + Vector(0,0,gB->dz())) - elem->eval(P);
      }
      else
      {
         fz = elem->eval(P) - elem->eval(P - Vector(0,0,gB->dz()));
      }
   }

   return Vector( fx/gB->dx(), fy/gB->dy(), fz/gB->dz() );
}

const Matrix FiniteDifferenceBoundedGradientVectorVolume::grad( const Vector& P ) const 
{
   int i,j,k, Nx=maxN, Ny=maxN, Nz=maxN;
   if( gB->isInside(P) )
   {
      gB->getGridIndex( P,i,j,k );
      if( i<Nx ) { Nx = i; }
      if( i+Nx>=(gB->nx()) ){ Nx = gB->nx() - i - 1; }
      if( j<Ny ) { Ny = j; }
      if( j+Ny>=(gB->ny()) ){ Ny = gB->ny() - j - 1; }
      if( k<Nz ) { Nz = k; }
      if( k+Nz>=(gB->nz()) ){ Nz = gB->nz() - k - 1; }
   }
   Vector fx;
   Vector fy;
   Vector fz;

   if( Nx > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Nx-1].size();i++)
      {
         ii++;
         Vector DX(ii*gB->dx(),0,0);
         fx += eval(P+DX)*coefficients[Nx-1][i];
         fx -= eval(P-DX)*coefficients[Nx-1][i];
      }
   }
   else
   {
      if( i == 0 )
      {
         fx = eval(P + Vector(gB->dx(),0,0)) - eval(P);
      }
      else
      {
         fx = eval(P) - eval(P - Vector(gB->dx(),0,0));
      }
   }
   if( Ny > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Ny-1].size();i++)
      {
         ii++;
         Vector DX(0,ii*gB->dy(),0);
         fy += eval(P+DX)*coefficients[Ny-1][i];
         fy -= eval(P-DX)*coefficients[Ny-1][i];
      }
   }
   else
   {
      if( j == 0 )
      {
         fy = eval(P + Vector(0,gB->dy(),0)) - eval(P);
      }
      else
      {
         fy = eval(P) - eval(P - Vector(0,gB->dy(),0));
      }
   }
   if( Nz > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Nz-1].size();i++)
      {
         ii++;
         Vector DX(0,0,ii*gB->dz());
         fz += eval(P+DX)*coefficients[Ny-1][i];
         fz -= eval(P-DX)*coefficients[Ny-1][i];
      }
   }
   else
   {
      if( k == 0 )
      {
         fz = eval(P + Vector(0,0,gB->dz())) - eval(P);
      }
      else
      {
         fz = eval(P) - eval(P - Vector(0,0,gB->dz()));
      }
   }
   Matrix result( fx/gB->dx(), fy/gB->dy(), fz/gB->dz() );
   return result;
}

















// N is the length of the difference kernel
FiniteDifferenceInteriorGradientVectorVolume::FiniteDifferenceInteriorGradientVectorVolume( Volume<float>* v, const int N, const double dx, const double dy, const double dz, Volume<float>* gb ) :
elem(v),
maxN(N),
dX(dx), dY(dy), dZ(dz),
gB(gb)
{
   for(int nn=1;nn<11;nn++)
   {
      vector<double> coeffs;
      if( nn==1 )
      {
         coeffs.push_back(0.5);
      }
      else if( nn==2 )
      {
         coeffs.push_back( 0.57142857142857151); 
         coeffs.push_back(-0.071428571428571438);
      }
      else if( nn==3 )
      {
         coeffs.push_back( 0.60810810810810767); 
         coeffs.push_back(-0.1188063063063063);
         coeffs.push_back( 0.010698198198198195);
      }
      else if( nn==4 )
      {
         coeffs.push_back( 0.63039399624765435);
         coeffs.push_back(-0.15064623723160311);
         coeffs.push_back( 0.02101573900354391);
         coeffs.push_back(-0.00076349801959558095);
      }
      else if( nn==5 )
      {
         coeffs.push_back( 0.64535955746772966);
         coeffs.push_back(-0.17321459058935865);
         coeffs.push_back( 0.029559879293859895); 
         coeffs.push_back(-0.0017378184115276803); 
         coeffs.push_back( 3.2972239295226059e-05);
      }
      else if( nn==6 )
      {
         coeffs.push_back( 0.65609985088633171);
         coeffs.push_back(-0.18996590374149147);
         coeffs.push_back( 0.036481469793162628);
         coeffs.push_back(-0.0026999821278442119);
         coeffs.push_back( 8.5528288191602389e-05);
         coeffs.push_back(-9.6309840124388086e-07);
      }
      else if( nn==7 )
      {
         coeffs.push_back( 0.66418180822321082);
         coeffs.push_back(-0.20286437997462345);
         coeffs.push_back( 0.042119156095069443);
         coeffs.push_back(-0.0035797086604024858);
         coeffs.push_back( 0.00014591197349313843);
         coeffs.push_back(-2.8080248121598451e-06);
         coeffs.push_back( 2.0367364328490342e-08);
      }
      else if( nn==8 )
      { 
         coeffs.push_back( 0.67048320765246205);
         coeffs.push_back(-0.213090576806653);
         coeffs.push_back( 0.046767457952830598);
         coeffs.push_back(-0.004362453117772256);
         coeffs.push_back( 0.00020752664200811052);
         coeffs.push_back(-5.2279046772863686e-06);
         coeffs.push_back( 6.598736368669243e-08);
         coeffs.push_back(-3.26811256996574e-10);
      }
      else if( nn==9 )
      {
         coeffs.push_back( 0.67553384734048427);
         coeffs.push_back(-0.22139161101119745);
         coeffs.push_back( 0.050651292408117839);
         coeffs.push_back(-0.0050528139760881441);
         coeffs.push_back( 0.00026710733705086263);
         coeffs.push_back(-7.9566998459801905e-06);
         coeffs.push_back( 1.3328639394368043e-07);
         coeffs.push_back(-1.165378723448066e-09);
         coeffs.push_back( 4.1187471834005603e-12);
      }
      else if( nn==10 )
      {
         coeffs.push_back( 0.67967266644170954);
         coeffs.push_back(-0.22826167257103971);
         coeffs.push_back( 0.053937736604997032);
         coeffs.push_back(-0.005661092015920628);
         coeffs.push_back( 0.00032320263464067625);
         coeffs.push_back(-1.0809865109373177e-05);
         coeffs.push_back( 2.1645208248117699e-07);
         coeffs.push_back(-2.5401252316582755e-09);
         coeffs.push_back( 1.6035414536489318e-11);
         coeffs.push_back(-4.1871266927217782e-14);
      }
      coefficients.push_back( coeffs );
   }
}





// N is the length of the difference kernel
FiniteDifferenceInteriorGradientVectorVolume::FiniteDifferenceInteriorGradientVectorVolume( const ScalarField& v, const int N, const double dx, const double dy, const double dz, const ScalarField& gb ) :
elem(v),
maxN(N),
dX(dx), dY(dy), dZ(dz),
gB(gb)
{
   for(int nn=1;nn<11;nn++)
   {
      vector<double> coeffs;
      if( nn==1 )
      {
         coeffs.push_back(0.5);
      }
      else if( nn==2 )
      {
         coeffs.push_back( 0.57142857142857151); 
         coeffs.push_back(-0.071428571428571438);
      }
      else if( nn==3 )
      {
         coeffs.push_back( 0.60810810810810767); 
         coeffs.push_back(-0.1188063063063063);
         coeffs.push_back( 0.010698198198198195);
      }
      else if( nn==4 )
      {
         coeffs.push_back( 0.63039399624765435);
         coeffs.push_back(-0.15064623723160311);
         coeffs.push_back( 0.02101573900354391);
         coeffs.push_back(-0.00076349801959558095);
      }
      else if( nn==5 )
      {
         coeffs.push_back( 0.64535955746772966);
         coeffs.push_back(-0.17321459058935865);
         coeffs.push_back( 0.029559879293859895); 
         coeffs.push_back(-0.0017378184115276803); 
         coeffs.push_back( 3.2972239295226059e-05);
      }
      else if( nn==6 )
      {
         coeffs.push_back( 0.65609985088633171);
         coeffs.push_back(-0.18996590374149147);
         coeffs.push_back( 0.036481469793162628);
         coeffs.push_back(-0.0026999821278442119);
         coeffs.push_back( 8.5528288191602389e-05);
         coeffs.push_back(-9.6309840124388086e-07);
      }
      else if( nn==7 )
      {
         coeffs.push_back( 0.66418180822321082);
         coeffs.push_back(-0.20286437997462345);
         coeffs.push_back( 0.042119156095069443);
         coeffs.push_back(-0.0035797086604024858);
         coeffs.push_back( 0.00014591197349313843);
         coeffs.push_back(-2.8080248121598451e-06);
         coeffs.push_back( 2.0367364328490342e-08);
      }
      else if( nn==8 )
      { 
         coeffs.push_back( 0.67048320765246205);
         coeffs.push_back(-0.213090576806653);
         coeffs.push_back( 0.046767457952830598);
         coeffs.push_back(-0.004362453117772256);
         coeffs.push_back( 0.00020752664200811052);
         coeffs.push_back(-5.2279046772863686e-06);
         coeffs.push_back( 6.598736368669243e-08);
         coeffs.push_back(-3.26811256996574e-10);
      }
      else if( nn==9 )
      {
         coeffs.push_back( 0.67553384734048427);
         coeffs.push_back(-0.22139161101119745);
         coeffs.push_back( 0.050651292408117839);
         coeffs.push_back(-0.0050528139760881441);
         coeffs.push_back( 0.00026710733705086263);
         coeffs.push_back(-7.9566998459801905e-06);
         coeffs.push_back( 1.3328639394368043e-07);
         coeffs.push_back(-1.165378723448066e-09);
         coeffs.push_back( 4.1187471834005603e-12);
      }
      else if( nn==10 )
      {
         coeffs.push_back( 0.67967266644170954);
         coeffs.push_back(-0.22826167257103971);
         coeffs.push_back( 0.053937736604997032);
         coeffs.push_back(-0.005661092015920628);
         coeffs.push_back( 0.00032320263464067625);
         coeffs.push_back(-1.0809865109373177e-05);
         coeffs.push_back( 2.1645208248117699e-07);
         coeffs.push_back(-2.5401252316582755e-09);
         coeffs.push_back( 1.6035414536489318e-11);
         coeffs.push_back(-4.1871266927217782e-14);
      }
      coefficients.push_back( coeffs );
   }
}



const Vector FiniteDifferenceInteriorGradientVectorVolume::eval( const Vector& P ) const
{
   int Nx=0, Ny = 0, Nz = 0;
   if( gB->eval(P) < 0.0 )
   {
      return Vector(0,0,0);
   }
   else
   {
      for( int l=1;l<=maxN;l++ )
      {
         Vector X = P + Vector( l*dX,0,0 );
         if( gB->eval(X) >= 0.0 )
         {
            Nx = l;
            break;
         }
      }
      for( int l=1;l<=maxN;l++ )
      {
         Vector X = P + Vector( 0,l*dY,0 );
         if( gB->eval(X) >= 0.0 )
         {
            Ny = l;
            break;
         }
      }
      for( int l=1;l<=maxN;l++ )
      {
         Vector X = P + Vector( 0,0,l*dZ );
         if( gB->eval(X) >= 0.0 )
         {
            Nz = l;
            break;
         }
      }
   }

   double fx = 0;
   double fy = 0;
   double fz = 0;
   if( Nx > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Nx-1].size();i++)
      {
         ii++;
         Vector DX(ii*dX,0,0);
         fx += elem->eval(P+DX)*coefficients[Nx-1][i];
         fx -= elem->eval(P-DX)*coefficients[Nx-1][i];
      }
   }
   else
   {
      if( gB->eval(P+Vector(dX,0,0)) >= 0.0 )
      {
         fx = elem->eval(P + Vector(dX,0,0)) - elem->eval(P);
      }
      else
      {
         fx = elem->eval(P) - elem->eval(P - Vector(dX,0,0));
      }
   }
   if( Ny > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Ny-1].size();i++)
      {
         ii++;
         Vector DX(0,ii*dY,0);
         fy += elem->eval(P+DX)*coefficients[Ny-1][i];
         fy -= elem->eval(P-DX)*coefficients[Ny-1][i];
      }
   }
   else
   {
      if( gB->eval(P+Vector(0,dY,0)) >= 0.0 )
      {
         fy = elem->eval(P + Vector(0,dY,0)) - elem->eval(P);
      }
      else
      {
         fy = elem->eval(P) - elem->eval(P - Vector(0,dY,0));
      }
   }
   if( Nz > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Nz-1].size();i++)
      {
         ii++;
         Vector DX(0,0,ii*dZ);
         fz += elem->eval(P+DX)*coefficients[Nz-1][i];
         fz -= elem->eval(P-DX)*coefficients[Nz-1][i];
      }
   }
   else
   {
      if( gB->eval(P+Vector(0,0,dZ)) >= 0.0 )
      {
         fz = elem->eval(P + Vector(0,0,dZ)) - elem->eval(P);
      }
      else
      {
         fz = elem->eval(P) - elem->eval(P - Vector(0,0,dZ));
      }
   }

   return Vector( fx/dX, fy/dY, fz/dZ );
}

const Matrix FiniteDifferenceInteriorGradientVectorVolume::grad( const Vector& P ) const 
{
  int Nx=maxN, Ny=maxN, Nz=maxN;
   if( gB->eval(P) < 0.0 )
   {
      Nx = Ny = Nz = 0;
   }
   else
   {
      for( int l=1;l<=maxN;l++ )
      {
         Vector X = P + Vector( l*dX,0,0 );
         if( gB->eval(X) < 0.0 )
         {
            Nx = l;
            break;
         }
      }
      for( int l=1;l<=maxN;l++ )
      {
         Vector X = P + Vector( 0,l*dY,0 );
         if( gB->eval(X) < 0.0 )
         {
            Ny = l;
            break;
         }
      }
      for( int l=1;l<=maxN;l++ )
      {
         Vector X = P + Vector( 0,0,l*dZ );
         if( gB->eval(X) < 0.0 )
         {
            Nz = l;
            break;
         }
      }
   }
   Vector fx;
   Vector fy;
   Vector fz;

   if( Nx > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Nx-1].size();i++)
      {
         ii++;
         Vector DX(ii*dX,0,0);
         fx += eval(P+DX)*coefficients[Nx-1][i];
         fx -= eval(P-DX)*coefficients[Nx-1][i];
      }
   }
   else
   {
      if( gB->eval(P+Vector(dX,0,0)) >= 0.0 )
      {
         fx = eval(P + Vector(dX,0,0)) - eval(P);
      }
      else
      {
         fx = eval(P) - eval(P - Vector(dX,0,0));
      }
   }
   if( Ny > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Ny-1].size();i++)
      {
         ii++;
         Vector DX(0,ii*dY,0);
         fy += eval(P+DX)*coefficients[Ny-1][i];
         fy -= eval(P-DX)*coefficients[Ny-1][i];
      }
   }
   else
   {
      if( gB->eval(P+Vector(0,dY,0)) >= 0.0 )
      {
         fy = eval(P + Vector(0,dY,0)) - eval(P);
      }
      else
      {
         fy = eval(P) - eval(P - Vector(0,dY,0));
      }
   }
   if( Nz > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Nz-1].size();i++)
      {
         ii++;
         Vector DX(0,0,ii*dZ);
         fz += eval(P+DX)*coefficients[Ny-1][i];
         fz -= eval(P-DX)*coefficients[Ny-1][i];
      }
   }
   else
   {
      if( gB->eval(P+Vector(0,0,dZ)) >= 0.0 )
      {
         fz = eval(P + Vector(0,0,dZ)) - eval(P);
      }
      else
      {
         fz = eval(P) - eval(P - Vector(0,0,dZ));
      }
   }
   Matrix result( fx/dX, fy/dY, fz/dZ );
   return result;
}











CrossProductVectorVolume::CrossProductVectorVolume( Volume<Vector>* v1, Volume<Vector>* v2 ) :
      elem1(v1),
      elem2(v2),
      step(0.001)
    {}

CrossProductVectorVolume::CrossProductVectorVolume( const VectorField& v1, const VectorField& v2 ) :
      elem1(v1),
      elem2(v2),
      step(0.001)
    {}


    const Vector CrossProductVectorVolume::eval( const Vector& P ) const { return (elem1->eval(P))^(elem2->eval(P)); }

    // evaluate gradient by finite difference
    // size of difference is controlled via the "dx" attribute
    const Matrix CrossProductVectorVolume::grad( const Vector& P ) const 
    {
       float dx = step;
       Vector e = eval(P);
       Vector ex = eval(P+Vector(dx,0,0));
       Vector ey = eval(P+Vector(0,dx,0));
       Vector ez = eval(P+Vector(0,0,dx));
       Matrix result( ex-e, ey-e, ez-e );
       result /= dx;
       return result.transpose();
    }








DivergenceVectorVolume::DivergenceVectorVolume( Volume<Vector>* v, const float dx ) :
      elem(v),
      step(dx)
    {}

DivergenceVectorVolume::DivergenceVectorVolume( const VectorField& v, const float dx ) :
      elem(v),
      step(dx)
    {}


    const float DivergenceVectorVolume::eval( const Vector& P ) const 
    {
       Matrix m = elem->grad(P);
       return trace(m);
       /*
       Vector e0 = elem->eval(P);
       float dx = attribute("dx");
       Vector g1 = (  elem->eval(P+Vector(dx,0,0)) - e0 )/dx;
       Vector g2 = (  elem->eval(P+Vector(0,dx,0)) - e0 )/dx;
       Vector g3 = (  elem->eval(P+Vector(0,0,dx)) - e0 )/dx;
       return g1[0] + g2[1] + g3[2];
       */
    }

    // evaluate gradient by finite difference
    // size of difference is controlled via the "dx" attribute
    const Vector DivergenceVectorVolume::grad( const Vector& P ) const 
    {
       float dx = step;
       float e = eval(P);
       float ex = eval(P+Vector(dx,0,0));
       float ey = eval(P+Vector(0,dx,0));
       float ez = eval(P+Vector(0,0,dx));
       Vector result( ex-e, ey-e, ez-e );
       result /= dx;
       return result;
    }




FiniteDifferenceDivergenceVectorVolume::FiniteDifferenceDivergenceVectorVolume( Volume<Vector>* v, const int N, const float dx, const float dy, const float dz ) :
elem(v),
dX(dx), dY(dy), dZ(dz)
{
   if( N <= 1 )
   {
      coefficients.push_back(0.5);
   }
   else if( N == 2 )
   {
      coefficients.push_back( 0.57142857142857151); 
      coefficients.push_back(-0.071428571428571438);
   }
   else if( N==3 )
   {
      coefficients.push_back( 0.60810810810810767); 
      coefficients.push_back(-0.1188063063063063);
      coefficients.push_back( 0.010698198198198195);
   }
   else if( N==4 )
   {
      coefficients.push_back( 0.63039399624765435);
      coefficients.push_back(-0.15064623723160311);
      coefficients.push_back( 0.02101573900354391);
      coefficients.push_back(-0.00076349801959558095);
   }
   else if( N==5 )
   {
      coefficients.push_back( 0.64535955746772966);
      coefficients.push_back(-0.17321459058935865);
      coefficients.push_back( 0.029559879293859895); 
      coefficients.push_back(-0.0017378184115276803); 
      coefficients.push_back( 3.2972239295226059e-05);
   }
   else if( N==6 )
   {
      coefficients.push_back( 0.65609985088633171);
      coefficients.push_back(-0.18996590374149147);
      coefficients.push_back( 0.036481469793162628);
      coefficients.push_back(-0.0026999821278442119);
      coefficients.push_back( 8.5528288191602389e-05);
      coefficients.push_back(-9.6309840124388086e-07);
   }
   else if( N==7 )
   {
      coefficients.push_back( 0.66418180822321082);
      coefficients.push_back(-0.20286437997462345);
      coefficients.push_back( 0.042119156095069443);
      coefficients.push_back(-0.0035797086604024858);
      coefficients.push_back( 0.00014591197349313843);
      coefficients.push_back(-2.8080248121598451e-06);
      coefficients.push_back( 2.0367364328490342e-08);
   }
   else if( N==8 )
   {
      coefficients.push_back( 0.67048320765246205);
      coefficients.push_back(-0.213090576806653);
      coefficients.push_back( 0.046767457952830598);
      coefficients.push_back(-0.004362453117772256);
      coefficients.push_back( 0.00020752664200811052);
      coefficients.push_back(-5.2279046772863686e-06);
      coefficients.push_back( 6.598736368669243e-08);
      coefficients.push_back(-3.26811256996574e-10);
   }
   else if( N==9 )
   {
      coefficients.push_back( 0.67553384734048427);
      coefficients.push_back(-0.22139161101119745);
      coefficients.push_back( 0.050651292408117839);
      coefficients.push_back(-0.0050528139760881441);
      coefficients.push_back( 0.00026710733705086263);
      coefficients.push_back(-7.9566998459801905e-06);
      coefficients.push_back( 1.3328639394368043e-07);
      coefficients.push_back(-1.165378723448066e-09);
      coefficients.push_back( 4.1187471834005603e-12);
   }
   else if( N>=10 )
   {
      coefficients.push_back( 0.67967266644170954);
      coefficients.push_back(-0.22826167257103971);
      coefficients.push_back( 0.053937736604997032);
      coefficients.push_back(-0.005661092015920628);
      coefficients.push_back( 0.00032320263464067625);
      coefficients.push_back(-1.0809865109373177e-05);
      coefficients.push_back( 2.1645208248117699e-07);
      coefficients.push_back(-2.5401252316582755e-09);
      coefficients.push_back( 1.6035414536489318e-11);
      coefficients.push_back(-4.1871266927217782e-14);
   }
}


FiniteDifferenceDivergenceVectorVolume::FiniteDifferenceDivergenceVectorVolume( const VectorField& v, const int N, const float dx, const float dy, const float dz ) :
elem(v),
dX(dx), dY(dy), dZ(dz)
{
   if( N <= 1 )
   {
      coefficients.push_back(0.5);
   }
   else if( N == 2 )
   {
      coefficients.push_back( 0.57142857142857151); 
      coefficients.push_back(-0.071428571428571438);
   }
   else if( N==3 )
   {
      coefficients.push_back( 0.60810810810810767); 
      coefficients.push_back(-0.1188063063063063);
      coefficients.push_back( 0.010698198198198195);
   }
   else if( N==4 )
   {
      coefficients.push_back( 0.63039399624765435);
      coefficients.push_back(-0.15064623723160311);
      coefficients.push_back( 0.02101573900354391);
      coefficients.push_back(-0.00076349801959558095);
   }
   else if( N==5 )
   {
      coefficients.push_back( 0.64535955746772966);
      coefficients.push_back(-0.17321459058935865);
      coefficients.push_back( 0.029559879293859895); 
      coefficients.push_back(-0.0017378184115276803); 
      coefficients.push_back( 3.2972239295226059e-05);
   }
   else if( N==6 )
   {
      coefficients.push_back( 0.65609985088633171);
      coefficients.push_back(-0.18996590374149147);
      coefficients.push_back( 0.036481469793162628);
      coefficients.push_back(-0.0026999821278442119);
      coefficients.push_back( 8.5528288191602389e-05);
      coefficients.push_back(-9.6309840124388086e-07);
   }
   else if( N==7 )
   {
      coefficients.push_back( 0.66418180822321082);
      coefficients.push_back(-0.20286437997462345);
      coefficients.push_back( 0.042119156095069443);
      coefficients.push_back(-0.0035797086604024858);
      coefficients.push_back( 0.00014591197349313843);
      coefficients.push_back(-2.8080248121598451e-06);
      coefficients.push_back( 2.0367364328490342e-08);
   }
   else if( N==8 )
   {
      coefficients.push_back( 0.67048320765246205);
      coefficients.push_back(-0.213090576806653);
      coefficients.push_back( 0.046767457952830598);
      coefficients.push_back(-0.004362453117772256);
      coefficients.push_back( 0.00020752664200811052);
      coefficients.push_back(-5.2279046772863686e-06);
      coefficients.push_back( 6.598736368669243e-08);
      coefficients.push_back(-3.26811256996574e-10);
   }
   else if( N==9 )
   {
      coefficients.push_back( 0.67553384734048427);
      coefficients.push_back(-0.22139161101119745);
      coefficients.push_back( 0.050651292408117839);
      coefficients.push_back(-0.0050528139760881441);
      coefficients.push_back( 0.00026710733705086263);
      coefficients.push_back(-7.9566998459801905e-06);
      coefficients.push_back( 1.3328639394368043e-07);
      coefficients.push_back(-1.165378723448066e-09);
      coefficients.push_back( 4.1187471834005603e-12);
   }
   else if( N>=10 )
   {
      coefficients.push_back( 0.67967266644170954);
      coefficients.push_back(-0.22826167257103971);
      coefficients.push_back( 0.053937736604997032);
      coefficients.push_back(-0.005661092015920628);
      coefficients.push_back( 0.00032320263464067625);
      coefficients.push_back(-1.0809865109373177e-05);
      coefficients.push_back( 2.1645208248117699e-07);
      coefficients.push_back(-2.5401252316582755e-09);
      coefficients.push_back( 1.6035414536489318e-11);
      coefficients.push_back(-4.1871266927217782e-14);
   }
}


const float FiniteDifferenceDivergenceVectorVolume::eval( const Vector& P ) const 
{
   Vector fx;
   Vector fy;
   Vector fz;
   int ii = 0;
   for( size_t i=0;i<coefficients.size();i++)
   {
      ii++;
      Vector DX(ii*dX,0,0);
      fx += elem->eval(P+DX)*coefficients[i];
      fx -= elem->eval(P-DX)*coefficients[i];
      Vector DY(0,ii*dY,0);
      fy += elem->eval(P+DY)*coefficients[i];
      fy -= elem->eval(P-DY)*coefficients[i];
      Vector DZ(0,0,ii*dZ);
      fz += elem->eval(P+DZ)*coefficients[i];
      fz -= elem->eval(P-DZ)*coefficients[i];
   }
   return (fx.X()/dX + fy.Y()/dY + fz.Z()/dZ);
}

const Vector FiniteDifferenceDivergenceVectorVolume::grad( const Vector& P ) const 
{
   double fx = 0;
   double fy = 0;
   double fz = 0;
   int ii = 0;
   for( size_t i=0;i<coefficients.size();i++)
   {
      ii++;
      Vector DX(ii*dX,0,0);
      fx += eval(P+DX)*coefficients[i];
      fx -= eval(P-DX)*coefficients[i];
      Vector DY(0,ii*dY,0);
      fy += eval(P+DY)*coefficients[i];
      fy -= eval(P-DY)*coefficients[i];
      Vector DZ(0,0,ii*dZ);
      fz += eval(P+DZ)*coefficients[i];
      fz -= eval(P-DZ)*coefficients[i];
   }
   return Vector(fx/dX, fy/dY, fz/dZ);
}








// N is the length of the difference kernel
FiniteDifferenceInteriorDivergenceVectorVolume::FiniteDifferenceInteriorDivergenceVectorVolume( Volume<Vector>* v, const int N, const double dx, const double dy, const double dz, Volume<float>* gb ) :
elem(v),
maxN(N),
dX(dx), dY(dy), dZ(dz),
gB(gb)
{
   for(int nn=1;nn<11;nn++)
   {
      vector<double> coeffs;
      if( nn==1 )
      {
         coeffs.push_back(0.5);
      }
      else if( nn==2 )
      {
         coeffs.push_back( 0.57142857142857151); 
         coeffs.push_back(-0.071428571428571438);
      }
      else if( nn==3 )
      {
         coeffs.push_back( 0.60810810810810767); 
         coeffs.push_back(-0.1188063063063063);
         coeffs.push_back( 0.010698198198198195);
      }
      else if( nn==4 )
      {
         coeffs.push_back( 0.63039399624765435);
         coeffs.push_back(-0.15064623723160311);
         coeffs.push_back( 0.02101573900354391);
         coeffs.push_back(-0.00076349801959558095);
      }
      else if( nn==5 )
      {
         coeffs.push_back( 0.64535955746772966);
         coeffs.push_back(-0.17321459058935865);
         coeffs.push_back( 0.029559879293859895); 
         coeffs.push_back(-0.0017378184115276803); 
         coeffs.push_back( 3.2972239295226059e-05);
      }
      else if( nn==6 )
      {
         coeffs.push_back( 0.65609985088633171);
         coeffs.push_back(-0.18996590374149147);
         coeffs.push_back( 0.036481469793162628);
         coeffs.push_back(-0.0026999821278442119);
         coeffs.push_back( 8.5528288191602389e-05);
         coeffs.push_back(-9.6309840124388086e-07);
      }
      else if( nn==7 )
      {
         coeffs.push_back( 0.66418180822321082);
         coeffs.push_back(-0.20286437997462345);
         coeffs.push_back( 0.042119156095069443);
         coeffs.push_back(-0.0035797086604024858);
         coeffs.push_back( 0.00014591197349313843);
         coeffs.push_back(-2.8080248121598451e-06);
         coeffs.push_back( 2.0367364328490342e-08);
      }
      else if( nn==8 )
      { 
         coeffs.push_back( 0.67048320765246205);
         coeffs.push_back(-0.213090576806653);
         coeffs.push_back( 0.046767457952830598);
         coeffs.push_back(-0.004362453117772256);
         coeffs.push_back( 0.00020752664200811052);
         coeffs.push_back(-5.2279046772863686e-06);
         coeffs.push_back( 6.598736368669243e-08);
         coeffs.push_back(-3.26811256996574e-10);
      }
      else if( nn==9 )
      {
         coeffs.push_back( 0.67553384734048427);
         coeffs.push_back(-0.22139161101119745);
         coeffs.push_back( 0.050651292408117839);
         coeffs.push_back(-0.0050528139760881441);
         coeffs.push_back( 0.00026710733705086263);
         coeffs.push_back(-7.9566998459801905e-06);
         coeffs.push_back( 1.3328639394368043e-07);
         coeffs.push_back(-1.165378723448066e-09);
         coeffs.push_back( 4.1187471834005603e-12);
      }
      else if( nn==10 )
      {
         coeffs.push_back( 0.67967266644170954);
         coeffs.push_back(-0.22826167257103971);
         coeffs.push_back( 0.053937736604997032);
         coeffs.push_back(-0.005661092015920628);
         coeffs.push_back( 0.00032320263464067625);
         coeffs.push_back(-1.0809865109373177e-05);
         coeffs.push_back( 2.1645208248117699e-07);
         coeffs.push_back(-2.5401252316582755e-09);
         coeffs.push_back( 1.6035414536489318e-11);
         coeffs.push_back(-4.1871266927217782e-14);
      }
      coefficients.push_back( coeffs );
   }
}





// N is the length of the difference kernel
FiniteDifferenceInteriorDivergenceVectorVolume::FiniteDifferenceInteriorDivergenceVectorVolume( const VectorField& v, const int N, const double dx, const double dy, const double dz, const ScalarField& gb ) :
elem(v),
maxN(N),
dX(dx), dY(dy), dZ(dz),
gB(gb)
{
   for(int nn=1;nn<11;nn++)
   {
      vector<double> coeffs;
      if( nn==1 )
      {
         coeffs.push_back(0.5);
      }
      else if( nn==2 )
      {
         coeffs.push_back( 0.57142857142857151); 
         coeffs.push_back(-0.071428571428571438);
      }
      else if( nn==3 )
      {
         coeffs.push_back( 0.60810810810810767); 
         coeffs.push_back(-0.1188063063063063);
         coeffs.push_back( 0.010698198198198195);
      }
      else if( nn==4 )
      {
         coeffs.push_back( 0.63039399624765435);
         coeffs.push_back(-0.15064623723160311);
         coeffs.push_back( 0.02101573900354391);
         coeffs.push_back(-0.00076349801959558095);
      }
      else if( nn==5 )
      {
         coeffs.push_back( 0.64535955746772966);
         coeffs.push_back(-0.17321459058935865);
         coeffs.push_back( 0.029559879293859895); 
         coeffs.push_back(-0.0017378184115276803); 
         coeffs.push_back( 3.2972239295226059e-05);
      }
      else if( nn==6 )
      {
         coeffs.push_back( 0.65609985088633171);
         coeffs.push_back(-0.18996590374149147);
         coeffs.push_back( 0.036481469793162628);
         coeffs.push_back(-0.0026999821278442119);
         coeffs.push_back( 8.5528288191602389e-05);
         coeffs.push_back(-9.6309840124388086e-07);
      }
      else if( nn==7 )
      {
         coeffs.push_back( 0.66418180822321082);
         coeffs.push_back(-0.20286437997462345);
         coeffs.push_back( 0.042119156095069443);
         coeffs.push_back(-0.0035797086604024858);
         coeffs.push_back( 0.00014591197349313843);
         coeffs.push_back(-2.8080248121598451e-06);
         coeffs.push_back( 2.0367364328490342e-08);
      }
      else if( nn==8 )
      { 
         coeffs.push_back( 0.67048320765246205);
         coeffs.push_back(-0.213090576806653);
         coeffs.push_back( 0.046767457952830598);
         coeffs.push_back(-0.004362453117772256);
         coeffs.push_back( 0.00020752664200811052);
         coeffs.push_back(-5.2279046772863686e-06);
         coeffs.push_back( 6.598736368669243e-08);
         coeffs.push_back(-3.26811256996574e-10);
      }
      else if( nn==9 )
      {
         coeffs.push_back( 0.67553384734048427);
         coeffs.push_back(-0.22139161101119745);
         coeffs.push_back( 0.050651292408117839);
         coeffs.push_back(-0.0050528139760881441);
         coeffs.push_back( 0.00026710733705086263);
         coeffs.push_back(-7.9566998459801905e-06);
         coeffs.push_back( 1.3328639394368043e-07);
         coeffs.push_back(-1.165378723448066e-09);
         coeffs.push_back( 4.1187471834005603e-12);
      }
      else if( nn==10 )
      {
         coeffs.push_back( 0.67967266644170954);
         coeffs.push_back(-0.22826167257103971);
         coeffs.push_back( 0.053937736604997032);
         coeffs.push_back(-0.005661092015920628);
         coeffs.push_back( 0.00032320263464067625);
         coeffs.push_back(-1.0809865109373177e-05);
         coeffs.push_back( 2.1645208248117699e-07);
         coeffs.push_back(-2.5401252316582755e-09);
         coeffs.push_back( 1.6035414536489318e-11);
         coeffs.push_back(-4.1871266927217782e-14);
      }
      coefficients.push_back( coeffs );
   }
}



const float FiniteDifferenceInteriorDivergenceVectorVolume::eval( const Vector& P ) const
{
  int Nx=maxN, Ny=maxN, Nz=maxN;
   if( gB->eval(P) < 0.0 )
   {
      Nx = Ny = Nz = 0;
   }
   else
   {
      for( int l=1;l<=maxN;l++ )
      {
         Vector X = P + Vector( l*dX,0,0 );
         if( gB->eval(X) < 0.0 )
         {
            Nx = l;
            break;
         }
      }
      for( int l=1;l<=maxN;l++ )
      {
         Vector X = P + Vector( 0,l*dY,0 );
         if( gB->eval(X) < 0.0 )
         {
            Ny = l;
            break;
         }
      }
      for( int l=1;l<=maxN;l++ )
      {
         Vector X = P + Vector( 0,0,l*dZ );
         if( gB->eval(X) < 0.0 )
         {
            Nz = l;
            break;
         }
      }
   }

   Vector fx;
   Vector fy;
   Vector fz;
   if( Nx > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Nx-1].size();i++)
      {
         ii++;
         Vector DX(ii*dX,0,0);
         fx += elem->eval(P+DX)*coefficients[Nx-1][i];
         fx -= elem->eval(P-DX)*coefficients[Nx-1][i];
      }
   }
   else
   {
      if( gB->eval(P+Vector(dX,0,0)) >= 0.0 )
      {
         fx = elem->eval(P + Vector(dX,0,0)) - elem->eval(P);
      }
      else
      {
         fx = elem->eval(P) - elem->eval(P - Vector(dX,0,0));
      }
   }
   if( Ny > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Ny-1].size();i++)
      {
         ii++;
         Vector DX(0,ii*dY,0);
         fy += elem->eval(P+DX)*coefficients[Ny-1][i];
         fy -= elem->eval(P-DX)*coefficients[Ny-1][i];
      }
   }
   else
   {
      if( gB->eval(P+Vector(0,dY,0)) >= 0.0 )
      {
         fy = elem->eval(P + Vector(0,dY,0)) - elem->eval(P);
      }
      else
      {
         fy = elem->eval(P) - elem->eval(P - Vector(0,dY,0));
      }
   }
   if( Nz > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Nz-1].size();i++)
      {
         ii++;
         Vector DX(0,0,ii*dZ);
         fz += elem->eval(P+DX)*coefficients[Nz-1][i];
         fz -= elem->eval(P-DX)*coefficients[Nz-1][i];
      }
   }
   else
   {
      if( gB->eval(P+Vector(0,0,dZ)) >= 0.0 )
      {
         fz = elem->eval(P + Vector(0,0,dZ)) - elem->eval(P);
      }
      else
      {
         fz = elem->eval(P) - elem->eval(P - Vector(0,0,dZ));
      }
   }

   return (fx.X()/dX + fy.Y()/dY + fz.Z()/dZ);
}

const Vector FiniteDifferenceInteriorDivergenceVectorVolume::grad( const Vector& P ) const 
{
  int Nx=maxN, Ny=maxN, Nz=maxN;
   if( gB->eval(P) < 0.0 )
   {
      Nx = Ny = Nz = 0;
   }
   else
   {
      for( int l=1;l<=maxN;l++ )
      {
         Vector X = P + Vector( l*dX,0,0 );
         if( gB->eval(X) < 0.0 )
         {
            Nx = l;
            break;
         }
      }
      for( int l=1;l<=maxN;l++ )
      {
         Vector X = P + Vector( 0,l*dY,0 );
         if( gB->eval(X) < 0.0 )
         {
            Ny = l;
            break;
         }
      }
      for( int l=1;l<=maxN;l++ )
      {
         Vector X = P + Vector( 0,0,l*dZ );
         if( gB->eval(X) < 0.0 )
         {
            Nz = l;
            break;
         }
      }
   }
   double fx = 0;
   double fy = 0;
   double fz = 0;
   if( Nx > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Nx-1].size();i++)
      {
         ii++;
         Vector DX(ii*dX,0,0);
         fx += eval(P+DX)*coefficients[Nx-1][i];
         fx -= eval(P-DX)*coefficients[Nx-1][i];
      }
   }
   else
   {
      if( gB->eval(P+Vector(dX,0,0)) >= 0.0 )
      {
         fx = eval(P + Vector(dX,0,0)) - eval(P);
      }
      else
      {
         fx = eval(P) - eval(P - Vector(dX,0,0));
      }
   }
   if( Ny > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Ny-1].size();i++)
      {
         ii++;
         Vector DX(0,ii*dY,0);
         fy += eval(P+DX)*coefficients[Ny-1][i];
         fy -= eval(P-DX)*coefficients[Ny-1][i];
      }
   }
   else
   {
      if( gB->eval(P+Vector(0,dY,0)) >= 0.0 )
      {
         fy = eval(P + Vector(0,dY,0)) - eval(P);
      }
      else
      {
         fy = eval(P) - eval(P - Vector(0,dY,0));
      }
   }
   if( Nz > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Nz-1].size();i++)
      {
         ii++;
         Vector DX(0,0,ii*dZ);
         fz += eval(P+DX)*coefficients[Ny-1][i];
         fz -= eval(P-DX)*coefficients[Ny-1][i];
      }
   }
   else
   {
      if( gB->eval(P+Vector(0,0,dZ)) >= 0.0 )
      {
         fz = eval(P + Vector(0,0,dZ)) - eval(P);
      }
      else
      {
         fz = eval(P) - eval(P - Vector(0,0,dZ));
      }
   }
   return Vector(fx/dX, fy/dY, fz/dZ);
}












































// N is the length of the difference kernel
FiniteDifferenceBoundedDivergenceVectorVolume::FiniteDifferenceBoundedDivergenceVectorVolume( Volume<Vector>* v, const int N, const GridBox& gb ) :
elem(v),
maxN(N),
gB(gb)
{
   for(int nn=1;nn<11;nn++)
   {
      vector<double> coeffs;
      if( nn==1 )
      {
         coeffs.push_back(0.5);
      }
      else if( nn==2 )
      {
         coeffs.push_back( 0.57142857142857151); 
         coeffs.push_back(-0.071428571428571438);
      }
      else if( nn==3 )
      {
         coeffs.push_back( 0.60810810810810767); 
         coeffs.push_back(-0.1188063063063063);
         coeffs.push_back( 0.010698198198198195);
      }
      else if( nn==4 )
      {
         coeffs.push_back( 0.63039399624765435);
         coeffs.push_back(-0.15064623723160311);
         coeffs.push_back( 0.02101573900354391);
         coeffs.push_back(-0.00076349801959558095);
      }
      else if( nn==5 )
      {
         coeffs.push_back( 0.64535955746772966);
         coeffs.push_back(-0.17321459058935865);
         coeffs.push_back( 0.029559879293859895); 
         coeffs.push_back(-0.0017378184115276803); 
         coeffs.push_back( 3.2972239295226059e-05);
      }
      else if( nn==6 )
      {
         coeffs.push_back( 0.65609985088633171);
         coeffs.push_back(-0.18996590374149147);
         coeffs.push_back( 0.036481469793162628);
         coeffs.push_back(-0.0026999821278442119);
         coeffs.push_back( 8.5528288191602389e-05);
         coeffs.push_back(-9.6309840124388086e-07);
      }
      else if( nn==7 )
      {
         coeffs.push_back( 0.66418180822321082);
         coeffs.push_back(-0.20286437997462345);
         coeffs.push_back( 0.042119156095069443);
         coeffs.push_back(-0.0035797086604024858);
         coeffs.push_back( 0.00014591197349313843);
         coeffs.push_back(-2.8080248121598451e-06);
         coeffs.push_back( 2.0367364328490342e-08);
      }
      else if( nn==8 )
      { 
         coeffs.push_back( 0.67048320765246205);
         coeffs.push_back(-0.213090576806653);
         coeffs.push_back( 0.046767457952830598);
         coeffs.push_back(-0.004362453117772256);
         coeffs.push_back( 0.00020752664200811052);
         coeffs.push_back(-5.2279046772863686e-06);
         coeffs.push_back( 6.598736368669243e-08);
         coeffs.push_back(-3.26811256996574e-10);
      }
      else if( nn==9 )
      {
         coeffs.push_back( 0.67553384734048427);
         coeffs.push_back(-0.22139161101119745);
         coeffs.push_back( 0.050651292408117839);
         coeffs.push_back(-0.0050528139760881441);
         coeffs.push_back( 0.00026710733705086263);
         coeffs.push_back(-7.9566998459801905e-06);
         coeffs.push_back( 1.3328639394368043e-07);
         coeffs.push_back(-1.165378723448066e-09);
         coeffs.push_back( 4.1187471834005603e-12);
      }
      else if( nn==10 )
      {
         coeffs.push_back( 0.67967266644170954);
         coeffs.push_back(-0.22826167257103971);
         coeffs.push_back( 0.053937736604997032);
         coeffs.push_back(-0.005661092015920628);
         coeffs.push_back( 0.00032320263464067625);
         coeffs.push_back(-1.0809865109373177e-05);
         coeffs.push_back( 2.1645208248117699e-07);
         coeffs.push_back(-2.5401252316582755e-09);
         coeffs.push_back( 1.6035414536489318e-11);
         coeffs.push_back(-4.1871266927217782e-14);
      }
      coefficients.push_back( coeffs );
   }
}

FiniteDifferenceBoundedDivergenceVectorVolume::FiniteDifferenceBoundedDivergenceVectorVolume( const VectorField& v, const int N, const GridBox& gb ) :
elem(v),
maxN(N),
gB(gb)
{
   for(int nn=1;nn<11;nn++)
   {
      vector<double> coeffs;
      if( nn==1 )
      {
         coeffs.push_back(0.5);
      }
      else if( nn==2 )
      {
         coeffs.push_back( 0.57142857142857151); 
         coeffs.push_back(-0.071428571428571438);
      }
      else if( nn==3 )
      {
         coeffs.push_back( 0.60810810810810767); 
         coeffs.push_back(-0.1188063063063063);
         coeffs.push_back( 0.010698198198198195);
      }
      else if( nn==4 )
      {
         coeffs.push_back( 0.63039399624765435);
         coeffs.push_back(-0.15064623723160311);
         coeffs.push_back( 0.02101573900354391);
         coeffs.push_back(-0.00076349801959558095);
      }
      else if( nn==5 )
      {
         coeffs.push_back( 0.64535955746772966);
         coeffs.push_back(-0.17321459058935865);
         coeffs.push_back( 0.029559879293859895); 
         coeffs.push_back(-0.0017378184115276803); 
         coeffs.push_back( 3.2972239295226059e-05);
      }
      else if( nn==6 )
      {
         coeffs.push_back( 0.65609985088633171);
         coeffs.push_back(-0.18996590374149147);
         coeffs.push_back( 0.036481469793162628);
         coeffs.push_back(-0.0026999821278442119);
         coeffs.push_back( 8.5528288191602389e-05);
         coeffs.push_back(-9.6309840124388086e-07);
      }
      else if( nn==7 )
      {
         coeffs.push_back( 0.66418180822321082);
         coeffs.push_back(-0.20286437997462345);
         coeffs.push_back( 0.042119156095069443);
         coeffs.push_back(-0.0035797086604024858);
         coeffs.push_back( 0.00014591197349313843);
         coeffs.push_back(-2.8080248121598451e-06);
         coeffs.push_back( 2.0367364328490342e-08);
      }
      else if( nn==8 )
      { 
         coeffs.push_back( 0.67048320765246205);
         coeffs.push_back(-0.213090576806653);
         coeffs.push_back( 0.046767457952830598);
         coeffs.push_back(-0.004362453117772256);
         coeffs.push_back( 0.00020752664200811052);
         coeffs.push_back(-5.2279046772863686e-06);
         coeffs.push_back( 6.598736368669243e-08);
         coeffs.push_back(-3.26811256996574e-10);
      }
      else if( nn==9 )
      {
         coeffs.push_back( 0.67553384734048427);
         coeffs.push_back(-0.22139161101119745);
         coeffs.push_back( 0.050651292408117839);
         coeffs.push_back(-0.0050528139760881441);
         coeffs.push_back( 0.00026710733705086263);
         coeffs.push_back(-7.9566998459801905e-06);
         coeffs.push_back( 1.3328639394368043e-07);
         coeffs.push_back(-1.165378723448066e-09);
         coeffs.push_back( 4.1187471834005603e-12);
      }
      else if( nn==10 )
      {
         coeffs.push_back( 0.67967266644170954);
         coeffs.push_back(-0.22826167257103971);
         coeffs.push_back( 0.053937736604997032);
         coeffs.push_back(-0.005661092015920628);
         coeffs.push_back( 0.00032320263464067625);
         coeffs.push_back(-1.0809865109373177e-05);
         coeffs.push_back( 2.1645208248117699e-07);
         coeffs.push_back(-2.5401252316582755e-09);
         coeffs.push_back( 1.6035414536489318e-11);
         coeffs.push_back(-4.1871266927217782e-14);
      }
      coefficients.push_back( coeffs );
   }
}


const float FiniteDifferenceBoundedDivergenceVectorVolume::eval( const Vector& P ) const
{
   int i,j,k, Nx=maxN, Ny=maxN, Nz=maxN;
   if( gB->isInside(P) )
   {
      gB->getGridIndex( P,i,j,k );
      if( i<Nx ) { Nx = i; }
      if( i+Nx>=(gB->nx()) ){ Nx = gB->nx() - i - 1; }
      if( j<Ny ) { Ny = j; }
      if( j+Ny>=(gB->ny()) ){ Ny = gB->ny() - j - 1; }
      if( k<Nz ) { Nz = k; }
      if( k+Nz>=(gB->nz()) ){ Nz = gB->nz() - k - 1; }
   }
   Vector fx;
   Vector fy;
   Vector fz;
   if( Nx > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Nx-1].size();i++)
      {
         ii++;
         Vector DX(ii*gB->dx(),0,0);
         fx += elem->eval(P+DX)*coefficients[Nx-1][i];
         fx -= elem->eval(P-DX)*coefficients[Nx-1][i];
      }
   }
   else
   {
      if( i == 0 )
      {
         fx = elem->eval(P + Vector(gB->dx(),0,0)) - elem->eval(P);
      }
      else
      {
         fx = elem->eval(P) - elem->eval(P - Vector(gB->dx(),0,0));
      }
   }
   if( Ny > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Ny-1].size();i++)
      {
         ii++;
         Vector DX(0,ii*gB->dy(),0);
         fy += elem->eval(P+DX)*coefficients[Ny-1][i];
         fy -= elem->eval(P-DX)*coefficients[Ny-1][i];
      }
   }
   else
   {
      if( j == 0 )
      {
         fy = elem->eval(P + Vector(0,gB->dy(),0)) - elem->eval(P);
      }
      else
      {
         fy = elem->eval(P) - elem->eval(P - Vector(0,gB->dy(),0));
      }
   }
   if( Nz > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Nz-1].size();i++)
      {
         ii++;
         Vector DX(0,0,ii*gB->dz());
         fz += elem->eval(P+DX)*coefficients[Nz-1][i];
         fz -= elem->eval(P-DX)*coefficients[Nz-1][i];
      }
   }
   else
   {
      if( k == 0 )
      {
         fz = elem->eval(P + Vector(0,0,gB->dz())) - elem->eval(P);
      }
      else
      {
         fz = elem->eval(P) - elem->eval(P - Vector(0,0,gB->dz()));
      }
   }

   return (fx.X()/gB->dx() + fy.Y()/gB->dy() + fz.Z()/gB->dz());
}

const Vector FiniteDifferenceBoundedDivergenceVectorVolume::grad( const Vector& P ) const 
{
   int i,j,k, Nx=maxN, Ny=maxN, Nz=maxN;
   if( gB->isInside(P) )
   {
      gB->getGridIndex( P,i,j,k );
      if( i<Nx ) { Nx = i; }
      if( i+Nx>=(gB->nx()) ){ Nx = gB->nx() - i - 1; }
      if( j<Ny ) { Ny = j; }
      if( j+Ny>=(gB->ny()) ){ Ny = gB->ny() - j - 1; }
      if( k<Nz ) { Nz = k; }
      if( k+Nz>=(gB->nz()) ){ Nz = gB->nz() - k - 1; }
   }
   double fx = 0;
   double fy = 0;
   double fz = 0;
   if( Nx > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Nx-1].size();i++)
      {
         ii++;
         Vector DX(ii*gB->dx(),0,0);
         fx += eval(P+DX)*coefficients[Nx-1][i];
         fx -= eval(P-DX)*coefficients[Nx-1][i];
      }
   }
   else
   {
      if( i == 0 )
      {
         fx = eval(P + Vector(gB->dx(),0,0)) - eval(P);
      }
      else
      {
         fx = eval(P) - eval(P - Vector(gB->dx(),0,0));
      }
   }
   if( Ny > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Ny-1].size();i++)
      {
         ii++;
         Vector DX(0,ii*gB->dy(),0);
         fy += eval(P+DX)*coefficients[Ny-1][i];
         fy -= eval(P-DX)*coefficients[Ny-1][i];
      }
   }
   else
   {
      if( j == 0 )
      {
         fy = eval(P + Vector(0,gB->dy(),0)) - eval(P);
      }
      else
      {
         fy = eval(P) - eval(P - Vector(0,gB->dy(),0));
      }
   }
   if( Nz > 0 )
   {
      int ii = 0;
      for( size_t i=0;i<coefficients[Nz-1].size();i++)
      {
         ii++;
         Vector DX(0,0,ii*gB->dz());
         fz += eval(P+DX)*coefficients[Ny-1][i];
         fz -= eval(P-DX)*coefficients[Ny-1][i];
      }
   }
   else
   {
      if( k == 0 )
      {
         fz = eval(P + Vector(0,0,gB->dz())) - eval(P);
      }
      else
      {
         fz = eval(P) - eval(P - Vector(0,0,gB->dz()));
      }
   }
   return Vector(fx/gB->dx(), fy/gB->dy(), fz/gB->dz());
}






































CurlVectorVolume::CurlVectorVolume( Volume<Vector>* v1, const float dx ) :
      elem(v1),
      step(dx)
    {}

CurlVectorVolume::CurlVectorVolume( const VectorField& v1, const float dx ) :
      elem(v1),
      step(dx)
    {}


    const Vector CurlVectorVolume::eval( const Vector& P ) const 
    {
       Matrix g = elem->grad(P);
       return Vector( g(1,2) - g(2,1), g(2,0) - g(0,2), g(0,1)-g(1,0) );
    }

    // evaluate gradient by finite difference
    // size of difference is controlled via the "dx" attribute
    const Matrix CurlVectorVolume::grad( const Vector& P ) const 
    {
       float dx = step;
       Vector e = eval(P);
       Vector ex = eval(P+Vector(dx,0,0));
       Vector ey = eval(P+Vector(0,dx,0));
       Vector ez = eval(P+Vector(0,0,dx));
       Matrix result( ex-e, ey-e, ez-e );
       result /= dx;
       return result.transpose();
    }







ReportVectorVolume::ReportVectorVolume( Volume<Vector>* v, const string tag ) :
      elem(v),
      label (tag)
    {}

ReportVectorVolume::ReportVectorVolume( const VectorField& v, const string tag ) :
      elem(v),
      label (tag)
    {}


    const Vector ReportVectorVolume::eval( const Vector& P ) const 
    {
       Vector value = elem->eval(P);
       cout << label << " eval: " << value[0] << " " << value[1] << " " << value[2] << endl;
       return value;
    }

    const Matrix ReportVectorVolume::grad( const Vector& P ) const { return elem->grad(P); }









AdvectVectorVolume::AdvectVectorVolume( Volume<Vector>* v, Volume<Vector>* u, const float delt ) :
       elem(v),
       velocity(u),
       dt (delt)
    {}

AdvectVectorVolume::AdvectVectorVolume( const VectorField& v, const VectorField& u, const float delt ) :
       elem(v),
       velocity(u),
       dt (delt)
    {}


    const Vector AdvectVectorVolume::eval( const Vector& P ) const 
    { 
       Vector X = P - ( velocity->eval(P) )*dt;
       return elem->eval(X); 
    }
    const Matrix AdvectVectorVolume::grad( const Vector& P ) const
    {
       Vector X = P - ( velocity->eval(P) )*dt;
       Matrix M = unitMatrix() - ( velocity->grad(P) )*dt;
       return M * (elem->grad(X)); 
    }

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

DisplaceVectorVolume::DisplaceVectorVolume( Volume<Vector>* v, Volume<Vector>* u ) :
      elem(v),
      disp(u)
    {}

DisplaceVectorVolume::DisplaceVectorVolume( const VectorField& v, const VectorField& u ) :
      elem(v),
      disp(u)
    {}


    const Vector DisplaceVectorVolume::eval( const Vector& P ) const 
    { 
       Vector X = P - ( disp->eval(P) );
       return elem->eval(X); 
    }
    const Matrix DisplaceVectorVolume::grad( const Vector& P ) const
    {
       Vector X = P - ( disp->eval(P) );
       Matrix M = unitMatrix() - ( disp->grad(P) );
       return M * (elem->grad(X)); 
    }



ContinuedFractionDisplacementVectorVolume::ContinuedFractionDisplacementVectorVolume( Volume<Vector>* u, int _iterations ) :
      elem(u),
      iterations(_iterations)
    {}

ContinuedFractionDisplacementVectorVolume::ContinuedFractionDisplacementVectorVolume( const VectorField& u, int _iterations ) :
      elem(u),
      iterations(_iterations)
    {}


    const Vector ContinuedFractionDisplacementVectorVolume::eval( const Vector& P ) const 
    { 
       Vector X = P;
       for( int i=0;i<iterations;i++)
       {
          X = P - elem->eval(X);
       }
       return X; 
    }

    const Matrix ContinuedFractionDisplacementVectorVolume::grad( const Vector& P ) const
    {
       float dx = 0.0001;
       Vector e = eval(P);
       Vector ex = eval(P+Vector(dx,0,0));
       Vector ey = eval(P+Vector(0,dx,0));
       Vector ez = eval(P+Vector(0,0,dx));
       Matrix result( ex-e, ey-e, ez-e );
       result /= dx;
       return result.transpose();
    }









WarpVectorVolume::WarpVectorVolume( Volume<Vector>* v, Volume<Vector>* u ) :
      elem(v),
      warp(u)
    {}

WarpVectorVolume::WarpVectorVolume( const VectorField& v, const VectorField& u ) :
      elem(v),
      warp(u)
    {}


    const Vector WarpVectorVolume::eval( const Vector& P ) const 
    { 
       Vector X = warp->eval(P);
       return elem->eval(X); 
    }
    const Matrix WarpVectorVolume::grad( const Vector& P ) const
    {
       Vector X = warp->eval(P);
       Matrix M = warp->grad(P);
       return M * (elem->grad(X)); 
    }





NoiseVectorVolume::NoiseVectorVolume( Noise* n, const float d ) : noise (n), dx (d) {}
NoiseVectorVolume::NoiseVectorVolume( NoiseMachine n, const float d ) : noise (n), dx (d) {}


NoiseSampleVectorVolume::NoiseSampleVectorVolume( Noise* n, const float d ) : noise (n), dx (d) {}
NoiseSampleVectorVolume::NoiseSampleVectorVolume( NoiseMachine n, const float d ) : noise (n), dx (d) {}






GriddedVectorVolume::GriddedVectorVolume( const VolumeGrid<Vector>* g ) :
       grid (g),
       sparsegrid(0)
    {
       dx = g->dx();
       dy = g->dy();
       dz = g->dz();
    }

GriddedVectorVolume::GriddedVectorVolume( const SparseVectorGrid* g ) :
       grid(0),
       sparsegrid (g)
    {
       dx = g->dx();
       dy = g->dy();
       dz = g->dz();
    }

GriddedSGridVectorVolume::GriddedSGridVectorVolume( const VectorGrid& g ) :
       elem(g)
    {
       dx = g->dx();
       dy = g->dy();
       dz = g->dz();
    }


GriddedFrustumVectorVolume::GriddedFrustumVectorVolume( const VectorFrustumGrid& g ) :
       felem( g )
    {
       dx = g->dx();
       dy = g->dy();
       dz = g->dz();
    }

/*
GriddedOGridVectorVolume::GriddedOGridVectorVolume( const VectorOGrid& g ) :
       elem(g)
    {
       dx = g->dx();
       dy = g->dy();
       dz = g->dz();
    }
*/


    const Vector GriddedVectorVolume::eval( const Vector& P ) const 
    {
       if( grid ){ return grid->eval(P); }
       if( sparsegrid ){ return sparsegrid->eval(P); }
       return Vector(0,0,0); 
    }

    const Matrix GriddedVectorVolume::grad( const Vector& P ) const
    {
       Matrix value( (eval(P+Vector(dx,0,0))-eval(P-Vector(dx,0,0)))/(2.0*dx),
                     (eval(P+Vector(0,dy,0))-eval(P-Vector(0,dy,0)))/(2.0*dy),
                     (eval(P+Vector(0,0,dz))-eval(P-Vector(0,0,dz)))/(2.0*dz)    );
       return value.transpose();
    }




    const Vector GriddedSGridVectorVolume::eval( const Vector& P ) const 
    {
       return elem->eval(P);
    }

    const Matrix GriddedSGridVectorVolume::grad( const Vector& P ) const
    {
       Matrix value( (eval(P+Vector(dx,0,0))-eval(P-Vector(dx,0,0)))/(2.0*dx),
                     (eval(P+Vector(0,dy,0))-eval(P-Vector(0,dy,0)))/(2.0*dy),
                     (eval(P+Vector(0,0,dz))-eval(P-Vector(0,0,dz)))/(2.0*dz)    );
       return value.transpose();
    }





    const Vector GriddedFrustumVectorVolume::eval( const Vector& P ) const 
    {
       return felem->eval(P);
    }

    const Matrix GriddedFrustumVectorVolume::grad( const Vector& P ) const
    {
       Matrix value( (eval(P+Vector(dx,0,0))-eval(P-Vector(dx,0,0)))/(2.0*dx),
                     (eval(P+Vector(0,dy,0))-eval(P-Vector(0,dy,0)))/(2.0*dy),
                     (eval(P+Vector(0,0,dz))-eval(P-Vector(0,0,dz)))/(2.0*dz)    );
       return value.transpose();
    }



/*
    const Vector GriddedOGridVectorVolume::eval( const Vector& P ) const 
    {
       return elem->eval(P);
    }

    const Matrix GriddedOGridVectorVolume::grad( const Vector& P ) const
    {
       Matrix value( (eval(P+Vector(dx,0,0))-eval(P-Vector(dx,0,0)))/(2.0*dx),
                     (eval(P+Vector(0,dy,0))-eval(P-Vector(0,dy,0)))/(2.0*dy),
                     (eval(P+Vector(0,0,dz))-eval(P-Vector(0,0,dz)))/(2.0*dz)    );
       return value.transpose();
    }

*/





PeriodicVectorVolume::PeriodicVectorVolume( Volume<Vector>* v, const Vector& o, const Vector& L ) :
      elem(v),
      origin (o),
      periods (L)
   {}

PeriodicVectorVolume::PeriodicVectorVolume( const VectorField& v, const Vector& o, const Vector& L ) :
      elem(v),
      origin (o),
      periods (L)
   {}


   const Vector PeriodicVectorVolume::eval( const Vector& P ) const
   {
      return elem->eval( periodP(P) );
   }
   
   const Matrix PeriodicVectorVolume::grad( const Vector& P ) const
   {
      return elem->grad( periodP(P) );
   }
   

   const Vector PeriodicVectorVolume::periodP( const Vector& P ) const
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


DivideVectorVolume::DivideVectorVolume( Volume<Vector> * v, const float a ) : 
       elem(v),
       factor( new ConstantVolume(1.0) ),
       constant (a) 
    {}

DivideVectorVolume::DivideVectorVolume( Volume<Vector> * v, Volume<float>* u ) : 
      elem(v),
      factor(u),
      constant (1.0) 
    {}

DivideVectorVolume::DivideVectorVolume( const VectorField& v, const float a ) : 
       elem(v),
       factor( new ConstantVolume(1.0) ),
       constant (a) 
    {}

DivideVectorVolume::DivideVectorVolume( const VectorField& v, const ScalarField& u ) : 
      elem(v),
      factor(u),
      constant (1.0) 
    {}

const Vector DivideVectorVolume::eval( const Vector& P ) const
{
   return (elem->eval(P))/(constant*(factor->eval(P)));
}



const Matrix DivideVectorVolume::grad( const Vector& P ) const
{
   Matrix eg = elem->grad(P);
   Vector e = elem->eval(P);
   float uv = factor->eval(P);
   Vector ug = factor->grad(P);
   return (eg*(1.0/(constant*uv)) - (ug&e)*(1.0/(constant*uv*uv)));
}


ComponentVectorVolume::ComponentVectorVolume( Volume<float>* x, Volume<float>* y, Volume<float>* z ) :
   X(x), Y(y), Z(z) {}
ComponentVectorVolume::ComponentVectorVolume( const ScalarField& x, const ScalarField& y, const ScalarField& z ) :
   X(x), Y(y), Z(z) {}


const Vector ComponentVectorVolume::eval( const Vector& P ) const
{
   return Vector( X->eval(P), Y->eval(P), Z->eval(P)  );
}

const Matrix ComponentVectorVolume::grad( const Vector& P ) const
{
   Matrix m( X->grad(P), Y->grad(P), Z->grad(P)  );
   return m.transpose();
}


SwitchVectorVolume::SwitchVectorVolume( Volume<Vector>* v1, Volume<Vector>* v2, Volume<float>* swtch ) :
   elem1(v1),
   elem2(v2),
   swtchelem(swtch)
{}

SwitchVectorVolume::SwitchVectorVolume( const VectorField& v1, const VectorField& v2, const ScalarField& swtch ) :
   elem1(v1),
   elem2(v2),
   swtchelem(swtch)
{}

const Vector SwitchVectorVolume::eval( const Vector& P ) const
{
   float swvalue = swtchelem->eval(P);
   if( swvalue <= 0 )
   {
      return elem1->eval(P);
   }
   return elem2->eval(P);
}
   
const Matrix SwitchVectorVolume::grad( const Vector& P ) const
{
   float swvalue = swtchelem->eval(P);
   if( swvalue <= 0 )
   {
      return elem1->grad(P);
   }
   return elem2->grad(P);
}







XVectorVolume::XVectorVolume( Volume<Vector>* v ) :
   elem(v)
   {}

XVectorVolume::XVectorVolume( const VectorField& v ) :
   elem(v)
   {}

const float XVectorVolume::eval( const Vector& P ) const
{
   Vector value = elem->eval(P);
   return value[0];
}

const Vector XVectorVolume::grad( const Vector& P ) const
{
   Matrix value = elem->grad(P);
   return (Vector(1,0,0)*value);
}





YVectorVolume::YVectorVolume( Volume<Vector>* v ) :
   elem(v)
   {}

YVectorVolume::YVectorVolume( const VectorField& v ) :
   elem(v)
   {}

const float YVectorVolume::eval( const Vector& P ) const
{
   Vector value = elem->eval(P);
   return value[1];
}

const Vector YVectorVolume::grad( const Vector& P ) const
{
   Matrix value = elem->grad(P);
   return (Vector(0,1,0)*value);
}





ZVectorVolume::ZVectorVolume( Volume<Vector>* v ) :
   elem(v)
   {}

ZVectorVolume::ZVectorVolume( const VectorField& v ) :
   elem(v)
   {}

const float ZVectorVolume::eval( const Vector& P ) const
{
   Vector value = elem->eval(P);
   return value[2];
}

const Vector ZVectorVolume::grad( const Vector& P ) const
{
   Matrix value = elem->grad(P);
   return (Vector(0,0,1)*value);
}

const Vector GradientStretchCMVolume::eval( const Vector& P ) const
{
    const double dt = totalTime/nbiterations;
    Vector u = elem->eval(P);
    Matrix gu = elem->grad(P);
    Matrix GX = exp(-dt*gu );
    Vector X = P;
    Matrix integral;
    Matrix expintegral = dt*sinch(dt*gu);
    //Matrix expintegral = dt*GX;
    for( int i=0;i<nbiterations;i++ )
    {
       integral +=  expintegral;
       X = P - u*integral;
       gu = elem->grad(X);
       GX = exp( -dt*gu );
       expintegral = expintegral * GX;
    }
    return X;
}

const Matrix GradientStretchCMVolume::grad( const Vector& P ) const 
{
    const double dt = totalTime/nbiterations;
    Vector u = elem->eval(P);
    Vector X = P;
    Matrix gu = elem->grad(P);
    Matrix GX = exp(-dt*gu );
    Matrix integral;
    Matrix expintegral = dt*sinch(dt*gu);
    Matrix gradvalue = GX; 
    //Matrix expintegral = dt*GX;
    for( int i=0;i<nbiterations;i++ )
    {
       integral +=  expintegral;
       //integral += dt * expintegral;
       X = P - u*integral;
       GX = exp( -dt*elem->grad(X) );
       expintegral = expintegral * GX;
       gradvalue = gradvalue * GX;
    }
    return gradvalue;
}





FFTDivFreeNoiseVectorVolume::FFTDivFreeNoiseVectorVolume( const float powerLaw, const float largeScale, const float smallScale, const float length, const int dim, const int seed ) 
{
   dx = length / (dim-1);
}



const Vector FFTDivFreeNoiseVectorVolume::eval( const Vector& P ) const { return data.eval(P); }
const Matrix FFTDivFreeNoiseVectorVolume::grad( const Vector& P ) const 
{
       Matrix value( (eval(P+Vector(dx,0,0))-eval(P-Vector(dx,0,0)))/(2.0*dx),
                     (eval(P+Vector(0,dx,0))-eval(P-Vector(0,dx,0)))/(2.0*dx),
                     (eval(P+Vector(0,0,dx))-eval(P-Vector(0,0,dx)))/(2.0*dx)    );
       return value.transpose();
}







GradDisplaceVectorVolume::GradDisplaceVectorVolume( Volume<Matrix> * m, Volume<float>* ls, const double step ) :
   elem (m),
   sdf  (ls),
   ds   (step)
{}

GradDisplaceVectorVolume::GradDisplaceVectorVolume( const MatrixField& m, const ScalarField& ls, const double step ) :
   elem (m),
   sdf  (ls),
   ds   (step)
{}


const Vector GradDisplaceVectorVolume::eval( const Vector& P ) const
{
   float sdfvalue = sdf->eval(P);
   if( sdfvalue >= 0.0 ){ return Vector(0,0,0); }  // no displacement inside sdf
   Vector sdfDir = sdf->grad(P);
   double sdfDirmag = sdfDir.magnitude();
   if( sdfDirmag > 1.01 || sdfDirmag < 0.99 ){ return Vector(0,0,0); } // allow some slop in the gradient magnitude, but not much
   Vector CPT = P - sdfvalue * (sdfDir);
   double L = (CPT-P).magnitude();
   int nsteps = (int)( L/ ds ) + 1;
   double dstep = L / nsteps;
   Vector Dir = (P-CPT).unitvector();
   Vector Y = CPT;
   for( int ii=0;ii<nsteps;ii++)
   {
      Matrix m = elem->eval(Y);
      Vector prod = Dir*m;
      Y += dstep * prod;
   }
   return Y-P;
}

const Matrix GradDisplaceVectorVolume::grad( const Vector& P ) const
{
   float sdfvalue = sdf->eval(P);
   if( sdfvalue >= 0.0 ){ return Matrix(); }  // no displacement inside sdf
   return elem->eval(P) - unitMatrix(); // subtract unitMatrix becauses this only returns displacement
}

/*
RayMarchDivFreeZeroNormalVectorVolume::RayMarchDivFreeZeroNormalVectorVolume( const VectorField& w, const ScalarField& LS, const float resolution, const VectorField& BC ) :
   W (w),
   sdf(LS),
   ds (resolution),
   boundary_vel (BC)
{
   divW = lux::div(W);
   normal = lux::grad(sdf);
   curvature = lux::div(normal);
   CPT = identity() - normal*sdf;
}
    

RayMarchDivFreeZeroNormalVectorVolume::RayMarchDivFreeZeroNormalVectorVolume( Volume<Vector> * w, Volume<float> * LS, const float resolution, Volume<Vector>* BC) :
   W (w),
   sdf(LS),
   ds (resolution),
   boundary_vel (BC)
{
   divW = lux::div(W);
   normal = lux::grad(sdf);
   curvature = lux::div(normal);
   CPT = identity() - normal*sdf;
}


const Vector RayMarchDivFreeZeroNormalVectorVolume::eval( const Vector& P ) const
{
   float inside = sdf->eval(P);
   // Only solve inside LS
   if( inside < 0.0 ){ return boundary_vel->eval(P); }
   if( inside == 0.0 )
   {
       Vector nhat = normal->eval(P);
       Vector value = W->eval(P);
       value = value - nhat*(nhat*value) + nhat*(nhat*boundary_vel->eval(P));
       return value;
   }
   Vector start = CPT->eval(P);
   Vector bvel = boundary_vel->eval(start);
   double result = 0.0;
   double extinction = 1.0;
   lux::EndPointsRayMarch( start, P, curvature, divW, ds, result, extinction );
   Vector nhat = (P-start);
   double nhatmag = nhat.magnitude();
   if( nhatmag != 0 )
   {
      nhat /= nhatmag;
   }
   else
   {
      nhat = normal->eval(P);
   }
   double F = ( (nhat * W->eval(start))  - nhat*bvel )*extinction + result;
   Vector value = W->eval(P) - nhat*F;
   return value;
}




RayMarchDivFreePlanarVectorVolume::RayMarchDivFreePlanarVectorVolume( const VectorField& vel, const Vector& planeP, const Vector& norm, const float resolution, const VectorField& BC ) :
   W            (vel),
   Pplane       (planeP),
   normal       (norm),
   ds           (resolution),
   boundary_vel (BC)
{
   xhat = Vector(1,0,0);
   yhat = Vector(0,1,0);
   zhat = Vector(0,0,1);

   xperphat = xhat - normal*(xhat*normal);
   yperphat = yhat - normal*(yhat*normal);
   zperphat = zhat - normal*(zhat*normal);
}




RayMarchDivFreePlanarVectorVolume::RayMarchDivFreePlanarVectorVolume( Volume<Vector>* vel, const Vector& planeP, const Vector& norm, const float resolution, Volume<Vector>* BC ):
   W            (vel),
   Pplane       (planeP),
   normal       (norm),
   ds           (resolution),
   boundary_vel (BC)
{
   xhat = Vector(1,0,0);
   yhat = Vector(0,1,0);
   zhat = Vector(0,0,1);

   xperphat = xhat - normal*(xhat*normal);
   yperphat = yhat - normal*(yhat*normal);
   zperphat = zhat - normal*(zhat*normal);
}


const Vector RayMarchDivFreePlanarVectorVolume::eval( const Vector& P ) const
{
   double S = (P-Pplane)*normal;
   if( S < 0.0 )
   {
      return boundary_vel->eval(P);
   }
   Vector start = P - normal*S;
   Vector result = W->eval(P);
   result += normal*(normal*boundary_vel->eval(start));

   //result -= normal*(normal*W->eval(start));
   result -= normal*(normal*result);


   if( S > 0 )
   {
      int nbsteps = (int)(S/ds) + 1;
      double actual_ds = S/nbsteps;
      double accum = 0.0;
      Vector X = start;
      for( int i=0;i<nbsteps;i++ )
      {
         X += normal*actual_ds;
         double divW = xperphat*(W->eval(X+xhat*ds) - W->eval(X-xhat*ds))
                     + yperphat*(W->eval(X+yhat*ds) - W->eval(X-yhat*ds))
                     + zperphat*(W->eval(X+zhat*ds) - W->eval(X-zhat*ds));
         //double divW = xhat*(W->eval(X+xhat*ds) - W->eval(X-xhat*ds))
         //            + yhat*(W->eval(X+yhat*ds) - W->eval(X-yhat*ds))
         //            + zhat*(W->eval(X+zhat*ds) - W->eval(X-zhat*ds));
         accum += divW;
      }
      //if(accum != 0.0 ) { cout << "accum " << accum << "  nbsteps " << nbsteps << endl; }
      accum *= actual_ds / (2.0*ds);
      result -= normal*accum;
   }
   return result;
}

*/
