
#include "ImplicitFormShapes.h"
#include <iostream>
using namespace std;

using namespace lux;




ReportFormVolume::ReportFormVolume( Volume<Form>* v, const string tag ) :
  elem (v),
  label (tag)
{}

ReportFormVolume::ReportFormVolume( const FormField& v, const string tag ) :
  elem (v),
  label (tag)
{}



const Form ReportFormVolume::eval( const Vector& P ) const 
{
   Form value =  elem->eval(P);
   cout << label << " eval: " << value.__str__() << " at position " << P[0] << " " << P[1] << " " << P[2] << endl; 
   return value; 
}

const Form ReportFormVolume::grad( const Vector& P ) const 
{ 
   Form value =  elem->grad(P);
   cout << label << " grad: " << value.__str__()  << " at position " << P[0] << " " << P[1] << " " << P[2]  << endl; 
   return value; 
}




NegateFormVolume::NegateFormVolume( Volume<Form> * v ) :
   elem(v)
   {}

NegateFormVolume::NegateFormVolume( const FormField& v ) :
   elem(v)
   {}

const Form NegateFormVolume::eval( const Vector& P ) const
{
   return -(elem->eval(P));
}

const Form NegateFormVolume::grad(  const Vector& P ) const
{
   return -(elem->grad(P));
}




GradientFormVolume::GradientFormVolume( Volume<Form> * v ) :
   elem(v)
   {}

GradientFormVolume::GradientFormVolume( const FormField& v ) :
   elem(v)
   {}

const Form GradientFormVolume::eval( const Vector& P ) const
{
   return elem->grad(P);
}

const Form GradientFormVolume::grad(  const Vector& P ) const
{
   return Form(0,Vector(0,0,0), Vector(0,0,0),0);
}





ConstantFormVolume::ConstantFormVolume( const Form v ) :
       value (v),
       gradvalue (Form(0,Vector(0,0,0), Vector(0,0,0), 0))
    {}


const Form ConstantFormVolume::eval( const Vector& P ) const 
{
   return value; 
}

const Form ConstantFormVolume::grad( const Vector& P ) const { return gradvalue; }




ScaleFormVolume::ScaleFormVolume( Volume<Form> * v, const Vector& s ) :
       elem (v),
       scale (s)
    {}

ScaleFormVolume::ScaleFormVolume( Volume<Form> * v, const float& s ) :
       elem (v),
       scale (Vector(s,s,s))
    {}


ScaleFormVolume::ScaleFormVolume( const FormField& v, const Vector& s ) :
       elem (v),
       scale (s)
    {}

ScaleFormVolume::ScaleFormVolume( const FormField& v, const float& s ) :
       elem (v),
       scale (Vector(s,s,s))
    {}


const Form ScaleFormVolume::eval( const Vector& P ) const
{
   Vector X = P;
   X[0] /= scale[0];
   X[1] /= scale[1];
   X[2] /= scale[2];
   return  (elem->eval(X));
}

const Form ScaleFormVolume::grad( const Vector& P ) const
{
   Vector X = P;
   X[0] /= scale[0];
   X[1] /= scale[1];
   X[2] /= scale[2];
   Form d =  (elem->grad(X));
   return d;
}





TranslateFormVolume::TranslateFormVolume( Volume<Form> * v, const Vector& s ) :
       elem (v),
       translate (s)
    {}

TranslateFormVolume::TranslateFormVolume( const FormField& v, const Vector& s ) :
       elem (v),
       translate (s)
    {}

const Form TranslateFormVolume::eval( const Vector& P ) const
{
   Vector X = P - translate;
   return  (elem->eval(X));
}

const Form TranslateFormVolume::grad( const Vector& P ) const
{
   Vector X = P - translate;
   return  (elem->grad(X));
}



RotateFormVolume::RotateFormVolume( Volume<Form> * v, const Vector& s ) :
        elem (v),
	axis( s.unitvector()),
	sina (cos(s.magnitude()*M_PI/180.0)),
	cosa (sin(s.magnitude()*M_PI/180.0))
    {
       R = inverse(rotation( s.unitvector() , s.magnitude() ));
    }

RotateFormVolume::RotateFormVolume( const FormField& v, const Vector& s ) :
        elem (v),
	axis( s.unitvector()),
	sina (cos(s.magnitude()*M_PI/180.0)),
	cosa (sin(s.magnitude()*M_PI/180.0))
    {
       R = inverse(rotation( s.unitvector() , s.magnitude() ));
    }

const Form RotateFormVolume::eval( const Vector& P ) const
{
   Vector X = P*cosa + (axis*P)*(1.0-cosa)*axis - (axis^P)*sina;
   return  (elem->eval(X));
}

const Form RotateFormVolume::grad( const Vector& P ) const
{
   Vector X = P*cosa + (axis*P)*(1.0-cosa)*axis - (axis^P)*sina;
   Form v = elem->grad(X);
   return v;
}






MultiplyFormVolume::MultiplyFormVolume( Volume<Form> * v, const float a ) : 
       elem (v),
       factor( new ConstantVolume(1.0) ),
       constant (a) 
    {}

MultiplyFormVolume::MultiplyFormVolume( Volume<Form> * v, Volume<float>* u ) : 
       elem (v),
       factor(u),
       constant (1.0) 
    {}

MultiplyFormVolume::MultiplyFormVolume( const FormField& v, const float a ) : 
       elem (v),
       factor( new ConstantVolume(1.0) ),
       constant (a) 
    {}

MultiplyFormVolume::MultiplyFormVolume( const FormField& v, const ScalarField& u ) : 
       elem (v),
       factor(u),
       constant (1.0) 
    {}

const Form MultiplyFormVolume::eval( const Vector& P ) const
{
   return  elem->eval(P) * factor->eval(P) * constant;
}

const Form MultiplyFormVolume::grad( const Vector& P ) const
{
   return elem->grad(P) * (factor->eval(P) * constant );
}




DivideFormVolume::DivideFormVolume( Volume<Form> * v, const float a ) : 
       elem (v),
       factor( new ConstantVolume(1.0) ),
       constant (a) 
    {}

DivideFormVolume::DivideFormVolume( Volume<Form> * v, Volume<float>* u ) : 
       elem (v),
       factor(u),
       constant (1.0) 
    {}

DivideFormVolume::DivideFormVolume( const FormField& v, const float a ) : 
       elem (v),
       factor( new ConstantVolume(1.0) ),
       constant (a) 
    {}

DivideFormVolume::DivideFormVolume( const FormField& v, const ScalarField& u ) : 
       elem (v),
       factor(u),
       constant (1.0) 
    {}

const Form DivideFormVolume::eval( const Vector& P ) const
{
   return  elem->eval(P) / (factor->eval(P) * constant);
}

const Form DivideFormVolume::grad( const Vector& P ) const
{
   float factv = factor->eval(P);
   return ( elem->grad(P) /factv ) / constant;
}





AddFormVolume::AddFormVolume( Volume<Form> * v1, Volume<Form> * v2 ) :
      e1 (v1),
      e2 (v2)
    {}

AddFormVolume::AddFormVolume( const FormField&  v1, const FormField& v2 ) :
      e1 (v1),
      e2 (v2)
    {}

const Form AddFormVolume::eval( const Vector& P ) const
{
   return  e1->eval(P) + e2->eval(P);
}

const Form AddFormVolume::grad( const Vector& P ) const
{
   return  e1->grad(P) + e2->grad(P);
}



SubtractFormVolume::SubtractFormVolume( Volume<Form> * v1, Volume<Form> * v2 ) :
      e1 (v1),
      e2 (v2)
    {}

SubtractFormVolume::SubtractFormVolume( const FormField& v1, const FormField& v2 ) :
      e1 (v1),
      e2 (v2)
    {}

const Form SubtractFormVolume::eval( const Vector& P ) const
{
   return  e1->eval(P) - e2->eval(P);
}

const Form SubtractFormVolume::grad( const Vector& P ) const
{
   return  e1->grad(P) - e2->grad(P);
}






WedgeFormVolume::WedgeFormVolume( Volume<Form> * v1, Volume<Form> * v2 ) :
       elem1(v1),
       elem2(v2)
    {}

WedgeFormVolume::WedgeFormVolume( const FormField& v1, const FormField& v2 ) :
       elem1(v1),
       elem2(v2)
    {}

const Form WedgeFormVolume::eval( const Vector& P ) const
{
       const Form e0 = elem1->eval(P);
       const Form e1 = elem2->eval(P);
       return ( e0^e1  );
}

const Form WedgeFormVolume::grad( const Vector& P ) const
{
   const Form e0 = elem1->eval(P);
   const Form e1 = elem2->eval(P);
   const Form g0 = elem1->grad(P);
   const Form g1 = elem2->grad(P);
   return ( g0^e1 - e0^g1  );
}

  






AdvectFormVolume::AdvectFormVolume( Volume<Form>* v, Volume<Vector>* u, const float delt ) : 
       elem(v),
       velocity(u),
       dt (delt)
    {}

AdvectFormVolume::AdvectFormVolume( const FormField& v, const VectorField& u, const float delt ) : 
       elem(v),
       velocity(u),
       dt (delt)
    {}


const Form AdvectFormVolume::eval( const Vector& P ) const 
{ 
   Vector X = P - ( velocity->eval(P) )*dt;
   return elem->eval(X); 
}

const Form AdvectFormVolume::grad( const Vector& P ) const
{
   Vector X = P - ( velocity->eval(P) )*dt;
   return  elem->grad(X); 
}





WarpFormVolume::WarpFormVolume( Volume<Form>* v, Volume<Vector>* map ):
      elem(v),
      mapX(map)
    {}

WarpFormVolume::WarpFormVolume( const FormField& v, VectorField& map ):
      elem(v),
      mapX(map)
    {}

const Form WarpFormVolume::eval( const Vector& P ) const 
{ 
   Vector X = mapX->eval(P);
   return elem->eval(X); 
}
const Form WarpFormVolume::grad( const Vector& P ) const
{
   Vector X = mapX->eval(P);
   return (elem->grad(X)); 
}

  


PeriodicFormVolume::PeriodicFormVolume( Volume<Form>* v, const Vector& o, const Vector& L ) :
      origin (o),
      periods (L),
      elem(v)
   {}

PeriodicFormVolume::PeriodicFormVolume( const FormField& v, const Vector& o, const Vector& L ) :
      origin (o),
      periods (L),
      elem(v)
   {}

const Form PeriodicFormVolume::eval( const Vector& P ) const
{
   return elem->eval( periodP(P) );
}
   
const Form PeriodicFormVolume::grad( const Vector& P ) const
{
   return elem->grad( periodP(P) );
}
   
const Vector PeriodicFormVolume::periodP( const Vector& P ) const
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




SwitchFormVolume::SwitchFormVolume( Volume<Form>* v1, Volume<Form>* v2, Volume<float>* swtch ) :
   elem1(v1),
   elem2(v2),
   swtchelem(swtch)
{}

SwitchFormVolume::SwitchFormVolume( const FormField& v1, const FormField& v2, const ScalarField& swtch ) :
   elem1(v1),
   elem2(v2),
   swtchelem(swtch)
{}

const Form SwitchFormVolume::eval( const Vector& P ) const
{
   float swvalue = swtchelem->eval(P);
   if( swvalue <= 0 )
   {
      return elem1->eval(P);
   }
   return elem2->eval(P);
}
   
const Form SwitchFormVolume::grad( const Vector& P ) const
{
   float swvalue = swtchelem->eval(P);
   if( swvalue <= 0 )
   {
      return elem1->grad(P);
   }
   return elem2->grad(P);
}

   




ZeroFormVolume::ZeroFormVolume( Volume<Form>* v1 ) :
   elem1(v1)
   {}

ZeroFormVolume::ZeroFormVolume( const FormField& v1 ) :
   elem1(v1)
   {}


const float ZeroFormVolume::eval( const Vector& P ) const
{
   return elem1->eval(P).zero();
}
   
const Vector ZeroFormVolume::grad( const Vector& P ) const
{
   return elem1->grad(P).one();
}
   



ThreeFormVolume::ThreeFormVolume( Volume<Form>* v1 ) :
   elem1(v1)
   {}

ThreeFormVolume::ThreeFormVolume( const FormField& v1 ) :
   elem1(v1)
   {}


const float ThreeFormVolume::eval( const Vector& P ) const
{
   return elem1->eval(P).three();
}
   
const Vector ThreeFormVolume::grad( const Vector& P ) const
{
   return Vector(0,0,0);
}



OneFormVolume::OneFormVolume( Volume<Form>* v1 ) :
   elem1(v1)
   {}

OneFormVolume::OneFormVolume( const FormField& v1 ) :
   elem1(v1)
   {}


const Vector OneFormVolume::eval( const Vector& P ) const
{
   return elem1->eval(P).one();
}
   
const Matrix OneFormVolume::grad( const Vector& P ) const
{
   return unitMatrix();
}
   



TwoFormVolume::TwoFormVolume( Volume<Form>* v1 ) :
   elem1(v1)
   {}

TwoFormVolume::TwoFormVolume( const FormField& v1 ) :
   elem1(v1)
   {}


const Vector TwoFormVolume::eval( const Vector& P ) const
{
   return elem1->eval(P).two();
}
   
const Matrix TwoFormVolume::grad( const Vector& P ) const
{
   return unitMatrix();
} 





ComponentFormVolume::ComponentFormVolume( Volume<float>* v0, Volume<Vector>* v1, Volume<Vector>* v2, Volume<float>* v3 ) :
   elem0(v0),
   elem1(v1),
   elem2(v2),
   elem3(v3),
   cur(new CurlVectorVolume(v1)),
   dive(new DivergenceVectorVolume(v2))
   {}

ComponentFormVolume::ComponentFormVolume( const ScalarField& v0, const VectorField& v1, const VectorField& v2, const ScalarField& v3 ) :
   elem0(v0),
   elem1(v1),
   elem2(v2),
   elem3(v3),
   cur(new CurlVectorVolume(v1)),
   dive(new DivergenceVectorVolume(v2))
   {}


const Form ComponentFormVolume::eval( const Vector& P ) const
{
   return Form( elem0->eval(P), elem1->eval(P), elem2->eval(P), elem3->eval(P) );
}
   
const Form ComponentFormVolume::grad( const Vector& P ) const
{
   return Form( 0.0, elem0->grad(P), cur->eval(P), dive->eval(P) );
}
   

StarFormVolume::StarFormVolume( Volume<Form>* v1 ) :
   elem(v1)
   {}

StarFormVolume::StarFormVolume( const FormField& v1 ) :
   elem(v1)
   {}

const Form StarFormVolume::eval( const Vector& P ) const { return star(elem->eval(P)); }

// This one is incorrect and is just a placeholder
const Form StarFormVolume::grad( const Vector& P ) const { return star(elem->grad(P)); }






ContractionFormVolume::ContractionFormVolume( Volume<Vector>* v, Volume<Form>* f ) :
   X(v),
   elem(f)
   {
      ScalarField zero = X*VectorField(new OneFormVolume(elem));
      VectorField one = VectorField(new TwoFormVolume(elem))^X;
      VectorField two = X;
      ScalarField three = ScalarField(new ConstantVolume(0.0));
      contraxion = FormField(new ComponentFormVolume(zero,one,two,three));
   }

ContractionFormVolume::ContractionFormVolume( const VectorField& v, const FormField& f ) :
   X(v),
   elem(f)
   {
      ScalarField zero = X*VectorField(new OneFormVolume(elem));
      VectorField one = VectorField(new TwoFormVolume(elem))^X;
      VectorField two = X;
      ScalarField three = ScalarField(new ConstantVolume(0.0));
      contraxion = FormField(new ComponentFormVolume(zero,one,two,three));
   }

const Form ContractionFormVolume::eval( const Vector& P ) const 
{
   return contraxion->eval(P);
}
   
const Form ContractionFormVolume::grad( const Vector& P ) const
{
   return contraxion->grad(P);
}
   

