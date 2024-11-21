
//#include "FFTNoise.h"
#include "Fields.h"
//#include "ImplicitVolumeShapes.h"
#include <iostream>

using namespace std;
using namespace lux;




int main()
{

   ScalarField a = Sphere(Vector( 0,0,0 ), 1.0 );
   ScalarField b = Sphere( Vector( 0.5,0,0 ), 1.0 );
   ScalarField c = a+b;
   ScalarField d = c*b/a;
   ScalarField e = constant(0.0);
   ScalarField f = a && b;
   ScalarField g = a || b;
   ScalarField h = a ^ b;

   VectorField I = identity();
   VectorField sgrad = grad(a);
   ScalarField isgrad = I*sgrad;

   int count = 0;
   for (float x=-1.1;x<=1.1;x+=0.01)
   {
      Vector P(x,0,0);
      float valuea = a->eval( P );
      float valueb = b->eval( P );
      float valuec = c->eval( P );
      float valued = d->eval( P );
      e = e + a;
      float valuee = e->eval( P );
      ++count;
      float valuef = f->eval( P );
      float valueg = g->eval( P );
      float valueh = h->eval( P );
      float valueisgrad = isgrad->eval( P );
      cout << x << " " << valuea << " " << valueb << " " << valuea+valueb << " " << valuec << " " << valuec*valueb/valuea << " " << valued << " " << count*valuea << " " << valuee << " " << valuef << " " << valueg << " " << valueh << " " << valueisgrad << endl;
      float teste = evaluate( a, P );
   }




}
