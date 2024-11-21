

#include <iostream>
#include "GeometryVolumeShapes.h"


using namespace std;
using namespace lux;

TriangleLevelSet::TriangleLevelSet( const Vector& p0, const Vector& p1, const Vector& p2, const float bb ) :
   //P0     (p0),
   //edge10 (p1-p0),
   //edge20 (p2-p0),
   //normal ( ((p1-p0)^(p2-p0)).unitvector() ),
   T (Triangle(p0,p1,p2)),
   useBB  (false),
   reject (false)
   {
     /* 
      double denom = (edge10*edge10) * ( edge20*edge20) - pow( (edge10*edge20), 2 );
      if( denom == 0 )
      { 
	 reject = true; 
      }
      else
      {
      U = ( edge10 *(edge20*edge20) - edge20 * (edge20*edge10) ) / denom;
      V = ( edge20 *(edge10*edge10) - edge10 * (edge20*edge10) ) / denom;

      llc = P0;
      urc = P0;

      for( int i=0;i<3;i++ )
      {
         llc[i] = ( llc[i] < p1[i] ) ? llc[i] : p1[i];
         llc[i] = ( llc[i] < p2[i] ) ? llc[i] : p2[i];
         urc[i] = ( urc[i] > p1[i] ) ? urc[i] : p1[i];
         urc[i] = ( urc[i] > p2[i] ) ? urc[i] : p2[i];
      }

      Vector center = (llc + urc)*0.5;
      Vector lengths = urc-llc;
      lengths *= (1.0+fabs(bb))*0.5;
      llc = center - lengths;
      urc = center + lengths;

      useBB = true;
      }
      */
      reject = !T.isGood();
   }


TriangleLevelSet::TriangleLevelSet( const Vector& p0, const Vector& p1, const Vector& p2 ) :
   //P0     (p0),
   //edge10 (p1-p0),
   //edge20 (p2-p0),
   //normal ( ((p1-p0)^(p2-p0)).unitvector() ),
   T (Triangle(p0,p1,p2)),
   useBB  (false),
   reject (false)
   {
     /*
      double denom = (edge10*edge10) * ( edge20*edge20) - pow( (edge10*edge20), 2 );
      if( denom == 0 )
      {
         reject = true;
      }
      else
      {
      U = ( edge10 *(edge20*edge20) - edge20 * (edge20*edge10) ) / denom;
      V = ( edge20 *(edge10*edge10) - edge10 * (edge20*edge10) ) / denom;

      llc = P0;
      urc = P0;

      for( int i=0;i<3;i++ )
      {
         llc[i] = ( llc[i] < p1[i] ) ? llc[i] : p1[i];
         llc[i] = ( llc[i] < p2[i] ) ? llc[i] : p2[i];
         urc[i] = ( urc[i] > p1[i] ) ? urc[i] : p1[i];
         urc[i] = ( urc[i] > p2[i] ) ? urc[i] : p2[i];
     }
     }
     */
     reject = !T.isGood();
 
   }


const float TriangleLevelSet::eval( const Vector& P ) const
{
   return ClosestDistance( T, P );
   /*
   // figure out if it is inside the triangle

   float closest = 1.0e10;
   if( reject ){ return closest; }

   Vector dP = P - P0;
   double side = (dP*normal);
   dP -= side*normal;
   
   double u = dP * U;
   double v = dP * V;
   if( u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0 && (u+v) <= 1.0 )
   {
         closest = fabs(side);
   }

   
   else
   {
      // check for closest edge point

      float t = (dP * edge10 ) / (edge10*edge10);
      if( t >= 0 && t <= 1.0 )
      {
         float d = ( dP-t*edge10).magnitude();
	 closest = ( d < closest ) ? d : closest;
      }
      else
      {
         // check endpoints
         float d = ( dP-edge10).magnitude();
	 closest = ( d < closest ) ? d : closest;
         d = dP.magnitude();
	 closest = ( d < closest ) ? d : closest;
      }

      dP = P - P1;
      t = (dP * edge21 ) / (edge21*edge21);
      if( t >= 0 && t <= 1.0 )
      {
         float d = ( dP-t*edge21).magnitude();
	 closest = ( d < closest ) ? d : closest;
      }
      else
      {
         // check endpoints
         float d = ( dP-edge21).magnitude();
	 closest = ( d < closest ) ? d : closest;
         d = dP.magnitude();
	 closest = ( d < closest ) ? d : closest;
      }

      dP =  P - P2;
      t = ( dP * edge02 ) / (edge02*edge02);
      if( t >= 0 && t <= 1.0 )
      {
         float d = ( dP-t*edge02).magnitude();
	 closest = ( d < closest ) ? d : closest;
      }
      else
      {
         // check endpoints
         float d = ( dP-edge02).magnitude();
	 closest = ( d < closest ) ? d : closest;
         d = dP.magnitude();
	 closest = ( d < closest ) ? d : closest;
      }
   }

   closest = ( side >= 0 ) ? -fabs(closest) : fabs(closest);
   return closest;
   */


}




const Vector TriangleLevelSet::grad( const Vector& P ) const
{

   // finite difference version
   float dx = T.Eu().magnitude();
   float dy = T.Ev().magnitude();
   if( dy<dx ) { dx = dy; }
   dx *= 0.1;
   float v0 = eval(P);
   float vx = eval(P+Vector(dx,0,0));
   float vy = eval(P+Vector(0,dx,0));
   float vz = eval(P+Vector(0,0,dx));

   Vector value = Vector( vx-v0, vy-v0, vz-v0 )/dx;
   return value;


   /*
   // figure out if it is inside the triangle

   float closest = 1.0e10;
   if( reject ){ return Vector(0,0,0); }


   Vector dP = P - P0;
   float side = (dP*normal);

   if( useBB )
   {
      // verify P is inside the BB
      if( P[0] < llc[0] ||  P[1] < llc[1] ||  P[2] < llc[2] )
      {
         return Vector(0,0,0);
      }

      if( P[0] > urc[0] ||  P[1] > urc[1] ||  P[2] > urc[2] )
      {
         return Vector(0,0,0);
      }
   }

   Vector g = normal;
   dP -= side*normal;
   double u = dP * U;
   if( u >= 0 && u <= 1 )
   {
      double v = dP * V;
      if( v >= 0 && v <= 1 && (u+v) <= 1 )
      {
          return normal * ( side > 0 ? 1.0 : -1.0 );
      }
   }

   else
   {
      // check for closest edge point

      float t = (( P - P0 ) * edge10 ) / (edge10*edge10);
      if( t >= 0 && t <= 1.0 )
      {
         float d = ( P-P0-t*edge10).magnitude();
	 if( d < closest ) { g =  (P-P0-t*edge10)/d; }
	 closest = ( d < closest ) ? d : closest;
      }
      else
      {
         // check endpoints
         float d = ( P-P0-edge10).magnitude();
	 closest = ( d < closest ) ? d : closest;
         d = (P-P0).magnitude();
	 if( d < closest ) { g =  (P-P0)/d; }
	 closest = ( d < closest ) ? d : closest;
      }


      t = (( P - P1 ) * edge21 ) / (edge21*edge21);
      if( t >= 0 && t <= 1.0 )
      {
         float d = ( P-P1-t*edge21).magnitude();
	 if( d < closest ) { g =  (P-P1-t*edge21)/d; }
	 closest = ( d < closest ) ? d : closest;
      }
      else
      {
         // check endpoints
         float d = ( P-P1-edge21).magnitude();
	 if( d < closest ) { g =  (P-P1-edge21)/d; }
	 closest = ( d < closest ) ? d : closest;
         d = (P-P1).magnitude();
	 if( d < closest ) { g =  (P-P1)/d; }
	 closest = ( d < closest ) ? d : closest;
      }


      t = (( P - P2 ) * edge02 ) / (edge02*edge02);
      if( t >= 0 && t <= 1.0 )
      {
         float d = ( P-P2-t*edge02).magnitude();
	 if( d < closest ) { g =  (P-P2-t*edge02)/d; }
	 closest = ( d < closest ) ? d : closest;
      }
      else
      {
         // check endpoints
         float d = ( P-P2-edge02).magnitude();
	 closest = ( d < closest ) ? d : closest;
         d = (P-P2).magnitude();
	 if( d < closest ) { g =  (P-P2)/d; }
	 closest = ( d < closest ) ? d : closest;
      }
   }

   return g;

   */

}



const Vector lux::Intersection( const Triangle& t, const Vector& P, const Vector& D )
{
	Vector Dhat = D.unitvector();
	Vector E1 = t.Eu() - Dhat*(t.Eu()*Dhat); 
	Vector E2 = t.Ev() - Dhat*(t.Ev()*Dhat); 
	Vector Y = t.P() - P;
	Y -= Dhat*(Dhat*Y);

	double E1E1 = E1*E1;
	double E2E2 = E2*E2;
	double E1E2 = E1*E2;
	double denom = E1E1*E2E2 - E1E2*E1E2;

	if( fabs(denom) < 0.00001 ){ return Vector( 0,0, 100000000.0 ); }

	double u = ( (E2*Y)*E1E2 - (E1*Y)*E2E2 ) / denom;
	double v = ( (E1*Y)*E1E2 - (E2*Y)*E1E1 ) / denom;
        double tt = ( t.P() - P - u*t.Eu() - v*t.Ev() ).magnitude()/D.magnitude();
	return Vector(u,v,tt);
}





const double lux::ClosestDistance( const Triangle& t, const Vector& P )
{
   double eueu = t.Eu()*t.Eu();
   double evev = t.Ev()*t.Ev();
   double euev = t.Eu()*t.Ev();
   double denom = eueu*evev - euev*euev;
   Vector R = t.P() - P;
   double tuR = t.Eu()*R;
   double tvR = t.Ev()*R;
   if( fabs(denom) > 0.0 )
   {
      // checking interior
      double u = ( euev*(tvR) - evev*(tuR) )/denom;
      if( u >= 0.0 && u <= 1.0 )
      {
         double v = ( euev*(tuR) - eueu*(tvR) )/denom;
	 if( v >=0.0 && v <= 1.0 && (u+v)<=1.0 )
	 {
	    Vector d = R  + u*t.Eu() + v*t.Ev();
	    double distance = d.magnitude();
	    return distance;
	 }
      }
   }
   // checking sides
   bool foundEdge = false;
   double distance = 1.0e20;
   double tu =  -tuR/eueu;
   if( tu >= 0.0 && tu <= 1.0 )
   {
      Vector X = R + tu*t.Eu();
      double xmag = X.magnitude();
      if( foundEdge )
      {
         distance = ( distance < xmag ) ? distance : xmag;
      }
      else
      {
         distance = xmag;
      }
      foundEdge = true;
   }
   double tv = -tvR/evev;
   if( tv >= 0.0 && tv <= 1.0 )
   {
      Vector X = R + tv*t.Ev();
      double xmag = X.magnitude();
      if( foundEdge )
      {
         distance = ( distance < xmag ) ? distance : xmag;
      }
      else
      {
         distance = xmag;
      }
      foundEdge = true;
   }
   Vector e3 = t.Ev() - t.Eu();
   double t3 = -(R*e3)/(e3*e3);
   if( t3 >= 0.0 && t3 <= 1.0 )
   {
      Vector X = R + t3*e3;
      double xmag = X.magnitude();
      if( foundEdge )
      {
         distance = ( distance < xmag ) ? distance : xmag;
      }
      else
      {
         distance = xmag;
      }
      foundEdge = true;
   }

   if( foundEdge )
   {
      return distance;
   }

   // search vertices
   double d1 = R.magnitude();
   double d2 = (R + t.Eu() ).magnitude();
   double d3 = (R + t.Ev() ).magnitude();
   distance = d1;
   distance = ( distance < d2 ) ? distance : d2;
   distance = ( distance < d3 ) ? distance : d3;

   return distance;

}









const Vector lux::ClosestDistance( const TexturedTriangle& t, const Vector& P )
{
   double eueu = t.Eu()*t.Eu();
   double evev = t.Ev()*t.Ev();
   double euev = t.Eu()*t.Ev();
   double denom = eueu*evev - euev*euev;
   Vector R = t.P() - P;
   if( fabs(denom) > 0.0000001 )
   {
      // checking interior
      double u = ( euev*(t.Ev()*R) - evev*(t.Eu()*R) )/denom;
      if( u >= 0 && u <= 1 )
      {
         double v = ( euev*(t.Eu()*R) - eueu*(t.Ev()*R) )/denom;
	 if( v >=0 && v <= 1 && (u+v)<=1 )
	 {
	    Vector d = R  + u*t.Eu() + v*t.Ev();
	    double distance = d.magnitude();
	    return Vector( u, v, d.magnitude() );
	 }
      }
   }
   // checking sides
   bool foundEdge = false;
   Vector distance(0,0, 1.0e20);
   double tu =  -(R*t.Eu())/eueu;
   if( tu >= 0 && tu <= 1 )
   {
      Vector X = R + tu*t.Eu();
      double xmag = X.magnitude();
      if( foundEdge )
      {
         if( distance[2] >= xmag )
	 {
	    distance = Vector( tu, 0.0, xmag );
	 }
      }
      else
      {
         distance = Vector( tu, 0.0, xmag );
      }
      foundEdge = true;
   }
   double tv = -(R*t.Ev())/evev;
   if( tv >= 0 && tv <= 1 )
   {
      Vector X = R + tv*t.Ev();
      double xmag = X.magnitude();
      if( foundEdge )
      {
         if( distance[2] >= xmag )
	 {
	    distance = Vector( 0, tv, xmag );
	 }
      }
      else
      {
	 distance = Vector( 0, tv, xmag );
      }
      foundEdge = true;
   }
   Vector e3 = t.Ev() - t.Eu();
   double t3 = -(R*e3)/(e3*e3);
   if( t3 >= 0 && t3 <= 1 )
   {
      Vector X = R + t3*e3;
      double tu = -t3*( t.Eu()*e3 )/( e3*e3 );
      double tv =  t3*( t.Ev()*e3 )/( e3*e3 );
      double xmag = X.magnitude();
      if( foundEdge )
      {
         if( distance[2] >= xmag )
	 {
	    distance = Vector( tu, tv, xmag );
	 }
      }
      else
      {
	 distance = Vector( tu, tv, xmag );
      }
      foundEdge = true;
   }

   if( foundEdge )
   {
      return distance;
   }

   // search vertices
   double d1 = R.magnitude();
   double d2 = (R + t.Eu() ).magnitude();
   double d3 = (R + t.Ev() ).magnitude();
   distance = Vector( 0,0,d1 );
   if( distance[2] >= d2 )
   {
      distance = Vector( 1, 0, d2 );
   }
   if( distance[2] >= d3 )
   {
      distance = Vector( 0, 1, d3 );
   }

   return distance;

}

