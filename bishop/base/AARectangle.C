

#include "AARectangle.h"
#include <cmath>

using namespace lux;

AARectangle::AARectangle() :
       _llc (Vector(0,0,0)),
       _urc (Vector(0,0,0))
{
   _length = _urc - _llc;
   _center = (_urc + _llc )/2.0;
}


AARectangle::AARectangle( const Vector& llc, const Vector& urc ) :
       _llc (llc),
       _urc (urc)
{
   _length = _urc - _llc;
}

AARectangle::~AARectangle(){}

const bool AARectangle::isInside( const Vector& P ) const
{
   Vector relativeP = P - _llc;
   if( relativeP[0] > _length[0] || relativeP[0] < 0 ) { return false; }
   if( relativeP[1] > _length[1] || relativeP[1] < 0) { return false; }
   if( relativeP[2] > _length[2] || relativeP[2] < 0) { return false; }
   return true;
}

const double AARectangle::signedDistance( const Vector& P ) const 
{ 
   double dist = fabs( P[0] - _llc[0] );
   double cand = fabs( P[0] - _urc[0] );
   dist = ( dist > cand ) ? cand : dist;

   cand = fabs( P[1] - _urc[1] );
   dist = ( dist > cand ) ? cand : dist;
   cand = fabs( P[1] - _llc[1] );
   dist = ( dist > cand ) ? cand : dist;

   cand = fabs( P[2] - _urc[2] );
   dist = ( dist > cand ) ? cand : dist;
   cand = fabs( P[2] - _llc[2] );
   dist = ( dist > cand ) ? cand : dist;

   if( isInside(P) ) { dist = -dist; }

   return dist; 
}

const bool AARectangle::intersection( const Vector& P, const Vector& D, double& close, double& far ) const
{
   double divx =  1 / D[0];
   double tmin, tmax, tymin, tymax, tzmin, tzmax; 
   if ( divx >= 0) 
   {
      tmin = (_llc[0] - P[0]) * divx; 
      tmax = (_urc[0] - P[0]) * divx;
   } 
   else 
   {
      tmin = (_urc[0] - P[0]) * divx; 
      tmax = (_llc[0] - P[0]) * divx;
   } 
      
   divx = 1 / D[1];
   if ( divx >= 0) 
   {
      tymin = (_llc[1] - P[1]) * divx; 
      tymax = (_urc[1] - P[1]) * divx;
   } 
   else 
   {
      tymin = (_urc[1] - P[1]) * divx; 
      tymax = (_llc[1] - P[1]) * divx;
   } 

   if( (tmin > tymax) || ( tymin > tmax )  ) { return false; }

   if (tymin > tmin) { tmin = tymin; }
   if (tymax < tmax) { tmax = tymax; }

   divx = 1 / D[2];
   if (divx >= 0) 
   {
      tzmin = (_llc[2] - P[2]) * divx;
      tzmax = (_urc[2] - P[2]) * divx;
   }
   else 
   {
      tzmin = (_urc[2] - P[2]) * divx;
      tzmax = (_llc[2] - P[2]) * divx;
   }

   if ( (tmin > tzmax) || (tzmin > tmax) ) { return false; }

   if (tzmin > tmin) { tmin = tzmin; }
   if (tzmax < tmax) { tmax = tzmax; }
       
   if( tmin < 0 && tmax > 0 ){ tmin = 0; }
   else if ( tmax <= 0 ){ return false; }

   if( std::isnan(tmin) || std::isnan(tmax) ){ return false; }

   close = tmin;
   far = tmax;

   return true;
}


const double AARectangle::farIntersection( const Vector& P, const Vector& D ) const
{
   double near = -1;
   double far = -1;
   if( intersection( P, D, near, far ) )
   {
      return far;
   }
   return -1.0;
}

const double AARectangle::nearIntersection( const Vector& P, const Vector& D ) const
{
   double near = -1;
   double far = -1;
   if( intersection( P, D, near, far ) )
   {
      return near;
   }
   return -1.0;
}


const Vector AARectangle::normal( const Vector& P ) const
{
   Vector localP  = P -  _llc - _length*0.5;
   if ( fabs( localP[0] - _length[0]*0.5 ) <  _length[0]*0.0001 )
   {
      return Vector(1,0,0);
   }
   if ( fabs( localP[0] + _length[0]*0.5 ) <  _length[0]*0.0001 )
   {
      return Vector(-1,0,0);
   }
   if ( fabs( localP[1] - _length[1]*0.5 ) <  _length[1]*0.0001 )
   {
      return Vector(0,1,0);
   }
   if ( fabs( localP[1] + _length[1]*0.5 ) <  _length[1]*0.0001 )
   {
      return Vector(0,-1,0);
   }
   if ( fabs( localP[2] - _length[2]*0.5 ) <  _length[2]*0.0001 )
   {
      return Vector(0,0,1);
   }
   if ( fabs( localP[2] + _length[2]*0.5 ) <  _length[2]*0.0001 )
   {
      return Vector(0,0,-1);
   }
   return Vector(0,0,0);
}


void AARectangle::split( const int component, AARectangle& aabb1, AARectangle& aabb2 ) const
{
   Vector split_plane = _llc;
   split_plane[component] = _center[component];
   aabb1._llc = split_plane;
   aabb1._urc = _urc;
   split_plane = _urc;
   split_plane[component] = _center[component];
   aabb2._llc = _llc;
   aabb2._urc = split_plane;

   aabb1._center = (aabb1._llc + aabb1._urc)/2.0;
   aabb1._length = (aabb1._urc - aabb1._llc);

   aabb2._center = (aabb2._llc + aabb2._urc)/2.0;
   aabb2._length = (aabb2._urc - aabb2._llc);
}

const bool AARectangle::intersects( const AARectangle& aabb ) const
{
   if( _llc[0] > aabb._urc[0] ){ return false; }
   if( _llc[1] > aabb._urc[1] ){ return false; }
   if( _llc[2] > aabb._urc[2] ){ return false; }
   if( aabb._llc[0] > _urc[0] ){ return false; }
   if( aabb._llc[1] > _urc[1] ){ return false; }
   if( aabb._llc[2] > _urc[2] ){ return false; }
   return true;
}

AARectangle AARectangle::Union( const AARectangle& aabb ) const
{
   Vector nllc = _llc;
   Vector nurc = _urc;
   for( int i=0;i<3;i++ )
   {
      if( nllc[i] > aabb._llc[i] ){ nllc[i] = aabb._llc[i]; }
      if( nurc[i] < aabb._urc[i] ){ nurc[i] = aabb._urc[i]; }
   }
   return AARectangle(nllc,nurc);
}


AABB::AABB( const Vector& llc, const Vector& urc ) : AABBBase( new AARectangle(llc, urc) ) {}
AABB::~AABB() {}
AABB::AABB( const AABB& aabb ) : AABBBase( aabb ) {}


AABB lux::makeAABB( const Vector& llc, const Vector& urc )
{
   return  AABB( llc, urc );
}

AABB lux::expandAABB( const AABB& aabb, const double l )
{
   Vector llc = aabb->llc();
   Vector urc = aabb->urc();
   llc[0] -= l;
   llc[1] -= l;
   llc[2] -= l;
   urc[0] += l;
   urc[1] += l;
   urc[2] += l;
   return makeAABB( llc, urc );
}

AABB lux::expandAABB( const AABB& aabb, const Vector& l )
{
   Vector llc = aabb->llc();
   Vector urc = aabb->urc();
   llc -= l;
   urc += l;
   return makeAABB( llc, urc );
}


AABB lux::shrinkAABB( const AABB& aabb, const double l )
{
   Vector llc = aabb->llc();
   Vector urc = aabb->urc();
   llc[0] += l;
   llc[1] += l;
   llc[2] += l;
   urc[0] -= l;
   urc[1] -= l;
   urc[2] -= l;
   return makeAABB( llc, urc );
}

AABB lux::shrinkAABB( const AABB& aabb, const Vector& l )
{
   Vector llc = aabb->llc();
   Vector urc = aabb->urc();
   llc += l;
   urc -= l;
   return makeAABB( llc, urc );
}

AABB lux::translateAABB( const AABB& aabb, const Vector& t )
{
   Vector llc = aabb->llc();
   Vector urc = aabb->urc();
   llc += t;
   urc += t;
   return makeAABB( llc, urc );
}

const bool lux::isInside( const AABB& aabb, const Vector& P )
{
   return aabb->isInside(P);
}




const Vector& lux::getAABBLLC( const AABB& aabb )
{
    return aabb->llc();
}
const Vector& lux::getAABBURC( const AABB& aabb )
{
    return aabb->urc();
}
const Vector& lux::getAABBCenter( const AABB& aabb )
{
    return aabb->center();
}
const Vector& lux::getAABBLength( const AABB& aabb )
{
    return aabb->length();
}
