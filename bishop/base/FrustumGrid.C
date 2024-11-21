
#include "FrustumGrid.h"
using namespace lux;

#include <iostream>
using namespace std;


FrustumBox::FrustumBox( FrustumGrid* f ) : FrustumBoxBase(f) {}

FrustumBox::~FrustumBox(){}

FrustumGrid::~FrustumGrid(){}

const Vector FrustumGrid::evalP( int i, int j, int k ) const
{
   double x = (double)i/(double)nX;
   double y = (double)j/(double)nY;
   double z = (double)k/(double)nZ;
   Vector P = view(x,y);
   P = eye() + P * (z * (farPlane()-nearPlane()) + nearPlane());
   return P;
}


const bool FrustumGrid::isInside( const Vector& P ) const
{
   Vector X = transform(P);
   float test = X[0];
   if( test < 0 || test > 1.0 ){ return false; }
   test = X[1];
   if( test < 0 || test > 1.0 ){ return false; }
   test = X[2];
   if( test < 0 || test > 1.0 ){ return false; }
   return true;
}




int FrustumGrid::index( int i, int j, int k ) const { return ( i + nX*( j + nY*k) ); }
int FrustumGrid::index( const Vector& P ) const
{
   Vector X = transform(P);
   int ii = (int) (  X[0]*nX );
   int jj = (int) (  X[1]*nY );
   int kk = (int) (  X[2]*nZ );
   return index( ii, jj, kk );
}

const float FrustumGrid::fade( const float t ) const
{
   return  ( t * t * t * ( t * ( t * 6.0 - 15.0 ) + 10.0 ) ); 
}



void FrustumGrid::linearInterpolation( float x, int nx, int& i, int& ii, float& w, float& ww ) const
{
    float r = x*nx;
    i  =  r;
    if( i >= 0 && i < nx )
    {
       ii = i + 1;
       w = r-i;
       ww = 1.0 - w;
       if(  ii >= nx )
       {
	  ii = nx-1;
       }
    }
    else
    {
       i = ii = 0;
       w = ww = 0;
    }
}

void FrustumGrid::getLinearInterpolation( const Vector& P,  
                                 int& ix, int& iix, float& wx, float& wwx, 
				 int& iy, int& iiy,  float& wy, float& wwy,
				 int& iz, int& iiz,  float& wz, float& wwz ) const
{
   Vector X = transform( P );
   linearInterpolation( X[0], nX, ix, iix, wx, wwx );
   linearInterpolation( X[1], nY, iy, iiy, wy, wwy );
   linearInterpolation( X[2], nZ, iz, iiz, wz, wwz );

   //wx = fade(wx);
   //wy = fade(wy);
   //wz = fade(wz);
   //wwx = 1.0-wx;
   //wwy = 1.0-wy;
   //wwz = 1.0-wz;
}

const bool FrustumGrid::getGridIndex( const Vector& P, int& ix, int& iy, int& iz ) const
{
   const Vector lower = transform( P );
   ix =  (int)(lower[0]*nX);
   iy =  (int)(lower[1]*nY);
   iz =  (int)(lower[2]*nZ);
   if( ix < 0 ){ return false; }
   if( iy < 0 ){ return false; }
   if( iz < 0 ){ return false; }
   if( ix >= nX ){ return false; }
   if( iy >= nY ){ return false; }
   if( iz >= nZ ){ return false; }
   return true;
}



const bool FrustumGrid::getBox( const Vector& ll, const Vector& uu, int& ixl, int& ixu, int& iyl, int& iyu, int& izl, int& izu ) const
{
   const Vector lower = transform( ll );
   const Vector upper = transform( uu );
   ixl =  (int)(lower[0]*nX - 0.0);
   ixu =  (int)(upper[0]*nX + 1.0);
   iyl =  (int)(lower[1]*nY - 0.0);
   iyu =  (int)(upper[1]*nY + 1.0);
   izl =  (int)(lower[2]*nZ - 0.0);
   izu =  (int)(upper[2]*nZ + 1.0);

   if( ixu < 0 ){ return false; }
   if( ixl >= nX ){ return false; }
   if( ixl < 0 ){ ixl = 0; }
   if( ixu >= nX ){ ixu = nX-1; }

   if( iyu < 0 ){ return false; }
   if( iyl >= nY ){ return false; }
   if( iyl < 0 ){ iyl = 0; }
   if( iyu >= nY ){ iyu = nY-1; }

   if( izu < 0 ){ return false; }
   if( izl >= nZ ){ return false; }
   if( izl < 0 ){ izl = 0; }
   if( izu >= nZ ){ izu = nZ-1; }
   return true;
}



const Camera FrustumGrid::camera() const
{
   Camera cam;
   cam.setEyeViewUp( eye(), view(), up() );
   cam.setNearPlane( nearPlane() );
   cam.setFarPlane( farPlane() );
   cam.setFov( fov() );
   cam.setAspectRatio( aspectRatio() );
   return cam;
}


FrustumBox lux::makeFrustumBox( const int nx, const int ny, const int nz, const Camera& cam )
{
   FrustumBox fb( new FrustumGrid() );
   fb->setEyeViewUp( cam.eye(), cam.view(), cam.up() );
   fb->setFov( cam.fov() );
   fb->setAspectRatio( cam.aspectRatio() );
   fb->setNearPlane( cam.nearPlane() );
   fb->setFarPlane( cam.farPlane() );
   fb->init( nx, ny, nz );
   return fb;
}


const Vector FrustumGrid::llc() const
{
   Vector p = evalP(0,0,0);
   double xmin = p[0];
   double ymin = p[1];
   double zmin = p[2];
   p = evalP(0,0,nZ);
   xmin = ( xmin <= p[0] ) ? xmin : p[0];
   ymin = ( ymin <= p[1] ) ? ymin : p[1];
   zmin = ( zmin <= p[2] ) ? zmin : p[2];
   p = evalP(0,nY,0);
   xmin = ( xmin <= p[0] ) ? xmin : p[0];
   ymin = ( ymin <= p[1] ) ? ymin : p[1];
   zmin = ( zmin <= p[2] ) ? zmin : p[2];
   p = evalP(nX,0,0);
   xmin = ( xmin <= p[0] ) ? xmin : p[0];
   ymin = ( ymin <= p[1] ) ? ymin : p[1];
   zmin = ( zmin <= p[2] ) ? zmin : p[2];
   p = evalP(nX,0,nZ);
   xmin = ( xmin <= p[0] ) ? xmin : p[0];
   ymin = ( ymin <= p[1] ) ? ymin : p[1];
   zmin = ( zmin <= p[2] ) ? zmin : p[2];
   p = evalP(nX,nY,0);
   xmin = ( xmin <= p[0] ) ? xmin : p[0];
   ymin = ( ymin <= p[1] ) ? ymin : p[1];
   zmin = ( zmin <= p[2] ) ? zmin : p[2];
   p = evalP(0,nY,nZ);
   xmin = ( xmin <= p[0] ) ? xmin : p[0];
   ymin = ( ymin <= p[1] ) ? ymin : p[1];
   zmin = ( zmin <= p[2] ) ? zmin : p[2];
   p = evalP(nX,nY,nZ);
   xmin = ( xmin <= p[0] ) ? xmin : p[0];
   ymin = ( ymin <= p[1] ) ? ymin : p[1];
   zmin = ( zmin <= p[2] ) ? zmin : p[2];
   return Vector( xmin,ymin,zmin );
}
const Vector FrustumGrid::urc() const
{
   Vector p = evalP(0,0,0);
   double xmax = p[0];
   double ymax = p[1];
   double zmax = p[2];
   p = evalP(0,0,nZ);
   xmax = ( xmax >= p[0] ) ? xmax : p[0];
   ymax = ( ymax >= p[1] ) ? ymax : p[1];
   zmax = ( zmax >= p[2] ) ? zmax : p[2];
   p = evalP(0,nY,0);
   xmax = ( xmax >= p[0] ) ? xmax : p[0];
   ymax = ( ymax >= p[1] ) ? ymax : p[1];
   zmax = ( zmax >= p[2] ) ? zmax : p[2];
   p = evalP(nX,0,0);
   xmax = ( xmax >= p[0] ) ? xmax : p[0];
   ymax = ( ymax >= p[1] ) ? ymax : p[1];
   zmax = ( zmax >= p[2] ) ? zmax : p[2];
   p = evalP(nX,0,nZ);
   xmax = ( xmax >= p[0] ) ? xmax : p[0];
   ymax = ( ymax >= p[1] ) ? ymax : p[1];
   zmax = ( zmax >= p[2] ) ? zmax : p[2];
   p = evalP(nX,nY,0);
   xmax = ( xmax >= p[0] ) ? xmax : p[0];
   ymax = ( ymax >= p[1] ) ? ymax : p[1];
   zmax = ( zmax >= p[2] ) ? zmax : p[2];
   p = evalP(0,nY,nZ);
   xmax = ( xmax >= p[0] ) ? xmax : p[0];
   ymax = ( ymax >= p[1] ) ? ymax : p[1];
   zmax = ( zmax >= p[2] ) ? zmax : p[2];
   p = evalP(nX,nY,nZ);
   xmax = ( xmax >= p[0] ) ? xmax : p[0];
   ymax = ( ymax >= p[1] ) ? ymax : p[1];
   zmax = ( zmax >= p[2] ) ? zmax : p[2];
   return Vector( xmax,ymax,zmax );
}


const int lux::nx( const FrustumBox& gb ){ return gb->nx(); }
const int lux::ny( const FrustumBox& gb ){ return gb->ny(); }
const int lux::nz( const FrustumBox& gb ){ return gb->nz(); }
const Camera lux::camera( const FrustumBox& gb ){ return gb->camera(); }
