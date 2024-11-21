
#ifndef __FRUSTUMGRID_H__
#define __FRUSTUMGRID_H__

#include "Vector.h"
#include "Camera.h"
#include <memory>

namespace lux
{


class FrustumGrid : public Camera
{
  public:

    FrustumGrid() :
       nX   (0),
       nY   (0),
       nZ   (0)
    {}

   ~FrustumGrid();


    void init( const int& nx, const int& ny, const int& nz )
    {
       nX = nx;
       nY = ny;
       nZ = nz;
    }

    int index( int i, int j, int k ) const;
    int index( const Vector& P ) const;

    const int& nx() const { return nX; }
    const int& ny() const { return nY; }
    const int& nz() const { return nZ; }
    const float dx() const { return htanfov/nX; }
    const float dy() const { return vtanfov/nY; }
    const float dz() const { return 1.0/nZ; }
    const Camera camera() const;

    const Vector evalP( int i, int j, int k ) const;

    const bool isInside( const Vector& P ) const;

    const bool getBox( const Vector& lower, const Vector& upper, int& ixl, int& ixu, int& iyl, int& iyu, int& izl, int& izu ) const;
 
    const bool getGridIndex( const Vector& P, int& ix, int& iy, int& iz ) const;
    void getLinearInterpolation( const Vector& P,
                                 int& ix, int& iix, float& wx, float& wwx,
				 int& iy, int& iiy,  float& wy, float& wwy,
				 int& iz, int& iiz,  float& wz, float& wwz ) const;

    virtual const Vector transform( const Vector& P ) const
    {
       double x,y,z;
       XYZ( P, x, y, z );
       return Vector(x,y,z);
    }

    const Vector llc() const;
    const Vector urc() const;

  protected:

    int nX, nY, nZ;

  private:

    const float fade( const float t ) const;
    void linearInterpolation( float x, int nx, int& i, int& ii, float& w, float& ww ) const;
};



typedef std::shared_ptr<FrustumGrid>  FrustumBoxBase;


class FrustumBox : public FrustumBoxBase
{
  public:

    FrustumBox( FrustumGrid* f );
   ~FrustumBox();

};

FrustumBox makeFrustumBox( const int nx, const int ny, const int nz, const Camera& cam );
const int nx( const FrustumBox& gb );
const int ny( const FrustumBox& gb );
const int nz( const FrustumBox& gb );
const Camera camera( const FrustumBox& gb );




}




#endif
