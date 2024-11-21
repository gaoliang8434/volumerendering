
#ifndef __RECTANGULARGRID_H__
#define __RECTANGULARGRID_H__

#include "Vector.h"
#include <vector>
#include <memory>
#include "HighOrderInterpolator.h"
#include "AARectangle.h"

namespace lux
{


class RectangularGrid
{
  public:

    RectangularGrid() :
       nX   (0),
       nY   (0),
       nZ   (0),
       dX   (0),
       dY   (0),
       dZ   (0),
       lX   (0),
       lY   (0),
       lZ   (0),
       origin   (Vector(0,0,0)),
       isPeriodicX (false),
       isPeriodicY (false),
       isPeriodicZ (false),
       ho_order (1)
    {}

   ~RectangularGrid(){}


    void init( int nx, int ny, int nz,
               double Lx, double Ly, double Lz,
               const Vector& Origin        );

    void init( const Vector&llc, const Vector& urc, const Vector& res );

    void reset( double Lx, double Ly, double Lz, const Vector& origin );

    int index( int i, int j, int k ) const;
    int index( const Vector& P ) const;
    void triple( const int indx, int& i, int& j, int& k ) const;

    const int& nx() const { return nX; }
    const int& ny() const { return nY; }
    const int& nz() const { return nZ; }

    const Vector& llc() const { return origin; }
    const Vector& urc() const { return topleft; }

    const double& dx() const { return dX; }
    const double& dy() const { return dY; }
    const double& dz() const { return dZ; }
    const double& Lx() const { return lX; }
    const double& Ly() const { return lY; }
    const double& Lz() const { return lZ; }

    const Vector evalP( int i, int j, int k ) const;

    const bool isInside( const Vector& P ) const;
    const bool isInside( int i, int j, int k ) const;

    const bool getBox( const Vector& lower, const Vector& upper, int& ixl, int& ixu, int& iyl, int& iyu, int& izl, int& izu ) const;
    const bool getGridIndex( const Vector& P, int& ix, int& iy, int& iz ) const;
    void getLinearInterpolation( const Vector& P,  
                                 int& ix, int& iix, float& wx, float& wwx, 
				 int& iy, int& iiy,  float& wy, float& wwy,
				 int& iz, int& iiz,  float& wz, float& wwz ) const;

    void getHighOrderInterpolation( const Vector& P, 
                                    std::vector<int>& indices_x, std::vector<int>& indices_y, std::vector<int>& indices_z, 
                                    std::vector<double>& weights_x,
                                    std::vector<double>& weights_y,
                                    std::vector<double>& weights_z ) const;

    virtual const Vector transform( const Vector& P ) const { return P-origin; }

    void setPeriodic() { isPeriodicX = isPeriodicY = isPeriodicZ = true; }
    void setPeriodicX() { isPeriodicX = true; }
    void setPeriodicY() { isPeriodicY = true; }
    void setPeriodicZ() { isPeriodicZ = true; }

    bool periodicX() const { return isPeriodicX; }
    bool periodicY() const { return isPeriodicY; }
    bool periodicZ() const { return isPeriodicZ; }

    void setInterpolationOrder( const int o ) { ho_order = o; }
    const int getInterpolationOrder() const { return ho_order; }


  protected:

    int nX, nY, nZ;
    double dX, dY, dZ;
    double lX, lY, lZ;
    Vector origin, topleft;
    bool isPeriodicX, isPeriodicY, isPeriodicZ;

  private:

    HighOrderInterpolator ho_interp;
    int ho_order;

};


RectangularGrid Union( const RectangularGrid& g1, const RectangularGrid& g2 );


typedef std::shared_ptr<RectangularGrid>  GridBoxBase;


class GridBox : public GridBoxBase
{
  public:

    GridBox( RectangularGrid* f );
   ~GridBox();


    GridBox operator+=( const GridBox& e2 ); // Union

    const Vector evalP( int i, int j, int k ) const;
    int index( const Vector& P ) const;
    char * __str__();
};

GridBox makeGridBox( const Vector& llc, const Vector& urc, const Vector& dx );
GridBox makeGridBox( const RectangularGrid& rg );

const double dx( const GridBox& gb );
const double dy( const GridBox& gb );
const double dz( const GridBox& gb );
const int nx( const GridBox& gb );
const int ny( const GridBox& gb );
const int nz( const GridBox& gb );
const double Lx( const GridBox& gb );
const double Ly( const GridBox& gb );
const double Lz( const GridBox& gb );
const Vector llc( const GridBox& gb );
const Vector urc( const GridBox& gb );
void setPeriodic( GridBox& gb );
void setInterpOrder( GridBox& gb, int order );
int getInterpOrder( const GridBox& gb );
AABB getAABB( const GridBox& gb );

class SparseMapRectangularGrid
{
  public:

    SparseMapRectangularGrid() :
       dX   (0),
       dY   (0),
       dZ   (0),
       origin   (Vector(0,0,0))
    {}

   ~SparseMapRectangularGrid(){}


    void init( float dx, float dy, float dz, const Vector& Origin );

    const Vector& llc() const { return origin; }

    const float& dx() const { return dX; }
    const float& dy() const { return dY; }
    const float& dz() const { return dZ; }

    const Vector evalP( int i, int j, int k ) const;

    const bool getBox( const Vector& lower, const Vector& upper, int& ixl, int& ixu, int& iyl, int& iyu, int& izl, int& izu ) const;
    const bool getGridIndex( const Vector& P, int& ix, int& iy, int& iz ) const;
    void getLinearInterpolation( const Vector& P,  
                                 int& ix, int& iix, float& wx, float& wwx, 
				 int& iy, int& iiy,  float& wy, float& wwy,
				 int& iz, int& iiz,  float& wz, float& wwz ) const;

    virtual const Vector transform( const Vector& P ) const { return P-origin; }


  protected:

    float dX, dY, dZ;
    Vector origin;
};




}




#endif


