

#ifndef ____OVDBGRID_H____
#define ____OVDBGRID_H____

#include "OpenVDBTypes.h"
#include "HighOrderInterpolator.h"
#include <openvdb/io/File.h>
#include <vector>
#include "Vector.h"
#include "Color.h"
#include <string>



namespace lux
{

// Setting up logic to be able to determine the data type of the openvdb data 
template <typename U>
struct OVDB_Grid_Style
{
   typedef int openvdb_type;
   typedef int data_type;
   typedef int openvdb_data_type;
};

template<>
struct OVDB_Grid_Style<float>
{
   typedef FloatGrid openvdb_type;
   typedef float data_type;
   typedef float openvdb_data_type;
};

template<>
struct OVDB_Grid_Style<Vector>
{
   typedef Vec3dGrid openvdb_type;
   typedef Vector data_type;
   typedef openvdb::Vec3d openvdb_data_type;
};


template<>
struct OVDB_Grid_Style<Color>
{
   typedef Vec4fGrid openvdb_type;
   typedef Color data_type;
   typedef openvdb::Vec4f openvdb_data_type;
};

float equals( const FloatGrid::ValueType& b ){ return (float)b; }

openvdb::Vec3d equals( const Vector& b )
{
   openvdb::Vec3d a;
   a[0] = b.X();
   a[1] = b.Y();
   a[2] = b.Z();
   return a;
}


openvdb::Vec4f equals( const Color& b )
{
   openvdb::Vec4f a;
   a[0] = b.X();
   a[1] = b.Y();
   a[2] = b.Z();
   a[3] = b.W();
   return a;
}


Vector equals( const openvdb::Vec3d&  b )
{
   return Vector( b[0], b[1], b[2] );
}


Color equals( const openvdb::Vec4f& b )
{
   return Color( b[0], b[1], b[2], b[3] );
}









template< typename T >
class OVDBGrid
{

  public:

    typedef typename T::data_type DataType;
    typedef typename T::openvdb_type OVDBType;
    typedef typename T::openvdb_type::ValueAllIter Iterator;


    OVDBGrid() : ho_order(1)
    {
       defVal = defVal - defVal;
       init( 1.0, 1.0, 1.0, Vector(0,0,0) );
    }

   ~OVDBGrid(){ grid->clear(); }

    void init( double dx, double dy, double dz, Vector& origin )
    {
       grid = T::openvdb_type::create( defVal );
       dX = dx;
       dY = dy;
       dZ = dz;
       Origin = origin;
    }

    const Vector evalP( int i, int j, int k )
    {
       Vector P = Vector( i*dX, j*dY, k*dZ ) + Origin;
       return P;
    }

    const bool getGridIndex( const Vector& P, int& ix, int& iy, int& iz ) const
    {
       Vector X = P-Origin;
       ix = (int)(X.X()/dX);
       iy = (int)(X.Y()/dY);
       iz = (int)(X.Z()/dZ);
       return true;
    }

    void set( int i, int j, int k, const DataType& value )
    {
       OVDBType temp = equals(value);
       typename T::openvdb_type::Accessor accessor = grid->getAccessor();
       Coord ijk(i,j,k);
       accessor.setValue( ijk, temp );
    }

    const DataType get(int i, int j, int k ) const
    {
       typename T::openvdb_type::ConstAccessor accessor = grid->getConstAccessor();
       Coord ijk(i,j,k);
       DataType temp = equals(accessor.getValue( ijk ));
       return temp;
    }

    void setDefVal( const DataType& def )
    {
       defVal = def;
       changeBackground( grid, defVal );
    }

    const DataType& getDefVal() const { return defVal; }

    const DataType eval( const Vector& P ) const
    {
       DataType accum = defVal * 0.0;
       std::vector<int> indices_x, indices_y, indices_z;
       std::vector<double> weights_x, weights_y, weights_z;
       getHighOrderInterpolation( P, indices_x, indices_y, indices_z, weights_x, weights_y, weights_z );
       if( indices_x.empty() || indices_y.empty() || indices_z.empty() ){ return defVal; }
       for( size_t k=0;k<indices_z.size();k++ )
       {
          for( size_t j=0;j<indices_y.size();j++ )
          {
             for( size_t i=0;i<indices_x.size();i++ )
             {
                accum = accum + get( indices_x[i], indices_y[j], indices_z[k] ) * weights_x[i] * weights_y[j] * weights_z[k];
             }
          }
       }
       return accum;
    }

    Iterator getIterator() { return grid->beginValueAll(); } 


    void setInterpolationOrder( const int o ) { ho_order = o; }
    const int getInterpolationOrder() const { return ho_order; }
    void setName( const std::string& gname )
    {
       grid->setName( gname );
    }

    void write( const std::string& fname ) const
    {
       // create a VDB file object.
       openvdb::io::File fileW(fname);
       // add the grid pointer to a container (since we can write multiple grids out).
       openvdb::GridPtrVec grids;
       grids.push_back(grid);
       // write out the container and close the file.
       fileW.write(grids);
       fileW.close();
    }


    void read( const std::string& fname ) 
    {
       openvdb::io::File fileR(fname);
       fileR.open();
 
       // we need to use the generic GridBase type to temporarily track 
       // the grid until we convert it over to its proper type.
       openvdb::GridBase::Ptr base;
 
       // Read the first grid
       openvdb::io::File::NameIterator iter = fileR.beginName(); 
       if( iter != fileR.endName() )
       {
          base = fileR.readGrid(iter.gridName());
       }
 
       // we are done with the file, so close it
       fileR.close();
 
       // convert the baseGrid over to the proper grid.
       grid = openvdb::gridPtrCast<T::openvdb_type>(base);
    }



    void read( const std::string& fname, const std::string& gname )
    {
       openvdb::io::File fileR(fname);
       fileR.open();
       // we need to use the generic GridBase type to temporarily track 
       // the grid until we convert it over to its proper type.
       openvdb::GridBase::Ptr base;
       // Read the first grid
       openvdb::io::File::NameIterator iter = fileR.beginName(); 
       while( iter != fileR.endName() )
       {
          if( iter.gridName() == gname )
          {
             base = fileR.readGrid(iter.gridName());
             // convert the baseGrid over to the proper grid.
             grid = openvdb::gridPtrCast<T::openvdb_type>(base);
          }
       }
       // we are done with the file, so close it
       fileR.close();
    }

    const float dx() const { return dX; }
    const float dy() const { return dY; }
    const float dz() const { return dZ; }
    const Vector origin() const { return Origin; }


  private:


    typename OVDBType::Ptr grid;
    DataType defVal;

    double dX, dY, dZ;
    Vector Origin;

    HighOrderInterpolator ho_interp;
    int ho_order;


    void getHighOrderInterpolation( const Vector& P, 
                                    std::vector<int>& indices_x, std::vector<int>& indices_y, std::vector<int>& indices_z, 
                                    std::vector<double>& weights_x, std::vector<double>&
                                    weights_y, std::vector<double>& weights_z )
                                    const
    {
       Vector X = P - Origin;
       indices_x.clear();
       indices_y.clear();
       indices_z.clear();

       weights_x.clear();
       weights_y.clear();
       weights_z.clear();

       double rx = X[0]/dX;
       int ix = (int)rx;
       double eps_x = rx -(double)ix;
       int order_x = ho_order;
       ho_interp.weights( eps_x, order_x, weights_x );
       for( int i = ix - order_x + 1; i<= ix+order_x; i++ ){ indices_x.push_back(i); }

       double ry = X[1]/dY;
       int iy = (int)ry;
       double eps_y = ry -(double)iy;
       int order_y = ho_order;
       ho_interp.weights( eps_y, order_y, weights_y );
       for( int i = iy - order_y + 1; i<= iy+order_y; i++ ){ indices_y.push_back(i); }

       double rz = X[2]/dZ;
       int iz = (int)rz;
       double eps_z = rz -(double)iz;
       int order_z = ho_order;
       ho_interp.weights( eps_z, order_z, weights_z );
       for( int i = iz - order_z + 1; i<= iz+order_z; i++ ){ indices_z.push_back(i); }
    }
};



 

class ScalarField;
class VectorField;
class ColorField;




typedef std::shared_ptr<OVDBGrid<OVDB_Grid_Style<float> > > ScalarOGridBase;
typedef std::shared_ptr<OVDBGrid<OVDB_Grid_Style<Vector> > > VectorOGridBase;
typedef std::shared_ptr<OVDBGrid<OVDB_Grid_Style<Color> > > ColorOGridBase;


class ScalarOGrid : public ScalarOGridBase
{
  public:

    ScalarOGrid( OVDBGrid<OVDB_Grid_Style<float> >* f );
   ~ScalarOGrid();


    ScalarOGrid operator+=( const ScalarField& e2 );
    ScalarOGrid operator-=( const ScalarField& e2 );
    ScalarOGrid operator-();
    ScalarOGrid operator*=( const ScalarField& e2 );
    ScalarOGrid operator/=( const ScalarField& e2 );
};


class VectorOGrid : public VectorOGridBase
{
  public:

    VectorOGrid( OVDBGrid<OVDB_Grid_Style<Vector> >* f );
   ~VectorOGrid();


    VectorOGrid operator+=( const VectorField& e2 );
    VectorOGrid operator-=( const VectorField& e2 );
    VectorOGrid operator-();
    VectorOGrid operator*=( const ScalarField& e2 );
    VectorOGrid operator/=( const ScalarField& e2 );


};





class ColorOGrid : public ColorOGridBase
{
  public:

    ColorOGrid( OVDBGrid<OVDB_Grid_Style<Color> >* f );
   ~ColorOGrid();


    ColorOGrid operator+=( const ColorField& e2 );
    ColorOGrid operator*=( const ColorField& e2 );
    ColorOGrid operator*=( const ScalarField& e2 );
    ColorOGrid operator/=( const ScalarField& e2 );


};






}

#endif
