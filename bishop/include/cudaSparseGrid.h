
#ifndef __SPARSEGRID_H__
#define __SPARSEGRID_H__

#include "RectangularGrid.h"
#include "FrustumGrid.h"
#include "ProgressMeter.h"
#include "Color.h"
#include "Matrix.h"
#include "FrustumGrid.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include <map>
#include <memory>

// Forward declarations of DataObject
#include "data_interface/include/DataObjectDefs.h"

using namespace std;


namespace lux
{

class SparseGrid: public RectangularGrid
{
public:
   SparseGrid();
   SparseGrid( int psize );
   void init(int,int,int,float,float,float,const Vector &);
   const float get(int,int,int) const;
   void setDefVal(float);
   const float getDefVal() const;
   void setPartitionSize(int);
   void set(float,int,int,int);
   const long& size() const { return nbBlocks; }
   void blockBounds( int block, int& i0, int& j0, int& k0, int& i1, int& j1, int& k1 ) const;
   const int blockSize() const { return partitionSize; }
   const float eval( const Vector& P ) const;
   const bool goodBlock( long i ) const { return ( data[i] != NULL ); }
   const bool goodBlock( const int i, const int j, const int k ) const { return ( data[index(i,j,k)] != NULL ); }
   const int index(int,int,int) const;
   ~SparseGrid();
private:
   float **data;
   float defVal;
   int partitionSize;
   long nbBlocks;
   int nnx, nny, nnz;
};

void WriteVolumeGrid(SparseGrid& grid, ofstream& out );
void ReadVolumeGrid(SparseGrid& grid, ifstream& out);
void Blur( SparseGrid& grid );



class SparseColorGrid: public RectangularGrid
{
public:
   SparseColorGrid();
   void init(int,int,int,float,float,float,const Vector &);
   const Color& get(int,int,int) const;
   void setDefVal( Color );
   const Color& getDefVal() const;
   void setPartitionSize(int);
   void set(const Color&,int,int,int);
   const Color eval( const Vector& P ) const;
   void normalize( const SparseGrid& g );
   ~SparseColorGrid();
private:
   Color **data;
   Color defVal;
   int partitionSize;
   const int index(int,int,int) const;
};

void WriteVolumeGrid(SparseColorGrid& grid, ofstream& out );
void ReadVolumeGrid(SparseColorGrid& grid, ifstream& out);
void Blur( SparseColorGrid& grid );




class SparseVectorGrid: public RectangularGrid
{
public:
   SparseVectorGrid();
   void init(int,int,int,float,float,float,const Vector &);
   const Vector& get(int,int,int) const;
   void setDefVal( Vector );
   const Vector& getDefVal() const;
   void setPartitionSize(int);
   void set(const Vector&,int,int,int);
   const Vector eval( const Vector& P ) const;
   void normalize( const SparseGrid& g );
   ~SparseVectorGrid();
private:
   Vector **data;
   Vector defVal;
   int partitionSize;
   const int index(int,int,int) const;
};

void WriteVolumeGrid(SparseVectorGrid& grid, ofstream& out );
void ReadVolumeGrid(SparseVectorGrid& grid, ifstream& out);
void Blur( SparseVectorGrid& grid );




/*
  -- Working on new CUDA-capable SGrid implementation;
  -- ~zshore, 2017-01-10
*/
template <typename T>
class SGrid: public RectangularGrid
{
public:
   SGrid( int psize = 4 ) :
      partitionSize (psize),
      map_size(0),
      data_size(0),
      nbOccupied (0),
      nbMapPartitionsAllocated (0),
      nbDataPartitionsAllocated (0),
      mapOfMap(nullptr),
      map(nullptr),
      data(nullptr),
      cu_mapOfMap(nullptr),
      cu_map(nullptr),
      cu_data(nullptr)
   {
       T empty;
       defVal = empty * 0.0;
   }


   void init(int Nx,int Ny,int Nz, float lx,float ly,float lz,const Vector & Origin)
   {
      int psize = partitionSize*partitionSize;
      nnx = Nx/psize;
      nny = Ny/psize;
      nnz = Nz/psize;
      if( nnx*psize < Nx ){ nnx++; }
      if( nny*psize < Ny ){ nny++; }
      if( nnz*psize < Nz ){ nnz++; }
      int NNx = nnx * psize;
      int NNy = nny * psize;
      int NNz = nnz * psize;
      RectangularGrid::init(NNx, NNy, NNz, lx, ly, lz, Origin);
      totalSize = nnx*nny*nnz;

      deallocate();

      mapOfMap = new long[totalSize];
      for(long i = 0; i < totalSize; ++i)
         mapOfMap[i] = -1;

      nbOccupied = 0;
      nbMapPartitionsAllocated = 0;
      nbDataPartitionsAllocated = 0;
   }

   const T& get(int i,int j,int k) const
   {
      if(!isInside(i,j,k)){ return(defVal); }

      int ii = i / (partitionSize * partitionSize);
      int jj = j / (partitionSize * partitionSize);
      int kk = k / (partitionSize * partitionSize);
      int mapIndex = ii + nnx * (jj + nny * kk);
      if(mapOfMap[mapIndex] == -1)
      {
         return(defVal);
      }
      else
      {
         long offset = mapOfMap[mapIndex];
         ii = (i / partitionSize) % partitionSize;
         jj = (j / partitionSize) % partitionSize;
         kk = (k / partitionSize) % partitionSize;
         mapIndex = ii + partitionSize * (jj + partitionSize * kk);
         if(map[mapIndex + offset] == -1)
         {
            return(defVal);
         }
         else
         {
            offset = map[mapIndex + offset];
            ii = i % partitionSize;
            jj = j % partitionSize;
            kk = k % partitionSize;
            mapIndex = ii + partitionSize * (jj + partitionSize * kk);
            return data[mapIndex + offset];
         }
      }
   }

   void setDefVal( const T& def ) { defVal = def; }

   const T& getDefVal() const { return defVal; }

   void set( int i,int j,int k, const T& val )
   {
      if(val != defVal)
      {
         int ii = i / (partitionSize * partitionSize);
         int jj = j / (partitionSize * partitionSize);
         int kk = k / (partitionSize * partitionSize);
         int mapIndex = ii + nnx * (jj + nny * kk);
         if(mapOfMap[mapIndex] == -1)
         {
            // Allocate an array one block larger than the map array
            int psize = partitionSize * partitionSize * partitionSize;
            long* new_map = new long[map_size + psize];

            for(long x = map_size; x < map_size + psize; ++x)
               new_map[x] = -1;

            // Memcpy over old stuff from map
            if( map )
            {
               std::memcpy(new_map, map, map_size * sizeof(long));
               delete [] map;
            }

            map = new_map;
            mapOfMap[mapIndex] = map_size;
            map_size += psize;

            nbMapPartitionsAllocated++;
         }

         long offset = mapOfMap[mapIndex];
         ii = (i / partitionSize) % partitionSize;
         jj = (j / partitionSize) % partitionSize;
         kk = (k / partitionSize) % partitionSize;
         mapIndex = ii + partitionSize * (jj + partitionSize * kk);
         if(map[mapIndex + offset] == -1)
         {
            // Allocate an array one block larger than the data array
            int psize = partitionSize * partitionSize * partitionSize;
            T* new_data = new T[data_size + psize];

            for(long x = data_size; x < data_size + psize; ++x)
               new_data[x] = defVal;

            // Memcpy over old stuff from data
            if( data )
            {
               std::memcpy(new_data, data, data_size * sizeof(T));
               delete [] data;
            }

            data = new_data;
            map[mapIndex + offset] = data_size;
            data_size += psize;

            nbOccupied++;
            nbDataPartitionsAllocated++;
         }

         offset = map[mapIndex + offset];
         ii = i % partitionSize;
         jj = j % partitionSize;
         kk = k % partitionSize;
         mapIndex = ii + partitionSize * (jj + partitionSize * kk);

         data[mapIndex + offset] = val;
      }
   }

   const T eval( const Vector& P ) const
   {
       if( !isInside(P) ){ return defVal; }


       T accum = defVal * 0.0;
       vector<int> indices_x, indices_y, indices_z;
       vector<double> weights_x, weights_y, weights_z;
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

   ~SGrid() {}

   void deallocate()
   {
      if( !map )
         delete [] map;

      if( !data )
         delete [] data;

      map = NULL;
      data = NULL;
      map_size = 0;
      data_size = 0;

      if( !mapOfMap )
         delete [] mapOfMap;
   }

   void clear()
   {
      deallocate();

      mapOfMap = new long[totalSize];
      for(long i = 0; i < totalSize; ++i)
         mapOfMap[i] = -1;

      nbOccupied = 0;
      nbMapPartitionsAllocated = 0;
      nbDataPartitionsAllocated = 0;
   }

   long Size() const { return nbOccupied; }
   long NbPartitions() const { return nnx*nny*nnz; }

   int blockSize() const { return partitionSize; }

   const bool goodBlock( long i ) const 
   {
      int Nx = nnx * partitionSize * partitionSize;
      int Ny = nny * partitionSize * partitionSize;
      int Nz = nnz * partitionSize * partitionSize;
      int ii = i / (Nx * Ny);
      int jj = (i % (Nx * Ny)) / Nz;
      int kk = (i % (Nx * Ny)) % Nz;
      return goodBlock(ii, jj, kk);
   }

   const bool goodBlock( const int i, const int j, const int k ) const
   {
      int ii = i / (partitionSize * partitionSize);
      int jj = j / (partitionSize * partitionSize);
      int kk = k / (partitionSize * partitionSize);
      int mapIndex = ii + nnx * (jj + nny * kk);
      if(mapOfMap[mapIndex] == -1)
      {
         return false;
      }
      else
      {
         int offset = mapOfMap[mapIndex];
         ii = (i % (partitionSize * partitionSize)) / partitionSize;
         jj = (j % (partitionSize * partitionSize)) / partitionSize;
         kk = (k % (partitionSize * partitionSize)) / partitionSize;
         mapIndex = ii + partitionSize * (jj + partitionSize * kk);
         if(map[mapIndex + offset] == -1)
            return false;
         else
            return true;
      }
   }

   void blockBounds( long block, int& i0, int& j0, int& k0, int& i1, int& j1, int& k1 ) const
   {
      if ( block <0 || block >= totalSize || nbOccupied == 0 )
      {
         i0 = j0 = k0 = i1 = j1 = k1 = -1;
         return;
      }

      k0 = block / (nnx*nny);
      long bblock = block - k0*nnx*nny;
      j0 = bblock / nnx;
      i0 = bblock - j0*nnx;

      i0 *= partitionSize;
      j0 *= partitionSize;
      k0 *= partitionSize;

      i1 = i0 + partitionSize;
      j1 = j0 + partitionSize;
      k1 = k0 + partitionSize;

      return;
   }

   const long index(int i,int j,int k) const
   {
      int ii = i/partitionSize;
      int jj = j/partitionSize;
      int kk = k/partitionSize;
      return ii + nnx*( jj + nny*kk );
   }

   void initCUDA() {}

   gilligan::DataObject<T>* cuData() const { return cu_data; }
   gilligan::DataObject<long>* cuMap() const { return cu_map; }
   gilligan::DataObject<long>* cuMapOfMap() const { return cu_mapOfMap; }

private:
   long* mapOfMap;
   long* map;
   T* data;
   T defVal;

   gilligan::DataObject<long>* cu_mapOfMap;
   gilligan::DataObject<long>* cu_map;
   gilligan::DataObject<T>* cu_data;

   int partitionSize;
   int nnx, nny, nnz;
   size_t map_size, data_size;

   long nbOccupied;
   long nbMapPartitionsAllocated;
   long nbDataPartitionsAllocated;

   long totalSize;

};

/*

template <typename T>
class SGrid: public RectangularGrid
{
  public:
   SGrid( int psize = 4 ) :
      partitionSize (psize),
      nbOccupied (0),
      data(NULL)
   {
       T empty;
       defVal = empty * 0.0;
   }


   void init(int Nx,int Ny,int Nz, float lx,float ly,float lz,const Vector & Origin)
   {
      nnx = Nx/partitionSize;
      nny = Ny/partitionSize;
      nnz = Nz/partitionSize;
      if( nnx*partitionSize > Nx ){ nnx++; }
      if( nny*partitionSize > Ny ){ nny++; }
      if( nnz*partitionSize > Nz ){ nnz++; }
      int NNx = nnx * partitionSize;
      int NNy = nny * partitionSize;
      int NNz = nnz * partitionSize;
      RectangularGrid::init(NNx, NNy, NNz, lx, ly, lz, Origin);
      totalSize = nnx*nny*nnz;
      deallocate();
      data = new T* [totalSize];
      for(int i = 0; i < totalSize; i++)
      {
         data[i] = NULL;
      }
      nbOccupied = 0;
   }

   void clear()
   {
      deallocate();
      data = new T* [totalSize];
      for(int i = 0; i < totalSize; i++)
      {
         data[i] = NULL;
      }
      nbOccupied = 0;
   }

   const T& get(int i,int j,int k) const
   {
      if(!isInside(i,j,k)){ return(defVal); }
      int myIndex = sindex(i, j, k);
      if(data[myIndex] == NULL)
      {
         return(defVal);
      }
      else
      {
         int ii = i % partitionSize;
         int jj = j % partitionSize;
         int kk = k % partitionSize;
         int partitionIndex = ii + partitionSize*( jj + partitionSize * kk );
         return data[myIndex][partitionIndex];
      }
   }

   void setDefVal( const T& def ) { defVal = def; }

   const T& getDefVal() const { return defVal; }

   void set( int i,int j,int k, const T& val )
   {
      if(val != defVal)
      {
         int myIndex = index(i, j, k);
         if(data[myIndex] == NULL)
	 {
            data[myIndex] = new T[partitionSize * partitionSize * partitionSize];
            for(int i = 0; i < partitionSize * partitionSize * partitionSize; i++)
	    {
               data[myIndex][i] = defVal;
            }
	    nbOccupied++;
         }
         int ii = i % partitionSize;
         int jj = j % partitionSize;
         int kk = k % partitionSize;
         int partitionIndex = ii + partitionSize*( jj + partitionSize * kk );
         data[myIndex][partitionIndex] = val;
      }
   }

   const T eval( const Vector& P ) const
   {
       if( !isInside(P) ){ return defVal; }


       T accum = defVal * 0.0;
       vector<int> indices_x, indices_y, indices_z;
       vector<double> weights_x, weights_y, weights_z;
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

       int ix, iix;
       float wx, wwx;
       int iy, iiy;
       float wy, wwy;
       int iz, iiz;
       float wz, wwz;

       getLinearInterpolation( P,  
                               ix, iix, wx,  wwx, 
			       iy, iiy, wy,  wwy,
			       iz, iiz, wz,  wwz );

       T accum = get( ix, iy, iz ) * wwx * wwy * wwz;
       accum  += get( iix, iy, iz ) *  wx * wwy * wwz;
       accum  += get( ix, iiy, iz ) * wwx *  wy * wwz;
       accum  += get( iix, iiy, iz ) *  wx *  wy * wwz;
       accum  += get( ix, iy, iiz ) * wwx * wwy *  wz;
       accum  += get( iix, iy, iiz ) *  wx * wwy *  wz;
       accum  += get( ix, iiy, iiz ) * wwx *  wy *  wz;
       accum  += get( iix, iiy, iiz )*  wx *  wy *  wz;

       return accum;
   }

   ~SGrid()
   {
     deallocate();
   }

   void deallocate()
   {
      // Already deallocated.
      if(data == NULL)
      {
        return;
      }

      long totalSize = nnx*nny*nnz;
      for(long i = 0; i < totalSize; i++)
      {
         if(data[i] != NULL){ delete [] data[i]; }
      }
      delete [] data;
      data = NULL;
      nbOccupied = 0;
   }

   long Size() const { return nbOccupied; }
   long NbPartitions() const { return nnx*nny*nnz; }

   int blockSize() const { return partitionSize; }
   const bool goodBlock( const int i, const int j, const int k ) const { return ( data[sindex(i,j,k)] != NULL ); }
   const bool goodBlock( const int i ) const { return ( data[i] != NULL ); }

   void blockBounds( long block, int& i0, int& j0, int& k0, int& i1, int& j1, int& k1 ) const
   {
      if ( block <0 || block >= totalSize || nbOccupied == 0 )
      {
         i0 = j0 = k0 = i1 = j1 = k1 = -1;
         return;
      }

      k0 = block / (nnx*nny);
      long bblock = block - k0*nnx*nny;
      j0 = bblock / nnx;
      i0 = bblock - j0*nnx;

      i0 *= partitionSize;
      j0 *= partitionSize;
      k0 *= partitionSize;

      i1 = i0 + partitionSize;
      j1 = j0 + partitionSize;
      k1 = k0 + partitionSize;

      return;
   }

   const long sindex(int i,int j,int k) const
   {
      int ii = i/partitionSize;
      int jj = j/partitionSize;
      int kk = k/partitionSize;
      return ii + nnx*( jj + nny*kk );
   }

  private:
      T **data;
      T defVal;
      int partitionSize;
      int nnx, nny, nnz;
      long nbOccupied;
      long totalSize;

};

*/


template<typename T>
void WriteVolumeGrid(SGrid<T>& grid, ofstream& out )
{
   int nx = grid.nx();
   int ny = grid.ny();
   int nz = grid.nz();

   Vector llc = grid.llc();

   float x = llc[0];
   float y = llc[1];
   float z = llc[2];

   float Lx = grid.Lx();
   float Ly = grid.Ly();
   float Lz = grid.Lz();

   out.write( (char*)&nx, sizeof(nx) );
   out.write( (char*)&ny, sizeof(ny) );
   out.write( (char*)&nz, sizeof(nz) );

   out.write( (char*)&Lx, sizeof(Lx) );
   out.write( (char*)&Ly, sizeof(Ly) );
   out.write( (char*)&Lz, sizeof(Lz) );
   out.write( (char*)&x, sizeof(x) );
   out.write( (char*)&y, sizeof(y) );
   out.write( (char*)&z, sizeof(z) );

   T temp = grid.getDefValue();
   out.write( (char*)&temp, sizeof(T));

   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
            temp = grid.get(i,j,k);
            if(temp != grid.getDefVal()){
               out.write((char*)&i, sizeof(int));
               out.write((char*)&j, sizeof(int));
               out.write((char*)&k, sizeof(int));
               out.write((char*)&temp, sizeof(T));
            }
         }

      }

   }

}

template<typename T>
void ReadVolumeGrid(SGrid<T>& grid, ifstream& out)
{
   int nx;
   int ny;
   int nz;

   int i, j, k;

   float x;
   float y;
   float z;

   float Lx;
   float Ly;
   float Lz;

   T value;

   out.read( (char*)&nx, sizeof(nx) );
   out.read( (char*)&ny, sizeof(ny) );
   out.read( (char*)&nz, sizeof(nz) );
   out.read( (char*)&Lx, sizeof(Lx) );
   out.read( (char*)&Ly, sizeof(Ly) );
   out.read( (char*)&Lz, sizeof(Lz) );

   out.read( (char*)&x, sizeof(x) );
   out.read( (char*)&y, sizeof(y) );
   out.read( (char*)&z, sizeof(z) );

   out.read((char *)&value, sizeof(T));
   grid.setDefValue( value );
   grid.init( nx, ny, nz, Lx, Ly, Lz, Vector(x,y,z) );

   while(out.read((char *)&i, sizeof(int)))
   {
       out.read((char *)&j, sizeof(int));
       out.read((char *)&k, sizeof(int));
       out.read((char *)&value, sizeof(T));
       grid.set(value, i, j, k);
   }
}

template<typename T>
void Blur( SGrid<T>& g )
{
   SGrid<T> temp;
   temp.setDefVal( g.getDefVal() );
   temp.init( g.nx(), g.ny(), g.nz(), g.Lx(), g.Ly(), g.Lz(), g.llc() );
   ProgressMeter meter( (g.nx())*(g.ny())*(g.nz()) * 2, "Blur" );
   for( int k=0;k<g.nz();k++ )
   {
      for( int j=0;j<g.ny();j++ )
      {
         for( int i=0;i<g.nx();i++ )
         {
	    temp.set( g.get(i,j,k), i,j,k);
	    meter.update();
         }
      }
   }


   for( int k=0;k<g.nz();k++ )
   {
      int kmin = k-1;
      if( kmin < 0 ){ kmin = k; }
      int kmax = k+1;
      if( kmax >= g.nz() ){ kmax = k;}
      for( int j=0;j<g.ny();j++ )
      {
         int jmin = j-1;
         if( jmin < 0 ){ jmin = j; }
         int jmax = j+1;
         if( jmax >= g.ny() ){ jmax = j;}
         for( int i=0;i<g.nx();i++ )
         {
            int imin = i-1;
            if( imin < 0 ){ imin = i; }
            int imax = i+1;
            if( imax >= g.nx() ){ imax = i;}
	    T sum;
	    sum *= 0.0;
	    int nb = 0;
	    for( int kk=kmin;kk<=kmax;kk++ )
	    {
	       for( int jj=jmin;jj<=jmax;jj++ )
	       {
	          for( int ii=imin;ii<=imax;ii++ )
	          {
		     sum += temp.get( ii, jj, kk );
		     ++nb;
	          }
	       }
	    }
	    sum = sum/nb;
	    g.set(sum, i,j,k);
	    meter.update();
         }
      }
   }
}


class ScalarField;
class VectorField;
class ColorField;
class MatrixField;



typedef std::shared_ptr<SGrid<float> > ScalarGridBase;
typedef std::shared_ptr<SGrid<Vector> > VectorGridBase;
typedef std::shared_ptr<SGrid<Color> > ColorGridBase;
typedef std::shared_ptr<SGrid<Matrix> > MatrixGridBase;

class VectorGrid;
class ColorGrid;
class MatrixGrid;

class ScalarGrid : public ScalarGridBase
{
  public:

    ScalarGrid() {};
    ScalarGrid( SGrid<float>* f );
   ~ScalarGrid();


    ScalarGrid operator+=( const ScalarField& e2 );
    ScalarGrid operator-=( const ScalarField& e2 );
    ScalarGrid operator-();
    ScalarGrid operator*=( const ScalarField& e2 );
    ScalarGrid operator/=( const ScalarField& e2 );

};


class VectorGrid : public VectorGridBase
{
  public:

    VectorGrid() {};
    VectorGrid( SGrid<Vector>* f );
   ~VectorGrid();


    VectorGrid operator+=( const VectorField& e2 );
    VectorGrid operator-=( const VectorField& e2 );
    VectorGrid operator-();
    VectorGrid operator*=( const ScalarField& e2 );
    VectorGrid operator/=( const ScalarField& e2 );


};





class ColorGrid : public ColorGridBase
{
  public:

    ColorGrid() {};
    ColorGrid( SGrid<Color>* f );
   ~ColorGrid();


    ColorGrid operator+=( const ColorField& e2 );
    ColorGrid operator-=( const ColorField& e2 );
    ColorGrid operator-();
    ColorGrid operator*=( const ScalarField& e2 );
    ColorGrid operator/=( const ScalarField& e2 );


};



class MatrixGrid : public MatrixGridBase
{
  public:

    MatrixGrid() {};
    MatrixGrid( SGrid<Matrix>* f );
   ~MatrixGrid();


    MatrixGrid operator+=( const MatrixField& e2 );
    MatrixGrid operator-=( const MatrixField& e2 );
    MatrixGrid operator-();
    MatrixGrid operator*=( const ScalarField& e2 );
    MatrixGrid operator/=( const ScalarField& e2 );


};


void initCUDA( ScalarGrid& grid);
void initCUDA( VectorGrid& grid);
void initCUDA( ColorGrid& grid);
void initCUDA( MatrixGrid& grid);
























template <typename T>
class FSGrid: public FrustumGrid
{
  public:
   FSGrid( int psize = 4 ) :
      partitionSize (psize),
      map_size(0),
      data_size(0),
      nbOccupied (0),
      nbMapPartitionsAllocated (0),
      nbDataPartitionsAllocated (0),
      mapOfMap(nullptr),
      map(nullptr),
      data(nullptr),
      cu_mapOfMap(nullptr),
      cu_map(nullptr),
      cu_data(nullptr)
   {
       T empty;
       defVal = empty * 0.0;
   }


   void init(int Nx,int Ny,int Nz, const Camera& cam)
   {
      int psize = partitionSize*partitionSize;
      nnx = Nx/psize;
      nny = Ny/psize;
      nnz = Nz/psize;
      if( nnx*psize < Nx ){ nnx++; }
      if( nny*psize < Ny ){ nny++; }
      if( nnz*psize < Nz ){ nnz++; }
      int NNx = nnx * psize;
      int NNy = nny * psize;
      int NNz = nnz * psize;
      setEyeViewUp( cam.eye(), cam.view(), cam.up() );
      setFov( cam.fov() );
      setAspectRatio( cam.aspectRatio() );
      setNearPlane( cam.nearPlane() );
      setFarPlane( cam.farPlane() );
      FrustumGrid::init(NNx, NNy, NNz);
      totalSize = nnx*nny*nnz;

      deallocate();

      mapOfMap = new long[totalSize];
      for(long i = 0; i < totalSize; ++i)
         mapOfMap[i] = -1;

      nbOccupied = 0;
      nbMapPartitionsAllocated = 0;
      nbDataPartitionsAllocated = 0;
   }

   const T& get(int i,int j,int k) const
   {
      if(!inside(i,j,k)){ return(defVal); }

      int ii = i / (partitionSize * partitionSize);
      int jj = j / (partitionSize * partitionSize);
      int kk = k / (partitionSize * partitionSize);
      int mapIndex = ii + nnx * (jj + nny * kk);
      if(mapOfMap[mapIndex] == -1)
      {
         return(defVal);
      }
      else
      {
         long offset = mapOfMap[mapIndex];
         ii = (i / partitionSize) % partitionSize;
         jj = (j / partitionSize) % partitionSize;
         kk = (k / partitionSize) % partitionSize;
         mapIndex = ii + partitionSize * (jj + partitionSize * kk);
         if(map[mapIndex + offset] == -1)
         {
            return(defVal);
         }
         else
         {
            offset = map[mapIndex + offset];
            ii = i % partitionSize;
            jj = j % partitionSize;
            kk = k % partitionSize;
            mapIndex = ii + partitionSize * (jj + partitionSize * kk);
            return data[mapIndex + offset];
         }
      }
   }

   void setDefVal( const T& def ) { defVal = def; }

   const T& getDefVal() const { return defVal; }

   void set( int i,int j,int k, const T& val )
   {
      if(val != defVal)
      {
         int ii = i / (partitionSize * partitionSize);
         int jj = j / (partitionSize * partitionSize);
         int kk = k / (partitionSize * partitionSize);
         int mapIndex = ii + nnx * (jj + nny * kk);
         if(mapOfMap[mapIndex] == -1)
         {
            // Allocate an array one block larger than the map array
            int psize = partitionSize * partitionSize * partitionSize;
            long* new_map = new long[map_size + psize];

            for(long x = map_size; x < map_size + psize; ++x)
               new_map[x] = -1;

            // Memcpy over old stuff from map
            if( map )
            {
               std::memcpy(new_map, map, map_size * sizeof(long));
               delete [] map;
            }

            map = new_map;
            mapOfMap[mapIndex] = map_size;
            map_size += psize;

            nbMapPartitionsAllocated++;
         }

         long offset = mapOfMap[mapIndex];
         ii = (i / partitionSize) % partitionSize;
         jj = (j / partitionSize) % partitionSize;
         kk = (k / partitionSize) % partitionSize;
         mapIndex = ii + partitionSize * (jj + partitionSize * kk);
         if(map[mapIndex + offset] == -1)
         {
            // Allocate an array one block larger than the data array
            int psize = partitionSize * partitionSize * partitionSize;
            T* new_data = new T[data_size + psize];

            for(long x = data_size; x < data_size + psize; ++x)
               new_data[x] = defVal;

            // Memcpy over old stuff from data
            if( data )
            {
               std::memcpy(new_data, data, data_size * sizeof(T));
               delete [] data;
            }

            data = new_data;
            map[mapIndex + offset] = data_size;
            data_size += psize;

            nbOccupied++;
            nbDataPartitionsAllocated++;
         }

         offset = map[mapIndex + offset];
         ii = i % partitionSize;
         jj = j % partitionSize;
         kk = k % partitionSize;
         mapIndex = ii + partitionSize * (jj + partitionSize * kk);

         data[mapIndex + offset] = val;
      }
   }

   const T eval( const Vector& P ) const
   {
       if( !isInside(P) ){ return defVal; }
       int ix, iix;
       float wx, wwx;
       int iy, iiy;
       float wy, wwy;
       int iz, iiz;
       float wz, wwz;

       getLinearInterpolation( P,
                               ix, iix, wx,  wwx,
			       iy, iiy, wy,  wwy,
			       iz, iiz, wz,  wwz );

       T accum = get( ix, iy, iz ) * wwx * wwy * wwz;
       accum  += get( iix, iy, iz ) *  wx * wwy * wwz;
       accum  += get( ix, iiy, iz ) * wwx *  wy * wwz;
       accum  += get( iix, iiy, iz ) *  wx *  wy * wwz;
       accum  += get( ix, iy, iiz ) * wwx * wwy *  wz;
       accum  += get( iix, iy, iiz ) *  wx * wwy *  wz;
       accum  += get( ix, iiy, iiz ) * wwx *  wy *  wz;
       accum  += get( iix, iiy, iiz )*  wx *  wy *  wz;

       return accum;
   }

   ~FSGrid() {}

   void deallocate()
   {
      if( !map )
         delete [] map;

      if( !data )
         delete [] data;

      map = NULL;
      data = NULL;
      map_size = 0;
      data_size = 0;

      if( !mapOfMap )
         delete [] mapOfMap;
   }

   void clear()
   {
      deallocate();

      mapOfMap = new long[totalSize];
      for(long i = 0; i < totalSize; ++i)
         mapOfMap[i] = -1;

      nbOccupied = 0;
      nbMapPartitionsAllocated = 0;
      nbDataPartitionsAllocated = 0;
   }

   long Size() const { return nbOccupied; }
   long NbPartitions() const { return nnx*nny*nnz; }

   int blockSize() const { return partitionSize; }

   const bool goodBlock( long i ) const 
   {
      int Nx = nnx * partitionSize * partitionSize;
      int Ny = nny * partitionSize * partitionSize;
      int Nz = nnz * partitionSize * partitionSize;
      int ii = i / (Nx * Ny);
      int jj = (i % (Nx * Ny)) / Nz;
      int kk = (i % (Nx * Ny)) % Nz;
      return goodBlock(ii, jj, kk);
   }

   const bool goodBlock( const int i, const int j, const int k ) const
   {
      int ii = i / (partitionSize * partitionSize);
      int jj = j / (partitionSize * partitionSize);
      int kk = k / (partitionSize * partitionSize);
      int mapIndex = ii + nnx * (jj + nny * kk);
      if(mapOfMap[mapIndex] == -1)
      {
         return false;
      }
      else
      {
         int offset = mapOfMap[mapIndex];
         ii = (i % (partitionSize * partitionSize)) / partitionSize;
         jj = (j % (partitionSize * partitionSize)) / partitionSize;
         kk = (k % (partitionSize * partitionSize)) / partitionSize;
         mapIndex = ii + partitionSize * (jj + partitionSize * kk);
         if(map[mapIndex + offset] == -1)
            return false;
         else
            return true;
      }
   }

   void blockBounds( long block, int& i0, int& j0, int& k0, int& i1, int& j1, int& k1 ) const
   {
      if ( block <0 || block >= totalSize || nbOccupied == 0 )
      {
         i0 = j0 = k0 = i1 = j1 = k1 = -1;
         return;
      }

      k0 = block / (nnx*nny);
      long bblock = block - k0*nnx*nny;
      j0 = bblock / nnx;
      i0 = bblock - j0*nnx;

      i0 *= partitionSize;
      j0 *= partitionSize;
      k0 *= partitionSize;

      i1 = i0 + partitionSize;
      j1 = j0 + partitionSize;
      k1 = k0 + partitionSize;

      return;
   }

   const long index(int i,int j,int k) const
   {
      int ii = i/partitionSize;
      int jj = j/partitionSize;
      int kk = k/partitionSize;
      return ii + nnx*( jj + nny*kk );
   }

   const bool inside(int i, int j, int k) const
   {
      return (
         (i>-1) && (i<nX) &&
         (j>-1) && (j<nY) &&
         (k>-1) && (k<nZ)
      );
   }

   void initCUDA() {}

   gilligan::DataObject<T>* cuData() const { return cu_data; }
   gilligan::DataObject<long>* cuMap() const { return cu_map; }
   gilligan::DataObject<long>* cuMapOfMap() const { return cu_mapOfMap; }

private:
   long* mapOfMap;
   long* map;
   T* data;
   T defVal;

   gilligan::DataObject<long>* cu_mapOfMap;
   gilligan::DataObject<long>* cu_map;
   gilligan::DataObject<T>* cu_data;

   int partitionSize;
   int nnx, nny, nnz;
   size_t map_size, data_size;

   long nbOccupied;
   long nbMapPartitionsAllocated;
   long nbDataPartitionsAllocated;

   long totalSize;

};




typedef std::shared_ptr<FSGrid<float> > ScalarFrustumGridBase;
typedef std::shared_ptr<FSGrid<Vector> > VectorFrustumGridBase;
typedef std::shared_ptr<FSGrid<Color> > ColorFrustumGridBase;
typedef std::shared_ptr<FSGrid<Matrix> > MatrixFrustumGridBase;



class ScalarFrustumGrid : public ScalarFrustumGridBase
{
  public:

    ScalarFrustumGrid( FSGrid<float>* f );
   ~ScalarFrustumGrid();


    ScalarFrustumGrid operator+=( const ScalarField& e2 );
    ScalarFrustumGrid operator-=( const ScalarField& e2 );
    ScalarFrustumGrid operator-();
    ScalarFrustumGrid operator*=( const ScalarField& e2 );
    ScalarFrustumGrid operator/=( const ScalarField& e2 );


};


class VectorFrustumGrid : public VectorFrustumGridBase
{
  public:

    VectorFrustumGrid( FSGrid<Vector>* f );
   ~VectorFrustumGrid();


    VectorFrustumGrid operator+=( const VectorField& e2 );
    VectorFrustumGrid operator-=( const VectorField& e2 );
    VectorFrustumGrid operator-();
    VectorFrustumGrid operator*=( const ScalarField& e2 );
    VectorFrustumGrid operator/=( const ScalarField& e2 );
};





class ColorFrustumGrid : public ColorFrustumGridBase
{
  public:

    ColorFrustumGrid( FSGrid<Color>* f );
   ~ColorFrustumGrid();


    ColorFrustumGrid operator+=( const ColorField& e2 );
    ColorFrustumGrid operator-=( const ColorField& e2 );
    ColorFrustumGrid operator-();
    ColorFrustumGrid operator*=( const ScalarField& e2 );
    ColorFrustumGrid operator/=( const ScalarField& e2 );


};



class MatrixFrustumGrid : public MatrixFrustumGridBase
{
  public:

    MatrixFrustumGrid( FSGrid<Matrix>* f );
   ~MatrixFrustumGrid();


    MatrixFrustumGrid operator+=( const MatrixField& e2 );
    MatrixFrustumGrid operator-=( const MatrixField& e2 );
    MatrixFrustumGrid operator-();
    MatrixFrustumGrid operator*=( const ScalarField& e2 );
    MatrixFrustumGrid operator/=( const ScalarField& e2 );


};

void initCUDA( ScalarFrustumGrid& grid);
void initCUDA( VectorFrustumGrid& grid);
void initCUDA( ColorFrustumGrid& grid);
void initCUDA( MatrixFrustumGrid& grid);




}




#endif
