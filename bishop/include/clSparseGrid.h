                    
#ifndef __CLSPARSEGRID_H__
#define __CLSPARSEGRID_H__
       
#include <iostream>
#include <fstream>
#include <map>
#include <sys/time.h>
#include "VectorMath.h"


 namespace lux
{
typedef struct sparseValue{
    float defVal;
    float defColorx;
    float defColory;
    float defColorz;
    float defColorw;
    int partitionSize;
    int nx;
    int ny;
    int nz;
    int nnx;
    int nny;
    int nnz;
   } sparseVal;    

template<typename U>
class DSGrid : public RectangularGrid
{
public: 
   DSGrid();
   void init(int&,int&,int&,float,float,float,float,float,float);
   void initColor(int&,int&,int&,float,float,float,float,float,float);
   float get(int,int,int) const;
          
   void setDefVal(float);
   const float getDefVal() const;
   void setPartitionSize(int);
   void set(float,int,int,int);
   void setColor(Vector4d val, int i, int j, int k );

   int *getMapMap();
   int *getMap();
   float* rawPtr() { return clData; }
   void sparseValInit(sparseVal *);
   int getNumMapBlocks();
   int getNumDataBlocks();
   int getPartSize();
   int getDefVal(); 
   
   const float& llcx() const { return x0; }
   const float& llcy() const { return y0; }
   const float& llcz() const { return z0; }
   
   const int& nx() const { return nX; }
   const int& ny() const { return nY; }
   const int& nz() const { return nZ; }
    
   int* nxp() { return &nX; }
   int* nyp() { return &nY; }
   int* nzp() { return &nZ; }
   
   const float& dx() const { return dX; }     
   const float& dy() const { return dY; }
   const float& dz() const { return dZ; }   
    
   int indexToMap(int,int,int) const;
   int index( float i, float j, float k ) const;
  
  Vector4d* rawColorPtr(){ return clColor; }
   void setPeriodic() { periodic = true; }

   
    ~cl_SparseGrid();
private:
   float defVal;
   Vector4d defColor;
   int partitionSize;
   int setUpSpace(int,int,int);
   int *clMap;
   int *clMapMap;
   float *clData;
   Vector4d* clColor;
   int numUsedMap;
   int numUsedMapMap;
   int nX,nY,nZ;
   float x0,y0,z0;
   float dX,dY,dZ;
   bool periodic;
}; 

        
}
//


#endif

