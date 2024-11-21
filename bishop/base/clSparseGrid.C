#include "clSparseGrid.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace lux;

cl_SparseGrid::cl_SparseGrid(){
   defVal = 0.0;
   partitionSize = 4;
   numUsedMapMap = 0;
   numUsedMap = 0;
   clData = (float *) malloc(sizeof(float) * partitionSize * partitionSize * partitionSize);
   clMap = (int *) malloc(sizeof(int) * partitionSize * partitionSize * partitionSize);
   clColor =  NULL;
}

void cl_SparseGrid::init(int &nx, int &ny, int &nz,float dx, float dy, float dz,
	   float xo, float yo, float zo){
   //odd shaped grids
  nX = nx; nY = ny; nZ = nz;  
  if(nx% (partitionSize * partitionSize) != 0){
      nX = nx + partitionSize * partitionSize - nx % ( partitionSize* partitionSize);
      nx = nX;
   }
   if(ny % (partitionSize* partitionSize ) != 0){
      nY = ny + partitionSize* partitionSize - ny% (partitionSize* partitionSize);
      ny = nY;
   }
   if(nz % (partitionSize* partitionSize) != 0){ 
      nZ = nz + partitionSize* partitionSize - nz% (partitionSize* partitionSize);
      nz = nZ;   
   }

   dX = dx; dY = dy; dZ = dz;
   x0 = xo; y0 = yo; z0 = zo;

   int mapX = int(ceil(1.0 *nX/partitionSize));
   int mapY = int(ceil(1.0 *ny/partitionSize));
   int mapZ = int(ceil(1.0 *nz/partitionSize));
   
   int totalSize  = int(ceil(1.0 *mapX/partitionSize) * ceil(1.0*mapY/partitionSize)*ceil(1.0*mapZ/partitionSize));
   clMapMap = new int [totalSize];
      
   for(int ndx = 0; ndx < totalSize; ndx++){
        clMapMap[ndx] = -1;
   }
}

void cl_SparseGrid::initColor(int &nx, int &ny, int &nz,
      float dx, float dy, float dz,
	   float xo, float yo, float zo){
  
  init(nx,ny,nz,dx,dy,dz,xo,yo,zo);
  clColor = (Vector4d *) malloc(sizeof(Vector4d) * partitionSize * partitionSize * partitionSize);
}

cl_SparseGrid::~cl_SparseGrid(){
   delete [] clMapMap;
   delete [] clMap;
   delete [] clData;
}

void cl_SparseGrid::setDefVal(float newVal){
   defVal = newVal;
}

const float cl_SparseGrid::getDefVal() const {
   return defVal;
}

void cl_SparseGrid::setPartitionSize(int newSize){
   partitionSize = newSize;
}

int cl_SparseGrid::indexToMap(int i, int j, int k) const{  
   int ii = i/(partitionSize * partitionSize);
   int jj = j/(partitionSize * partitionSize);
   int kk = k/(partitionSize * partitionSize);
   int nnx = nX/(partitionSize * partitionSize);
   int nny = nY/(partitionSize * partitionSize);
   int nnz = nZ/(partitionSize * partitionSize);
   return ii + nnx*( jj + nny*kk );
}

float cl_SparseGrid::get(int i, int j, int k) const {
   int myMapMapIndex = indexToMap(i, j, k);
   
   if(clMapMap[myMapMapIndex] == -1){
      return(defVal);
   }
   else{
      int mapIndex = clMapMap[myMapMapIndex];
      mapIndex *= partitionSize * partitionSize * partitionSize;
     
      int ii = (i/partitionSize) % partitionSize;
      int jj = (j/partitionSize) % partitionSize;
      int kk = (k/partitionSize) % partitionSize;
      
      int mapPartIndex = ii + partitionSize*( jj + partitionSize * kk );

      if(clMap[mapIndex + mapPartIndex] == -1){
         return(defVal);
      }
      else{
         int dataIndex = clMap[mapIndex + mapPartIndex];
         dataIndex *= partitionSize * partitionSize * partitionSize;
         ii = i % partitionSize;
         jj = j % partitionSize;
         kk = k % partitionSize;
         int dataPartIndex = ii + partitionSize*( jj + partitionSize * kk );
         return clData[dataIndex + dataPartIndex];
     }
   }	
}

Vector4d cl_SparseGrid::getColor( int i, int j, int k ) const{
   int myMapIndex = indexToMap(i, j, k);
   if(clMapMap[myMapIndex] == -1){
      return(defVal);
   }
   else{
   
      int mapIndex = clMapMap[myMapIndex];
      mapIndex *= partitionSize * partitionSize * partitionSize;
     
      int ii = (i/partitionSize) % partitionSize;
      int jj = (j/partitionSize) % partitionSize;
      int kk = (k/partitionSize) % partitionSize;
      
      int mapPartIndex = ii + partitionSize*( jj + partitionSize * kk );

      if(clMap[mapIndex + mapPartIndex] == -1){
         return(defVal);
      }
      else{
         int dataIndex = clMap[mapIndex + mapPartIndex];
         dataIndex *= partitionSize * partitionSize * partitionSize;
         ii = i % partitionSize;
         jj = j % partitionSize;
         kk = k % partitionSize;
         int dataPartIndex = ii + partitionSize*( jj + partitionSize * kk );
         return clColor[dataIndex + dataPartIndex];
     }
   }	
}

int cl_SparseGrid::setUpSpace(int i,int j,int k){
   int pSizeCube = partitionSize * partitionSize * partitionSize;
   int myMapMapIndex = indexToMap(i, j,k); 
   int mapIndex = clMapMap[myMapMapIndex];
   mapIndex *= partitionSize * partitionSize * partitionSize;
     
   if(clMapMap[myMapMapIndex] == -1){
      clMapMap[myMapMapIndex] = numUsedMapMap++;
      mapIndex = clMapMap[myMapMapIndex];
      mapIndex *= partitionSize * partitionSize * partitionSize;
      if(clMap == NULL){
         clMap = (int *) malloc(sizeof(int) * pSizeCube);
      }
      else{
        clMap = (int*) realloc(clMap,numUsedMapMap *  sizeof(int) * pSizeCube);
      }
      for(int i = 0; i < pSizeCube; i++){
         clMap[mapIndex + i] = -1; 
      }
   }  
   int ii = (i/partitionSize) % partitionSize;
   int jj = (j/partitionSize) % partitionSize;
   int kk = (k/partitionSize) % partitionSize;
      
   int mapPartIndex = ii + partitionSize*( jj + partitionSize * kk );
   int dataIndex = clMap[mapIndex + mapPartIndex];
   dataIndex *= pSizeCube;

   if(clMap[mapIndex + mapPartIndex] == -1){
      clMap[mapIndex + mapPartIndex] = numUsedMap++;
      dataIndex = clMap[mapIndex + mapPartIndex];
      dataIndex *= pSizeCube;
      if(clData == NULL){
         clData = (float *) malloc(sizeof(float) * pSizeCube);
      }
         
      else{  
         clData = (float*) realloc(clData,numUsedMap  * sizeof(float) * pSizeCube);
         
      }
      if(clColor != NULL){
         clColor = (Vector4d *) realloc(clColor,numUsedMap*sizeof(Vector4d)*pSizeCube);
      }    
      for(int i = 0; i < pSizeCube; i++){
         clData[dataIndex  + i] = defVal;
       	if(clColor != NULL){
           clColor[dataIndex + i] = defColor;
         }   
      }  
   }    
   ii = i % partitionSize;
   jj = j % partitionSize;
   kk = k % partitionSize;
   int partitionIndex = ii + partitionSize*( jj + partitionSize * kk );
   return dataIndex + partitionIndex;
}

void cl_SparseGrid::set(float val,int i,int j,int k){     
   if(val != defVal && i < nX && j < nY && k < nZ ){   
      int dataIndex = setUpSpace(i,j,k);
      clData[dataIndex] = val;
   }
}

void cl_SparseGrid::setColor(Vector4d val, int i, int j, int k ){
  if(val != defVal && i < nX && j < nY && k < nZ ){   
      int dataIndex = setUpSpace(i,j,k);
      clColor[dataIndex] = val;
   }
}

int * cl_SparseGrid::getMap(){
   return clMap;
}

int * cl_SparseGrid::getMapMap(){
   return clMapMap;
}

int cl_SparseGrid::getNumMapBlocks(){
   if (numUsedMapMap == 0)
      return 1;
   return numUsedMapMap;
}
int cl_SparseGrid::getNumDataBlocks(){
   if (numUsedMap == 0)
      return 1;
   return numUsedMap;
}

int cl_SparseGrid::getPartSize(){
   return partitionSize;
}

int cl_SparseGrid::getDefVal(){
   return defVal;
}

void cl_SparseGrid::sparseValInit(sparseVal *vals){
   if(vals == NULL){
      std::cerr << "NULL SPARSEVAL PTR" << std::endl;
      return;
   }
   vals->defVal = defVal;
   vals->partitionSize = partitionSize;
   vals->nx = nX;
   vals->ny = nY;
   vals->nz = nZ;
   vals->nnx = nX/(partitionSize * partitionSize);
   vals->nny = nY/(partitionSize * partitionSize);
   vals->nnz = nZ/(partitionSize * partitionSize);
   vals->defColorx = defColor[0];
   vals->defColory = defColor[1];
   vals->defColorz = defColor[2];
   vals->defColorw = defColor[3];

}


