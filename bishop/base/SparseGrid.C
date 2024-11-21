#include "SparseGrid.h"
#include "Volume.h"

#include <stdio.h>
#include <stdlib.h>


using namespace lux;

SparseGrid::SparseGrid(){
   defVal = 0.0;
   partitionSize = 16;
   nbBlocks = 0;
   data = 0;
}


SparseGrid::SparseGrid( int psize ){
   defVal = 0.0;
   partitionSize = psize;
   nbBlocks = 0;
   data = 0;
}


SparseGrid::~SparseGrid(){
   for(int i = 0; i < nbBlocks; i++){
      if(data[i] != NULL){
         delete [] data[i];
      }
   }
   delete [] data;
}

void SparseGrid::setDefVal(float newVal){
   defVal = newVal;
}

const float SparseGrid::getDefVal() const {
   return defVal;
}

void SparseGrid::setPartitionSize(int newSize){
   partitionSize = newSize;
}

void SparseGrid::init(int i, int j, int k, float dx, float dy, float dz, const Vector &Origin){
   RectangularGrid::init(i, j, k, dx, dy, dz, Origin);

   nnx = nx()/partitionSize;
   if( nx() > nnx*partitionSize ){ nnx++; }
   nny = ny()/partitionSize;
   if( ny() > nny*partitionSize ){ nny++; }
   nnz = nz()/partitionSize;
   if( nz() > nnz*partitionSize ){ nnz++; }
   nbBlocks = nnx*nny*nnz;

   data = new float* [nbBlocks];
   for(int i = 0; i < nbBlocks; i++){
      data[i] = NULL;
   }
}

const int SparseGrid::index(int i, int j, int k) const{  
   int ii = i/partitionSize;
   int jj = j/partitionSize;
   int kk = k/partitionSize;
   return ii + nnx*( jj + nny*kk );
}

void SparseGrid::blockBounds( int block, int& i0, int& j0, int& k0, int& i1, int& j1, int& k1 ) const
{
   if ( block < 0 || block >= nbBlocks )
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

const float SparseGrid::get(int i, int j, int k) const {
   int myIndex = index(i, j, k);
   if(data[myIndex] == NULL){
      return(defVal);
   }
   else{
      int ii = i % partitionSize;
      int jj = j % partitionSize;
      int kk = k % partitionSize;
      int partitionIndex = ii + partitionSize*( jj + partitionSize * kk );
      return data[myIndex][partitionIndex];
   }	
}

void SparseGrid::set(float val,int i,int j,int k){
   if(val != defVal){
      int myIndex = index(i, j, k);  
      if(data[myIndex] == NULL){
         data[myIndex] = new float[partitionSize * partitionSize * partitionSize];
         for(int i = 0; i < partitionSize * partitionSize * partitionSize; i++){
            data[myIndex][i] = defVal;
         }
      }
      int ii = i % partitionSize;
      int jj = j % partitionSize;
      int kk = k % partitionSize;
      int partitionIndex = ii + partitionSize*( jj + partitionSize * kk );
      data[myIndex][partitionIndex] = val;
   }
}


const float SparseGrid::eval( const Vector& P ) const
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

       float accum = get( ix, iy, iz ) * wwx * wwy * wwz;
       accum += get( iix, iy, iz ) *  wx * wwy * wwz;
       accum += get( ix, iiy, iz ) * wwx *  wy * wwz;
       accum += get( iix, iiy, iz ) *  wx *  wy * wwz;
       accum += get( ix, iy, iiz ) * wwx * wwy *  wz;
       accum += get( iix, iy, iiz ) *  wx * wwy *  wz;
       accum += get( ix, iiy, iiz ) * wwx *  wy *  wz;
       accum += get( iix, iiy, iiz )*  wx *  wy *  wz;

       return accum;
}
        
void lux::ReadVolumeGrid(SparseGrid& grid, ifstream& out)
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
   
   float value;
   
   out.read( (char*)&nx, sizeof(nx) );
   out.read( (char*)&ny, sizeof(ny) );
   out.read( (char*)&nz, sizeof(nz) );
   out.read( (char*)&Lx, sizeof(Lx) );
   out.read( (char*)&Ly, sizeof(Ly) );
   out.read( (char*)&Lz, sizeof(Lz) );
       
   out.read( (char*)&x, sizeof(x) );
   out.read( (char*)&y, sizeof(y) );
   out.read( (char*)&z, sizeof(z) );
            
   grid.init( nx, ny, nz, Lx, Ly, Lz, Vector(x,y,z) );
               
   //while(out.read((char *)&i, sizeof(int))!=NULL)
   while(out.read((char *)&i, sizeof(int)))
   {
       out.read((char *)&j, sizeof(int));
       out.read((char *)&k, sizeof(int));
       out.read((char *)&value, sizeof(float));
       grid.set(value, i, j, k);
   }
}   
 
void lux::WriteVolumeGrid(SparseGrid& grid, ofstream& out )
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
   
   float temp;
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
               out.write((char*)&temp, sizeof(float));
            }
         }

      }
   
   }
}









SparseColorGrid::SparseColorGrid(){
   defVal = Color(0.0,0.0,0.0,0.0);
   partitionSize = 16;
}

SparseColorGrid::~SparseColorGrid(){
   int totalSize = int(ceil(1.0 *nx()/partitionSize)*ceil(1.0*ny()/partitionSize)*ceil(1.0*nz()/partitionSize));

   for(int i = 0; i < totalSize; i++){
      if(data[i] != NULL){
         delete [] data[i];
      }
   }
   delete [] data;
}

void SparseColorGrid::setDefVal(Color newVal){
   defVal = newVal;
}

const Color& SparseColorGrid::getDefVal() const {
   return defVal;
}

void SparseColorGrid::setPartitionSize(int newSize){
   partitionSize = newSize;
}

void SparseColorGrid::init(int i, int j, int k, float dx, float dy, float dz, const Vector &Origin){
   RectangularGrid::init(i, j, k, dx, dy, dz, Origin);
   int totalSize = int(ceil(1.0 *nx()/partitionSize)*ceil(1.0*ny()/partitionSize)*ceil(1.0*nz()/partitionSize));

   data = new Color* [totalSize];
   for(int i = 0; i < totalSize; i++){
      data[i] = NULL;
   }
}

const int SparseColorGrid::index(int i, int j, int k) const{  
   int ii = i/partitionSize;
   int jj = j/partitionSize;
   int kk = k/partitionSize;
   int nnx = nx()/partitionSize;
   int nny = ny()/partitionSize;
   return ii + nnx*( jj + nny*kk );
}


const Color& SparseColorGrid::get(int i, int j, int k) const {
   int myIndex = index(i, j, k);
   if(data[myIndex] == NULL){
      return(defVal);
   }
   else{
      int ii = i % partitionSize;
      int jj = j % partitionSize;
      int kk = k % partitionSize;
      int partitionIndex = ii + partitionSize*( jj + partitionSize * kk );
      return data[myIndex][partitionIndex];
   }	
}

void SparseColorGrid::set(const Color& val,int i,int j,int k){
   if(val != defVal){
      int myIndex = index(i, j, k);  
      if(data[myIndex] == NULL){
         data[myIndex] = new Color[partitionSize * partitionSize * partitionSize];
         for(int i = 0; i < partitionSize * partitionSize * partitionSize; i++){
            data[myIndex][i] = defVal;
         }
      }
      int ii = i % partitionSize;
      int jj = j % partitionSize;
      int kk = k % partitionSize;
      int partitionIndex = ii + partitionSize*( jj + partitionSize * kk );
      data[myIndex][partitionIndex] = val;
   }
}


const Color SparseColorGrid::eval( const Vector& P ) const
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

       Color accum = get( ix, iy, iz ) * wwx * wwy * wwz;
       accum += get( iix, iy, iz ) *  wx * wwy * wwz;
       accum += get( ix, iiy, iz ) * wwx *  wy * wwz;
       accum += get( iix, iiy, iz ) *  wx *  wy * wwz;
       accum += get( ix, iy, iiz ) * wwx * wwy *  wz;
       accum += get( iix, iy, iiz ) *  wx * wwy *  wz;
       accum += get( ix, iiy, iiz ) * wwx *  wy *  wz;
       accum += get( iix, iiy, iiz )*  wx *  wy *  wz;

       return accum;
}

void SparseColorGrid::normalize( const SparseGrid& g )
{
   // Need to do this directly in the partition blocks
   ProgressMeter meter( nx()*ny()*nz(), "Sparse Grid Normalize");
   if( g.nx() != nx() || g.ny() != ny() || g.nz() != nz() )
   {
   for( int k=0;k<nz();k++ )
   {
      for( int j=0;j<ny();j++ )
      {
         for( int i=0;i<nx();i++ )
	 {
	    float gvalue = g.eval( evalP(i,j,k) );
	    if( gvalue > 0 )
	    {
	       Color cd = get(i,j,k);
	       cd /= gvalue;
	       set( cd, i,j,k );
	    }
	    meter.update();
	 }
      }
   }
   }
   else
   {
   for( int k=0;k<nz();k++ )
   {
      for( int j=0;j<ny();j++ )
      {
         for( int i=0;i<nx();i++ )
	 {
	    float gvalue = g.get( i,j,k );
	    if( gvalue > 0 )
	    {
	       Color cd = get(i,j,k);
	       cd /= gvalue;
	       set( cd, i,j,k );
	    }
	    meter.update();
	 }
      }
   }
   }
}

       
void lux::ReadVolumeGrid(SparseColorGrid& grid, ifstream& out)
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
   
   Color value;
   
   out.read( (char*)&nx, sizeof(nx) );
   out.read( (char*)&ny, sizeof(ny) );
   out.read( (char*)&nz, sizeof(nz) );
   out.read( (char*)&Lx, sizeof(Lx) );
   out.read( (char*)&Ly, sizeof(Ly) );
   out.read( (char*)&Lz, sizeof(Lz) );
       
   out.read( (char*)&x, sizeof(x) );
   out.read( (char*)&y, sizeof(y) );
   out.read( (char*)&z, sizeof(z) );
            
   grid.init( nx, ny, nz, Lx, Ly, Lz, Vector(x,y,z) );
               
   float tempr, tempg, tempb, tempa;
   //while(out.read((char *)&i, sizeof(int))!=NULL)
   while(out.read((char *)&i, sizeof(int)))
   {
       out.read((char *)&j, sizeof(int));
       out.read((char *)&k, sizeof(int));
       out.read((char *)&tempr, sizeof(float));
       out.read((char *)&tempg, sizeof(float));
       out.read((char *)&tempb, sizeof(float));
       out.read((char *)&tempa, sizeof(float));
       value[0] = tempr;
       value[1] = tempg;
       value[2] = tempb;
       value[3] = tempa;
       grid.set(value, i, j, k);
   }
}   
 
void lux::WriteVolumeGrid(SparseColorGrid& grid, ofstream& out )
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
   
   Color temp;
   float tempr, tempg, tempb, tempa;
   for( int k=0;k<nz;k++ )   
   {    
      for( int j=0;j<ny;j++ )  
      {
         for( int i=0;i<nx;i++ )
         {
            temp = grid.get(i,j,k);
            if(temp != grid.getDefVal()){
	       tempr = temp[0];
	       tempg = temp[1];
	       tempb = temp[2];
	       tempa = temp[3];
               out.write((char*)&i, sizeof(int));
               out.write((char*)&j, sizeof(int));
               out.write((char*)&k, sizeof(int));
               out.write((char*)&tempr, sizeof(float));
               out.write((char*)&tempg, sizeof(float));
               out.write((char*)&tempb, sizeof(float));
               out.write((char*)&tempa, sizeof(float));
            }
         }

      }
   
   }
}


void lux::Blur( SparseGrid& g )
{
   SparseGrid temp;
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
	    float sum = 0;
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





void lux::Blur( SparseColorGrid& g )
{
   SparseColorGrid temp;
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
	    Color sum;
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







SparseVectorGrid::SparseVectorGrid(){
   defVal = Vector(0.0,0.0,0.0);
   partitionSize = 16;
}

SparseVectorGrid::~SparseVectorGrid(){
   int totalSize = int(ceil(1.0 *nx()/partitionSize)*ceil(1.0*ny()/partitionSize)*ceil(1.0*nz()/partitionSize));

   for(int i = 0; i < totalSize; i++){
      if(data[i] != NULL){
         delete [] data[i];
      }
   }
   delete [] data;
}

void SparseVectorGrid::setDefVal(Vector newVal){
   defVal = newVal;
}

const Vector& SparseVectorGrid::getDefVal() const {
   return defVal;
}

void SparseVectorGrid::setPartitionSize(int newSize){
   partitionSize = newSize;
}

void SparseVectorGrid::init(int i, int j, int k, float dx, float dy, float dz, const Vector &Origin){
   RectangularGrid::init(i, j, k, dx, dy, dz, Origin);
   int totalSize = int(ceil(1.0 *nx()/partitionSize)*ceil(1.0*ny()/partitionSize)*ceil(1.0*nz()/partitionSize));

   data = new Vector* [totalSize];
   for(int i = 0; i < totalSize; i++){
      data[i] = NULL;
   }
}

const int SparseVectorGrid::index(int i, int j, int k) const{  
   int ii = i/partitionSize;
   int jj = j/partitionSize;
   int kk = k/partitionSize;
   int nnx = nx()/partitionSize;
   int nny = ny()/partitionSize;
   return ii + nnx*( jj + nny*kk );
}


const Vector& SparseVectorGrid::get(int i, int j, int k) const {
   int myIndex = index(i, j, k);
   if(data[myIndex] == NULL){
      return(defVal);
   }
   else{
      int ii = i % partitionSize;
      int jj = j % partitionSize;
      int kk = k % partitionSize;
      int partitionIndex = ii + partitionSize*( jj + partitionSize * kk );
      return data[myIndex][partitionIndex];
   }	
}

void SparseVectorGrid::set(const Vector& val,int i,int j,int k){
   if(val != defVal){
      int myIndex = index(i, j, k);  
      if(data[myIndex] == NULL){
         data[myIndex] = new Vector[partitionSize * partitionSize * partitionSize];
         for(int i = 0; i < partitionSize * partitionSize * partitionSize; i++){
            data[myIndex][i] = defVal;
         }
      }
      int ii = i % partitionSize;
      int jj = j % partitionSize;
      int kk = k % partitionSize;
      int partitionIndex = ii + partitionSize*( jj + partitionSize * kk );
      data[myIndex][partitionIndex] = val;
   }
}


const Vector SparseVectorGrid::eval( const Vector& P ) const
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

       Vector accum = get( ix, iy, iz ) * wwx * wwy * wwz;
       accum += get( iix, iy, iz ) *  wx * wwy * wwz;
       accum += get( ix, iiy, iz ) * wwx *  wy * wwz;
       accum += get( iix, iiy, iz ) *  wx *  wy * wwz;
       accum += get( ix, iy, iiz ) * wwx * wwy *  wz;
       accum += get( iix, iy, iiz ) *  wx * wwy *  wz;
       accum += get( ix, iiy, iiz ) * wwx *  wy *  wz;
       accum += get( iix, iiy, iiz )*  wx *  wy *  wz;

       return accum;
}

void SparseVectorGrid::normalize( const SparseGrid& g )
{
   // Need to do this directly in the partition blocks
   ProgressMeter meter( nx()*ny()*nz(), "Sparse Grid Normalize");
   if( g.nx() != nx() || g.ny() != ny() || g.nz() != nz() )
   {
   for( int k=0;k<nz();k++ )
   {
      for( int j=0;j<ny();j++ )
      {
         for( int i=0;i<nx();i++ )
	 {
	    float gvalue = g.eval( evalP(i,j,k) );
	    if( gvalue > 0 )
	    {
	       Vector cd = get(i,j,k);
	       cd /= gvalue;
	       set( cd, i,j,k );
	    }
	    meter.update();
	 }
      }
   }
   }
   else
   {
   for( int k=0;k<nz();k++ )
   {
      for( int j=0;j<ny();j++ )
      {
         for( int i=0;i<nx();i++ )
	 {
	    float gvalue = g.get( i,j,k );
	    if( gvalue > 0 )
	    {
	       Vector cd = get(i,j,k);
	       cd /= gvalue;
	       set( cd, i,j,k );
	    }
	    meter.update();
	 }
      }
   }
   }
}

       
void lux::ReadVolumeGrid(SparseVectorGrid& grid, ifstream& out)
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
   
   Vector value;
   
   out.read( (char*)&nx, sizeof(nx) );
   out.read( (char*)&ny, sizeof(ny) );
   out.read( (char*)&nz, sizeof(nz) );
   out.read( (char*)&Lx, sizeof(Lx) );
   out.read( (char*)&Ly, sizeof(Ly) );
   out.read( (char*)&Lz, sizeof(Lz) );
       
   out.read( (char*)&x, sizeof(x) );
   out.read( (char*)&y, sizeof(y) );
   out.read( (char*)&z, sizeof(z) );
            
   grid.init( nx, ny, nz, Lx, Ly, Lz, Vector(x,y,z) );
               
   float tempr, tempg, tempb;
   //while(out.read((char *)&i, sizeof(int))!=NULL)
   while(out.read((char *)&i, sizeof(int)))
   {
       out.read((char *)&j, sizeof(int));
       out.read((char *)&k, sizeof(int));
       out.read((char *)&tempr, sizeof(float));
       out.read((char *)&tempg, sizeof(float));
       out.read((char *)&tempb, sizeof(float));
       value[0] = tempr;
       value[1] = tempg;
       value[2] = tempb;
       grid.set(value, i, j, k);
   }
}   
 
void lux::WriteVolumeGrid(SparseVectorGrid& grid, ofstream& out )
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
   
   Vector temp;
   float tempr, tempg, tempb;
   for( int k=0;k<nz;k++ )   
   {    
      for( int j=0;j<ny;j++ )  
      {
         for( int i=0;i<nx;i++ )
         {
            temp = grid.get(i,j,k);
            if(temp != grid.getDefVal()){
	       tempr = temp[0];
	       tempg = temp[1];
	       tempb = temp[2];
               out.write((char*)&i, sizeof(int));
               out.write((char*)&j, sizeof(int));
               out.write((char*)&k, sizeof(int));
               out.write((char*)&tempr, sizeof(float));
               out.write((char*)&tempg, sizeof(float));
               out.write((char*)&tempb, sizeof(float));
            }
         }

      }
   
   }
}



void lux::Blur( SparseVectorGrid& g )
{
   SparseVectorGrid temp;
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
	    Vector sum;
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



