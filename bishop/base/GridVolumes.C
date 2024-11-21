


#include "GridVolumes.h"
using namespace lux;


void lux::Union( VolumeGrid<float>& grid, const Volume<float> & v )
{
   float* ptr = grid.rawPtr();
   for( int k=0;k<grid.nz();k++ )
   {
      for( int j=0;j<grid.ny();j++ )
      {
         for( int i=0;i<grid.nx();i++ )
         {
	    Vector P = grid.evalP(i,j,k);
	    float val = v.eval(P);
	    int index = grid.index(i,j,k);
	    if( val > ptr[index] ){ ptr[index] = val; }
         }
      }
   }
}





void lux::Intersection( VolumeGrid<float>& grid, const Volume<float> & v )
{
   float* ptr = grid.rawPtr();
   for( int k=0;k<grid.nz();k++ )
   {
      for( int j=0;j<grid.ny();j++ )
      {
         for( int i=0;i<grid.nx();i++ )
         {
	    Vector P = grid.evalP(i,j,k);
	    float val = v.eval(P);
	    int index = grid.index(i,j,k);
	    if( val < ptr[index] ){ ptr[index] = val; }
         }
      }
   }
}


void lux::Union( VolumeGrid<float>& grid, const Volume<float> & v, const Vector& llc, const Vector& urc )
{
   float* ptr = grid.rawPtr();

   int imin, imax, jmin, jmax, kmin, kmax;
   if( grid.getBox( llc, urc, imin, imax, jmin, jmax, kmin, kmax ) )
   {

   for( int k=kmin;k<=kmax;k++ )
   {
      for( int j=jmin;j<=jmax;j++ )
      {
         for( int i=imin;i<=imax;i++ )
         {
	    Vector P = grid.evalP(i,j,k);
	    float val = v.eval(P);
	    int index = grid.index(i,j,k);
	    if( val > ptr[index] ){ ptr[index] = val; }
         }
      }
   }

   }
}



void lux::Intersection( VolumeGrid<float>& grid, const Volume<float> & v, const Vector& llc, const Vector& urc )
{
   float* ptr = grid.rawPtr();

   int imin, imax, jmin, jmax, kmin, kmax;
   if( grid.getBox( llc, urc, imin, imax, jmin, jmax, kmin, kmax ) )
   {

   for( int k=kmin;k<=kmax;k++ )
   {
      for( int j=jmin;j<=jmax;j++ )
      {
         for( int i=imin;i<=imax;i++ )
         {
	    Vector P = grid.evalP(i,j,k);
	    float val = v.eval(P);
	    int index = grid.index(i,j,k);
	    if( val < ptr[index] ){ ptr[index] = val; }
         }
      }
   }


   }
}









void lux::Sample( SparseGrid* grid, const Volume<float>* field )
{
   int nx = grid->nx();
   int ny = grid->ny();
   int nz = grid->nz();
   ProgressMeter meter( nx*ny*nz, "Sample");
   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
	    float value = field->eval( grid->evalP(i,j,k) );
	    grid->set(value, i,j,k);
	    meter.update();
         }
      }
   }
}


void lux::Sample( SparseColorGrid* grid, const Volume<Color>* field )
{
   int nx = grid->nx();
   int ny = grid->ny();
   int nz = grid->nz();
   ProgressMeter meter( nx*ny*nz, "Sample");
   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
	    Color value = field->eval( grid->evalP(i,j,k) );
	    grid->set(value, i,j,k);
	    meter.update();
         }
      }
   }
}



void lux::Sample( SparseVectorGrid* grid, const Volume<Vector>* field )
{
   int nx = grid->nx();
   int ny = grid->ny();
   int nz = grid->nz();
   ProgressMeter meter( nx*ny*nz, "Sample");
   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
	    Vector value = field->eval( grid->evalP(i,j,k) );
	    grid->set(value, i,j,k);
	    meter.update();
         }
      }
   }
}



