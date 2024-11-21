

#ifndef __GRIDVOLUME_H__
#define __GRIDVOLUME_H__


#include "VolumeGrid.h"
#include "Volume.h"
#include "ProgressMeter.h"
#include "SparseGrid.h"

namespace lux
{

void Union( VolumeGrid<float>& grid, const Volume<float> & v );
void Intersection( VolumeGrid<float>& grid, const Volume<float> & v );
void Union( VolumeGrid<float>& grid, const Volume<float> & v, const Vector& llc, const Vector& urc );
void Intersection( VolumeGrid<float>& grid, const Volume<float> & v, const Vector& llc, const Vector& urc );


template<class U>
void Sample( VolumeGrid<U>* grid, const Volume<U>* field )
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
	    U value = field->eval( grid->evalP(i,j,k) );
	    grid->value(i,j,k) = value;
	    meter.update();
         }
      }
   }
}



void Sample( SparseGrid* grid, const Volume<float>* field );
void Sample( SparseColorGrid* grid, const Volume<Color>* field );
void Sample( SparseVectorGrid* grid, const Volume<Vector>* field );


}

#endif
