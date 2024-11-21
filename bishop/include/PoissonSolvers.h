

#ifndef __POISSONSOLVERS_H__
#define __POISSONSOLVERS_H__

#include "Volume.h"
#include "VolumeGrid.h"
#include "SparseGrid.h"


namespace lux
{

   void GaussSeidelPoissonSolver( VolumeGrid<float>& p, Volume<float>* source, int nbiterations, float tolerance ); 
   void GaussSeidelDivFree( VolumeGrid<Vector>& divfreeU, Volume<Vector>* U, int nbiterations, float tolerance );
   void FFTDivFree( VolumeGrid<Vector>& divfreeU, Volume<Vector>* U );
   void FFTDivFree( VectorGrid& divfreeU, const VectorField& U );
   void FFTDivFree( VectorGrid& divfreeU );
   VectorField FFTDivFree( const GridBox& b, const VectorField& U );
   VectorField FFTVolumePreserve( const GridBox& b, const VectorField& U );



//   VectorField RayMarchDivFreeZeroNormal( const VectorField& W, const ScalarField& LS, const float resolution, const VectorField& BC );
//   VectorField RayMarchDivFreePlanar( const VectorField& W, const Vector& PlaneP, const Vector& normal, const float resolution, const VectorField& BC );
}


#endif
