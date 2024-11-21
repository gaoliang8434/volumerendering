#ifndef __INCOMPRESSIBLE_H__
#define __INCOMPRESSIBLE_H__

#include "Volume.h"
#include "SparseGrid.h"



namespace lux{


// Ray march version based on a implicit geometry, starting from each voxel.
const VectorField RMIncompressibleGradient( const VectorField& gf, const GridBox& g, const ScalarField& S, const VectorField& nS, const VectorField& base_nS, const float threshold, const int div_width, const int nbblur );





}







#endif
