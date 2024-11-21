


#ifndef __LEVELSETGEOMETRY_H__
#define __LEVELSETGEOMETRY_H__


#include "VolumeGeometry.h"
#include "SparseGrid.h"



namespace lux
{

const bool ClosestPerpDistance( const Triangle& t, const Vector& P, const double slop,  double& dist );

void GenerateLevelSet( const TriangleGeometry& geom, const int nbIterations, const double slop, const double narrow_band_expansion, ScalarGrid& ls );



}

#endif
