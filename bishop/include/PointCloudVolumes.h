#ifndef ____POINTCLOUDVOLUMES_H____
#define ____POINTCLOUDVOLUMES_H____

#include "PointCloud.h"
#include "Volume.h"

namespace lux
{

void warp( PointCloud& pc, const VectorField& X ); 
void advect( PointCloud& pc, const VectorField& u, const float dt ); 
void self_advect( PointCloud& pc, const float dt ); 


void stamp( const PointCloud& pc, ScalarGrid& g, const string attr );
void stamp( const PointCloud& pc, VectorGrid& g, const string attr );
void stamp( const PointCloud& pc, ColorGrid& g, const string attr );

void stampDensityVelocityColor( const PointCloud& pc, ScalarGrid& gs, VectorGrid& gv, ColorGrid& gc );
void stampDensityColor( const PointCloud& pc, ScalarGrid& gs, ColorGrid& gc );

//void stamp( const PointCloud& pc, ScalarGrid& g, const string attr, const string blendmethod );
//void stamp( const PointCloud& pc, VectorGrid& g, const string attr, const string blendmethod );
//void stamp( const PointCloud& pc, ColorGrid& g, const string attr, const string blendmethod );

void stampNoise( const PointCloud& pc, ScalarGrid& g, const float fade );

}
#endif
