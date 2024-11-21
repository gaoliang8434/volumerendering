#ifndef __GRIDS_H__
#define __GRIDS_H__

#include "Volume.h"
#include "VolumeGrid.h"
#include "SparseGrid.h"
//#include "OVDBGrid.h"
#include "Noise.h"
#include "UniformPRN.h"
#include "PerlinNoise.h"
#include "Color.h"
#include "Particle.h"
#include "IntervalTree.h"

namespace lux
{

ScalarGrid makeGrid( const RectangularGrid& rg, float defValue );
VectorGrid makeGrid( const RectangularGrid& rg, Vector defValue );
ColorGrid makeGrid( const RectangularGrid& rg, Color defValue );
MatrixGrid makeGrid( const RectangularGrid& rg, Matrix defValue );

ScalarGrid makeGrid( const GridBox& rg, float defValue );
VectorGrid makeGrid( const GridBox& rg, Vector defValue );
ColorGrid makeGrid( const GridBox& rg, Color defValue );
MatrixGrid makeGrid( const GridBox& rg, Matrix defValue );

ScalarGrid makeGrid( const GridBox& rg, float defValue, int partitionSize );
VectorGrid makeGrid( const GridBox& rg, Vector defValue, int partitionSize );
ColorGrid makeGrid( const GridBox& rg, Color defValue, int partitionSize );
MatrixGrid makeGrid( const GridBox& rg, Matrix defValue, int partitionSize );


ScalarGrid makeScalarGrid( const RectangularGrid& rg, float defValue );
VectorGrid makeVectorGrid( const RectangularGrid& rg, Vector defValue );
ColorGrid makeColorGrid( const RectangularGrid& rg, Color defValue );
MatrixGrid makeMatrixGrid( const RectangularGrid& rg, Matrix defValue );


void Blur( ScalarGrid& g );
void Blur( VectorGrid& g );
void Blur( ColorGrid& g );



void writeGrid( ScalarGrid& grid, string fname );
void writeScalarGrid( ScalarGrid& grid, string fname );
void readGrid( ScalarGrid& grid, string fname );
void readScalarGrid( ScalarGrid& grid, string fname );

void writeGrid( VectorGrid& grid, string fname );
void writeVectorGrid( VectorGrid& grid, string fname );
void readGrid( VectorGrid& grid, string fname );
void readVectorGrid( VectorGrid& grid, string fname );

void writeGrid( ColorGrid& grid, string fname );
void writeColorGrid( ColorGrid& grid, string fname );
void readGrid( ColorGrid& grid, string fname );
void readColorGrid( ColorGrid& grid, string fname );



VectorGrid optimumVelocityFromGrad( const MatrixField& m, const GridBox& gb, const int nb_iterations );

VectorField gradientStretchCMFromGrad( const MatrixField& m, const GridBox& gb, const int nb_iterations, const float T, const int nbgs  );



void stampNoise( ScalarGrid& g, Noise* noise, const Vector& center, const float radius );


void stampNoise( ScalarGrid& g, Noise* noise, const Vector& center, const float radius, const Vector& vel, const Vector& accel, const float timestep, const int seed = 757584 );

//void StampNoiseAndColor( ScalarGrid& grid, ColorGrid& cgrid,  Noise* noise, const Vector& center, const float radius, const float fade, const Color& Cd );

//void StampNoiseAndColor( VolumeGrid<float>& grid, VolumeGrid<Color >& cgrid,  Noise* noise, const Vector& center, const float radius, const float fade, const Color& Cd, const Vector& vel, const Vector& accel, const float timestep, const int seed = 757584  );

void stampNoise( ScalarGrid& grid, Noise* noise, const ParticleGroupA& particles );
//void StampNoiseAndColor( VolumeGrid<float>& grid, VolumeGrid<Color>& cgrid,  Noise* noise, const ParticleGroupA& particles );


//Volume<Vector>* StampSELMA( VolumeGrid<Vector>& Xgrid, Volume<Vector>* X, Volume<Vector>* velocity, float dt );

//void stampPointWisps( >& grid, VolumeGrid<Color>& cgrid, const ParticleGroupA& particles );


//void StampPointWisps( SparseGrid& grid, SparseColorGrid& cgrid, const vector<string> files, const vector<Color> typeColors, vector<float>& opacity, float blurScale  );

void stampPointWisps( ScalarGrid& grid, const ParticleGroupA& particles );
//void StampPointWisps( SparseGrid& grid, SparseColorGrid& cgrid, const ParticleGroupA& particles );
//void StampFlameWisps( SparseGrid& grid, SparseColorGrid& cgrid, const ParticleGroupA& particles, const Volume<Color>& baseColor );

void stampBlurredWisps( ScalarGrid& grid, const Vector& P, const float timestep, const Vector& velocity, const Vector& acceleration, const float opacity, int seed );
void stampBlurredWisps( ScalarGrid& grid, ColorGrid& cgrid, const Vector& P, const float timestep, const Vector& velocity, const Vector& acceleration, const float opacity, const Color& cd, int seed );


void stamp( ScalarGrid& grid, const ScalarField& field, const int nbsamples );
void stamp( VectorGrid& grid, const VectorField& field, const int nbsamples );
void stamp( ColorGrid& grid, const ColorField& field, const int nbsamples );
void stamp( MatrixGrid& grid, const MatrixField& field, const int nbsamples );


// stamp a value into a spherical volume of a grid.  This is an additive stamp
void stamp( ScalarGrid& grid, const Vector& pos, const float pscale, const float value );
void stamp( VectorGrid& grid, const Vector& pos, const float pscale, const Vector& value );
void stamp( ColorGrid& grid, const Vector& pos, const float pscale, const Color& value );








//void StampPyro(VolumeGrid<float>& grid, VolumeGrid<Color>& cgrid, const Particle& part );
//void StampPyro( VolumeGrid<float>& grid, VolumeGrid<Color>& cgrid, const ParticleGroupA& particles );

void stampParticles( ScalarGrid& grid, const ParticleGroupA& particles, const float timestep );
//void StampParticles( VolumeGrid<float>& grid, VolumeGrid<Color>& cgrid, const ParticleGroupA& particles, const float timestep );

//void StampGridPattern( SparseGrid& grid, SparseColorGrid& cgrid, int spacing, const Color& col );
//void StampGridPattern( SparseGrid& grid, int spacing );

void gridStatistics( const ScalarGrid& g, double& mean, double& stddev, double& max, double& min );
void gridStatistics( const VectorGrid& g, Vector& mean, double& stddev, Vector& max, Vector& min );
double gridMean( const ScalarGrid& g );
double gridStdDev( const ScalarGrid& g );
double gridMax( const ScalarGrid& g );
double gridMin( const ScalarGrid& g );


int getNbAvailablePartitions( const ScalarGrid& g );
int getNbUsedPartitions( const ScalarGrid& g );
int getPartitionSize( const ScalarGrid& g );



void GreenConvolve( const MatrixField& m, VectorGrid& result );
void GreenSurfaceConvolve( const VectorField& X, VectorGrid& result );


// makeGrid() aliases for Frustum grids
ScalarFrustumGrid makeGrid( const FrustumBox& rg, float defValue, int partitionSize=4 );
VectorFrustumGrid makeGrid( const FrustumBox& rg, Vector defValue, int partitionSize=4 );
ColorFrustumGrid makeGrid( const FrustumBox& rg, Color defValue, int partitionSize=4 );
MatrixFrustumGrid makeGrid( const FrustumBox& rg, Matrix defValue, int partitionSize=4 );

ScalarFrustumGrid makeFrustumGrid( const FrustumBox& rg, float defValue );
VectorFrustumGrid makeFrustumGrid( const FrustumBox& rg, Vector defValue );
ColorFrustumGrid makeFrustumGrid( const FrustumBox& rg, Color defValue );
MatrixFrustumGrid makeFrustumGrid( const FrustumBox& rg, Matrix defValue );

ScalarFrustumGrid makeFrustumGrid( const FrustumGrid& rg, float defValue );
VectorFrustumGrid makeFrustumGrid( const FrustumGrid& rg, Vector defValue );
ColorFrustumGrid makeFrustumGrid( const FrustumGrid& rg, Color defValue );
MatrixFrustumGrid makeFrustumGrid( const FrustumGrid& rg, Matrix defValue );

ScalarFrustumGrid makeFrustumGrid( const FrustumBox& rg, float defValue );
VectorFrustumGrid makeFrustumGrid( const FrustumBox& rg, Vector defValue );
ColorFrustumGrid makeFrustumGrid( const FrustumBox& rg, Color defValue );
MatrixFrustumGrid makeFrustumGrid( const FrustumBox& rg, Matrix defValue );

ScalarFrustumGrid makeFrustumGrid( const FrustumBox& rg, float defValue, int partitionSize );
VectorFrustumGrid makeFrustumGrid( const FrustumBox& rg, Vector defValue, int partitionSize );
ColorFrustumGrid makeFrustumGrid( const FrustumBox& rg, Color defValue, int partitionSize );
MatrixFrustumGrid makeFrustumGrid( const FrustumBox& rg, Matrix defValue, int partitionSize );

ScalarFrustumGrid makeScalarFrustumGrid( const FrustumGrid& rg, float defValue );
VectorFrustumGrid makeVectorFrustumGrid( const FrustumGrid& rg, Vector defValue );
ColorFrustumGrid makeColorFrustumGrid( const FrustumGrid& rg, Color defValue );
MatrixFrustumGrid makeMatrixFrustumGrid( const FrustumGrid& rg, Matrix defValue );



















/*
void readScalarOGrid( const string& fname, ScalarOGrid& g );
void readVectorOGrid( const string& fname, VectorOGrid& g );
void readColorOGrid( const string& fname, ColorOGrid& g );

void writeScalarOGrid( const string& fname, const ScalarOGrid& g );
void writeVectorOGrid( const string& fname, const VectorOGrid& g );
void writeColorOGrid( const string& fname, const ColorOGrid& g );

ScalarOGrid makeOGrid( const GridBox& rg, float defValue );
VectorOGrid makeOGrid( const GridBox& rg, Vector defValue );
ColorOGrid makeOGrid( const GridBox& rg, Color defValue );
MatrixOGrid makeOGrid( const GridBox& rg, Matrix defValue );


void stamp( ScalarOGrid& grid, const ScalarField& field, const int nbsamples );
void stamp( VectorOGrid& grid, const VectorField& field, const int nbsamples );
void stamp( ColorOGrid& grid, const ColorField& field, const int nbsamples );
void stamp( MatrixOGrid& grid, const MatrixField& field, const int nbsamples );
*/




IntervalSet getIntervalSet( const ScalarGrid& g, int maxlvl, int minobj );
IntervalSet getIntervalSet( const ScalarFrustumGrid& g, int maxlvl, int minobj );




}

#endif
