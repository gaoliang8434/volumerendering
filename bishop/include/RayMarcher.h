

#ifndef __RAYMARCHER_H__
#define __RAYMARCHER_H__

#include <vector>

#include "Vector.h"
#include "Color.h"
#include "Volume.h"
#include "VolumeGrid.h"
#include "AARectangle.h"
#include "SparseGrid.h"
#include "PhaseFunction.h"
#include "IntervalTree.h"
#include "VolumeLight.h"

namespace lux
{

using namespace std;


class VShader;

class RenderData
{
  public:
   RenderData() :
      use_dsm_range      (false),
      dsm_range          (1.0),
      scatterCoefficient (Color(1,1,1,1)),
      ds                 (0.01),
      maxPathlength      (100000000.0),
      sparseGrid         (0),
      phaseFunction      (PhaseFunction( new UniformPhaseFunction(1.0) )),
      useHoldOut         (false),
      intervalTree       (makeIntervalSet( Vector(-1,-1,-1), Vector(1,1,1),0,0,0 ) ),
      vlf_samples        (4)
   {}

   ~RenderData(){}

   ScalarField densityField;
   ScalarField ambientDensityField;
   vector<ScalarField> dsmField;
   bool use_dsm_range;
   double dsm_range;
   vector<Color> lightColor;
   vector<Vector> lightPosition;
   ColorField colorField;  
   ColorField ambientColorField;  
   vector<Vector> startPosition;
   vector<Vector> startDirection;
   Color scatterCoefficient;
   float ds;
   float maxPathlength;

   vector<AABB> boundingBoxes;
   ScalarGrid sparseGrid;

   vector<VShader*> shaders;

   PhaseFunction phaseFunction;

   ScalarField holdOut;
   bool useHoldOut;

   IntervalSet intervalTree;

   vector<VolumeLightField> volumeLights;
   int vlf_samples; 
};

void AddDSM( RenderData* d, const ScalarField& dsm );
const ScalarField& GetDSM( RenderData* d, int i );

void SetDensityField(  RenderData* d, const ScalarField& field );
void SetHoldOut(  RenderData* d, const ScalarField& field );
void SetAmbientDensityField(  RenderData* d, const ScalarField& field );
void SetColorField(  RenderData* d, const ColorField& field );
void SetAmbientColorField(  RenderData* d, const ColorField& field );
void SetSparseGrid(  RenderData* d, ScalarGrid& field );
void SetUniformPhaseFunction( RenderData*d, float value = 1.0/(4.0*M_PI) );
void SetHenyeyGreensteinPhaseFunction( RenderData*d, float g );
void SetDoubleHenyeyGreensteinPhaseFunction( RenderData*d, float g0, float g1, float mix );
void SetFournierForandPhaseFunction( RenderData*d, float en, float mu );
void AddBoundingBox( RenderData *d, const Vector llc, const Vector urc );
void AddBoundingBox( RenderData *d, const AABB& aabb );
void AddBoundingBoxes( RenderData *d, const std::vector<AABB>& boxes );
void SetIntervalTree( RenderData *d, IntervalSet t );
void SetDSMRange( RenderData *d, const double value );
void AddVolumeLight( RenderData *d, VolumeLightField& vlf );
void SetVolumeLightSamples( RenderData *d, int samples );

void RayMarchDSMAccumulation( const ScalarField& densityField, 
                              const Vector& lightPosition, 
                              float ds,
                              ScalarGrid& dsmField
                            );

void RayMarchDSMAccumulation( const RenderData& input, 
                              const Vector& lightPosition, 
                              ScalarGrid& dsmField
                            );

void RayMarchDSMAccumulation( const RenderData& input, 
                              const Vector& lightPosition, 
                              ScalarFrustumGrid& dsmField
                            );


void RayMarchVisibleDSMAccumulation( const RenderData& input, 
                                     const Vector& lightPosition, 
                                     const Camera& cam,
                                     ScalarGrid& dsmField
                                   );


void RayMarchDSMAccumulation( const ScalarField& densityField, 
                              const Vector& lightPosition, 
                              float ds,
                              ScalarGrid& dsmField,
                              ScalarField& holdout
                            );

void RayMarchDSMAccumulation( const ScalarField& densityField, 
                              ScalarFrustumGrid& dsmField,
			      int nbsamples
                            );




void ssRayMarchAccumulation( const RenderData& input, vector<Color>& output_color, vector<Color>& output_opacity );
void ssPointToPointRayMarchAccumulation( const RenderData& input, const Vector& start, const Vector& end, Color& output_color, Color& output_opacity );

//void renderPointLight( const RenderData& input, vector<Color>& output, const float psf_width );


const Color gatherLight( const RenderData& data, const Vector P, const Vector Dir );

void setNbCores( int nb );


//void EndPointsRayMarch( const Vector& start, const Vector& end, const ColorField& attenuation, const ScalarField& integrand, const float ds,  double& result, double& extinction );

const Color PointToPointTransmissivity( const RenderData& d, const Vector& start, const Vector& end );
const float PointToPointLWP( const RenderData& d, const Vector& start, const Vector& end );

}
#endif
