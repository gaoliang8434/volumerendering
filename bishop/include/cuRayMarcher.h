
#ifndef __CURAYMARCHER_H__
#define __CURAYMARCHER_H__

#ifndef __PHASEFUNCTIONTYPE__
#define __PHASEFUNCTIONTYPE__

struct cuPhaseFunction
{
   enum PhaseFunctionType
   {
      UNIFORM = 0,
      HENYEYGREENSTEIN = 1,
      DOUBLE_HENYEYGREENSTEIN = 2,
      FOURNIERFORAND = 3
   };

   int type;
   float val;
   float g0, g1, mix;
   float en, mu;
   
   // For Fournier Forand
   float delta180;
   float nu;
   double p1factor;
};

#endif

#include "RayMarcher.h"
#include "Color.h"

// #include "data_interface/include/DataObjectDefs.h"
// #include "cuda/cuFieldStructs.cuh"

using namespace std;

namespace lux
{

class cu_RenderData
{
public:
   cu_RenderData()
   {
      phaseF.type = cuPhaseFunction::PhaseFunctionType::UNIFORM;
      phaseF.val = 1.0;
   }

   ~cu_RenderData(){}

   ScalarGrid densityField;
   ScalarGrid ambientDensityField;
   vector<ScalarFrustumGrid> dsmField;
   vector<Color> lightColor;
   vector<Vector> lightPosition;
   ColorGrid colorField;
   ColorGrid ambientColorField;
   vector<Vector> startPosition;
   vector<Vector> startDirection;
   float scatterCoefficient;
   float ds;
   float maxPathlength;
   float clampv;

   vector<AARectangle> boundingBoxes;

   vector<VShader*> shaders;

   cuPhaseFunction phaseF;

   // cu_ScalarField density;
   // cu_ScalarField ambientDensity;
   // std::unique_ptr< gilligan::DataObject<cu_DSM> > dsmList;
   // cu_ColorField color;
   // cu_ColorField ambientColor;

};

void cu_AddDSM( cu_RenderData* d, const ScalarFrustumGrid& dsm );
const ScalarFrustumGrid& cu_GetDSM( cu_RenderData* d, int i );

void cu_InitFields( cu_RenderData* d );
void cu_SetDensityField(  cu_RenderData* d, const ScalarGrid& field );
void cu_SetAmbientDensityField(  cu_RenderData* d, const ScalarGrid& field );
void cu_SetColorField(  cu_RenderData* d, const ColorGrid& field );
void cu_SetAmbientColorField(  cu_RenderData* d, const ColorGrid& field );
void cu_SetUniformPhaseFunction( cu_RenderData*d, float value = 1.0/(4.0*M_PI) );
void cu_SetHenyeyGreensteinPhaseFunction( cu_RenderData*d, float g );
void cu_SetDoubleHenyeyGreensteinPhaseFunction( cu_RenderData*d, float g0, float g1, float mix );
void cu_SetFournierForandPhaseFunction( cu_RenderData*d, float en, float mu );
void cu_AddBoundingBox( cu_RenderData *d, const Vector llc, const Vector urc );

void cu_ssRayMarchAccumulation( const cu_RenderData& input, vector<Color>& Cd );

}

#endif
