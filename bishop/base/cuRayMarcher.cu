
#include "cuRayMarcher.h"

#include "data_interface/include/DataObject.h"
#include "data_interface/include/DataObjectBishop.h"

#include "data_interface/include/ExecutionPolicy.cuh"

#include "cuda/ssRaymarchKernels.cu"

/*
__global__ void ssRayMarchAccumulation(
   const float* densityField, const long* densityMap, const long* densityMapOfMap,
   int densityWidth, int densityHeight, int densityDepth, int densityPartitionSize,
   float4 densityLLC, float4 densityRes, float densityDefVal,

   const float* dsmField, const long* dsmMap, const long* dsmMapOfMap,
   int dsmWidth, int dsmHeight, int dsmDepth, int dsmPartitionSize,
   float4 dsmLLC, float4 dsmRes, float dsmDefVal,

   const float4* colorField, const long* colorMap, const long* colorMapOfMap,
   int colorWidth, int colorHeight, int colorDepth, int colorPartitionSize,
   float4 colorLLC, float4 colorRes, float4 colorDefVal,

   const float4* startPosition, const float4* startDirection, int nbRays,
   float scatterCoefficient, float ds, float maxPathlength,
   float4* Cd
)
*/

using namespace std;

namespace lux
{

void cu_AddDSM( cu_RenderData* d, const ScalarFrustumGrid& dsm )
{
   d->dsmField.push_back( dsm );
}

const ScalarFrustumGrid& cu_GetDSM( cu_RenderData* d, int i )
{
   return d->dsmField[i];
}

void cu_SetDensityField(  cu_RenderData* d, const ScalarGrid& field )
{
    d->densityField = field;
}
void cu_SetAmbientDensityField(  cu_RenderData* d, const ScalarGrid& field )
{
    d->ambientDensityField = field;
}
void cu_SetColorField(  cu_RenderData* d, const ColorGrid& field )
{
    d->colorField = field;
}
void cu_SetAmbientColorField(  cu_RenderData* d, const ColorGrid& field )
{
    d->ambientColorField = field;
}
void cu_AddBoundingBox( cu_RenderData *d, const Vector llc, const Vector urc )
{
   d->boundingBoxes.push_back( AARectangle(llc, urc) );
}

void cu_InitFields( cu_RenderData* d )
{
/*
   // Setup density field struct
   d->density.data = d->densityField->cuData()->getDevicePointer();
   d->density.map = d->densityField->cuMap()->getDevicePointer();
   d->density.mapOfMap = d->densityField->cuMapOfMap()->getDevicePointer();
   d->density.nx = d->densityField->nx();
   d->density.ny = d->densityField->ny();
   d->density.nz = d->densityField->nz();
   d->density.partitionSize = d->densityField->blockSize();

   int psize = d->densityField->blockSize();
   psize *= psize;

   d->density.nnx = d->density.nx / psize;
   d->density.nnx += (d->density.nnx*psize < d->density.nx) ? 1 : 0;

   d->density.nny = d->density.ny / psize;
   d->density.nny += (d->density.nny*psize < d->density.ny) ? 1 : 0;

   d->density.LLC.x = d->densityField->llc()[0];
   d->density.LLC.y = d->densityField->llc()[1];
   d->density.LLC.z = d->densityField->llc()[2];
   d->density.LLC.w = 0.0f;
   d->density.Res.x = d->densityField->dx();
   d->density.Res.y = d->densityField->dy();
   d->density.Res.z = d->densityField->dz();
   d->density.Res.w = 0.0f;
   d->density.defVal = d->densityField->getDefVal();

   // Setup ambient density field struct
   d->ambientDensity.data = d->ambientDensityField->cuData()->getDevicePointer();
   d->ambientDensity.map = d->ambientDensityField->cuMap()->getDevicePointer();
   d->ambientDensity.mapOfMap = d->ambientDensityField->cuMapOfMap()->getDevicePointer();
   d->ambientDensity.nx = d->ambientDensityField->nx();
   d->ambientDensity.ny = d->ambientDensityField->ny();
   d->ambientDensity.nz = d->ambientDensityField->nz();
   d->ambientDensity.partitionSize = d->ambientDensityField->blockSize();

   psize = d->ambientDensityField->blockSize();
   psize *= psize;

   d->ambientDensity.nnx = d->ambientDensity.nx / psize;
   d->ambientDensity.nnx += (d->ambientDensity.nnx*psize < d->ambientDensity.nx) ? 1 : 0;

   d->ambientDensity.nny = d->ambientDensity.ny / psize;
   d->ambientDensity.nny += (d->ambientDensity.nny*psize < d->ambientDensity.ny) ? 1 : 0;

   d->ambientDensity.LLC.x = d->ambientDensityField->llc()[0];
   d->ambientDensity.LLC.y = d->ambientDensityField->llc()[1];
   d->ambientDensity.LLC.z = d->ambientDensityField->llc()[2];
   d->ambientDensity.LLC.w = 0.0f;
   d->ambientDensity.Res.x = d->ambientDensityField->dx();
   d->ambientDensity.Res.y = d->ambientDensityField->dy();
   d->ambientDensity.Res.z = d->ambientDensityField->dz();
   d->ambientDensity.Res.w = 0.0f;
   d->ambientDensity.defVal = d->ambientDensityField->getDefVal();

   // Setup DSM struct array
   d->dsmList.reset( new gilligan::DataObject<cu_DSM>(d->dsmField.size()) );
   cu_DSM* dsm_h = d->dsmList->getHostPointer();
   for(int i = 0; i < d->dsmField.size(); ++i)
   {
      // Camera properties
      float fov = d->dsmField[i]->camera().fov();
      float ar = d->dsmField[i]->camera().aspectRatio();
      Vector eye = d->dsmField[i]->camera().eye();
      Vector view = d->dsmField[i]->camera().view();
      Vector up = d->dsmField[i]->camera().up();
      Vector right = (view ^ up).unitvector();

      dsm_h[i].htanfov = tan(fov * 0.5 * CUDART_PI_F / 180.0);
      dsm_h[i].vtanfov = dsm_h[i].htanfov / ar;
      dsm_h[i].near = d->dsmField[i]->camera().nearPlane();
      dsm_h[i].far = d->dsmField[i]->camera().farPlane();

      dsm_h[i].camera_position.x = eye[0];
      dsm_h[i].camera_position.y = eye[1];
      dsm_h[i].camera_position.z = eye[2];
      dsm_h[i].camera_position.w = 0.0f;

      dsm_h[i].camera_right.x = right[0];
      dsm_h[i].camera_right.y = right[1];
      dsm_h[i].camera_right.z = right[2];
      dsm_h[i].camera_right.w = 0.0f;

      dsm_h[i].camera_up.x = up[0];
      dsm_h[i].camera_up.y = up[1];
      dsm_h[i].camera_up.z = up[2];
      dsm_h[i].camera_up.w = 0.0f;

      dsm_h[i].camera_view.x = view[0];
      dsm_h[i].camera_view.y = view[1];
      dsm_h[i].camera_view.z = view[2];
      dsm_h[i].camera_view.w = 0.0f;

      // Field properties
      dsm_h[i].field.data = d->dsmField[i]->cuData()->getDevicePointer();
      dsm_h[i].field.map = d->dsmField[i]->cuMap()->getDevicePointer();
      dsm_h[i].field.mapOfMap = d->dsmField[i]->cuMapOfMap()->getDevicePointer();

      dsm_h[i].field.nx = d->dsmField[i]->nx();
      dsm_h[i].field.ny = d->dsmField[i]->ny();
      dsm_h[i].field.nz = d->dsmField[i]->nz();

      dsm_h[i].field.partitionSize = d->dsmField[i]->blockSize();

      int psize = d->dsmField[i]->blockSize();
      psize *= psize;

      dsm_h[i].field.nnx = dsm_h[i].field.nx / psize;
      dsm_h[i].field.nnx += (dsm_h[i].field.nnx*psize < dsm_h[i].field.nx) ? 1 : 0;

      dsm_h[i].field.nny = dsm_h[i].field.ny / psize;
      dsm_h[i].field.nny += (dsm_h[i].field.nny*psize < dsm_h[i].field.ny) ? 1 : 0;

      Vector llc = d->dsmField[i]->evalP(0, 0, 0);
      dsm_h[i].field.LLC.x = llc[0];
      dsm_h[i].field.LLC.y = llc[1];
      dsm_h[i].field.LLC.z = llc[2];
      dsm_h[i].field.LLC.w = 0.0f;

      dsm_h[i].field.Res.x = d->dsmField[i]->dx();
      dsm_h[i].field.Res.y = d->dsmField[i]->dy();
      dsm_h[i].field.Res.z = d->dsmField[i]->dz();
      dsm_h[i].field.Res.w = 0.0f;

      dsm_h[i].field.defVal = d->dsmField[i]->getDefVal();
   }
   d->dsmList->updateDevice();

   // Setup color field struct
   d->color.data = d->colorField->cuData()->getDevicePointer();
   d->color.map = d->colorField->cuMap()->getDevicePointer();
   d->color.mapOfMap = d->colorField->cuMapOfMap()->getDevicePointer();
   d->color.nx = d->colorField->nx();
   d->color.ny = d->colorField->ny();
   d->color.nz = d->colorField->nz();
   d->color.partitionSize = d->colorField->blockSize();

   psize = d->colorField->blockSize();
   psize *= psize;

   d->color.nnx = d->color.nx / psize;
   d->color.nnx += (d->color.nnx*psize < d->color.nx) ? 1 : 0;

   d->color.nny = d->color.ny / psize;
   d->color.nny += (d->color.nny*psize < d->color.ny) ? 1 : 0;

   d->color.LLC.x = d->colorField->llc()[0];
   d->color.LLC.y = d->colorField->llc()[1];
   d->color.LLC.z = d->colorField->llc()[2];
   d->color.LLC.w = 0.0f;
   d->color.Res.x = d->colorField->dx();
   d->color.Res.y = d->colorField->dy();
   d->color.Res.z = d->colorField->dz();
   d->color.Res.w = 0.0f;
   d->color.defVal.x = d->colorField->getDefVal()[0];
   d->color.defVal.y = d->colorField->getDefVal()[1];
   d->color.defVal.z = d->colorField->getDefVal()[2];
   d->color.defVal.w = d->colorField->getDefVal()[3];

   // Setup color field struct
   d->ambientColor.data = d->ambientColorField->cuData()->getDevicePointer();
   d->ambientColor.map = d->ambientColorField->cuMap()->getDevicePointer();
   d->ambientColor.mapOfMap = d->ambientColorField->cuMapOfMap()->getDevicePointer();
   d->ambientColor.nx = d->ambientColorField->nx();
   d->ambientColor.ny = d->ambientColorField->ny();
   d->ambientColor.nz = d->ambientColorField->nz();
   d->ambientColor.partitionSize = d->ambientColorField->blockSize();

   psize = d->ambientColorField->blockSize();
   psize *= psize;

   d->ambientColor.nnx = d->ambientColor.nx / psize;
   d->ambientColor.nnx += (d->ambientColor.nnx*psize < d->ambientColor.nx) ? 1 : 0;

   d->ambientColor.nny = d->ambientColor.ny / psize;
   d->ambientColor.nny += (d->ambientColor.nny*psize < d->ambientColor.ny) ? 1 : 0;

   d->ambientColor.LLC.x = d->ambientColorField->llc()[0];
   d->ambientColor.LLC.y = d->ambientColorField->llc()[1];
   d->ambientColor.LLC.z = d->ambientColorField->llc()[2];
   d->ambientColor.LLC.w = 0.0f;
   d->ambientColor.Res.x = d->ambientColorField->dx();
   d->ambientColor.Res.y = d->ambientColorField->dy();
   d->ambientColor.Res.z = d->ambientColorField->dz();
   d->ambientColor.Res.w = 0.0f;
   d->ambientColor.defVal.x = d->ambientColorField->getDefVal()[0];
   d->ambientColor.defVal.y = d->ambientColorField->getDefVal()[1];
   d->ambientColor.defVal.z = d->ambientColorField->getDefVal()[2];
   d->ambientColor.defVal.w = d->ambientColorField->getDefVal()[3];
*/
}

void cu_ssRayMarchAccumulation( const cu_RenderData& input, vector<Color>& Cd )
{
   int psize;

   Cd.clear();
   Cd.resize( input.startPosition.size() );

   /* Allocate device memory for std::vector<T>s:
    * Cd
    * startPosition
    * startDirection
    * lightColor
    * lightPosition
    */

   gilligan::DataObject<Color> cu_Cd(Cd.size());
   for(int i = 0; i < Cd.size(); ++i)
   {
      (*cu_Cd)[i].x = Cd[i][0];
      (*cu_Cd)[i].y = Cd[i][1];
      (*cu_Cd)[i].z = Cd[i][2];
      (*cu_Cd)[i].w = Cd[i][3];
   }
   cu_Cd.updateDevice();

   gilligan::DataObject<Vector> cu_StartPosition(input.startPosition.size());
   for(int i = 0; i < input.startPosition.size(); ++i)
   {
      (*cu_StartPosition)[i].x = input.startPosition[i][0];
      (*cu_StartPosition)[i].y = input.startPosition[i][1];
      (*cu_StartPosition)[i].z = input.startPosition[i][2];
      (*cu_StartPosition)[i].w = 0.0f;
   }
   cu_StartPosition.updateDevice();

   gilligan::DataObject<Vector> cu_StartDirection(input.startDirection.size());
   for(int i = 0; i < input.startDirection.size(); ++i)
   {
      (*cu_StartDirection)[i].x = input.startDirection[i][0];
      (*cu_StartDirection)[i].y = input.startDirection[i][1];
      (*cu_StartDirection)[i].z = input.startDirection[i][2];
      (*cu_StartDirection)[i].w = 0.0f;
   }
   cu_StartDirection.updateDevice();

   gilligan::DataObject<Color> cu_LightColor(input.lightColor.size());
   for(int i = 0; i < input.lightColor.size(); ++i)
   {
      (*cu_LightColor)[i].x = input.lightColor[i][0];
      (*cu_LightColor)[i].y = input.lightColor[i][1];
      (*cu_LightColor)[i].z = input.lightColor[i][2];
      (*cu_LightColor)[i].w = input.lightColor[i][3];
   }
   cu_LightColor.updateDevice();

   gilligan::DataObject<Vector> cu_LightPosition(input.lightPosition.size());
   for(int i = 0; i < input.lightPosition.size(); ++i)
   {
      (*cu_LightPosition)[i].x = input.lightPosition[i][0];
      (*cu_LightPosition)[i].y = input.lightPosition[i][1];
      (*cu_LightPosition)[i].z = input.lightPosition[i][2];
      (*cu_LightPosition)[i].w = 0.0f;
   }
   cu_LightPosition.updateDevice();

   // Setup density field struct
   cu_ScalarField density;
   density.data = input.densityField->cuData()->getDevicePointer();
   density.map = input.densityField->cuMap()->getDevicePointer();
   density.mapOfMap = input.densityField->cuMapOfMap()->getDevicePointer();
   density.nx = input.densityField->nx();
   density.ny = input.densityField->ny();
   density.nz = input.densityField->nz();
   density.partitionSize = input.densityField->blockSize();

   psize = input.densityField->blockSize();
   psize *= psize;

   density.nnx = density.nx / psize;
   density.nnx += (density.nnx*psize < density.nx) ? 1 : 0;

   density.nny = density.ny / psize;
   density.nny += (density.nny*psize < density.ny) ? 1 : 0;

   density.LLC.x = input.densityField->llc()[0];
   density.LLC.y = input.densityField->llc()[1];
   density.LLC.z = input.densityField->llc()[2];
   density.LLC.w = 0.0f;
   density.Res.x = input.densityField->dx();
   density.Res.y = input.densityField->dy();
   density.Res.z = input.densityField->dz();
   density.Res.w = 0.0f;
   density.defVal = input.densityField->getDefVal();

   // Setup ambient density field struct
   cu_ScalarField ambientDensity;
   ambientDensity.data = input.ambientDensityField->cuData()->getDevicePointer();
   ambientDensity.map = input.ambientDensityField->cuMap()->getDevicePointer();
   ambientDensity.mapOfMap = input.ambientDensityField->cuMapOfMap()->getDevicePointer();
   ambientDensity.nx = input.ambientDensityField->nx();
   ambientDensity.ny = input.ambientDensityField->ny();
   ambientDensity.nz = input.ambientDensityField->nz();
   ambientDensity.partitionSize = input.ambientDensityField->blockSize();

   psize = input.ambientDensityField->blockSize();
   psize *= psize;

   ambientDensity.nnx = ambientDensity.nx / psize;
   ambientDensity.nnx += (ambientDensity.nnx*psize < ambientDensity.nx) ? 1 : 0;

   ambientDensity.nny = ambientDensity.ny / psize;
   ambientDensity.nny += (ambientDensity.nny*psize < ambientDensity.ny) ? 1 : 0;

   ambientDensity.LLC.x = input.ambientDensityField->llc()[0];
   ambientDensity.LLC.y = input.ambientDensityField->llc()[1];
   ambientDensity.LLC.z = input.ambientDensityField->llc()[2];
   ambientDensity.LLC.w = 0.0f;
   ambientDensity.Res.x = input.ambientDensityField->dx();
   ambientDensity.Res.y = input.ambientDensityField->dy();
   ambientDensity.Res.z = input.ambientDensityField->dz();
   ambientDensity.Res.w = 0.0f;
   ambientDensity.defVal = input.ambientDensityField->getDefVal();

   // Setup DSM struct array
   gilligan::DataObject<cu_DSM> dsmList(input.dsmField.size());
   cu_DSM* dsm_h = dsmList.getHostPointer();
   for(int i = 0; i < input.dsmField.size(); ++i)
   {
      // Camera properties
      float fov = input.dsmField[i]->camera().fov();
      float ar = input.dsmField[i]->camera().aspectRatio();
      Vector eye = input.dsmField[i]->camera().eye();
      Vector view = input.dsmField[i]->camera().view();
      Vector up = input.dsmField[i]->camera().up();
      Vector right = (view ^ up).unitvector();

      dsm_h[i].htanfov = tan(fov * 0.5 * CUDART_PI_F / 180.0);
      dsm_h[i].vtanfov = dsm_h[i].htanfov / ar;
      dsm_h[i].near = input.dsmField[i]->camera().nearPlane();
      dsm_h[i].far = input.dsmField[i]->camera().farPlane();

      dsm_h[i].camera_position.x = eye[0];
      dsm_h[i].camera_position.y = eye[1];
      dsm_h[i].camera_position.z = eye[2];
      dsm_h[i].camera_position.w = 0.0f;

      dsm_h[i].camera_right.x = right[0];
      dsm_h[i].camera_right.y = right[1];
      dsm_h[i].camera_right.z = right[2];
      dsm_h[i].camera_right.w = 0.0f;

      dsm_h[i].camera_up.x = up[0];
      dsm_h[i].camera_up.y = up[1];
      dsm_h[i].camera_up.z = up[2];
      dsm_h[i].camera_up.w = 0.0f;

      dsm_h[i].camera_view.x = view[0];
      dsm_h[i].camera_view.y = view[1];
      dsm_h[i].camera_view.z = view[2];
      dsm_h[i].camera_view.w = 0.0f;

      // Field properties
      dsm_h[i].field.data = input.dsmField[i]->cuData()->getDevicePointer();
      dsm_h[i].field.map = input.dsmField[i]->cuMap()->getDevicePointer();
      dsm_h[i].field.mapOfMap = input.dsmField[i]->cuMapOfMap()->getDevicePointer();

      dsm_h[i].field.nx = input.dsmField[i]->nx();
      dsm_h[i].field.ny = input.dsmField[i]->ny();
      dsm_h[i].field.nz = input.dsmField[i]->nz();

      dsm_h[i].field.partitionSize = input.dsmField[i]->blockSize();

      int psize = input.dsmField[i]->blockSize();
      psize *= psize;

      dsm_h[i].field.nnx = dsm_h[i].field.nx / psize;
      dsm_h[i].field.nnx += (dsm_h[i].field.nnx*psize < dsm_h[i].field.nx) ? 1 : 0;

      dsm_h[i].field.nny = dsm_h[i].field.ny / psize;
      dsm_h[i].field.nny += (dsm_h[i].field.nny*psize < dsm_h[i].field.ny) ? 1 : 0;

      Vector llc = input.dsmField[i]->evalP(0, 0, 0);
      dsm_h[i].field.LLC.x = llc[0];
      dsm_h[i].field.LLC.y = llc[1];
      dsm_h[i].field.LLC.z = llc[2];
      dsm_h[i].field.LLC.w = 0.0f;

      dsm_h[i].field.Res.x = input.dsmField[i]->dx();
      dsm_h[i].field.Res.y = input.dsmField[i]->dy();
      dsm_h[i].field.Res.z = input.dsmField[i]->dz();
      dsm_h[i].field.Res.w = 0.0f;

      dsm_h[i].field.defVal = input.dsmField[i]->getDefVal();
   }
   dsmList.updateDevice();

   // Setup color field struct
   cu_ColorField color;
   color.data = input.colorField->cuData()->getDevicePointer();
   color.map = input.colorField->cuMap()->getDevicePointer();
   color.mapOfMap = input.colorField->cuMapOfMap()->getDevicePointer();
   color.nx = input.colorField->nx();
   color.ny = input.colorField->ny();
   color.nz = input.colorField->nz();
   color.partitionSize = input.colorField->blockSize();

   psize = input.colorField->blockSize();
   psize *= psize;

   color.nnx = color.nx / psize;
   color.nnx += (color.nnx*psize < color.nx) ? 1 : 0;

   color.nny = color.ny / psize;
   color.nny += (color.nny*psize < color.ny) ? 1 : 0;

   color.LLC.x = input.colorField->llc()[0];
   color.LLC.y = input.colorField->llc()[1];
   color.LLC.z = input.colorField->llc()[2];
   color.LLC.w = 0.0f;
   color.Res.x = input.colorField->dx();
   color.Res.y = input.colorField->dy();
   color.Res.z = input.colorField->dz();
   color.Res.w = 0.0f;
   color.defVal.x = input.colorField->getDefVal()[0];
   color.defVal.y = input.colorField->getDefVal()[1];
   color.defVal.z = input.colorField->getDefVal()[2];
   color.defVal.w = input.colorField->getDefVal()[3];

   // Setup color field struct
   cu_ColorField ambientColor;
   ambientColor.data = input.ambientColorField->cuData()->getDevicePointer();
   ambientColor.map = input.ambientColorField->cuMap()->getDevicePointer();
   ambientColor.mapOfMap = input.ambientColorField->cuMapOfMap()->getDevicePointer();
   ambientColor.nx = input.ambientColorField->nx();
   ambientColor.ny = input.ambientColorField->ny();
   ambientColor.nz = input.ambientColorField->nz();
   ambientColor.partitionSize = input.ambientColorField->blockSize();

   psize = input.ambientColorField->blockSize();
   psize *= psize;

   ambientColor.nnx = ambientColor.nx / psize;
   ambientColor.nnx += (ambientColor.nnx*psize < ambientColor.nx) ? 1 : 0;

   ambientColor.nny = ambientColor.ny / psize;
   ambientColor.nny += (ambientColor.nny*psize < ambientColor.ny) ? 1 : 0;

   ambientColor.LLC.x = input.ambientColorField->llc()[0];
   ambientColor.LLC.y = input.ambientColorField->llc()[1];
   ambientColor.LLC.z = input.ambientColorField->llc()[2];
   ambientColor.LLC.w = 0.0f;
   ambientColor.Res.x = input.ambientColorField->dx();
   ambientColor.Res.y = input.ambientColorField->dy();
   ambientColor.Res.z = input.ambientColorField->dz();
   ambientColor.Res.w = 0.0f;
   ambientColor.defVal.x = input.ambientColorField->getDefVal()[0];
   ambientColor.defVal.y = input.ambientColorField->getDefVal()[1];
   ambientColor.defVal.z = input.ambientColorField->getDefVal()[2];
   ambientColor.defVal.w = input.ambientColorField->getDefVal()[3];

   size_t threads_per_block = 256;
   size_t grid_size = input.startPosition.size() / threads_per_block;
   if( grid_size * threads_per_block < input.startPosition.size() ) grid_size++;
   gilligan::util::ExecutionPolicy policy(threads_per_block, grid_size);

   int nbRays = Cd.size();

   // Launch CUDA kernel with ExecutionPolicy
   kernel_ssRayMarchAccumulation<<< policy.gridSize(), policy.blockSize() >>>(
      density, dsmList.getDevicePointer(), color, ambientDensity, ambientColor,

      cu_LightColor.getDevicePointer(), cu_LightPosition.getDevicePointer(), input.dsmField.size(),

      cu_StartPosition.getDevicePointer(), cu_StartDirection.getDevicePointer(), nbRays,
      input.scatterCoefficient, input.ds, input.maxPathlength, input.clampv, input.phaseF,

      cu_Cd.getDevicePointer()
   );

   cudaThreadSynchronize();
   cudaError_t res = cudaPeekAtLastError();
   if (res != cudaSuccess)
   {
      std::cout << "CUDA Error: " << cudaGetErrorString(res) << std::endl;
   }

   // Transfer back cu_Cd from device
   cu_Cd.updateHost();
   for(int i = 0; i < Cd.size(); ++i)
   {
      Cd[i][0] = (*cu_Cd)[i].x;
      Cd[i][1] = (*cu_Cd)[i].y;
      Cd[i][2] = (*cu_Cd)[i].z;
      Cd[i][3] = (*cu_Cd)[i].w;
   }
}


void cu_SetUniformPhaseFunction( cu_RenderData *d, float value )
{
   d->phaseF.type = cuPhaseFunction::PhaseFunctionType::UNIFORM;
   d->phaseF.val = value;
}

void cu_SetHenyeyGreensteinPhaseFunction( cu_RenderData *d, float g )
{
   d->phaseF.type = cuPhaseFunction::PhaseFunctionType::HENYEYGREENSTEIN;
   d->phaseF.g0 = g;
}

void cu_SetDoubleHenyeyGreensteinPhaseFunction( cu_RenderData *d, float g0, float g1, float mix )
{
   d->phaseF.type = cuPhaseFunction::PhaseFunctionType::DOUBLE_HENYEYGREENSTEIN;
   d->phaseF.g0 = g0;
   d->phaseF.g1 = g1;
   d->phaseF.mix = mix;
}

void cu_SetFournierForandPhaseFunction( cu_RenderData *d, float en, float mu )
{
   d->phaseF.type = cuPhaseFunction::PhaseFunctionType::FOURNIERFORAND;
   d->phaseF.en = en;
   d->phaseF.mu = mu;

   d->phaseF.nu = (3.0 - mu) / 2.0;
   d->phaseF.delta180 = 4.0 / ( 3.0 * (en-1.0)*(en-1.0) );
   double delta180nupower = std::pow( d->phaseF.delta180, d->phaseF.nu );
   d->phaseF.p1factor =  (1.0 - delta180nupower) / ( 16.0 * M_PI * delta180nupower * (d->phaseF.delta180-1.0) );
}

}
