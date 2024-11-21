#ifndef __CUDAFIELDSTRUCTS_H__
#define __CUDAFIELDSTRUCTS_H__

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

#ifndef __CUDACC__
typedef struct
{
   float x;
   float y;
   float z;
   float w;
} float4;
#endif

typedef struct
{
   const float* data;
   const long* map;
   const long* mapOfMap;
   int nx;
   int ny;
   int nz;
   int partitionSize;
   int nnx;
   int nny;
   float4 LLC;
   float4 Res;
   float defVal;
} cu_ScalarField;

typedef struct 
{
   float htanfov;
   float vtanfov;
   float near;
   float far;
   float4 camera_position;
   float4 camera_right;
   float4 camera_up;
   float4 camera_view;

   cu_ScalarField field;
} cu_DSM;

typedef struct
{
   const float4* data;
   const long* map;
   const long* mapOfMap;
   int nx;
   int ny;
   int nz;
   int partitionSize;
   int nnx;
   int nny;
   float4 LLC;
   float4 Res;
   float4 defVal;
} cu_ColorField;

#endif