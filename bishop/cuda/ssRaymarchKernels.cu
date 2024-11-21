#include <float.h>
#include <helper_math.h>
#include "cuda/cuFieldStructs.cuh"

#define CUDART_PI_F (3.141592654f)

inline __device__ float4 floor(float4 a){
   return ( (float4){floor(a.x), floor(a.y), floor(a.z), floor(a.w)} );
}

__device__ float4 unitvector(float4 v)
{
   float vmag = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
   return v/vmag;
}

__device__ bool isInside(int i, int j, int k, int nx, int ny, int nz )
{
   return (
      (i > -1) && (i < nx) &&
      (j > -1) && (j < ny) &&
      (k > -1) && (k < nz)
   );
}

__device__ float4 cameraTransform(float4 X, cu_DSM& dsm)
{
   X = X - dsm.camera_position;
   float z = dot(X, dsm.camera_view);
   if( z < 0 )
      return make_float4(0.0f);

   float x = dot(X, dsm.camera_right) / z;
   float y = dot(X, dsm.camera_up) / z;
   x = (x / dsm.htanfov + 1.0) / 2.0;
   y = (y / dsm.vtanfov + 1.0) / 2.0;
   z = (length(X) - dsm.near) / (dsm.far - dsm.near);
   return make_float4(x, y, z, 0.0f);
}

template <typename T>
__device__ T get(
      const T* data, const long* map, const long* mapOfMap,
      int i, int j, int k,
      int nx, int ny, int nz,
      int nnx, int nny, int partitionSize,
      T defVal
)
{
   if(!isInside( i, j, k, nx, ny, nz )){ return (defVal); }

   int ii = i / (partitionSize * partitionSize);
   int jj = j / (partitionSize * partitionSize);
   int kk = k / (partitionSize * partitionSize);
   int mapIndex = ii + nnx * (jj + nny * kk);
   if(mapOfMap[mapIndex] == -1)
   {
      return (defVal);
   }
   else
   {
      long offset = mapOfMap[mapIndex];
      ii = (i / partitionSize) % partitionSize;
      jj = (j / partitionSize) % partitionSize;
      kk = (k / partitionSize) % partitionSize;
      mapIndex = ii + partitionSize * (jj + partitionSize * kk);
      if(map[mapIndex + offset] == -1)
      {
         return (defVal);
      }
      else
      {
         offset = map[mapIndex + offset];
         ii = i % partitionSize;
         jj = j % partitionSize;
         kk = k % partitionSize;
         mapIndex = ii + partitionSize * (jj + partitionSize * kk);
         return data[mapIndex + offset];
      }
   }
}

__device__ float valueAtLinearInterpolation( float4 X, cu_ScalarField& field )
{
   X = X - field.LLC;

   // Get integer coordinates
   float4 a = make_float4( X.x/field.Res.x, X.y/field.Res.y, X.z/field.Res.z, 0.0f );

   if(
      a.x < 0 || a.x >= field.nx ||
      a.y < 0 || a.y >= field.ny ||
      a.z < 0 || a.z >= field.nz
   )
   { return field.defVal; }

   // Take the floor of each
   int4 b;
   b.x = (int)(floor(a.x));
   b.y = (int)(floor(a.y));
   b.z = (int)(floor(a.z));
   b.w = 0;

   // Get the weights
   float4 s = a - floor(a);
   float4 t = make_float4( 1.f - s.x, 1.f - s.y, 1.f - s.z, 0.0f );

   float psize = field.partitionSize;
   return
      get(field.data, field.map, field.mapOfMap, b.x,   b.y,   b.z,   field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (t.x * t.y * t.z) +
      get(field.data, field.map, field.mapOfMap, b.x,   b.y,   b.z+1, field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (t.x * t.y * s.z) +
      get(field.data, field.map, field.mapOfMap, b.x,   b.y+1, b.z,   field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (t.x * s.y * t.z) +
      get(field.data, field.map, field.mapOfMap, b.x,   b.y+1, b.z+1, field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (t.x * s.y * s.z) +
      get(field.data, field.map, field.mapOfMap, b.x+1, b.y,   b.z,   field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (s.x * t.y * t.z) +
      get(field.data, field.map, field.mapOfMap, b.x+1, b.y,   b.z+1, field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (s.x * t.y * s.z) +
      get(field.data, field.map, field.mapOfMap, b.x+1, b.y+1, b.z,   field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (s.x * s.y * t.z) +
      get(field.data, field.map, field.mapOfMap, b.x+1, b.y+1, b.z+1, field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (s.x * s.y * s.z);
}

__device__ float dsmValueAtLinearInterpolation( float4 X, cu_DSM& dsm )
{
   X = cameraTransform(X, dsm);

   if(
      X.x < 0 || X.x > 1.0 ||
      X.y < 0 || X.y > 1.0 ||
      X.z < 0 || X.z > 1.0
   )
   { return dsm.field.defVal; }

   // Get integer coordinates
   float4 a = make_float4( X.x * dsm.field.nx, X.y * dsm.field.ny, X.z * dsm.field.nz, 0.0f );

   // Take the floor of each
   int4 b;
   b.x = (int)(floor(a.x));
   b.y = (int)(floor(a.y));
   b.z = (int)(floor(a.z));
   b.w = 0;

   // Get the weights
   float4 s = a - floor(a);
   float4 t = make_float4( 1.f - s.x, 1.f - s.y, 1.f - s.z, 0.0f );

   float psize = dsm.field.partitionSize;
   return
      get(dsm.field.data, dsm.field.map, dsm.field.mapOfMap, b.x,   b.y,   b.z,   dsm.field.nx, dsm.field.ny, dsm.field.nz, dsm.field.nnx, dsm.field.nny, psize, dsm.field.defVal) * (t.x * t.y * t.z) +
      get(dsm.field.data, dsm.field.map, dsm.field.mapOfMap, b.x,   b.y,   b.z+1, dsm.field.nx, dsm.field.ny, dsm.field.nz, dsm.field.nnx, dsm.field.nny, psize, dsm.field.defVal) * (t.x * t.y * s.z) +
      get(dsm.field.data, dsm.field.map, dsm.field.mapOfMap, b.x,   b.y+1, b.z,   dsm.field.nx, dsm.field.ny, dsm.field.nz, dsm.field.nnx, dsm.field.nny, psize, dsm.field.defVal) * (t.x * s.y * t.z) +
      get(dsm.field.data, dsm.field.map, dsm.field.mapOfMap, b.x,   b.y+1, b.z+1, dsm.field.nx, dsm.field.ny, dsm.field.nz, dsm.field.nnx, dsm.field.nny, psize, dsm.field.defVal) * (t.x * s.y * s.z) +
      get(dsm.field.data, dsm.field.map, dsm.field.mapOfMap, b.x+1, b.y,   b.z,   dsm.field.nx, dsm.field.ny, dsm.field.nz, dsm.field.nnx, dsm.field.nny, psize, dsm.field.defVal) * (s.x * t.y * t.z) +
      get(dsm.field.data, dsm.field.map, dsm.field.mapOfMap, b.x+1, b.y,   b.z+1, dsm.field.nx, dsm.field.ny, dsm.field.nz, dsm.field.nnx, dsm.field.nny, psize, dsm.field.defVal) * (s.x * t.y * s.z) +
      get(dsm.field.data, dsm.field.map, dsm.field.mapOfMap, b.x+1, b.y+1, b.z,   dsm.field.nx, dsm.field.ny, dsm.field.nz, dsm.field.nnx, dsm.field.nny, psize, dsm.field.defVal) * (s.x * s.y * t.z) +
      get(dsm.field.data, dsm.field.map, dsm.field.mapOfMap, b.x+1, b.y+1, b.z+1, dsm.field.nx, dsm.field.ny, dsm.field.nz, dsm.field.nnx, dsm.field.nny, psize, dsm.field.defVal) * (s.x * s.y * s.z);
}

__device__ float4 valueAtLinearInterpolation4( float4 X, cu_ColorField& field )
{
   X = X - field.LLC;

   // Get integer coordinates
   float4 a = make_float4( X.x/field.Res.x, X.y/field.Res.y, X.z/field.Res.z, 0.0f );

   if(
      a.x < 0 || a.x >= field.nx ||
      a.y < 0 || a.y >= field.ny ||
      a.z < 0 || a.z >= field.nz
   )
   { return field.defVal; }

   // Take the floor of each
   int4 b;
   b.x = (int)(floor(a.x));
   b.y = (int)(floor(a.y));
   b.z = (int)(floor(a.z));
   b.w = 0;

   // Get the weights
   float4 s = a - floor(a);
   float4 t = make_float4( 1.f - s.x, 1.f - s.y, 1.f - s.z, 0.0f );

   float psize = field.partitionSize;
   return
      get(field.data, field.map, field.mapOfMap, b.x,   b.y,   b.z,   field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (t.x * t.y * t.z) +
      get(field.data, field.map, field.mapOfMap, b.x,   b.y,   b.z+1, field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (t.x * t.y * s.z) +
      get(field.data, field.map, field.mapOfMap, b.x,   b.y+1, b.z,   field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (t.x * s.y * t.z) +
      get(field.data, field.map, field.mapOfMap, b.x,   b.y+1, b.z+1, field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (t.x * s.y * s.z) +
      get(field.data, field.map, field.mapOfMap, b.x+1, b.y,   b.z,   field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (s.x * t.y * t.z) +
      get(field.data, field.map, field.mapOfMap, b.x+1, b.y,   b.z+1, field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (s.x * t.y * s.z) +
      get(field.data, field.map, field.mapOfMap, b.x+1, b.y+1, b.z,   field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (s.x * s.y * t.z) +
      get(field.data, field.map, field.mapOfMap, b.x+1, b.y+1, b.z+1, field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (s.x * s.y * s.z);
}

// Like valueAtLinearInterpolation4, but with periodic boundary conditions (maybe)
__device__ float4 valueAtLinearInterpolation4Periodic( float4 X, cu_ColorField& field )
{
   X = X - field.LLC;

   // Get integer coordinates
   float4 a = make_float4( X.x/field.Res.x, X.y/field.Res.y, X.z/field.Res.z, 0.0f );

   // Periodic boundary
   a.x += (a.x < 0.0f) ? field.nx : (a.x > field.nx) ? -field.nx : 0.0f;
   a.y += (a.y < 0.0f) ? field.ny : (a.y > field.ny) ? -field.ny : 0.0f;
   a.z += (a.z < 0.0f) ? field.nz : (a.z > field.nz) ? -field.nz : 0.0f;

   // Take the floor of each
   int4 b;
   b.x = (int)(floor(a.x));
   b.y = (int)(floor(a.y));
   b.z = (int)(floor(a.z));
   b.w = 0;

   // Get the weights
   float4 s = a - floor(a);
   float4 t = make_float4( 1.f - s.x, 1.f - s.y, 1.f - s.z, 0.0f );

   float psize = field.partitionSize;
   return
      get(field.data, field.map, field.mapOfMap, b.x,   b.y,   b.z,   field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (t.x * t.y * t.z) +
      get(field.data, field.map, field.mapOfMap, b.x,   b.y,   b.z+1, field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (t.x * t.y * s.z) +
      get(field.data, field.map, field.mapOfMap, b.x,   b.y+1, b.z,   field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (t.x * s.y * t.z) +
      get(field.data, field.map, field.mapOfMap, b.x,   b.y+1, b.z+1, field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (t.x * s.y * s.z) +
      get(field.data, field.map, field.mapOfMap, b.x+1, b.y,   b.z,   field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (s.x * t.y * t.z) +
      get(field.data, field.map, field.mapOfMap, b.x+1, b.y,   b.z+1, field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (s.x * t.y * s.z) +
      get(field.data, field.map, field.mapOfMap, b.x+1, b.y+1, b.z,   field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (s.x * s.y * t.z) +
      get(field.data, field.map, field.mapOfMap, b.x+1, b.y+1, b.z+1, field.nx, field.ny, field.nz, field.nnx, field.nny, psize, field.defVal) * (s.x * s.y * s.z);
}

__device__ float Theta( float4 nhat, float4 nhatp )
{
   return acos( dot(nhat, nhatp) );
}

__device__ float Ftau( float theta, float K, float tau, float taup )
{
   return( (theta + 2.f * CUDART_PI_F * K) * (taup / tau) );
}

__device__ float M( float theta, float K, float tau )
{
   return( (theta + (2.f * CUDART_PI_F * K)) / (sqrt(2.0f) * tau) );
}

__device__ float N( float mu )
{
   return( sqrt(CUDART_PI_F * 0.5f * mu) / (1.f - exp(-2.f / mu)) );
}

__device__ float Delmo( float mu, float M, float N, float R, float tau, float taup )
{
   float lhs = mu * N / R;
   float rhs = sin(M * taup) * sin( M * (tau - taup) );
   rhs /= ( M * sin(M * tau) );
   return ( lhs * rhs );
}

__device__ float Deloo( float mu, float N, float R, float tau, float taup )
{
   return ( (mu * N / R) * (taup * (tau - taup) / tau) );
}

__device__ float4 Bo( float4 nhatp, float4 rhat, float f )
{
   return( nhatp * cos(f) + rhat * sin(f) );
}

__device__ float4 DXDS( float4 nhat, float4 nhatp, float4 rhat, float tau, float taup, float theta, float K, float mu, float n, float R )
{
//   float theta = Theta(nhat, nhatp);
   float f = Ftau(theta, K, tau, taup);
   float4 bnaught = Bo(nhatp, rhat, f);
//   float n = N(mu);
   float m = M(theta, K, tau);

   float delmo = Delmo(mu, m, n, R, tau, taup);
   float deloo = Deloo(mu, n, R, tau, taup);

   return( bnaught * (1.f + (0.5f * (delmo + deloo))) );
}

__device__ float getR( const float* lookup, int size, float eta )
{
   int minv = 0;
   int maxv = size - 1;

   while( minv <= maxv )
   {
      int ndx = (minv + maxv) / 2;
      if( lookup[ndx] == eta ){
         float R = (float)ndx / (float)(size - 1);
         return ( R );
      }
      else if( lookup[ndx] > eta ){
         minv = ndx + 1;
      }
      else{
         maxv = ndx - 1;
      }
   }

   if( lookup[minv] == eta )
      return ( (float)minv / (float)(size - 1) );
   else if( lookup[maxv] == eta )
      return ( (float)maxv / (float)(size - 1) );
   else if( lookup[minv] > eta && lookup[maxv] > eta )
   {
      while( lookup[minv] > eta && lookup[maxv] > eta )
      {
         minv += 1;
         maxv += 1;
      }
   }
   else if( lookup[minv] < eta && lookup[maxv] < eta )
   {
      while( lookup[minv] < eta && lookup[maxv] < eta )
      {
         minv -= 1;
         maxv -= 1;
      }
   }

   float Rn = (float)minv / (float)(size - 1);
   float Rnm = (float)maxv / (float)(size - 1);
   float peta = (eta - lookup[minv]) / (lookup[maxv] - lookup[minv]);
   return ( peta * (Rn - Rnm) + Rnm );
}

__device__ float getC( const float* clookup, int wsize, float wdtau, float tau )
{
   int index = (int)(floor(tau / wdtau));
   return ( clookup[index] );
}

__device__ float evalPhaseFunction(cuPhaseFunction& phaseF, float theta)
{
   float ret = 0.0f;
   if( phaseF.type == cuPhaseFunction::PhaseFunctionType::HENYEYGREENSTEIN )
   {
      float costheta = cos(theta);
      float denom = 1.0 + phaseF.g0 * phaseF.g0 - 2.0 * phaseF.g0 * costheta;
      ret = (1.0 - phaseF.g0 * phaseF.g0)/(4.0*3.14159265);
      ret /= denom;
   } else if( phaseF.type == cuPhaseFunction::PhaseFunctionType::DOUBLE_HENYEYGREENSTEIN ) {
      float costheta = cos(theta);
      float denom = 1.0 + phaseF.g0 * phaseF.g0 - 2.0 * phaseF.g0 * costheta;
      float pf0 = (1.0 - phaseF.g0 * phaseF.g0)/(4.0*3.14159265);
      pf0 /= denom;

      denom = 1.0 + phaseF.g1 * phaseF.g1 - 2.0 * phaseF.g1 * costheta;
      float pf1 = (1.0 - phaseF.g1 * phaseF.g1)/(4.0*3.14159265);
      pf1 /= denom;

      ret = phaseF.mix * pf0 + (1.0 - phaseF.mix) * pf1;
   } else if( phaseF.type == cuPhaseFunction::PhaseFunctionType::FOURNIERFORAND ) {
      double s = sin( theta/2.0 );
      if( fabs(theta) < 0.00001 ){ s = sin(0.000001/2.0); }
      double s2 = s*s;
      double delta = phaseF.delta180*s2;
      double deltanupower = pow( (double)delta, (double)phaseF.nu );
      double c = cos(theta);
      if( fabs(theta) < 0.00001 ){ c = cos(0.000001); }
      double c2 = c*c;
      double p1 = phaseF.p1factor * ( 3.0*c2 - 1.0 );
      double p0 = phaseF.nu*(1.0-delta) - (1.0-deltanupower) + ( delta*(1.0-deltanupower) - phaseF.nu*(1.0-delta) )/s2;
      p0 /= 4.0*3.14159265*(1.0-delta)*(1.0-delta)*deltanupower;
      float value = p0 + p1;

      ret = value;
   } else {
      ret = phaseF.val;
   }

   return ret;
}

__device__ float4 gatherLight(float4 X, float4 D, float k, cu_DSM *dsmFields, int nbLights, cuPhaseFunction& phaseF, const float4 *lightColor, const float4 *lightPosition)
{
   float4 accum = make_float4(0.0f);
   for(int i = 0; i < nbLights; ++i)
   {
      float dsm = dsmValueAtLinearInterpolation(X, dsmFields[i]);
      float4 scatterDirection = lightPosition[i] - X;
      float cosd = dot(unitvector(scatterDirection), unitvector(D));
      float theta = acos(cosd);

      float pf = evalPhaseFunction(phaseF, theta);
      accum += pf * lightColor[i] * exp( -dsm * k );
   }
   return accum;
}

__global__ void kernel_ssRayMarchAccumulation(
   cu_ScalarField densityField, cu_DSM* dsmFields, cu_ColorField colorField,
   cu_ScalarField ambientDensityField, cu_ColorField ambientColorField,

   const float4* lightColor, const float4* lightPosition, int nbLights,

   const float4* startPosition, const float4* startDirection, int nbRays,
   float scatterCoefficient, float ds, float maxPathlength, float clampv,
   cuPhaseFunction phaseFunction,

   float4* Cd
)
{
   uint item = blockDim.x * blockIdx.x + threadIdx.x;
   if( item >= nbRays ) return;

   float4 accum = make_float4(0.0f);

   float T = 1;
   float s = 0;
   float density, ambDensity, totalDensity, dT;
   float4 color, amb, light;

   float4 X;
   float4 Y = startPosition[item];
   float4 D = startDirection[item];

   float4 boxdem;
   boxdem.x = densityField.Res.x * (float)(densityField.nx);
   boxdem.y = densityField.Res.y * (float)(densityField.ny);
   boxdem.z = densityField.Res.z * (float)(densityField.nz);
   boxdem.w = 0.f;  

   float4 boxmax = densityField.LLC + boxdem;

   float t1, t2, tnear, tfar;
   float4 hits;

   tnear = -FLT_MAX;
   tfar = FLT_MAX;

   /* AABB Test */

   // x-axis
   if( D.x != 0.f || (D.x == 0.f && Y.x > densityField.LLC.x && Y.x < boxmax.x) ){
      t1 = (densityField.LLC.x - Y.x) / D.x;
      t2 = (boxmax.x - Y.x) / D.x;
      
      if( t1 > t2 )
      {
         float temp = t1;
         t1 = t2;
         t2 = temp;
      }
      if( t1 > tnear )
         tnear = t1;
      if( t2 < tfar )
         tfar = t2;
      hits.x = tnear;
   }

   tnear = -FLT_MAX;
   tfar = FLT_MAX;

   // y-axis
   if( D.y != 0.f || (D.y == 0.f && Y.y > densityField.LLC.y && Y.y < boxmax.y) ){
      t1 = (densityField.LLC.y - Y.y) / D.y;
      t2 = (boxmax.y - Y.y) / D.y;
      
      if( t1 > t2 )
      {
         float temp = t1;
         t1 = t2;
         t2 = temp;
      }
      if( t1 > tnear )
         tnear = t1;
      if( t2 < tfar )
         tfar = t2;
      hits.y = tnear;
   }

   tnear = -FLT_MAX;
   tfar = FLT_MAX;

   // z-axis
   if( D.z != 0.f || (D.z == 0.f && Y.z > densityField.LLC.z && Y.z < boxmax.z) ){
      t1 = (densityField.LLC.z - Y.z) / D.z;
      t2 = (boxmax.z - Y.z) / D.z;
      if( t1 > t2 )
      {
         float temp = t1;
         t1 = t2;
         t2 = temp;
      }
      if( t1 > tnear )
         tnear = t1;
      if( t2 < tfar )
         tfar = t2;
      hits.z = tnear;
   }

   tnear = hits.x;
   if( tnear < hits.y )
      tnear = hits.y;
   if( tnear < hits.z )
      tnear = hits.z;

   // tnear = 0;
   // tfar = maxPathlength;
   X = Y + tnear * D;
   s = tnear;
   maxPathlength = tfar - tnear;
   while( s <= maxPathlength && T > 1.0e-6 )
   {
      density = valueAtLinearInterpolation(X, densityField);

      // P = X - ambientDensityField.LLC;
      ambDensity = valueAtLinearInterpolation(X, ambientDensityField);

      density = clamp(density / clampv, 0.0f, 1.0f);
      ambDensity = clamp(ambDensity / clampv, 0.0f, 1.0f);
      if( density + ambDensity > 0 )
      {
         totalDensity = density + ambDensity;

         // P = X - colorField.LLC;
         color = valueAtLinearInterpolation4(X, colorField);

         // P = X - ambientColorField.LLC;
         amb = valueAtLinearInterpolation4(X, ambientColorField);

         light = gatherLight(X, D, scatterCoefficient, dsmFields, nbLights, phaseFunction, lightColor, lightPosition);

         dT = exp( -scatterCoefficient * ds * totalDensity );
         color = (amb * ambDensity + color * light * density) / totalDensity;
         color = color * T * (1.0 - dT);

         accum += color;
         T *= dT;
      }
      s += ds;
      X += D * ds;

   }

   accum.w = 1.0 - T;
   Cd[item] = accum;
}

__global__ void kernel_ssWarpRayMarchAccumulation(
   cu_ScalarField densityField, cu_DSM* dsmFields, cu_ColorField warpField, cu_ColorField colorField,
   cu_ScalarField ambientDensityField, cu_ColorField ambientColorField,

   const float4* lightColor, const float4* lightPosition, int nbLights,

   const float4* startPosition, const float4* startDirection, int nbRays, int nbWarps,
   float scatterCoefficient, float ds, float maxPathlength, float clampv,
   cuPhaseFunction phaseFunction,

   float4* Cd
)
{
   uint item = blockDim.x * blockIdx.x + threadIdx.x;
   if( item >= nbRays ) return;
   {
      float4 accum = make_float4(0.0f);
      float T = 1;
      float s = 0;
      float4 X = startPosition[item];
      float4 D = startDirection[item];
      float density, ambDensity, totalDensity, dT;
      float4 color, amb, light, W, dW;
      int iw;
      while( s < maxPathlength && T > 1.0e-6 )
      {
         W = X - warpField.LLC;
         iw = 0;
         while( iw < nbWarps )
         {
            dW = valueAtLinearInterpolation4Periodic(W, warpField);
            W = W - dW;
            iw++;
         }
         W = W + warpField.LLC;
         // P = W - densityField.LLC;
         density = valueAtLinearInterpolation(W, densityField);

         // P = W - ambientDensityField.LLC;
         ambDensity = valueAtLinearInterpolation(W, ambientDensityField);

         density = clamp(density / clampv, 0.0f, 1.0f);
         ambDensity = clamp(ambDensity / clampv, 0.0f, 1.0f);
         if( density + ambDensity > 0.0f )
         {
            totalDensity = density + ambDensity;

            // P = W - colorField.LLC;
            color = valueAtLinearInterpolation4(W, colorField);
            
            // P = W - ambientColorField.LLC;
            amb = valueAtLinearInterpolation4(W, ambientColorField);

            light = gatherLight(W, D, scatterCoefficient, dsmFields, nbLights, phaseFunction, lightColor, lightPosition);

            dT = exp( -scatterCoefficient * ds * totalDensity );
            color = (amb * ambDensity + color * light * density) / totalDensity;
            color = color * T * (1.0f - dT);

            accum += color;
            T *= dT;
         }
         s += ds;
         X += D * ds;
      }
      accum.w = 1.0 - T;
      Cd[item] = accum;
   }
}


// DSM generator

/* Comment out for now; come back later after figuring out host-side
__global__ void ssRayMarchDSMAccumulation(
   float* densityField,
   int densityWidth, int densityHeight, int densityDepth,
   float4 densityLLC, float4 densityRes,
   float* dsmField,
   int dsmWidth, int dsmHeight, int dsmDepth,
   float4 dsmLLC, float4 dsmRes,
   float ds, float4 lightP
)
{
   int4 i;
   i.x = get_global_id(0);
   i.y = get_global_id(1);
   i.z = get_global_id(2);
   i.w = 0;

   float4 tempdist = dsmRes;
   tempdist.x *= dsmWidth;
   tempdist.y *= dsmHeight;
   tempdist.z *= dsmDepth;

   float accum = 0.f;
   float s = 0.f;

   float4 X;
   X.x = (i.x + 0.5f) * dsmRes.x + dsmLLC.x;
   X.y = (i.y + 0.5f) * dsmRes.y + dsmLLC.y;
   X.z = (i.z + 0.5f) * dsmRes.z + dsmLLC.z;
   X.w = 0.f;

   float4 D = lightP - X;
   float maxPathLength = length(D);
   D /= maxPathLength;

   while( s < maxPathLength ){
      float4 Y = X - densityLLC;
      float density = valueAtLinearInterpolation(
         densityField, Y.x, Y.y, Y.z,
         densityRes.x, densityRes.y, densityRes.z,
         densityWidth, densityHeight, densityDepth
      );
      if( density > 0.f && !isnan(density) ){
         accum += density * ds;
      }
      s += ds;
      X += (D * ds);
      Y = X - dsmLLC;

      if( Y.x < 0.f || Y.x > tempdist.x ) { s = maxPathLength; }
      if( Y.y < 0.f || Y.y > tempdist.y ) { s = maxPathLength; }
      if( Y.z < 0.f || Y.z > tempdist.z ) { s = maxPathLength; }
   }

   int index = i.x + dsmWidth * i.y + dsmWidth * dsmHeight * i.z;
   dsmField[index] = accum;
}

__global__ void ssRayMarchDSMWarpAccumulation(
   __read_only image3d_t densityField, float4 densityLLC, float4 densityRes,
   float* dsmField, int dsmWidth, int dsmHeight, int dsmDepth, float4 dsmLLC, float4 dsmRes,
   float ds, float4 lightP,
   __read_only image3d_t warpField, float4 warpLLC, float4 warpRes, int nbWarps
)
{

   int4 i;
   i.x = get_global_id(0);
   i.y = get_global_id(1);
   i.z = get_global_id(2);
   i.w = 0;

   float4 tempdist = dsmRes;
   tempdist.x *= dsmWidth;
   tempdist.y *= dsmHeight;
   tempdist.z *= dsmDepth;
   float maxDistance = length(tempdist);

   float accum = 0.0;
   float s = 0.0;

   float4 X;
   X.x = (i.x + 0.5f) * dsmRes.x + dsmLLC.x;
   X.y = (i.y + 0.5f) * dsmRes.y + dsmLLC.y;
   X.z = (i.z + 0.5f) * dsmRes.z + dsmLLC.z;
   X.w = 0.f;

   float4 D = lightP - X;
   float maxPathLength = length(D);
   D /= maxPathLength;

   while( s < maxPathLength ){
      float4 W;
      float4 Y = X - warpLLC;
      float4 SampY;
      for( int iw = 0; iw < nbWarps; iw++ ){
         W = read_imagef( warpField, nearest, SampY );
         Y -= W;
      }
      Y += warpLLC;
      float density = ( read_imagef(densityField, linear, Y) ).w;
      if( density > 0 && !isnan(density) ){
         accum += density * ds;
      }
      s += ds;
      X += (D * ds);
      Y = X - dsmLLC;
      if( Y.x < 0 || Y.x > tempdist.x ){ s = maxPathLength; }
      if( Y.y < 0 || Y.y > tempdist.y ){ s = maxPathLength; }
      if( Y.z < 0 || Y.z > tempdist.z ){ s = maxPathLength; }
   }

   int index = i.x + dsmWidth * i.y + dsmWidth * dsmHeight * i.z;

   dsmField[index] = accum;
}
*/
