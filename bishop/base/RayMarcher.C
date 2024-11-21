

//=======================================================================
//
//  Single Scatter straight line ray marcher.
//  Integrates along a single straight line path until exiting
//  the density or accumulating sufficient opacity.
//
//  Returns accumulated color and transmissivity as a float array
//
//=======================================================================

#include "RayMarcher.h"
#include "VShader.h"
#include "ProgressMeter.h"
#include "Tracking.h"
#include "UniformPRN.h"
#include <iostream>
#include <algorithm>
#include <omp.h>
#include <cmath>

using namespace std;
using namespace lux;
// Runs a ray march along a direction from a start position.  Samples density along the way, accumulating color/transmittance.
// Volume color is assumed to be white
// Returns color and opacity
//
// Basic accumulation is
//
//  dT = exp(-b*density*ds)
//  Cd += color * (1-dT) * T/b
//  T *= dT
//
//  Density is assumed to be axis aligned, with a lower left hand corner and a resolution.
//


void lux::AddDSM( RenderData* d, const ScalarField& dsm )
{
	 d->dsmField.push_back( dsm );
}

const ScalarField& lux::GetDSM( RenderData* d, int i )
{
	 return d->dsmField[i];
}

void lux::SetDensityField(  RenderData* d, const ScalarField& field )
{
		d->densityField = field;
}
void lux::SetAmbientDensityField(  RenderData* d, const ScalarField& field )
{
		d->ambientDensityField = field;
}
void lux::SetColorField(  RenderData* d, const ColorField& field )
{
		d->colorField = field;
}
void lux::SetAmbientColorField(  RenderData* d, const ColorField& field )
{
		d->ambientColorField = field;
}
void lux::SetSparseGrid(  RenderData* d, ScalarGrid& field )
{
		d->sparseGrid= field;
}
void lux::AddBoundingBox( RenderData *d, const Vector llc, const Vector urc )
{
   d->boundingBoxes.push_back( AABB(llc, urc) );
}
void lux::AddBoundingBox( RenderData *d, const AABB& aabb )
{
   d->boundingBoxes.push_back( aabb );
}
void lux::AddBoundingBoxes( RenderData *d, const std::vector<AABB>& boxes )
{
   d->boundingBoxes = boxes;
}
void lux::SetHoldOut(  RenderData* d, const ScalarField& field )
{
   d->holdOut = field;
   d->useHoldOut = true;
}
void lux::SetIntervalTree( RenderData *d, IntervalSet t )
{
   d->intervalTree = t;
   d->boundingBoxes = t->objects();
   cout << "\n\nNumber of bounding boxes: " << d->boundingBoxes.size() << endl << endl;
}

void lux::SetDSMRange( RenderData *d, const double value )
{
   if( value <= 0.0 )
   {
      d->use_dsm_range = false;
   }
   else
   {
      d->use_dsm_range = true;
      d->dsm_range = value;
   }
}


void lux::AddVolumeLight( RenderData* d, VolumeLightField& vlf ) { d->volumeLights.push_back(vlf); }
void lux::SetVolumeLightSamples( RenderData *d, int samples ) { d->vlf_samples = samples; }

// calculate DSM the long hard way
void lux::RayMarchDSMAccumulation( const ScalarField& densityField,
			           const Vector& lightPosition,
			           float ds,
				   ScalarGrid& dsmField
				   //VolumeGrid<float>& dsmField
				 )
{
    size_t dsmWidth = dsmField->nx();
    size_t dsmHeight = dsmField->ny();
    size_t dsmDepth = dsmField->nz();

    ProgressMeter meter( dsmHeight*dsmDepth, "DSM" );

    for( int idsmz=0;idsmz<dsmDepth;idsmz++)
    {
       for( int idsmy=0;idsmy<dsmHeight;idsmy++)
       {
          for( int idsmx=0;idsmx<dsmWidth;idsmx++)
          {
             float accum = 0;
             Vector X = dsmField->evalP( idsmx, idsmy, idsmz );
	     float density = densityField->eval( X );
	     bool doAccum = density > 0;
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx+1,idsmy, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx-1,idsmy, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx,idsmy+1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx,idsmy-1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx,idsmy, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx,idsmy, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx+1,idsmy+1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx+1,idsmy-1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx+1,idsmy, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx+1,idsmy, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx-1,idsmy+1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx-1,idsmy-1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx-1,idsmy, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx-1,idsmy, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx,idsmy+1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx,idsmy+1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx,idsmy-1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx,idsmy-1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx+1,idsmy+1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx+1,idsmy+1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx+1,idsmy-1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx+1,idsmy-1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx-1,idsmy+1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx-1,idsmy+1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx-1,idsmy-1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx-1,idsmy-1, idsmz-1 ) ) > 0; }

         if( doAccum )
             {
                Vector D = lightPosition-X;
                float maxPathlength = D.magnitude();
                D /= maxPathlength;
                float s = 0;
                while( s < maxPathlength )
                {
	               float density = densityField->eval( X );
                   if( density > 0 )
                   {
                      accum += density * ds;
                   }
                   s += ds;
                   X += D * ds;
	           if( !dsmField->isInside( X ) ) { s = maxPathlength; }
                }
             }
          dsmField->set(idsmx,idsmy,idsmz, accum);
          }
          meter.update();
       }
    }
}




void lux::RayMarchDSMAccumulation( const RenderData& input, 
                                   const Vector& lightPosition, 
                                   ScalarGrid& dsmField
                                 )
{
    size_t dsmWidth = dsmField->nx();
    size_t dsmHeight = dsmField->ny();
    size_t dsmDepth = dsmField->nz();


    //if( input.boundingBoxes.size() < 2 )
    //{
       ProgressMeter meter( dsmHeight*dsmDepth, "DSM Interval" );
    for( int idsmz=0;idsmz<dsmDepth;idsmz++)
    {
       for( int idsmy=0;idsmy<dsmHeight;idsmy++)
       {
#pragma omp parallel for
          for( int idsmx=0;idsmx<dsmWidth;idsmx++)
          {
             float accum = 0;
             Vector X = dsmField->evalP( idsmx, idsmy, idsmz );
             float density = input.densityField->eval( X );
	     bool doAccum = density > 0;
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy+1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy-1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy+1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy-1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy+1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy-1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy+1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy+1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy-1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy-1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy+1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy+1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy-1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy-1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy+1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy+1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy-1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy-1, idsmz-1 ) ) > 0; }

             if( doAccum )
             {
                
                Vector D = lightPosition-X;
                float maxPathlength = D.magnitude();
                bool go = true;
                if( input.use_dsm_range )
                {
                   if( maxPathlength >= input.dsm_range ){ go = false; }
                }
                if(go)
                {
                D /= maxPathlength;
                float s = 0;
                float ds = input.ds;
                while( s < maxPathlength )
                {
	           float density = input.densityField->eval( X );
                   if( density > 0 )
                   {
                      accum += density * ds;
                   }
                   s += ds;
                   X += D * ds;
	           if( !dsmField->isInside( X ) ) { s = maxPathlength; }
                }
                
                //accum = PointToPointLWP( input, X, lightPosition );
                if( input.use_dsm_range )
                {
                   accum = accum * std::pow( 1.0 - (maxPathlength/input.dsm_range), 2.0 );
                }
                dsmField->set(idsmx,idsmy,idsmz, accum);
                }
             }
          }
          meter.update();
       }
    }
/*
    }
    else
    {
       ProgressMeter meter( input.boundingBoxes.size(), "DSM AABB" );
       for( size_t ib=0;ib<input.boundingBoxes.size();ib++ )
       {
          int i0,j0,k0,i1,j1,k1;
          const AABB& aabb = input.boundingBoxes[ib];
          if( dsmField->getBox( aabb->llc(), aabb->urc(), i0, j0, k0, i1, j1, k1 ) )
          {
             for( int k=k0;k<k1;k++ )
             {
                for( int j=j0;j<j1;j++ )
                {
#pragma omp parallel for
                   for( int i=i0;i<i1;i++ )
                   {
                      float test = dsmField->get(i,j,k);
                      if( test == 0 )
                      {
                         Vector X = dsmField->evalP( i,j,k );
        
                float accum = 0.0;
                Vector D = lightPosition-X;
                float maxPathlength = D.magnitude();
                D /= maxPathlength;
                float s = 0;
                float ds = input.ds;
                while( s < maxPathlength )
                {
	           float density = input.densityField->eval( X );
                   if( density > 0 )
                   {
                      accum += density * ds;
                   }
                   s += ds;
                   X += D * ds;
	           if( !dsmField->isInside( X ) ) { s = maxPathlength; }
                }
                

                         //float accum = PointToPointLWP( input, X, lightPosition );
                         dsmField->set(i,j,k, accum);
                      }
                   }
                }
             }
          }
          meter.update(); 
       }
    }
*/
}




void lux::RayMarchDSMAccumulation( const RenderData& input, 
                                   const Vector& lightPosition, 
                                   ScalarFrustumGrid& dsmField
                                 )
{
    size_t dsmWidth = dsmField->nx();
    size_t dsmHeight = dsmField->ny();
    size_t dsmDepth = dsmField->nz();


    if( input.boundingBoxes.size() < 2 )
    {
       ProgressMeter meter( dsmHeight*dsmDepth, "DSM Interval" );
    for( int idsmz=0;idsmz<dsmDepth;idsmz++)
    {
       for( int idsmy=0;idsmy<dsmHeight;idsmy++)
       {
#pragma omp parallel for
          for( int idsmx=0;idsmx<dsmWidth;idsmx++)
          {
             float accum = 0;
             Vector X = dsmField->evalP( idsmx, idsmy, idsmz );
             float density = input.densityField->eval( X );
	     bool doAccum = density > 0;
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy+1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy-1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy+1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy-1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy+1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy-1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy+1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy+1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy-1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy-1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy+1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy+1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy-1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy-1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy+1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy+1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy-1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy-1, idsmz-1 ) ) > 0; }

             if( doAccum )
             {
               /* 
                Vector D = lightPosition-X;
                float maxPathlength = D.magnitude();
                D /= maxPathlength;
                float s = 0;
                float ds = input.ds;
                while( s < maxPathlength )
                {
	           float density = input.densityField->eval( X );
                   if( density > 0 )
                   {
                      accum += density * ds;
                   }
                   s += ds;
                   X += D * ds;
	           if( !dsmField->isInside( X ) ) { s = maxPathlength; }
                }
                */
                accum = PointToPointLWP( input, X, lightPosition );
                dsmField->set(idsmx,idsmy,idsmz, accum);
             }
          }
          meter.update();
       }
    }
    }
    else
    {
       ProgressMeter meter( input.boundingBoxes.size(), "DSM AABB" );
       for( size_t ib=0;ib<input.boundingBoxes.size();ib++ )
       {
          int i0,j0,k0,i1,j1,k1;
          const AABB& aabb = input.boundingBoxes[ib];
          if( dsmField->getBox( aabb->llc(), aabb->urc(), i0, j0, k0, i1, j1, k1 ) )
          {
             for( int k=k0;k<k1;k++ )
             {
                for( int j=j0;j<j1;j++ )
                {
#pragma omp parallel for
                   for( int i=i0;i<i1;i++ )
                   {
                      float test = dsmField->get(i,j,k);
                      if( test == 0 )
                      {
                         Vector X = dsmField->evalP( i,j,k );
        
                float accum = 0.0;
                Vector D = lightPosition-X;
                float maxPathlength = D.magnitude();
                D /= maxPathlength;
                float s = 0;
                float ds = input.ds;
                while( s < maxPathlength )
                {
	           float density = input.densityField->eval( X );
                   if( density > 0 )
                   {
                      accum += density * ds;
                   }
                   s += ds;
                   X += D * ds;
	           if( !dsmField->isInside( X ) ) { s = maxPathlength; }
                }
                

                         //float accum = PointToPointLWP( input, X, lightPosition );
                         dsmField->set(i,j,k, accum);
                      }
                   }
                }
             }
          }
          meter.update(); 
       }
    }
}






void lux::RayMarchVisibleDSMAccumulation( const RenderData& input, 
                                          const Vector& lightPosition, 
                                          const Camera& cam,
                                          ScalarGrid& dsmField
                                        )
{
    size_t dsmWidth = dsmField->nx();
    size_t dsmHeight = dsmField->ny();
    size_t dsmDepth = dsmField->nz();


    //if( input.boundingBoxes.size() < 2 )
    //{
       ProgressMeter meter( dsmHeight*dsmDepth, "DSM Interval" );
    for( int idsmz=0;idsmz<dsmDepth;idsmz++)
    {
       for( int idsmy=0;idsmy<dsmHeight;idsmy++)
       {
#pragma omp parallel for
          for( int idsmx=0;idsmx<dsmWidth;idsmx++)
          {
             float accum = 0;
             Vector X = dsmField->evalP( idsmx, idsmy, idsmz );
             if( cam.isVisible(X) )
             {
             float density = input.densityField->eval( X );
	     bool doAccum = density > 0;
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy+1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy-1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy+1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy-1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy+1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy-1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy+1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy+1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy-1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx,idsmy-1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy+1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy+1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy-1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx+1,idsmy-1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy+1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy+1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy-1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = input.densityField->eval( dsmField->evalP(idsmx-1,idsmy-1, idsmz-1 ) ) > 0; }

             if( doAccum )
             {
                
                Vector D = lightPosition-X;
                float maxPathlength = D.magnitude();
                bool go = true;
                if( input.use_dsm_range )
                {
                   if( maxPathlength >= input.dsm_range ){ go = false; }
                }
                if(go)
                {
                D /= maxPathlength;
                float s = 0;
                float ds = input.ds;
                while( s < maxPathlength )
                {
	           float density = input.densityField->eval( X );
                   if( density > 0 )
                   {
                      accum += density * ds;
                   }
                   s += ds;
                   X += D * ds;
	           if( !dsmField->isInside( X ) ) { s = maxPathlength; }
                }
                
                //accum = PointToPointLWP( input, X, lightPosition );
                dsmField->set(idsmx,idsmy,idsmz, accum);
                }
             }
             }
          }
          meter.update();
       }
    }
/*
    }
    else
    {
       ProgressMeter meter( input.boundingBoxes.size(), "DSM AABB" );
       for( size_t ib=0;ib<input.boundingBoxes.size();ib++ )
       {
          int i0,j0,k0,i1,j1,k1;
          const AABB& aabb = input.boundingBoxes[ib];
          if( dsmField->getBox( aabb->llc(), aabb->urc(), i0, j0, k0, i1, j1, k1 ) )
          {
             for( int k=k0;k<k1;k++ )
             {
                for( int j=j0;j<j1;j++ )
                {
#pragma omp parallel for
                   for( int i=i0;i<i1;i++ )
                   {
                      float test = dsmField->get(i,j,k);
                      if( test == 0 )
                      {
                         Vector X = dsmField->evalP( i,j,k );
        
                float accum = 0.0;
                Vector D = lightPosition-X;
                float maxPathlength = D.magnitude();
                D /= maxPathlength;
                float s = 0;
                float ds = input.ds;
                while( s < maxPathlength )
                {
	           float density = input.densityField->eval( X );
                   if( density > 0 )
                   {
                      accum += density * ds;
                   }
                   s += ds;
                   X += D * ds;
	           if( !dsmField->isInside( X ) ) { s = maxPathlength; }
                }
                

                         //float accum = PointToPointLWP( input, X, lightPosition );
                         dsmField->set(i,j,k, accum);
                      }
                   }
                }
             }
          }
          meter.update(); 
       }
    }
*/
}















// calculate DSM the long hard way, including a holdout object
void lux::RayMarchDSMAccumulation( const ScalarField& densityField, 
                                   const Vector& lightPosition, 
                                   float ds,
                                   ScalarGrid& dsmField,
                                   ScalarField& holdout
                                   //VolumeGrid<float>& dsmField
                                 )
{
    size_t dsmWidth = dsmField->nx();
    size_t dsmHeight = dsmField->ny();
    size_t dsmDepth = dsmField->nz();

    ProgressMeter meter( dsmHeight*dsmDepth, "DSM" );

    for( int idsmz=0;idsmz<dsmDepth;idsmz++)
    {
       for( int idsmy=0;idsmy<dsmHeight;idsmy++)
       {
          for( int idsmx=0;idsmx<dsmWidth;idsmx++)
          {
             float accum = 0;
             Vector X = dsmField->evalP( idsmx, idsmy, idsmz );
	     float density = densityField->eval( X );
	     bool doAccum = density > 0;
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx+1,idsmy, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx-1,idsmy, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx,idsmy+1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx,idsmy-1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx,idsmy, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx,idsmy, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx+1,idsmy+1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx+1,idsmy-1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx+1,idsmy, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx+1,idsmy, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx-1,idsmy+1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx-1,idsmy-1, idsmz ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx-1,idsmy, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx-1,idsmy, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx,idsmy+1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx,idsmy+1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx,idsmy-1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx,idsmy-1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx+1,idsmy+1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx+1,idsmy+1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx+1,idsmy-1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx+1,idsmy-1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx-1,idsmy+1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx-1,idsmy+1, idsmz-1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx-1,idsmy-1, idsmz+1 ) ) > 0; }
	     if (!doAccum) { doAccum = densityField->eval( dsmField->evalP(idsmx-1,idsmy-1, idsmz-1 ) ) > 0; }

             if( doAccum )
             {
                Vector D = lightPosition-X;
                float maxPathlength = D.magnitude();
                D /= maxPathlength;
                float s = 0;
                while( s < maxPathlength )
                {
                   if( holdout->eval(X) >= 0.0 )
                   {
                      s = maxPathlength;
                      accum = 1.0e10;
                      break;
                   }
	               float density = densityField->eval( X );
                   if( density > 0 )
                   {
                      accum += density * ds;
                   }
                   s += ds;
                   X += D * ds;
	               if( !dsmField->isInside( X ) ) { s = maxPathlength; }
                }
             }
             dsmField->set(idsmx,idsmy,idsmz, accum);
          }
          meter.update();
       }
    }
}





void lux::RayMarchDSMAccumulation( const ScalarField& densityField,
																	 ScalarFrustumGrid& dsmField,
					 int nbsamples
																 )
{
		size_t dsmWidth = dsmField->nx();
		size_t dsmHeight = dsmField->ny();
		size_t dsmDepth = dsmField->nz();

		ProgressMeter meter( dsmHeight, "Frustum DSM" );
		UniformPRN rng;
		Noise_t p;
		p.seed = 47457459;
		rng.setParameters(p);

		for( int idsmy=0;idsmy<dsmHeight;idsmy++)
		{
			 for( int idsmx=0;idsmx<dsmWidth;idsmx++)
			 {
					float accum = 0;
					Vector X0 = dsmField->evalP( idsmx, idsmy, 0 );
					for( int idsmz=0;idsmz<dsmDepth;idsmz++)
					{
						 Vector X = dsmField->evalP( idsmx, idsmy, idsmz );
						float density = densityField->eval( X );
						float ds = (X-X0).magnitude();
						if((nbsamples > 1) &&
							 (idsmz > 0) &&
							(idsmz < dsmDepth-1) &&
							(idsmx > 0) &&
							(idsmx < dsmWidth-1) &&
							(idsmy > 0) &&
							(idsmy < dsmHeight-1)   )
						{
							Vector dX = (dsmField->evalP( idsmx+1, idsmy, idsmz ) - dsmField->evalP( idsmx-1, idsmy, idsmz ));
							Vector dY = (dsmField->evalP( idsmx, idsmy+1, idsmz ) - dsmField->evalP( idsmx, idsmy-1, idsmz ));
							Vector dZ = (dsmField->evalP( idsmx, idsmy, idsmz+1 ) - dsmField->evalP( idsmx, idsmy, idsmz-1 ));
							for( int s=1;s<nbsamples;s++ )
							{
								 Vector P = X + dX * (rng.eval()-0.5) + dY * (rng.eval()-0.5) + dZ * (rng.eval()-0.5);
								 density += densityField->eval(P);
							}
							density /= nbsamples;
						}
						accum += density * ds;
						dsmField->set(idsmx,idsmy,idsmz, accum);
						X0 = X;
					}
			 }
			 meter.update();
		}
}











void lux::ssRayMarchAccumulation( const RenderData& input, vector<Color>& Cd, vector<Color>& Oi )
{
    Cd.clear();
    Cd.resize( input.startPosition.size() );

    Oi.clear();
    Oi.resize( input.startPosition.size() );

    Color onec(1.0,1.0,1.0,1.0);
    Color zeroc(0.0,0.0,0.0,0.0);
   //cout << "Nb IntervalTree objects : " << input.intervalTree->nbObjects() << endl;
   if( input.intervalTree->nbObjects() > 0 )
   {
       double maxpathlength = input.maxPathlength;
#pragma omp parallel for
      for( int item=0;item<(int)input.startPosition.size();item++)
      {
         Color accum, CC;
         Color ci, amb, light;
         Color T = onec;
         Color dT = onec;
         float density = 0.0;
         float ambDensity = 0.0;
         Vector X = input.startPosition[item];
         Vector D = input.startDirection[item];
         double s = 0.0;

         IntervalData id = input.intervalTree->interval_intersect( X, D );
         while( id.status )
         {
            //cout << item << "   tmin, tmax: " << id.tmin << " " << id.tmax << endl;
            X += D * id.tmin;
            s += id.tmin;
            if( s >= maxpathlength || (T[0] < 1.0e-6 && T[1] < 1.0e-6 && T[2] < 1.0e-6 && T[3] < 1.0e-6) ){ break; }
            double sinterval = id.tmin;
            int nbsteps = (int)((id.tmax - id.tmin)/input.ds) + 1;
            double ds = (id.tmax-id.tmin)/nbsteps;
            while( sinterval < id.tmax )
            {
               if( input.useHoldOut )
               {
                  if( input.holdOut->eval(X) >= 0.0 )
                  {
                     s = maxpathlength;
                     T = zeroc;
                     break;
                  }
               }

               density = input.densityField->eval( X );
	       ambDensity = input.ambientDensityField->eval( X );
	       if( density < 0 ){ density = 0; }
	       if( ambDensity < 0 ){ ambDensity = 0; }
	       if( density+ambDensity > 0 )
	       {
                  float totaldensity = density+ambDensity;
                  dT = exp( -input.scatterCoefficient * ds * totaldensity );
	          ci = input.colorField->eval( X );
                  amb = input.ambientColorField->eval(X);
	          light = gatherLight( input, X, D );
	          CC =  ( amb*ambDensity + ci * light*density )/totaldensity;
                  CC = CC * T * (onec-dT);
		  accum += CC;
                  T *= dT;
	       }
               X += D * ds;
               sinterval += ds;
            }
            id = input.intervalTree->interval_intersect( X, D );
         }
         Cd[item] = accum;
         Oi[item] = onec - T;
      }
   }
   else if( input.sparseGrid.get() == 0 )
   {
    // Consider OMP here
#pragma omp parallel for
    for( int item=0;item<(int)input.startPosition.size();item++)
    {
       Color accum, CC;
       Color ci, amb, light;
       Color T = onec;
       Color dT = onec;
       float ds = input.ds;
       float density = 0;
       float ambDensity = 0;
       Vector X = input.startPosition[item];
       Vector D = input.startDirection[item];
       Vector Y = X;

       double totals = 0;
       double maxpathlength = input.maxPathlength;

	if(input.boundingBoxes.size()>0)
	{
	   double closeMin = input.maxPathlength;
	   double farMax = 0;
	   bool foundBox = false;
	   for( size_t b=0;b<input.boundingBoxes.size();b++)
	   {
	      double close, far;
	      if( input.boundingBoxes[b]->intersection( X, D, close, far ) )
	      {
	         closeMin = (closeMin > close ) ? close : closeMin;
		 farMax = (farMax < far ) ? far : farMax;
		 foundBox = true;
	      }
	   }
	   if( foundBox )
	   {
	      X = X + D*closeMin;
	      maxpathlength = farMax - closeMin;
	   }
	   else
	   {
	      totals = maxpathlength + 1.0;
	   }

		}

    while( totals <= maxpathlength && (T[0] > 1.0e-6 || T[1] > 1.0e-6 || T[2] > 1.0e-6 || T[3] > 1.0e-6) )
    {
        if( input.useHoldOut )
        {
           if( input.holdOut->eval(X) >= 0.0 )
           {
              totals = maxpathlength;
              T = zeroc;
              break;
           }
        }
        density = input.densityField->eval( X );
		ambDensity = input.ambientDensityField->eval( X );
	    if( density < 0 ){ density = 0; }
	    if( ambDensity < 0 ){ ambDensity = 0; }
		if( density+ambDensity > 0 )
		{
           float totaldensity = density+ambDensity;
           dT = exp( -input.scatterCoefficient * ds * totaldensity );
	       ci = input.colorField->eval( X );
           amb = input.ambientColorField->eval(X);
	       light = gatherLight( input, X, D );
	       CC =  ( amb*ambDensity + ci * light*density )/totaldensity;
           CC = CC * T * (onec-dT);
		   accum += CC;
           T *= dT;
		}
	    totals += ds;
        X += D * ds;
	}
    Cd[item] = accum;
    Oi[item] = onec - T;
    }
   }
   else
   {
#pragma omp parallel for
    for( int item=0;item<(int)input.startPosition.size();item++)
    {
       Color accum, CC;
       Color ci, amb, light;
       Color T = onec;
       Color dT = onec;
       float ds = input.ds;
       float density = 0;
       float ambDensity = 0;
       Vector X = input.startPosition[item];
       Vector D = input.startDirection[item];
       Vector Y = X;
       vector<TraceSegment> sortedBoxes;
	   findRayMarchBoxes( input.sparseGrid, X, D, sortedBoxes, input.sparseGrid->dx() );
       for( size_t b=0;b<sortedBoxes.size();b++ )
       {
           X = sortedBoxes[b].pos();
           double s = 0;
	       double pathLength = sortedBoxes[b].thickness();
           while( s <= pathLength && (T[0] > 1.0e-6 || T[1] > 1.0e-6 || T[2] > 1.0e-6 || T[3] > 1.0e-6) )
           {
              if( input.useHoldOut )
              {
                 if( input.holdOut->eval(X) >= 0.0 )
                 {
                    s = pathLength;
                    T = zeroc;
                    break;
                 }
              }

              density = input.densityField->eval( X );
		      ambDensity = input.ambientDensityField->eval( X );
	          if( density < 0 ){ density = 0; }
	          if( ambDensity < 0 ){ ambDensity = 0; }
		      if( density+ambDensity > 0 )
		      {
                 float totaldensity = density+ambDensity;
                 dT = exp( -input.scatterCoefficient * ds * totaldensity );
	             ci = input.colorField->eval( X );
                 amb = input.ambientColorField->eval(X);
	             light = gatherLight( input, X, D );
	             CC =  ( amb*ambDensity + ci * light*density )/totaldensity;
                 CC = CC * T * (onec-dT);
		         accum += CC;
                 T *= dT;
		      }
	          s += ds;
              X += D * ds;
	       }
          }
	      sortedBoxes.clear();
          Cd[item] = accum;
          Oi[item] = onec - T;
       }
    }
}




void lux::ssPointToPointRayMarchAccumulation( const RenderData& input, const Vector& start, const Vector& end, Color& Cd, Color& Oi )
{
    Color onec(1.0,1.0,1.0,1.0);
    Color zeroc(0.0,0.0,0.0,0.0);
    Cd = zeroc;

         double maxpathlength = (start-end).magnitude();
         Color accum, CC;
         Color ci, amb, light;
         Color T = onec;
         Color dT = onec;
         float density = 0.0;
         float ambDensity = 0.0;
         Vector X = start;
         Vector D = (end-start).unitvector();
         double s = 0.0;

         int nbsteps = (int)(maxpathlength/input.ds) + 1;
         double ds = maxpathlength/nbsteps;



         while( s<=maxpathlength && s<=input.maxPathlength && (T[0] > 1.0e-6 || T[1] > 1.0e-6 || T[2] > 1.0e-6 || T[3] > 1.0e-6) )
         {
            density = input.densityField->eval( X );
	        ambDensity = input.ambientDensityField->eval( X );
	        if( density < 0 ){ density = 0; }
	        if( ambDensity < 0 ){ ambDensity = 0; }
	        if( density+ambDensity > 0 )
	        {
               float totaldensity = density+ambDensity;
               dT = exp( -input.scatterCoefficient * ds * totaldensity );
	           ci = input.colorField->eval( X );
               amb = input.ambientColorField->eval(X);
	           light = gatherLight( input, X, D );
	           CC =  ( amb*ambDensity + ci * light*density )/totaldensity;
               CC = CC * T * (onec-dT);
		       accum += CC;
               T *= dT;
	        }
            X += D * ds;
            s += ds;
         }
         Cd = accum;
         Oi = onec - T;



/*
   if( input.intervalTree->nbObjects() > 0 )
   {
         IntervalData id = input.intervalTree->interval_intersect( X, D );
         while( id.status )
         {
            X += D * id.tmin;
            s += id.tmin;
            if( s >= maxpathlength || T < 1.0e-6 ){ break; }
            double sinterval = id.tmin;
            int nbsteps = (int)((id.tmax - id.tmin)/input.ds) + 1;
            double ds = (id.tmax-id.tmin)/nbsteps;
            while( sinterval < id.tmax )
            {
               if( input.useHoldOut )
               {
                  if( input.holdOut->eval(X) >= 0.0 )
                  {
                     s = maxpathlength;
                     T = 0.0;
                     break;
                  }
               }
               density = input.densityField->eval( X );
	           ambDensity = input.ambientDensityField->eval( X );
	           if( density < 0 ){ density = 0; }
	           if( ambDensity < 0 ){ ambDensity = 0; }
	           if( density+ambDensity > 0 )
	           {
                  float totaldensity = density+ambDensity;
                  dT = exp( -input.scatterCoefficient * ds * totaldensity );
	              ci = input.colorField->eval( X );
                  amb = input.ambientColorField->eval(X);
	              light = gatherLight( input, X, D );
	              CC =  ( amb*ambDensity + ci * light*density )/totaldensity;
                  CC = CC * T * (1.0-dT);
		          accum += CC;
                  T *= dT;
	           }
               X += D * ds;
               sinterval += ds;
            }
            id = input.intervalTree->interval_intersect( X, D );
         }
         accum[3] = 1.0 - T;
         Cd = accum;
   }
*/
}






//void renderPointLight( const RenderData& input, vector<Color>& output, const float psf_width )
















const Color lux::gatherLight( const RenderData& input, const Vector P, const Vector raymarchDirection )
{
   Color accum;
   Vector raymarchUnitDirection = raymarchDirection.unitvector();
   if( input.dsmField.size() > 0 )
   {
      for( size_t i=0;i<input.dsmField.size();i++ )
      {
         float dsm = input.dsmField[i]->eval( P );
         Vector scatterDirection = input.lightPosition[i] - P;
         double range_factor = 1.0;
         if( input.use_dsm_range )
         {
	    double range = scatterDirection.magnitude();
            if( range >= input.dsm_range )
            {
               range_factor = 0.0;
            }
            else
            {
               range_factor = std::pow( 1.0 - (range/input.dsm_range), 2.0 );
            }
         }
         float cosinedirection = scatterDirection.unitvector() * raymarchUnitDirection;
         float theta = acos( cosinedirection );
         float phasefunction = input.phaseFunction->eval(theta);
         accum += range_factor * phasefunction * input.lightColor[i] * exp( - dsm * input.scatterCoefficient );
      }
   }
   else
   {
      for( size_t i=0;i<input.lightPosition.size();i++ )
      {
         Color T = Color(1.0,1.0,1.0,1.0);
         Vector scatterDirection = input.lightPosition[i] - P;
         double range_factor = 1.0;
         if( input.use_dsm_range )
         {
	    double range = scatterDirection.magnitude();
            if( range >= input.dsm_range )
            {
               range_factor = 0.0;
            }
            else
            {
               range_factor = std::pow( 1.0 - (range/input.dsm_range), 2.0 );
               T = PointToPointTransmissivity( input, P, input.lightPosition[i] );
            }
         }
         else
         {
            T = PointToPointTransmissivity( input, P, input.lightPosition[i] );
         }

         float cosinedirection = scatterDirection.unitvector() * raymarchUnitDirection;
         float theta = acos( cosinedirection );
         float phasefunction = input.phaseFunction->eval(theta);
         accum += phasefunction * input.lightColor[i] * T;
      }
   }


   if( !input.volumeLights.empty() )
   {
      for( size_t il=0;il<input.volumeLights.size(); il++ )
      {
         int count = 0;
         Color volume_accum;
         for( int is=0;is<input.vlf_samples;is++ )
         {
            Vector X;
            Color cd;
            if( input.volumeLights[il]->getLight( X, cd ) )
            {
                Color T = PointToPointTransmissivity( input, P, X );
                Vector scatterDirection = X-P;
                float cosinedirection = scatterDirection.unitvector() * raymarchUnitDirection;
                float theta = acos( cosinedirection );
                float phasefunction = input.phaseFunction->eval(theta);
                volume_accum += phasefunction * cd * T;
                count++;

                //cout << "Light volume hit is, il, distance, T, accum " << is << " " << il << " " << " " << (P-X).magnitude() << "   " << T << "  " << cd.red() << " " << cd.green() << " " << cd.blue() << endl;
            }
         }
         if( count > 0 )
         {
            accum += volume_accum/count;
         } 
      }
   }




   return accum;
}


void lux::setNbCores( int nb )
{
   int nbmax = omp_get_max_threads();
   if( nb <= nbmax ){ omp_set_num_threads( nb ); }
   else { omp_set_num_threads( nbmax ); }
}


void lux::SetUniformPhaseFunction( RenderData*d, float value )
{
	 d->phaseFunction = PhaseFunction( new UniformPhaseFunction( value ) );
}

void lux::SetHenyeyGreensteinPhaseFunction( RenderData*d, float g )
{
	 d->phaseFunction = PhaseFunction( new HenyeyGreensteinPhaseFunction( g ) );
}

void lux::SetDoubleHenyeyGreensteinPhaseFunction( RenderData*d, float g0, float g1, float mix )
{
	 d->phaseFunction = PhaseFunction( new DoubleHenyeyGreensteinPhaseFunction( g0, g1, mix ) );
}

void lux::SetFournierForandPhaseFunction( RenderData*d, float en, float mu )
{
	 d->phaseFunction = PhaseFunction( new FournierForandPhaseFunction( en, mu ) );
}



/*
void lux::EndPointsRayMarch( const Vector& start, const Vector& end, const
ColorField& attenuation, const ScalarField& integrand, const float nominal_ds,  double& result, Color& extinction )
{
// Switched to an algorithm that accounts for variability of attenuation within
// a cell, essentially behaving as approaching an integrable singularity. 


   result = 0.0;
   extinction = 1.0;

   Vector direction = end-start;
   double distance = direction.magnitude();
   if( distance == 0.0 ){ return; }

   direction /= distance;

   Vector X = start;
   Vector X0 = X;
   double steps = distance / nominal_ds;
   int nbsteps = (int)steps;
   if( steps - (double)nbsteps != 0 ){ nbsteps++; }
   double ds = distance / nbsteps;
   double kappa0 = attenuation->eval(X);
   double inter0 = integrand->eval(X);

   for( int i=0;i<nbsteps;i++ )
   {
      X += direction * ds;
      // Evaluate integrand at midpoint
      double inter = integrand->eval(X);
      double kappa = attenuation->eval(X);
      if( kappa==kappa0 )
      {
         if( kappa != 0.0 )
         {
            Color deltaT = exp( - kappa * ds );
            result = result * deltaT + (inter+inter0)*0.5*(1.0 - deltaT)/kappa;
            extinction *= deltaT;
         }
         else
         {
            result += (inter+inter0)*0.5*ds;
         }
      }
      else
      {
         double ratio = kappa0/kappa;
         if( !std::isnan(ratio) && ratio > 1.0e-2 && ratio < 1.0e2 )
         {
            double nX0 = direction*X0;
            double X0mag = X0.magnitude();
            double radical = nX0*nX0 - (X0*X0)*(1.0 - ratio*ratio);
            radical = sqrt(radical);
            if(nX0 < 0.0 ){ radical = -radical; }
            double D = (-nX0 + radical)/ds;
            double pwr = kappa0*X0mag/D;
            double Dtest = D*D*( X0mag*X0mag - nX0*nX0 );
            if( Dtest < 1.0e-4 )
            {
               double A = ( D*ds + nX0 )/nX0;
               double deltaT = pow(A,-pwr);
               //double interfactor = ds / (-log(deltaT));
               result = result * deltaT + (inter+inter0)*0.5*ds;
               //result = result * deltaT + (inter+inter0)*0.5*interfactor*(1.0-deltaT);
               extinction *= deltaT;
                           }
            else
            {
               double A = ( (X0 + direction*D*ds).magnitude() + D*ds + nX0 )/( X0mag + nX0 );
               double deltaT = pow(A,-pwr);
               //double interfactor = ds / (-log(deltaT));
               if(std::isnan(deltaT))
               { 
                  cout << "deltaT nan\n"; 
                  cout << "\tA " << A << endl; 
                  cout << "\tpwr " << pwr << endl; 
                  cout << "\tradical " << radical << endl; 
                  cout << "\tnX0 " << nX0 << endl; 
                  cout << "\tX0mag " << X0mag << endl; 
                  cout << "\tratio " << ratio << endl; 
                  cout << "\tkappa0 " << kappa0 << endl; 
                  cout << "\tD " << D << endl; 
                  cout << "\tDtest " << Dtest << endl; 
               }
               result = result * deltaT + (inter+inter0)*0.5*ds;
               //result = result * deltaT + (inter+inter0)*0.5*interfactor*(1.0-deltaT);
               extinction *= deltaT;
            }
         }
         else
         {
            cout << "EndPointsRayMarch TROUBLE: ratio ill conditioned:  ratio = " << ratio << endl;
         }
      }
      kappa0 = kappa;
      inter0 = inter;
      X0 = X;
   }

}
*/


const Color lux::PointToPointTransmissivity( const RenderData& input, const Vector& start, const Vector& end )
{
   Color onec(1,1,1,1);
   Color zeroc(0,0,0,0);
   if( input.intervalTree->nbObjects() > 0 )
   {
       double maxpathlength = (end-start).magnitude();
       Color T = onec;
       Color dT = onec;
       float ds = input.ds;
       int nbsteps = (int)( maxpathlength / ds ) + 1;
       ds = maxpathlength/nbsteps;
       float density = 0.0;
       float ambDensity = 0.0;
       Vector X = start;
       Vector D = (end-start).unitvector();
       double s = 0.0;

       IntervalData id = input.intervalTree->interval_intersect( X, D );
       while( id.status )
       {
          X += D * id.tmin;
          s += id.tmin;
          if( s >= maxpathlength ||  (T[0] < 1.0e-6 && T[1] < 1.0e-6 && T[2] < 1.0e-6 && T[3] < 1.0e-6) ){ break; }
          double sinterval = id.tmin;
          int nbsteps = (int)((id.tmax - id.tmin)/input.ds) + 1;
          double ds = (id.tmax-id.tmin)/nbsteps;
          while( sinterval < id.tmax )
          {
             if( input.useHoldOut )
             {
                if( input.holdOut->eval(X) >= 0.0 )
                {
                   s = maxpathlength;
                   T = zeroc;
                   break;
                }
             }

             density = input.densityField->eval( X );
	     ambDensity = input.ambientDensityField->eval( X );
	     if( density < 0 ){ density = 0; }
	     if( ambDensity < 0 ){ ambDensity = 0; }
	     if( density+ambDensity > 0 )
	     {
                float totaldensity = density+ambDensity;
                dT = exp( -input.scatterCoefficient * ds * totaldensity );
                T *= dT;
	     }
             X += D * ds;
             sinterval += ds;
          }
          id = input.intervalTree->interval_intersect( X, D );
      }
      return T;
   }
   else
   {
       Color T = onec;
       Color dT = onec;
       float ds = input.ds;
       float density = 0;
       float ambDensity = 0;
       Vector X = start; 
       Vector D = (end-start).unitvector();
       double s = 0;
       double maxpathlength = (end-start).magnitude();
       int nbsteps = (int)( maxpathlength / ds ) + 1;
       ds = maxpathlength/nbsteps;
       while( s <= maxpathlength && (T[0] > 1.0e-6 || T[1] > 1.0e-6 || T[2] > 1.0e-6 || T[3] > 1.0e-6) )
       {
          if( input.useHoldOut )
          {
             if( input.holdOut->eval(X) >= 0.0 )
             {
                s = maxpathlength;
                T = zeroc;
                    break;
             }
          }

          density = input.densityField->eval( X );
	  ambDensity = input.ambientDensityField->eval( X );
	  if( density < 0 ){ density = 0; }
	  if( ambDensity < 0 ){ ambDensity = 0; }
	  if( density+ambDensity > 0 )
	  {
             float totaldensity = density+ambDensity;
             dT = exp( -input.scatterCoefficient * ds * totaldensity );
             T *= dT;
	  }
	  s += ds;
          X += D * ds;
       }
       return T;
    }
}



const float lux::PointToPointLWP( const RenderData& input, const Vector& start, const Vector& end )
{
   if( input.intervalTree->nbObjects() > 0 )
   {
       double maxpathlength = (end-start).magnitude();
       float lwp = 0.0;
       float ds = input.ds;
       int nbsteps = (int)( maxpathlength / ds ) + 1;
       ds = maxpathlength/nbsteps;
       float density = 0.0;
       float ambDensity = 0.0;
       Vector X = start;
       Vector D = (end-start).unitvector();
       double s = 0.0;

       IntervalData id = input.intervalTree->interval_intersect( X, D );
       while( id.status )
       {
          X += D * id.tmin;
          s += id.tmin;
          if( s >= maxpathlength ){ break; }
          double sinterval = id.tmin;
          int nbsteps = (int)((id.tmax - id.tmin)/input.ds) + 1;
          double ds = (id.tmax-id.tmin)/nbsteps;
          while( sinterval < id.tmax )
          {
             if( input.useHoldOut )
             {
                if( input.holdOut->eval(X) >= 0.0 )
                {
                   s = maxpathlength;
                   lwp = 1.0e6;
                   break;
                }
             }

             density = input.densityField->eval( X );
	     ambDensity = input.ambientDensityField->eval( X );
	     if( density < 0 ){ density = 0; }
	     if( ambDensity < 0 ){ ambDensity = 0; }
             lwp += ds * (density+ambDensity);
             X += D * ds;
             sinterval += ds;
          }
          id = input.intervalTree->interval_intersect( X, D );
      }
      return lwp;
   }
   else
   {
       float lwp = 0.0;
       float ds = input.ds;
       float density = 0;
       float ambDensity = 0;
       Vector X = start; 
       Vector D = (end-start).unitvector();
       double s = 0;
       double maxpathlength = (end-start).magnitude();
       int nbsteps = (int)( maxpathlength / ds ) + 1;
       ds = maxpathlength/nbsteps;
       while( s <= maxpathlength )
       {
          if( input.useHoldOut )
          {
             if( input.holdOut->eval(X) >= 0.0 )
             {
                s = maxpathlength;
                lwp = 1.0e6;
                break;
             }
          }

          density = input.densityField->eval( X );
	  ambDensity = input.ambientDensityField->eval( X );
	  if( density < 0 ){ density = 0; }
	  if( ambDensity < 0 ){ ambDensity = 0; }
          lwp += ds * (density+ambDensity);
	  s += ds;
          X += D * ds;
       }
       return lwp;
    }
}

