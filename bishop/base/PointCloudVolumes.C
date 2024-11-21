
#include "PointCloudVolumes.h"
#include "Fields.h"
#include "SparseGrid.h"
#include "Grids.h"
//#include "PointCloudBlends.h"
#include "Noise.h"
#include "NoiseMachine.h"
#include "Stamp.h"


using namespace lux;

void lux::warp( PointCloud& pc, const VectorField& X )
{
   for(size_t i=0;i<pc->nb_particles();i++ )
   {
      const Vector& pos = evaluate( X, pc->pos(i) );
      pc->set_pos( i, pos );
   }   
}


void lux::self_advect( PointCloud& pc, const float dt )
{
   for(size_t i=0;i<pc->nb_particles();i++ )
   {
      Vector pos = pc->pos(i);
      const Vector& vel = pc->vel(i);
      pos += vel * dt;
      pc->set_pos( i, pos );
   }   
}


void lux::advect( PointCloud& pc, const VectorField& u, const float dt )
{
   for(size_t i=0;i<pc->nb_particles();i++ )
   {
      Vector pos = pc->pos(i);
      const Vector& vel = evaluate( u, pos );
      pos += vel * dt;
      pc->set_pos( i, pos );
      pc->set_vel( i, vel );
   }   
}

void lux::stamp( const PointCloud& pc, ScalarGrid& g, const string attr )
{
   if( pc->float_attr_exists( attr ) == false ){ return; }   
   GridBox gb = makeGridBox( *g );
   ScalarGrid count = makeGrid( gb, 0.0 );
   ProgressMeter pm( pc->nb_particles(), "PointCloud Stamp");
   for( size_t i=0;i<pc->nb_particles();i++ )
   {
      float pscale = pc->pscale(i);
      Vector pos = pc->pos(i);
      float value = pc->get_float_attr( attr, i );
      stamp( g, pos, pscale, value );
      stamp( count, pos, pscale, 1.0 );
      pm.update();
   }
   int nx = g->nx();
   int ny = g->ny();
   int nz = g->nz();
   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
            float cvalue = count->get(i,j,k);
            if( cvalue > 0.0 )
            {
               float value = g->get(i,j,k);
               g->set(i,j,k, value/cvalue );
            }
         }
      }
   }
}

void lux::stamp( const PointCloud& pc, VectorGrid& g, const string attr )
{
   if( pc->vector_attr_exists( attr ) == false ){ return; }   
   GridBox gb = makeGridBox( *g );
   ScalarGrid count = makeGrid( gb, 0.0 );
   ProgressMeter pm( pc->nb_particles(), "PointCloud Stamp");
   for( size_t i=0;i<pc->nb_particles();i++ )
   {
      float pscale = pc->pscale(i);
      Vector pos = pc->pos(i);
      const Vector value = pc->get_vector_attr( attr, i );
      stamp( g, pos, pscale, value );
      stamp( count, pos, pscale, 1.0 );
   }
   int nx = g->nx();
   int ny = g->ny();
   int nz = g->nz();
   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
            float cvalue = count->get(i,j,k);
            if( cvalue > 0.0 )
            {
               Vector value = g->get(i,j,k);
               g->set(i,j,k, value/cvalue );
            }
         }
      }
   }
}

void lux::stamp( const PointCloud& pc, ColorGrid& g, const string attr )
{
   if( pc->color_attr_exists( attr ) == false ){ return; }   
   GridBox gb = makeGridBox( *g );
   ScalarGrid count = makeGrid( gb, 0.0 );
   ProgressMeter pm( pc->nb_particles(), "PointCloud Stamp");
   for( size_t i=0;i<pc->nb_particles();i++ )
   {
      float pscale = pc->pscale(i);
      Vector pos = pc->pos(i);
      const Color value = pc->get_color_attr( attr, i );
      stamp( g, pos, pscale, value );
      stamp( count, pos, pscale, 1.0 );
   }
   int nx = g->nx();
   int ny = g->ny();
   int nz = g->nz();
   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
            float cvalue = count->get(i,j,k);
            if( cvalue > 0.0 )
            {
               Color value = g->get(i,j,k);
               g->set(i,j,k, value/cvalue );
            }
         }
      }
   }
}

void lux::stampDensityVelocityColor( const PointCloud& pc, ScalarGrid& gs, VectorGrid& gv, ColorGrid& gc )
{
   GridBox gb = makeGridBox( *gs );
   ScalarGrid count = makeGrid( gb, 0.0 );
   ProgressMeter pm( pc->nb_particles(), "PointCloud Stamp");
   for( size_t i=0;i<pc->nb_particles();i++ )
   {
      float pscale = pc->pscale(i);
      const Vector pos = pc->pos(i);
      const Vector vel = pc->vel(i);
      float density = pc->density(i);
      const Color ci = pc->ci(i);
      stamp( gs, pos, pscale, density );
      stamp( gv, pos, pscale, vel );
      stamp( gc, pos, pscale, ci );
      stamp( count, pos, pscale, 1.0 );
   }
   int nx = gs->nx();
   int ny = gs->ny();
   int nz = gs->nz();
   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
            float cvalue = count->get(i,j,k);
            if( cvalue > 0.0 )
            {
               Color civalue = gc->get(i,j,k);
               gc->set(i,j,k, civalue/cvalue );
               float fvalue = gs->get(i,j,k);
               gs->set(i,j,k, fvalue/cvalue );
               Vector vvalue = gv->get(i,j,k);
               gv->set(i,j,k, vvalue/cvalue );
            }
         }
      }
   }
}

void lux::stampDensityColor( const PointCloud& pc, ScalarGrid& gs, ColorGrid& gc )
{
   GridBox gb = makeGridBox( *gs );
   ScalarGrid count = makeGrid( gb, 0.0 );
   ProgressMeter pm( pc->nb_particles(), "PointCloud Stamp");
   for( size_t i=0;i<pc->nb_particles();i++ )
   {
      float pscale = pc->pscale(i);
      const Vector pos = pc->pos(i);
      float density = pc->density(i);
      const Color ci = pc->ci(i);
      stamp( gs, pos, pscale, density );
      stamp( gc, pos, pscale, ci );
      stamp( count, pos, pscale, 1.0 );
   }
   int nx = gs->nx();
   int ny = gs->ny();
   int nz = gs->nz();
   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
            float cvalue = count->get(i,j,k);
            if( cvalue > 0.0 )
            {
               Color civalue = gc->get(i,j,k);
               gc->set(i,j,k, civalue/cvalue );
               float fvalue = gs->get(i,j,k);
               gs->set(i,j,k, fvalue/cvalue );
            }
         }
      }
   }
}


/*
void stamp( const PointCloud& pc, ScalarGrid& g, const string attr, const string blendmethod )
{
   if( pc->float_attr_exists( attr ) == false ){ return; }   
   PointCloudBlendMethod* blender = 0;
   if( blendmethod == "add" ){ blender = new PointCloudAddBlend(); }
   if( blendmethod == "max" ){ blender = new PointCloudMaxBlend(); }

   if( blender == 0 ){ return; }
   ProgressMeter pm( pc->nb_particles(), "PointCloud Stamp");
   for( size_t i=0;i<pc->nb_particles();i++ )
   {
      float pscale = pc->pscale(i);
      Vector pos = pc->pos(i);
      float value = pc->get_float_attr( attr, i );
      stamp( g, pos, pscale, value );
      pm.update();
   }
   delete blender;
}




void stamp( const PointCloud& pc, VectorGrid& g, const string attr, const string blendmethod );
void stamp( const PointCloud& pc, ColorGrid& g, const string attr, const string blendmethod );

*/



void lux::stampNoise( const PointCloud& pc, ScalarGrid& g, const float fade)
{
   Noise_t noise_parameters;
   Noise* n = new FractalSum<PerlinNoiseGustavson>();
   n->setParameters( noise_parameters );
   ProgressMeter pm( pc->nb_particles(), "PointCloud Stamp Noise");
   for( size_t i=0;i<pc->nb_particles();i++ )
   {
      float pscale = pc->pscale(i);
      Vector pos = pc->pos(i);
      noise_parameters.translate = pc->translate(i);
      noise_parameters.roughness = pc->roughness(i);
      noise_parameters.octaves = pc->octaves(i);
      noise_parameters.fjump = pc->fjump(i);
      noise_parameters.frequency = pc->frequency(i);
      noise_parameters.amplitude = pc->density(i);
      n->setParameters(noise_parameters);
      StampNoise( g, n, pos, pscale, fade );
      pm.update();
   }
   delete n;
}
