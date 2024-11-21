
#include "cfd.h"
#include "GridVolumes.h"
#include "PoissonSolvers.h"

using namespace lux;



   GasSystem::GasSystem( const int nx, const int ny, const int nz, const float Lx, const float Ly, const float Lz, const Vector& origin, 
              const Vector& g, const float couple, const float refdensity) :
      gravity(g),
      coupling(couple),
      referenceDensity(refdensity),
      densitySource (0),
      externalForce(0),
      stepsPerUpdate(1),
      nbDivFreeIterations (1000), 
      toleranceDivFree (0.00001)
      {
         velocity.init(nx,ny,nz,Lx,Ly,Lz,origin);
         previousVelocity.init(nx,ny,nz,Lx,Ly,Lz,origin);
         density.init(nx,ny,nz,Lx,Ly,Lz,origin);
         previousDensity.init(nx,ny,nz,Lx,Ly,Lz,origin);

	 velocity.setClearValue( Vector(0,0,0) );
	 previousVelocity.setClearValue( Vector(0,0,0) );
	 density.setClearValue( 0 );
	 previousDensity.setClearValue( 0 );

	 nbDivFreeIterations = nx*ny*nz;
      }




   void GasSystem::update( const float dt, const string advectionMethod )
   {
      float ddt = dt / stepsPerUpdate;
      for( int spu=0;spu<stepsPerUpdate;spu++ )
      {
      // copy fields to previous values
      for( int k=0;k<velocity.nz();k++ )
      {
         for( int j=0;j<velocity.ny();j++ )
         {
            for( int i=0;i<velocity.nx();i++ )
            {
	       previousDensity.value(i,j,k) = density.value(i,j,k);
	       previousVelocity.value(i,j,k) = velocity.value(i,j,k);
            }
         }
      }


      Volume<Vector>* velocityfield = new GriddedVectorVolume( &previousVelocity );
      Volume<float>*  densityfield = new GriddedVolume( &previousDensity );

      // advect density and add source
      Volume<float>* advectedDensity = new  AdvectVolume( densityfield, velocityfield, ddt );
      VolumeGrid<float> forwardGrid;
      if( advectionMethod == "bfecc" )
      {
      // should try sampling to grid between forward and backward advections
	 forwardGrid.init( density.nx(),density.ny(),density.nz(),density.Lx(),density.Ly(),density.Lz(),density.llc() );
         Sample( &forwardGrid, advectedDensity );
	 Volume<float>* backwardDensity = new AdvectVolume( new GriddedVolume( &forwardGrid ), velocityfield, -ddt );
	 advectedDensity = new AddVolume( advectedDensity, new MultiplyVolume(  new SubtractVolume( densityfield, backwardDensity ) , 0.5 ) );
      }
      if( densitySource != 0 )
      {
         advectedDensity = new AddVolume( advectedDensity, new MultiplyVolume( densitySource, ddt ) );
      }
      Sample( &density, advectedDensity );

      // advect velocity and add source
     
      Volume<Vector>* advectedVelocity = new  AdvectVectorVolume( velocityfield, velocityfield, ddt );
      VolumeGrid<Vector> forwardVelGrid;
      if( advectionMethod == "bfecc" )
      {
	 forwardVelGrid.init( density.nx(),density.ny(),density.nz(),density.Lx(),density.Ly(),density.Lz(),density.llc() );
         Sample( &forwardVelGrid, advectedVelocity );
	 Volume<Vector>* backwardVelocity = new AdvectVectorVolume( new GriddedVectorVolume( &forwardVelGrid ), velocityfield, -ddt );
	 advectedVelocity = new AddVectorVolume( advectedVelocity, new MultiplyVectorVolume(  new SubtractVectorVolume( velocityfield, backwardVelocity ) , 0.5 ) );
      }
      Volume<float>* relativeDensity = new SubtractVolume( densityfield, new ConstantVolume( referenceDensity) );
      Volume<Vector>* bouyancyForce = new MultiplyVectorVolume( new ConstantVectorVolume( -gravity*coupling*ddt ), relativeDensity );
      advectedVelocity = new AddVectorVolume( advectedVelocity, bouyancyForce );

      if( externalForce != 0 )
      {
         advectedVelocity = new AddVectorVolume( advectedVelocity, new MultiplyVectorVolume( externalForce, ddt ) );
      }

      // project velocity to divergence-free
      //GaussSeidelDivFree( velocity, advectedVelocity, nbDivFreeIterations, toleranceDivFree );
      FFTDivFree( velocity, advectedVelocity );

      }
   }

   void GasSystem::InitializeDensity(  Volume<float>* d )
   {
      Sample( &density, d );
      Sample( &previousDensity, d );
   }
   void GasSystem::InitializeVelocity(  Volume<Vector>* v )
   {
      Sample( &velocity, v );
      Sample( &previousVelocity, v );
   }



