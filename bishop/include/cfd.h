
#ifndef __CFD_H__
#define __CFD_H__


#include "Volume.h"
#include "VolumeGrid.h"
#include "Vector.h"
#include "ImplicitVolumeShapes.h"
#include "ImplicitVectorShapes.h"

using namespace std;


namespace lux
{


class GasSystem
{
  public:

   GasSystem( const int nx, const int ny, const int nz, const float Lx, const float Ly, const float Lz, const Vector& origin, 
              const Vector& g, const float couple, const float refdensity );
  ~GasSystem(){}


   const VolumeGrid<Vector>& getVelocityGrid() const { return velocity; }
   const VolumeGrid<float>& getDensityGrid() const { return density; }

   void setDensitySource( Volume<float>* s ) { densitySource = s; }
   void setExternalForce( Volume<Vector>* f ) { externalForce = f; }

   void update( const float dt,  const string advectionMethod = "semilagrangian" );

   void InitializeDensity(  Volume<float>* d );
   void InitializeVelocity(  Volume<Vector>* v );

   void setStepsPerUpdate( const int nb ){ stepsPerUpdate = nb; }
   void setProjectionTolerance( const float tol ){ toleranceDivFree = tol; }
   void setNbDivFreeIterations( const int nb ){ nbDivFreeIterations = nb; }


  private:


     VolumeGrid<Vector> velocity, previousVelocity;
     VolumeGrid<float> density, previousDensity;

     Vector gravity;
     float coupling;
     float referenceDensity;

     Volume<float>* densitySource;
     Volume<Vector>* externalForce;

     int stepsPerUpdate;
     int nbDivFreeIterations;
     float toleranceDivFree;

};


}


#endif
