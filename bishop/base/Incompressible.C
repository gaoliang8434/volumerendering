
#include "Incompressible.h"
#include "Fields.h"
#include "Grids.h"
#include <cmath>
#include <cfloat>

using namespace std;
namespace lux
{


const VectorField RMIncompressibleGradient( const VectorField& gf, const GridBox& g, const ScalarField& S, const VectorField& nSu, const VectorField& base_nS, const float threshold, const int div_width, const int nbblur )
{
   int nx = g->nx();
   int ny = g->ny();
   int nz = g->nx();
   float dx = g->dx();
   float dy = g->dy();
   float dz = g->dx();

   double basestep = dx;
   if( basestep > dy ){ basestep = dy; }
   if( basestep > dz ){ basestep = dz; }

   // Setting up the necessary fields: 
   //   Tangential gradient \nabla_t \phi
   //   Curvature 
   //   Divergence of the tangential gradient \nabla\cdot\nabla_t\phi
   //   Several of them are sampled to a grid to (1) speed up the computation,
   //      and (2) allow blurring to reduce noise and instability
   VectorField nS = unitvector( nSu );
   VectorField bnS = unitvector( base_nS );
   ScalarField curvature = fdboundeddiv(nS, div_width, g);
   ScalarGrid curvg = makeGrid( g, 0.0 );
   stamp( curvg, curvature, 1 );
   curvature = gridded( curvg );
   ScalarField ngf = dot(nS,gf);
   VectorField tgrad = subtract(gf, multiply(nS, ngf));
   VectorGrid tgradg = makeGrid( g, Vector(0,0,0) );
   stamp( tgradg, tgrad, 1 );
   tgrad = gridded( tgradg );
   ScalarField divtgrad = fdboundeddiv(tgrad, div_width, g);
   ScalarGrid divtgradg = makeGrid( g, 0.0 );
   stamp( divtgradg, divtgrad, 1 );
   for( int i=0;i<nbblur;i++ )
   {
      Blur(divtgradg);
   }
   divtgrad = gridded( divtgradg );

   // Grid to store the incompressible gradient
   VectorGrid igrad = makeGrid( g, Vector(0,0,0) );
   // Initialize by filling the IG grid with the tangential gradient
   stamp( igrad, tgrad, 1 );

   const double lsresolution = 0.0;

   ProgressMeter pm( ny*nz, "RMIG" );
   double angle_change_threshold = std::cos( 20.0*M_PI/180.0 );

   // This can be parallelized because the sparsegrid has already been filled out.
   for(int k=0;k<nz;k++)
   {
      for(int j=0;j<ny;j++)
      {
#pragma omp parallel for
          for(int i=0;i<nx;i++)
          {
             std::vector<float> pathvalues;
             std::vector<float> pathcurvatures;
             Vector X0 = g->evalP(i,j,k);
             Vector X = X0 ;
             Vector Xold = X0;
             float lsvalue = S->eval(X);
             float lsvalueold = lsvalue;
             Vector dirold = nS->eval(X0);
             Vector dir = dirold;
             double step = basestep; 
             // build collection of values along the ray march
             pathvalues.push_back( step*divtgrad->eval(X0) );
             pathcurvatures.push_back( step*curvature->eval(X0) );
             int pass = 0;
             int nbvalues = 0;
             if( lsvalue >= 0 )
             {
                X -= dir * step;
                lsvalue = S->eval(X);
                while(lsvalue >= -lsresolution && g->isInside(X) && nbvalues < 100 )
                {
                   nbvalues = nbvalues + 1;

                   //This test is relevant when the input surface is dynamics,
                   //because the surface can become very irregular and the ray march
                   //experiences sudden turns.
                   if( dir*dirold < angle_change_threshold && pass < 10 )
                   {
                      // back up and step smaller
                      dir = dirold;
                      step *= 0.5;
                      X = Xold;
                      lsvalue = lsvalueold;
                      pass = pass + 1;
                   }
                   else
                   {
                      if( pass >= 10 )
                      {
                         // Having trouble stepping accurately
                         // Successive refinements have not worked,
                         // so we bail on it and
                         // use the base surface for a full step
                         dir = bnS->eval(Xold);
                         X = Xold - dir*basestep;
                      }
                      float kappa = step*curvature->eval(X);
                      float dd = step*divtgrad->eval(X);
                      pathcurvatures.push_back(kappa);
                      pathvalues.push_back(dd);
                      pass = 0;
                      step = basestep;
                   }
                   Xold = X;
                   dirold = dir;
                   lsvalueold = lsvalue;
                   dir = nS->eval(X);
                   X -= dir * step;
                   lsvalue = S->eval(X);
                }
             }
             else
             {
                X += dir * step;
                lsvalue = S->eval(X);
                while(lsvalue <= -lsresolution && g->isInside(X) && nbvalues < 100 )
                {
                   nbvalues = nbvalues + 1;
                   if( dir*dirold < angle_change_threshold && pass < 10 )
                   {
                      // back up and step smaller
                      dir = dirold;
                      step *= 0.5;
                      X = Xold;
                      lsvalue = lsvalueold;
                      pass = pass + 1;
                   }
                   else
                   {
                      if( pass >= 10 )
                      {
                         // Having trouble stepping accurately
                         // Successive refinements have not worked,
                         // so we bail on it and
                         // use the base surface for a full step
                         dir = bnS->eval(Xold);
                         X = Xold - dir*basestep;
                      }
                      float kappa = step*curvature->eval(X);
                      float dd = step*divtgrad->eval(X);
                      pathcurvatures.push_back(kappa);
                      pathvalues.push_back(dd);
                      pass = 0;
                      step = basestep;
                   }
                   Xold = X;
                   dirold = dir;
                   lsvalueold = lsvalue;
                   dir = nS->eval(X);
                   X += dir * step;
                   lsvalue = S->eval(X);
                }
             }
             if( !pathvalues.empty() && g->isInside(X) && nbvalues < 100)
             {
                size_t pathsize = pathvalues.size();
                double accum = 0;
                for( size_t l=0;l<pathsize;l++ )
                {
                   size_t ll = pathsize-l-1;
                   double kappa = pathcurvatures[ll];
                   double v = pathvalues[ll];
                   if( v < -threshold ){ v = -threshold; }
                   if( v > threshold ){ v = threshold; }
                   if( fabs(kappa) > FLT_EPSILON )
                   {
                      double ek = std::exp(-kappa);
                      // threshold this because it can get too big
                      double vvv = accum*ek + v * (1.0-ek)/(kappa);
                      if( fabs(kappa) >= 1.5 ){ vvv = 0.0; }
                      accum = vvv;
                   }
                   else
                   {
                      accum += v;
                   }
                }
                // threshold final result
                Vector result = igrad->get(i,j,k) + nS->eval(X0) * accum;
                if( result.magnitude() > threshold ){ result *= threshold/result.magnitude(); }
                igrad->set(i,j,k, result);
             }
             else
             {
                Vector result = igrad->get(i,j,k);
                // threshold final result
                if( result.magnitude() > threshold )
                {
                   result *= threshold/result.magnitude();
                   igrad->set(i,j,k, result);
                }
             }
          }
          pm.update();
       }
    }
    return gridded( igrad ); 
}
}
