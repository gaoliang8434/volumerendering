

#include "PoissonSolvers.h"
#include "GridVolumes.h"
#include "ImplicitVectorShapes.h"
#include "ImplicitVolumeShapes.h"
#include "Grids.h"
#include "Fields.h"
#include <fftw3.h>
#include <cmath>


using namespace std;
using namespace lux;




void lux::GaussSeidelPoissonSolver( VolumeGrid<float>& p, Volume<float>* source, int nbiterations, float tolerance )
{
   // The boundary condition is that p = 0 outside of the grid.
   // make a grid with the source sampled onto to
   VolumeGrid<float> gsource;
   gsource.init(   p.nx(), p.ny(), p.nz(), p.Lx(), p.Ly(), p.Lz(), p.llc()  );

   lux::Sample( &gsource, source );
   float dx = p.dx();
   dx *= dx;  // square voxels assumption

   VolumeGrid<float> temp;
   temp.init(   p.nx(), p.ny(), p.nz(), p.Lx(), p.Ly(), p.Lz(), p.llc()  );

   int nx = p.nx();
   int ny = p.ny();
   int nz = p.nz();
   p.setClearValue(0.0);
   ProgressMeter meter( nbiterations, "Gauss-Seidel" );
   for( int iter=0;iter<nbiterations;iter++ )
   {
      CopyVolumeGrid( p, temp );
      double error = 0;
      for( int k=1;k<nz-1;k++ )
      {
         for( int j=1;j<ny-1;j++ )
	 {
	    for( int i=1;i<nx-1;i++ )
	    {
	       float newvalue  = 0;
	       newvalue += temp.value(i+1,j,k);
	       newvalue += temp.value(i-1,j,k);
	       newvalue += temp.value(i,j+1,k);
	       newvalue += temp.value(i,j-1,k);
               newvalue += temp.value(i,j,k+1);
	       newvalue += temp.value(i,j,k-1);
	       newvalue = newvalue/6.0 - dx*gsource.value(i,j,k);
	       error += std::pow( newvalue - temp.value(i,j,k), 2 );
	       p.value(i,j,k) = newvalue;
	    }
	 }
      }


     // boundary conditions
         for( int j=0;j<ny;j++ )
	 {
	    for( int i=0;i<nx;i++ )
	    {
	       p.value(i,j,0) = p.value(i,j,1);
	       p.value(i,j,nz-1) = p.value(i,j,nz-2);
	    }
	 }

        for( int k=0;k<nz;k++ )
	 {
	    for( int i=0;i<nx;i++ )
	    {
	       p.value(i,0,k) = p.value(i,1,k);
	       p.value(i,ny-1,k) = p.value(i,ny-2,k);
	    }
	 }

        for( int k=0;k<nz;k++ )
	 {
	    for( int j=0;j<ny;j++ )
	    {
	       p.value(0,j,k) = p.value(1,j,k);
	       p.value(nx-1,j,k) = p.value(nx-2,j,k);
	    }
	 }


      error /= (nx-2)*(ny-2)*(nz-2);
      cout << "\rerror " << sqrt(error) << "         " << flush;
      meter.update();
      if( sqrt(error) <= tolerance )
      {
         cout << "\n\nGauss-Seidel: Passed tolerance test\n";
         return; 
      }
   }
}




void lux::GaussSeidelDivFree( VolumeGrid<Vector>& divfreeU, Volume<Vector>* U, int nbiterations, float tolerance )
{
   Volume<float>* source = new DivergenceVectorVolume( U, divfreeU.dx() );
   VolumeGrid<float> pressure;
   pressure.init( divfreeU.nx(), divfreeU.ny(), divfreeU.nz(), divfreeU.Lx(), divfreeU.Ly(), divfreeU.Lz(), divfreeU.llc() );
   pressure.setClearValue(0.0);
   lux::GaussSeidelPoissonSolver( pressure, source, nbiterations, tolerance );
   Volume<float>* pressurefield = new GriddedVolume( &pressure );
   Volume<Vector>* gradpressure = new GradientVectorVolume( pressurefield, divfreeU.dx() );
   Volume<Vector>* divfree = new SubtractVectorVolume( U, gradpressure );
   lux::Sample( &divfreeU, divfree );
}




void lux::FFTDivFree( VolumeGrid<Vector>& divfreeU, Volume<Vector>* U )
{
   fftw_complex *datax, *datay, *dataz;
   fftw_plan pforwardx, pbackwardx;
   fftw_plan pforwardy, pbackwardy;
   fftw_plan pforwardz, pbackwardz;
   datax = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * divfreeU.nx()*divfreeU.ny()*divfreeU.nz());
   datay = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * divfreeU.nx()*divfreeU.ny()*divfreeU.nz());
   dataz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * divfreeU.nx()*divfreeU.ny()*divfreeU.nz());
   pforwardx = fftw_plan_dft_3d( divfreeU.nx(), divfreeU.ny(), divfreeU.nz(), datax, datax, FFTW_FORWARD, FFTW_ESTIMATE);
   pbackwardx = fftw_plan_dft_3d( divfreeU.nx(), divfreeU.ny(), divfreeU.nz(), datax, datax, FFTW_BACKWARD, FFTW_ESTIMATE);
   pforwardy = fftw_plan_dft_3d( divfreeU.nx(), divfreeU.ny(), divfreeU.nz(), datay, datay, FFTW_FORWARD, FFTW_ESTIMATE);
   pbackwardy = fftw_plan_dft_3d( divfreeU.nx(), divfreeU.ny(), divfreeU.nz(), datay, datay, FFTW_BACKWARD, FFTW_ESTIMATE);
   pforwardz = fftw_plan_dft_3d( divfreeU.nx(), divfreeU.ny(), divfreeU.nz(), dataz, dataz, FFTW_FORWARD, FFTW_ESTIMATE);
   pbackwardz = fftw_plan_dft_3d( divfreeU.nx(), divfreeU.ny(), divfreeU.nz(), dataz, dataz, FFTW_BACKWARD, FFTW_ESTIMATE);



   Sample( &divfreeU, U );


   long index = 0;
   for( int k=0;k<divfreeU.nz();k++ )
   {
      for( int j=0;j<divfreeU.ny();j++ )
      {
         for( int i=0;i<divfreeU.nx();i++ )
	 {
            datax[index][0] = divfreeU.value(i,j,k)[0];
            datax[index][1] = 0;
            datay[index][0] = divfreeU.value(i,j,k)[1];
            datay[index][1] = 0;
            dataz[index][0] = divfreeU.value(i,j,k)[2];
            dataz[index][1] = 0;
            ++index;
	 }
      }
   }


   fftw_execute(pforwardx);
   fftw_execute(pforwardy);
   fftw_execute(pforwardz);


   float dKx = 2.0 * M_PI  / divfreeU.nx();
   float dKy = 2.0 * M_PI  / divfreeU.ny();
   float dKz = 2.0 * M_PI  / divfreeU.nz();
   float nyqK = 2.0 * M_PI;
   float scaler = 1.0/(divfreeU.nx() * divfreeU.ny() * divfreeU.nz() );
   index = 0;
   for( int k=0;k<divfreeU.nz();k++ )
   {
      float kz = k * dKz;
      if( k > divfreeU.nz()/2 ){ kz -= nyqK; }
      for( int j=0;j<divfreeU.ny();j++ )
      {
         float ky = j * dKy;
         if( j > divfreeU.ny()/2 ){ ky -= nyqK; }
         for( int i=0;i<divfreeU.nx();i++ )
         {
	    float kx = i * dKx;
            if( i > divfreeU.nx()/2 ){ kx -= nyqK; }
	    float ksq = kx*kx + ky*ky + kz*kz;

            if( ksq == 0.0 )
	    {
	       datax[index][0] = datax[index][1] = 0.0;
	       datay[index][0] = datay[index][1] = 0.0;
	       dataz[index][0] = dataz[index][1] = 0.0;
	    }
            else
	    {
            float vdotkr = kx*datax[index][0] + ky*datay[index][0] + kz*dataz[index][0];
            float vdotki = kx*datax[index][1] + ky*datay[index][1] + kz*dataz[index][1];

            datax[index][0] -= vdotkr * kx / ksq;
            datay[index][0] -= vdotkr * ky / ksq;
            dataz[index][0] -= vdotkr * kz / ksq;

            datax[index][1] -= vdotki * kx / ksq;
            datay[index][1] -= vdotki * ky / ksq;
            dataz[index][1] -= vdotki * kz / ksq;

	    datax[index][0] *= scaler;
	    datay[index][0] *= scaler;
	    dataz[index][0] *= scaler;

	    datax[index][1] *= scaler;
	    datay[index][1] *= scaler;
	    dataz[index][1] *= scaler;
            }
	    ++index;
         }
      }
   }


   fftw_execute(pbackwardx);
   fftw_execute(pbackwardy);
   fftw_execute(pbackwardz);


   index = 0;
   for( int k=0;k<divfreeU.nz();k++ )
   {
      for( int j=0;j<divfreeU.ny();j++ )
      {
         for( int i=0;i<divfreeU.nx();i++ )
	 {
	    divfreeU.value(i,j,k) = Vector( datax[index][0], datay[index][0], dataz[index][0] );
	    ++index;
         }
      }
   }


   fftw_destroy_plan(pforwardx);
   fftw_destroy_plan(pbackwardx);
   fftw_free(datax);

   fftw_destroy_plan(pforwardy);
   fftw_destroy_plan(pbackwardy);
   fftw_free(datay);
   
   fftw_destroy_plan(pforwardz);
   fftw_destroy_plan(pbackwardz);
   fftw_free(dataz);
}

void lux::FFTDivFree( VectorGrid& divfreeU, const VectorField& U )
{
   stamp( divfreeU, U, 1 );
   FFTDivFree(divfreeU);
}

void lux::FFTDivFree( VectorGrid& divfreeU )
{
   fftw_complex *datax, *datay, *dataz;
   fftw_plan pforwardx, pbackwardx;
   fftw_plan pforwardy, pbackwardy;
   fftw_plan pforwardz, pbackwardz;
   datax = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * divfreeU->nx()*divfreeU->ny()*divfreeU->nz());
   datay = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * divfreeU->nx()*divfreeU->ny()*divfreeU->nz());
   dataz = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * divfreeU->nx()*divfreeU->ny()*divfreeU->nz());
   pforwardx = fftw_plan_dft_3d( divfreeU->nx(), divfreeU->ny(), divfreeU->nz(), datax, datax, FFTW_FORWARD, FFTW_ESTIMATE);
   pbackwardx = fftw_plan_dft_3d( divfreeU->nx(), divfreeU->ny(), divfreeU->nz(), datax, datax, FFTW_BACKWARD, FFTW_ESTIMATE);
   pforwardy = fftw_plan_dft_3d( divfreeU->nx(), divfreeU->ny(), divfreeU->nz(), datay, datay, FFTW_FORWARD, FFTW_ESTIMATE);
   pbackwardy = fftw_plan_dft_3d( divfreeU->nx(), divfreeU->ny(), divfreeU->nz(), datay, datay, FFTW_BACKWARD, FFTW_ESTIMATE);
   pforwardz = fftw_plan_dft_3d( divfreeU->nx(), divfreeU->ny(), divfreeU->nz(), dataz, dataz, FFTW_FORWARD, FFTW_ESTIMATE);
   pbackwardz = fftw_plan_dft_3d( divfreeU->nx(), divfreeU->ny(), divfreeU->nz(), dataz, dataz, FFTW_BACKWARD, FFTW_ESTIMATE);

   long index = 0;
   for( int k=0;k<divfreeU->nz();k++ )
   {
      for( int j=0;j<divfreeU->ny();j++ )
      {
         for( int i=0;i<divfreeU->nx();i++ )
	 {
            datax[index][0] = divfreeU->get(i,j,k)[0];
            datax[index][1] = 0;
            datay[index][0] = divfreeU->get(i,j,k)[1];
            datay[index][1] = 0;
            dataz[index][0] = divfreeU->get(i,j,k)[2];
            dataz[index][1] = 0;
            ++index;
	 }
      }
   }


   fftw_execute(pforwardx);
   fftw_execute(pforwardy);
   fftw_execute(pforwardz);


   float dKx = 2.0 * M_PI  / divfreeU->nx();
   float dKy = 2.0 * M_PI  / divfreeU->ny();
   float dKz = 2.0 * M_PI  / divfreeU->nz();
   float nyqK = 2.0 * M_PI;
   float scaler = 1.0/(divfreeU->nx() * divfreeU->ny() * divfreeU->nz() );
   index = 0;
   for( int k=0;k<divfreeU->nz();k++ )
   {
      float kz = k * dKz;
      if( k > divfreeU->nz()/2 ){ kz -= nyqK; }
      for( int j=0;j<divfreeU->ny();j++ )
      {
         float ky = j * dKy;
         if( j > divfreeU->ny()/2 ){ ky -= nyqK; }
         for( int i=0;i<divfreeU->nx();i++ )
         {
	    float kx = i * dKx;
            if( i > divfreeU->nx()/2 ){ kx -= nyqK; }
	    float ksq = kx*kx + ky*ky + kz*kz;

            if( ksq == 0.0 )
	    {
	       datax[index][0] = datax[index][1] = 0.0;
	       datay[index][0] = datay[index][1] = 0.0;
	       dataz[index][0] = dataz[index][1] = 0.0;
	    }
            else
	    {
            float vdotkr = kx*datax[index][0] + ky*datay[index][0] + kz*dataz[index][0];
            float vdotki = kx*datax[index][1] + ky*datay[index][1] + kz*dataz[index][1];

            datax[index][0] -= vdotkr * kx / ksq;
            datay[index][0] -= vdotkr * ky / ksq;
            dataz[index][0] -= vdotkr * kz / ksq;

            datax[index][1] -= vdotki * kx / ksq;
            datay[index][1] -= vdotki * ky / ksq;
            dataz[index][1] -= vdotki * kz / ksq;

	    datax[index][0] *= scaler;
	    datay[index][0] *= scaler;
	    dataz[index][0] *= scaler;

	    datax[index][1] *= scaler;
	    datay[index][1] *= scaler;
	    dataz[index][1] *= scaler;
            }
	    ++index;
         }
      }
   }


   fftw_execute(pbackwardx);
   fftw_execute(pbackwardy);
   fftw_execute(pbackwardz);


   index = 0;
   for( int k=0;k<divfreeU->nz();k++ )
   {
      for( int j=0;j<divfreeU->ny();j++ )
      {
         for( int i=0;i<divfreeU->nx();i++ )
	 {
	    divfreeU->set(i,j,k, Vector( datax[index][0], datay[index][0], dataz[index][0] ) );
	    ++index;
         }
      }
   }


   fftw_destroy_plan(pforwardx);
   fftw_destroy_plan(pbackwardx);
   fftw_free(datax);

   fftw_destroy_plan(pforwardy);
   fftw_destroy_plan(pbackwardy);
   fftw_free(datay);
   
   fftw_destroy_plan(pforwardz);
   fftw_destroy_plan(pbackwardz);
   fftw_free(dataz);
}



VectorField lux::FFTDivFree( const GridBox& b, const VectorField& U )
{
   VectorGrid grid = makeGrid( b, Vector(0,0,0) );
   FFTDivFree( grid, U );
   return gridded(grid);
}


VectorField lux::FFTVolumePreserve( const GridBox& bb, const VectorField& U )
{

   MatrixField M = exp( grad(U) );
   VectorGrid Xgrid = makeGrid( bb, Vector(0,0,0) );

   fftw_complex *Mfft[3];
   fftw_complex *Xfft;
   fftw_plan pforward[3];
   fftw_plan pbackward;

   for( int i=0;i<3;i++)
   {
      Mfft[i] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * bb->nx()*bb->ny()*bb->nz());
      pforward[i] = fftw_plan_dft_3d( bb->nx(), bb->ny(), bb->nz(), Mfft[i], Mfft[i], FFTW_FORWARD, FFTW_ESTIMATE);
   }
   Xfft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * bb->nx()*bb->ny()*bb->nz());
   pbackward = fftw_plan_dft_3d( bb->nx(), bb->ny(), bb->nz(), Xfft, Xfft, FFTW_BACKWARD, FFTW_ESTIMATE);
   

   ///////////////////////////////////////////////////////////
   //
   //   dD/dx
   //
   ///////////////////////////////////////////////////////////

   float dKx = 2.0 * M_PI  / (bb->nx()*bb->dx());
   float dKy = 2.0 * M_PI  / (bb->ny()*bb->dy());
   float dKz = 2.0 * M_PI  / (bb->nz()*bb->dz());
   float nyqKx = 2.0 * M_PI/bb->dx();
   float nyqKy = 2.0 * M_PI/bb->dy();
   float nyqKz = 2.0 * M_PI/bb->dz();
   float scaler = 1.0/(bb->nx() * bb->ny() * bb->nz() );
   long index = 0;
   for( int col=0;col<3;col++ )
   {
   for( int k=0;k<bb->nz();k++ )
   {
      for( int j=0;j<bb->ny();j++ )
      {
         for( int i=0;i<bb->nx();i++ )
	 {
	    // Compute Displacement = Map - identity
	    Vector P = bb->evalP(i,j,k);
	    Matrix Mvalue = M->eval( P ) - unitMatrix() ;
            Mfft[0][index][0] = Mvalue(0,col);
            Mfft[0][index][1] = 0;
            Mfft[1][index][0] = Mvalue(1,col);
            Mfft[1][index][1] = 0;
            Mfft[2][index][0] = Mvalue(2,col);
            Mfft[2][index][1] = 0;
            ++index;
	 }
      }
   }


   fftw_execute(pforward[0]);
   fftw_execute(pforward[1]);
   fftw_execute(pforward[2]);


   // x component of gradient
   index = 0;
   for( int k=0;k<bb->nz();k++ )
   {
      float kz = k * dKz;
      if( k > bb->nz()/2 ){ kz -= nyqKz; }
      for( int j=0;j<bb->ny();j++ )
      {
         float ky = j * dKy;
         if( j > bb->ny()/2 ){ ky -= nyqKy; }
         for( int i=0;i<bb->nx();i++ )
         {
	    float kx = i * dKx;
            if( i > bb->nx()/2 ){ kx -= nyqKx; }

            double ksq = kx*kx + ky*ky + kz*kz;

	    if( ksq > 0.0 )
	    {
               Xfft[index][0] = -Mfft[0][index][1] * kx * scaler/ksq;
               Xfft[index][0] -= Mfft[1][index][1] * ky * scaler/ksq;
               Xfft[index][0] -= Mfft[2][index][1] * kz * scaler/ksq;

               Xfft[index][1] =  Mfft[0][index][0] * kx * scaler/ksq;
               Xfft[index][1] =  Mfft[1][index][0] * ky * scaler/ksq;
               Xfft[index][1] =  Mfft[2][index][0] * kz * scaler/ksq;
            }
	    else
	    {
	       Xfft[index][0] = Xfft[index][1] = 0.0;
	    }
	    ++index;
         }
      }
   }



   fftw_execute(pbackward);

   index = 0;
   for( int k=0;k<bb->nz();k++ )
   {
      for( int j=0;j<bb->ny();j++ )
      {
         for( int i=0;i<bb->nx();i++ )
	 {
	    Vector val = Xgrid->get(i,j,k);
	    Vector P = Xgrid->evalP(i,j,k);
	    val[col] = Xfft[index][0];
	    Xgrid->set(i,j,k,val);
	    ++index;
         }
      }
   }


   }


   fftw_destroy_plan(pforward[0]);
   fftw_destroy_plan(pforward[1]);
   fftw_destroy_plan(pforward[2]);
   fftw_destroy_plan(pbackward);
   fftw_free(Mfft[0]);
   fftw_free(Mfft[1]);
   fftw_free(Mfft[2]);
   fftw_free(Xfft);

   VectorField Map = identity() + gridded(Xgrid);

   return Map;
}





/*
VectorField lux::RayMarchDivFreeZeroNormal( const VectorField& W, const ScalarField& LS, const float resolution, const VectorField& BC )
{
   return VectorField( new RayMarchDivFreeZeroNormalVectorVolume( W, LS, resolution, BC ) ); 
}




VectorField lux::RayMarchDivFreePlanar( const VectorField& W, const Vector& PlaneP, const Vector& normal, const float resolution, const VectorField& BC )
{
   return VectorField( new RayMarchDivFreePlanarVectorVolume( W, PlaneP, normal, resolution, BC ) ); 
}
*/
