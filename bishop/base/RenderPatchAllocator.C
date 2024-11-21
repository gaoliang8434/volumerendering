
#include "RenderPatchAllocator.h"

#include <cstdlib>

using namespace std;
using namespace lux;

RenderPatchAllocator::RenderPatchAllocator( int width, int height, const Camera& cam ) :
   imgWidth  (width),
   imgHeight (height),
   camera    (cam)
{}

const int RenderPatchAllocator::setPatches( const int patchsize, const string pattern, int x, int y, int nbRaysPP )
{
   raysPerPixel = nbRaysPP;
   uncompletedPatches.clear();
   completedPatches.clear();
   // Right now the simplest possible pattern, more later

   int nbWholePatches = imgWidth*imgHeight/patchsize;
   PixelPatch patch;
   if( nbWholePatches > 0 )
   {
   size_t ps = (size_t)patchsize;
   patch.resize(ps);
   for( size_t ip=0;ip<ps;ip++ ){ patch[ip] = ip; } 
   for( int i=0;i<nbWholePatches;i++ )
   {
      uncompletedPatches.push_back( patch );
      for( size_t ip=0;ip<ps;ip++ ){ patch[ip] += ps; } 
   }
   }
   // Last patch is small
   int lastpatchsize = imgWidth*imgHeight - nbWholePatches*patchsize;

   if( lastpatchsize > 0 )
   {
      patch.clear();
      patch.resize( (size_t)lastpatchsize );
      size_t start = (size_t) (nbWholePatches*patchsize);
      for( size_t ip=0;ip<(size_t)lastpatchsize;ip++ ){ patch[ip] = start + ip; } 
      uncompletedPatches.push_back( patch );
   }
   
   Noise_t jitterparms;
   jitterparms.seed = raysPerPixel * imgWidth*imgHeight;
   jitter.setParameters( jitterparms );

   return (int)uncompletedPatches.size();
}


void RenderPatchAllocator::reset( int seed )
{
   completedPatches.insert( completedPatches.end(), uncompletedPatches.begin(), uncompletedPatches.end() );
   uncompletedPatches.clear();
   uncompletedPatches.swap( completedPatches );
   Noise_t jitterparms;
   jitterparms.seed = raysPerPixel * imgWidth*imgHeight + seed;
   jitter.setParameters( jitterparms );
}

void RenderPatchAllocator::popPatch( vector<Vector>& startP, vector<Vector>& startD, PixelPatch& pixels )
{
   startP.clear();
   startD.clear();
   pixels.clear();
   if( uncompletedPatches.empty() )
   {
      return;
   }


   pixels.resize( uncompletedPatches.front().size() * raysPerPixel );
   size_t index = 0;
   for( size_t p=0;p<uncompletedPatches.front().size();p++ )
   {
      for( int ii=0;ii<raysPerPixel;ii++)
      {
         pixels[index++] = uncompletedPatches.front()[p];
      }
   }
   //startP.resize( pixels.size() );
   //startD.resize( pixels.size() );
   for( size_t ip=0;ip<pixels.size();ip++ )
   {
      int iy = ((int)pixels[ip]) / imgWidth;
      int ix = ((int)pixels[ip]) - iy*imgWidth;

      float x = ((float)ix + (drand48()-0.5))/(float)imgWidth;
      float y = ((float)iy + (drand48()-0.5))/(float)imgHeight;

      Vector D = camera.view( x, y );
      startD.push_back( D );
      startP.push_back( camera.eye() + D*camera.nearPlane() );
   }

   uncompletedPatches.erase( uncompletedPatches.begin() );
   completedPatches.push_back( pixels );
}


void RenderPatchAllocator::selectPatch( const size_t patch, vector<Vector>& startP, vector<Vector>& startD, PixelPatch& pixels )
{
   startP.clear();
   startD.clear();
   if( uncompletedPatches.empty() )
   {
      pixels.clear();
      return;
   }

   pixels = uncompletedPatches[patch];
   startP.resize( pixels.size() );
   startD.resize( pixels.size() );
   for( size_t ip=0;ip<pixels.size();ip++ )
   {
      int iy = ((int)pixels[ip]) / imgWidth;
      int ix = ((int)pixels[ip]) - iy*imgWidth;

      float x = ((float)ix + jitter.eval())/(float)imgWidth;
      float y = ((float)iy + jitter.eval())/(float)imgHeight;

      startD[ip] = camera.view( x, y );
      startP[ip] = camera.eye() + camera.view(x,y)*camera.nearPlane();
   }
}

