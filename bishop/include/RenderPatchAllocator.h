
#ifndef __RENDERPATCHALLOCATOR_H__
#define __RENDERPATCHALLOCATOR_H__

#include <vector>
#include <string>

#include "Vector.h"
#include "Camera.h"

#include "UniformPRN.h"


class RenderPatchAllocator
{
  public:

    RenderPatchAllocator( int width, int height, const Camera& cam );
   ~RenderPatchAllocator(){}

    typedef std::vector<size_t> PixelPatch;

    const int setPatches( const int patchsize, const std::string pattern, int x, int y, int nbRaysPP );

    void reset( int seed = 0 );

    void popPatch( std::vector<lux::Vector>& startP, std::vector<lux::Vector>& startD, PixelPatch& pixels );
    void selectPatch( size_t patch, std::vector<lux::Vector>& startP, std::vector<lux::Vector>& startD, PixelPatch& pixels );

    size_t nbUncompletedPatches() const { return uncompletedPatches.size(); }
    size_t nbCompletedPatches() const { return completedPatches.size(); }

    const int getRaysPerPixel() const { return raysPerPixel; }

  private:

    int imgWidth, imgHeight;
    int patchSize;
    int raysPerPixel;
    Camera camera;

    std::vector<PixelPatch> uncompletedPatches;
    std::vector<PixelPatch> completedPatches;

    lux::UniformPRN jitter;
};


#endif

