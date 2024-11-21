
#ifndef __FFTNOISE_H__
#define __FFTNOISE_H__

#include "Vector.h"
#include "Volume.h"
#include "VolumeGrid.h"
#include "Noise.h"

namespace lux
{

class FFTNoise : public Noise
{
  public:

    FFTNoise(){}
   ~FFTNoise(){}

    const float eval( const Vector& P ) const { return grid.eval(P); }


    void setParameters( const Noise_t& parameters );
    void getParameters( Noise_t& parameters ) const {}


  private:
 
   VolumeGrid<float> grid;

};


}
#endif
