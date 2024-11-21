
#ifndef __NOISEMACHINE_H__
#define __NOISEMACHINE_H__

#include "Noise.h"
#include "PerlinNoise.h"
#include <memory>

namespace lux
{

typedef std::shared_ptr<Noise>  NoiseMachineBase;


class NoiseMachine : public NoiseMachineBase
{
  public:

    NoiseMachine() :  std::shared_ptr<Noise>() {}
    NoiseMachine( Noise* f ) :  std::shared_ptr<Noise>( f ) {}
   ~NoiseMachine(){}

};


NoiseMachine perlin( Noise_t n );






}

#endif
