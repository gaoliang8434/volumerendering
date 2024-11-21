

#include "NoiseMachine.h"

using namespace std;
using namespace lux;





NoiseMachine lux::perlin( Noise_t nparm )
{
   Noise* n = new FractalSum<PerlinNoiseGustavson>();
   n->setParameters( nparm );
   return NoiseMachine(n);
}




