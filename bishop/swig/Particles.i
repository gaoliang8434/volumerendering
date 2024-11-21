
%module bishop
%{
#include "Particle.h"
%}

%include "Particle.h"
%template(ParticleGroupA) std::vector<lux::Particle>;

