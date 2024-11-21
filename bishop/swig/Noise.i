
%module  bishop
%{
#include "Noise.h"
#include "PerlinNoise.h"
%}


/*
namespace lux
{

struct Noise_t 
{

  Noise_t() :
  wavelength   (1),
  translate    (Vector(0,0,0)),
  octaves      (1.0),
  amplitude    (1),
  offset       (0),
  fjump        (2),
  roughness    (0.5),
  time         (0.0),
  fftLowCutoff (0.01),
  fftHighCutoff (1.0),
  fftPower      (3.5),
  fftNbGridPoints (128),
  fftLength     (10.0),
  lognormalmean (1.0),
  gaussianstandarddeviation (1.0),
  seed         (12345)
  {}

	float wavelength;
	Vector translate;
	float octaves;
	float amplitude;
	float offset;
	float fjump;
        float roughness;
	float time;
	float fftLowCutoff;
	float fftHighCutoff;
	float fftPower;
	int   fftNbGridPoints;
	float fftLength;
	float lognormalmean;
	float gaussianstandarddeviation;
	int   seed;
};


class Noise
{
  public:

    const float eval( const float x ) const;
    const float eval( const Vector& x ) const;

    void setParameters( const Noise_t& parameters );
    void getParameters( Noise_t& parameters ) const;
};



}

*/

%include "Noise.h"
%include "PerlinNoise.h"
%template(PerlinFractalSum) lux::FractalSum<lux::PerlinNoise>;
%template(AnchorChain) vector<lux::Noise_t>;





