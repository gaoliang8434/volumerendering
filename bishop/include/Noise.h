

//
//
// noise.h
//
// General purpose noise code
//
// Integrates multiple noise sources into a single framework
//
//
//


#ifndef __NOISE_H__
#define __NOISE_H__

#include "Vector.h"
#include "Color.h"
#include <vector>

namespace lux
{

struct Noise_t 
{

  Noise_t() :
  frequency    (1),
  translate    (Vector(0,0,0)),
  octaves      (1.0),
  amplitude    (1),
  offset       (0),
  fjump        (2),
  roughness    (0.5),
  radius       (1.0),
  capradius    (0.1),
  pscale       (1.0),
  gamma        (1.0),
  time         (0.0),
  fftLowCutoff (0.01),
  fftHighCutoff (1.0),
  fftPower      (3.5),
  fftNbGridPoints (128),
  fftLength     (10.0),
  lognormalmean (1.0),
  gaussianstandarddeviation (1.0),
  seed         (12345),
  tangent      (Vector(0,0,1)),
  normal       (Vector(0,1,0)),
  binormal     (Vector(1,0,0)),
  axis         (Vector(0,1,0)),
  angle        (0.0),
  P            (Vector(0,0,0)),
  v            (Vector(0,0,0)),
  A            (Vector(0,0,0)),
  age          (0.0),
  lifeTime     (1.0),
  shutter      (0.5),
  frameRate    (1.0/24.0),
  falloff      (1.0),
  Cd           (Color(1,1,1,0)),
  nbWisps      (0),
  wispOctaves  (1.0),
  wispFreq     (1.0),
  wispTranslate(Vector(0,0,0)),
  wispOffset   (0.0),
  wispFjump    (2.13),
  wispRoughness(0.5),
  wispCorrelation(0.0),
  wispRadialGroup(1.0/3.0),
  wispDisplacementScale(1.0)
  {}

	float frequency;
	Vector translate;
	float octaves;
	float amplitude;
	float offset;
	float fjump;
    float roughness;
	float radius;
	float capradius;
	float pscale;
	float gamma;
	float time;
	float fftLowCutoff;
	float fftHighCutoff;
	float fftPower;
	int   fftNbGridPoints;
	float fftLength;
	float lognormalmean;
	float gaussianstandarddeviation;
	int   seed;
	Vector tangent;
	Vector normal;
	Vector binormal;
	Vector axis;
	float angle;
	Vector P;
	Vector v;
	Vector A;
	float  age;
	float lifeTime;
	float shutter;
	float frameRate;
	float falloff;
    Color Cd;
    // wisp parameters
    int nbWisps;
    float wispOctaves;
    float wispFreq;
    Vector wispTranslate;
    float wispOffset;
    float wispFjump;
    float wispRoughness;
    float wispCorrelation;
    float wispRadialGroup;
    float wispDisplacementScale;


};


class Noise
{
  public:

    Noise(){}
    virtual ~Noise(){}

    virtual const float eval( const float x ) const { return 0; }
    virtual const float eval( const Vector& x ) const { return 0; }

    virtual void setParameters( const Noise_t& parameters ){}
    virtual void getParameters( Noise_t& parameters ) const {}
};


typedef Noise_t Anchor;
typedef std::vector<Anchor> AnchorChain;
void setAnchor( AnchorChain& a, int i, const Anchor& c );
Anchor getAnchor( AnchorChain& a, int i );


const Anchor evaluateAnchorChain( const AnchorChain& alist, const double x );
const Anchor interpolateAnchors( const Anchor& a1, const Anchor& a2, const double x );


template< typename BaseNoise>
class FractalSum : public Noise
{
  public:

    FractalSum() : 
       octaves    (3.0), 
       fjump      (2.0), 
       roughness  (0.5),
       frequency  (0.666666),
       translate  (Vector(0,0,0)),
       offset     (0.0),
       axis       (Vector(0,0,1)),
       angle      (0.0),
       scaling    (1.0)
    {}

   ~FractalSum(){}

    const float eval( const float x ) const
    {
       float exponent = 1;
       float amplitude = 1;
       float accum = 0;
       int ioct = (int)octaves;
       for( int oc=0;oc<ioct;oc++ )
       {
          const float X = (x - translate[0]) * frequency * exponent;
          accum += amplitude * _noise.eval( X );
          exponent *= fjump;
	  amplitude *= roughness;
       }
       const float X = (x - translate[0]) * frequency * exponent;
       float val = amplitude * _noise.eval( X );
       accum += (octaves - (int)octaves) * val;
       return (accum + offset)*scaling;
    }

    const float eval( const Vector& x ) const
    {
       float exponent = 1;
       float amplitude = 1;
       float accum = 0;
       int ioct = (int)octaves;
       Vector X = (x - translate);
       if( angle != 0.0 )
       {
          float ca = cos(angle);
	  float sa = sin(angle);
	  X = X*ca + axis*(axis*X)*(1.0-ca) + (axis^X)*sa;
       }
       X *= frequency*exponent;
       for( int oc=0;oc<ioct;oc++ )
       {
          accum += amplitude * _noise.eval( X );
          X *= fjump;
	  amplitude *= roughness;
       }
       float val = amplitude * _noise.eval( X );
       accum += (octaves - (int)octaves) * val;
       return (accum+offset)*scaling;
    }


    void setParameters( const Noise_t& parameters )
    {
       octaves = parameters.octaves;
       fjump = parameters.fjump;
       roughness = parameters.roughness;
       frequency = parameters.frequency; 
       translate = parameters.translate; 
       offset = parameters.offset; 
       axis = parameters.axis;
       angle = parameters.angle;
       scaling = parameters.amplitude;
       _noise.setTime( parameters.time );
    }

   
    void getParameters( Noise_t& parameters ) const
    {
       parameters.octaves = octaves;
       parameters.fjump = fjump;
       parameters.roughness = roughness;
       parameters.frequency = frequency; 
       parameters.translate = translate; 
       parameters.offset = offset; 
       parameters.axis = axis;
       parameters.angle = angle;
       parameters.amplitude = scaling;
    }
 


  private:

    BaseNoise _noise;

    float octaves;
    float fjump;
    float roughness;
    float frequency;
    Vector translate;
    float offset;
    Vector axis;
    float angle;
    float scaling;
};



class PRN
{
  public:

    PRN(){}
    virtual ~PRN(){}

    virtual const float eval(){ return 0; }

    virtual void setParameters( const Noise_t& parameters ){}
};






}

#endif



