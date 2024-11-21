
#ifndef __PHASEFUNCTION_H__
#define __PHASEFUNCTION_H__

#include <memory>

#include <cmath>

namespace lux
{


class PhaseFunctionBase
{
  public:


    PhaseFunctionBase(){}
   ~PhaseFunctionBase(){}

    virtual const float eval( const float theta ) const { return 1.0; }

};


typedef std::shared_ptr<PhaseFunctionBase> PhaseFunction;


class UniformPhaseFunction : public PhaseFunctionBase
{
  public:

    UniformPhaseFunction( const float val = 1.0/(4.0*3.14159265) ) : value (val) {}
   ~UniformPhaseFunction(){}

    const float eval( const float theta ) const { return value; }


  private:

    float value;

};

class HenyeyGreensteinPhaseFunction : public PhaseFunctionBase
{
  public:

    HenyeyGreensteinPhaseFunction( const float val = 0.9 ) : g (val) {}
   ~HenyeyGreensteinPhaseFunction(){}

    const float eval( const float theta ) const
    {
       float costheta = std::cos(theta);
       float denom = 1.0 + g*g - 2.0*g*costheta;
       float value = (1.0 - g*g)/(4.0*3.14159265);
       value /= denom;
       return value;
    }


  private:

    float g;

};


class DoubleHenyeyGreensteinPhaseFunction : public PhaseFunctionBase
{
  public:

    DoubleHenyeyGreensteinPhaseFunction( const float g0, const float g1, const float mix ) : 
       mixture (mix),
       pf0 ( HenyeyGreensteinPhaseFunction(g0) ),
       pf1 ( HenyeyGreensteinPhaseFunction(g1) )
    {}

   ~DoubleHenyeyGreensteinPhaseFunction(){}

    const float eval( const float theta ) const
    {
       return mixture * pf0.eval(theta) + (1.0-mixture)*pf1.eval(theta);
    }


  private:

    float mixture;
    HenyeyGreensteinPhaseFunction pf0, pf1;
};




class FournierForandPhaseFunction : public PhaseFunctionBase
{
  public:

    FournierForandPhaseFunction( const float en, const float mu );
   ~FournierForandPhaseFunction(){}

    const float eval( const float theta ) const;

  private:

    float delta180;
    float nu;
    double p1factor;
};


}
#endif
