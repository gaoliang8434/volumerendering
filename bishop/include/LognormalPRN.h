
#ifndef __LOGNORMALPRN_H__
#define __LOGNORMALPRN_H__

#include "GaussianPRN.h"


namespace lux
{


class LognormalPRN : public PRN
{
  public:

    LognormalPRN();

   ~LognormalPRN();

    void setParameters( const Noise_t& n ); 
    const float eval();

  private:

    GaussianPRN generator;
    float mean;
    float devFactor, devScale;


};



}

#endif

