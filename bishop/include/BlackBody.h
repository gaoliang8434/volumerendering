

#ifndef __BLACKBODY_H__
#define __BLACKBODY_H__

#include <utility>
#include <vector>
#include <Volume.h>

namespace lux
{

class BlackBodyEmission
{
  public:

    BlackBodyEmission(){}
   ~BlackBodyEmission(){}

    //void setSpectralBands(  const std::vector< std::pair< float,float > >& bands );
    //void setNbLevels( const int n ) { nbLevels = n; }

    //void emission( const float T, std::vector<float>& value ) const;
    //std::vector<float> emission( const float T ) const;


    void emission( const float T, Color& value ) const;
    const Color emission( const float T ) const;

  private:

    //std::vector< std::pair< float, float >  > spectralBands;
    //std::vector< std::pair< double, double >  > spectralBandMinMaxs;
    //int nbLevels;
    //double hbar, speedOfLight, boltzmann, chbaroverk;


    //double sumup( const double x ) const;


};



class BlackBodyVolume : public Volume<Color>, public BlackBodyEmission
{
  public:

    BlackBodyVolume( Volume<float>* temp ) :
      elem(temp)
    {}

    BlackBodyVolume( const ScalarField& temp ) :
      elem(temp)
    {}

   ~BlackBodyVolume(){}
    

    const Color eval( const Vector& P ) const
    {
       const float T = elem->eval(P);
       if ( T <= 0.0 )
       {
          return Color( 0,0,0,1);
       }
       return emission( T );
    }

  private:

    ScalarField elem;
};



}
#endif
