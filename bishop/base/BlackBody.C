
#include "BlackBody.h"
#include <cmath>

using namespace std;
using namespace lux;

/*
BlackBodyEmission::BlackBodyEmission() :
   nbLevels     (3),
   hbar         (1.05457148e-34),
   speedOfLight (299792458.0),
   boltzmann    (1.3806503e-23)
{
   chbaroverk = speedOfLight * hbar / boltzmann;
   // standard bands for rgb
   vector< pair<float,float> > rgb;
   //rgb.push_back( make_pair( 400.0e-9,500.0e-9 )  );
   //rgb.push_back( make_pair( 450.0e-9,630.0e-9 )  );
   //rgb.push_back( make_pair( 500.0e-9,700.0e-9 )  );

   rgb.push_back( make_pair( 400.0e-9,525.0e-9 )  );
   rgb.push_back( make_pair( 525.0e-9,600.0e-9 )  );
   rgb.push_back( make_pair( 600.0e-9,700.0e-9 )  );

   setSpectralBands( rgb );

}



BlackBodyEmission::~BlackBodyEmission() {}


void BlackBodyEmission::emission( const float T, vector<float>& value ) const
{
   if( value.size() != spectralBands.size() )
   {
      value.clear();
      value.resize( spectralBands.size() );
   }

   float scale = pow( (double)T, (double)4.0 ) / 1.0e11 ;
   for( size_t i=0;i<spectralBands.size();i++ )
   {
      float val = scale*( sumup( spectralBandMinMaxs[i].first/T ) - sumup( spectralBandMinMaxs[i].second/T ) );
      value[i] = val;
   }
}


vector<float> BlackBodyEmission::emission( const float T ) const
{
   vector<float> result;
   emission( T, result );
   return result;
}

double BlackBodyEmission::sumup( const double x ) const
{
   double sum = 0.0;
   double x2 = x*x;
   double x3 = x2*x;

   for( int i=1;i<=nbLevels;i++)
   {
      sum += exp( -i*x ) * ( (6.0/(i*i*i*i)) + 6.0*x/(i*i*i) + 3.0*x2/(i*i)  + x3/i );   
   }
   return sum;
}


void BlackBodyEmission::setSpectralBands(  const vector<pair<float,float> >& bands )
{
   spectralBands = bands;
   spectralBandMinMaxs.clear();
   for( size_t i=0;i<spectralBands.size();i++ )
   {
      pair<float,float>& band = spectralBands[i];
      pair<double,double> minMax;
      minMax.first = chbaroverk/max( band.first, band.second );
      minMax.second = chbaroverk/min( band.first, band.second );
      spectralBandMinMaxs.push_back( minMax );
   }
}
*/


void BlackBodyEmission::emission( const float T, Color& value ) const
{
   if( T < 1667.0 || T > 25000.0 )
   {
      value[0] = value[1] = value[2] = value[3] = 0.0;
      return;
   }

   // Plankian locus approximation based on 
   // http://en.wikipedia.org/wiki/Planckian_locus
   float T1 = 1000.0/T;
   float T2 = T1*T1;
   float T3 = T2*T1;

   float CIEx = 0;
   float CIEy = 0;

   if( T <= 4000.0 )
   {
      CIEx = 0.179910 + 0.8776965*T1 - 0.2343580*T2 - 0.2661239*T3;
   }
   else
   {
      CIEx = 0.24039 + 0.2226347*T1 + 2.1070379*T2 - 3.0258469*T3;
   }

   if( T <= 2222.0 )
   {
      CIEy = -0.20219683 + 2.18555832*CIEx - 1.3481102*CIEx*CIEx - 1.1063814*CIEx*CIEx*CIEx;
   }
   else if( T <= 4000.0 )
   {
      CIEy = -0.16748867 + 2.09137015*CIEx - 1.37418593*CIEx*CIEx - 0.9549476*CIEx*CIEx*CIEx;
   }
   else if( T <= 25000.0 )
   {
      CIEy = -0.37001483 + 3.75112997*CIEx - 5.8733867*CIEx*CIEx + 3.081758*CIEx*CIEx*CIEx;
   }
   float CIEz = 1.0 - CIEx - CIEy;

   // Conversion to sRGB according to 
   // http://en.wikipedia.org/wiki/SRGB
   if( CIEy > 0.0 )
   {
      float Y = 1.0;
      float X = Y*CIEx/CIEy;
      float Z = Y*CIEz/CIEy;
      value[0] = 3.2406*X - 1.5372*Y - 0.4986*Z;
      value[1] = -0.9689*X + 1.8758*Y + 0.0415*Z;
      value[2] = 0.0557*X - 0.2040*Y + 1.0507*Z;
   }
   

   if( value[0] < 0 ){ value[0] = 0; }
   if( value[1] < 0 ){ value[1] = 0; }
   if( value[2] < 0 ){ value[2] = 0; }
}



const Color BlackBodyEmission::emission( const float T ) const
{
   Color result;
   emission( T, result );
   return result;
}
