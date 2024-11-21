


#include "VolumeLight.h"
#include "Noise.h"


using namespace lux;
using namespace std;




VolumeLight::VolumeLight( const AABB& bnds, const ScalarField& prob, const ColorField& cd, const float threshold ) :
   bounds (bnds),
   probabilityField (prob),
   cdField (cd),
   probability_threshold (threshold),
   sample_limit (100000)
   {
      Noise_t parm;
      parm.seed = 7758;
      prn.setParameters( parm );
   }

VolumeLight::~VolumeLight(){}

void VolumeLight::setSeed( const int seed ) 
{
   Noise_t parm;
   parm.seed = seed; 
   prn.setParameters( parm ); 
}

void VolumeLight::setSampleLimit( const int lim ) { sample_limit = lim; }
void VolumeLight::setThreshold( const float thresh ) { probability_threshold = thresh; }
//void VolumeLight::setBounds( const AABB& bnds ) { bounds = bnds; }

const int& lux::VolumeLight::getSampleLimit() const { return sample_limit; }
const float& lux::VolumeLight::getThreshold() const { return probability_threshold; }
const AABB& lux::VolumeLight::getBounds() const { return bounds; }






bool lux::VolumeLight::getLight( Vector& pos, Color& cd ) 
{
   const Vector dim = bounds->urc() - bounds->llc();
   for( int i=0;i<sample_limit;i++ )
   {
      pos = bounds->llc() + Vector( dim.X() * prn.eval(), dim.Y() * prn.eval(), dim.Z() * prn.eval() );
      float prob = probabilityField->eval(pos);
      if( prob >= prn.eval() * probability_threshold )
      {
         cd = cdField->eval(pos);
         return true;
      }
   }
   cout << "No volume lights found\n";
   return false;
}

VolumeLightField::VolumeLightField() :  std::shared_ptr<VolumeLight>() {}
VolumeLightField::VolumeLightField( VolumeLight* f ) :  std::shared_ptr<VolumeLight>( f ) {}
VolumeLightField::~VolumeLightField(){}


VolumeLightField lux::makeVolumeLight( const AABB& bnds, const ScalarField& prob, const ColorField& cd, const float threshold )
{
   return VolumeLightField( new VolumeLight( bnds, prob, cd, threshold ) );
}

void lux::SetSampleLimit( VolumeLightField& vlf, int value ){ vlf->setSampleLimit(value); }
void lux::SetThreshold( VolumeLightField& vlf, float value ){ vlf->setThreshold(value); }
void lux::SetSeed( VolumeLightField& vlf, int value ){ vlf->setSeed(value); }

