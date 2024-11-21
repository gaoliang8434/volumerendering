
#ifndef __VOLUMELIGHT_H__
#define __VOLUMELIGHT_H__


#include "Vector.h"
#include "Color.h"
#include "AARectangle.h"
#include "UniformPRN.h"
#include "Volume.h"


/*
 *
 *  This is a stochastically sampled volumetric light.  getLight() randomly pickes a point within the bounds
 *  and monte carlo accepts/rejects based on probability field and threshold. If rejects, it iterates up to
 *  the sample limit searching for a point.
 *
*/


namespace lux
{

using namespace std;

class VolumeLight
{
  public:

    VolumeLight( const AABB& bnds, const ScalarField& prob, const ColorField& cd, const float threshold );

   ~VolumeLight();

    void setSeed( const int seed );
    void setSampleLimit( const int lim );
    void setThreshold( const float thresh );
//    void setBounds( const AABB& bnds );

    const int& getSampleLimit() const;
    const float& getThreshold() const;
    const AABB& getBounds() const;

    bool getLight( Vector& pos, Color& cd );

    std::string typelabel() { return "Volume Light"; }
    std::string documentation() { return "No documentation"; }

  private:

    const AABB bounds;
    const ScalarField probabilityField;
    const ColorField cdField;
    float probability_threshold;
    int sample_limit;
    UniformPRN prn;

};

typedef std::shared_ptr<VolumeLight> VLBase;

class VolumeLightField : public VLBase
{
  public:

    VolumeLightField();
    VolumeLightField( VolumeLight* f );
   ~VolumeLightField();

     char* __str__() 
     {
       static char typeLabel[2048];
       std::string lbl = (*this)->typelabel();
       size_t lbllength = lbl.size();
       if( lbllength > 2047 ){ lbllength = 2047; }
       lbllength = lbl.copy( typeLabel, lbllength);
       typeLabel[lbllength] = '\0';
       return typeLabel;
    }


     char* __doc__() 
     {
       static char docLabel[2048];
       std::string lbl = (*this)->documentation();
       size_t lbllength = lbl.size();
       if( lbllength > 2047 ){ lbllength = 2047; }
       lbllength = lbl.copy( docLabel, lbllength);
       docLabel[lbllength] = '\0';
       return docLabel;
    }

};

VolumeLightField makeVolumeLight( const AABB& bnds, const ScalarField& prob, const ColorField& cd, const float threshold );
void SetSampleLimit( VolumeLightField& vlf, int value );
void SetThreshold( VolumeLightField& vlf, float value );
void SetSeed( VolumeLightField& vlf, int value );


}

#endif
