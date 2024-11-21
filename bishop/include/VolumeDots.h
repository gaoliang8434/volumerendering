

#ifndef __VOLUMEDOTS_H__
#define __VOLUMEDOTS_H__

#include "Vector.h"
#include "Volume.h"
#include <map>
#include <string>
#include "VolumeGrid.h"
#include "ImplicitVolumeShapes.h"



namespace lux
{

class VolumeDots
{
   public:

     VolumeDots( VolumeGrid<float>* g, VolumeGrid<Color>* gc ) : volume (g), color(gc) {}
     ~VolumeDots(){}

    typedef std::map<std::string,Vector> DotAttributes;



    void setAtt( const std::string& name, const Vector& value )
    {
       dotAttributes[name] = value;
    }

    void setAtt( const std::string& name, const float& value )
    {
       dotAttributes[name] = Vector(value, value, value);
    }

    void emit();

  private:

    DotAttributes dotAttributes;
    VolumeGrid<float>* volume;
    VolumeGrid<Color>* color;

};


class VolumeBlobs
{
   public:

     VolumeBlobs( Volume<float> * pat, float t ) : 
        volume ( new Volume<float> () ), 
	color  ( new Volume<Color>() ), 
	pattern(pat),
	thick (t)
     {}
     ~VolumeBlobs(){}

    typedef std::map<std::string,Vector> DotAttributes;

    void setAtt( const std::string& name, const Vector& value )
    {
       dotAttributes[name] = value;
    }

    void setAtt( const std::string& name, const float& value )
    {
       dotAttributes[name] = Vector(value, value, value);
    }
    

    void emit();

    VolumeFloatPtr EmittedDensity(){ return VolumeFloatPtr( new ClampVolume( new MultiplyVolume( volume, 1.0/thick ), 0.0, 1.0 ) ); }
    VolumeColorPtr  EmittedColor(){ return VolumeColorPtr(color); }

  private:

    DotAttributes dotAttributes;
    Volume<float> * volume;
    Volume<Color>* color;
    Volume<float> * pattern;
    float thick;

};





/*  cauliflower distribution  */
//float cauliflowerscale;
//int nbcauliflowerclumps, nbcauliflowerrecursions;
//The parameter cauliflowerscale controls the relative size of spheres in each recursion level. At each level, nbcauliflowerclumps spheres are placed on the surface of each sphere from the previous level. Finally, the total number of recursion levels is nbcauliflowerrecursions. These three parameters are used to recursively create the cauliflower. The approach is:


void putbumps(Vector gpos, float gradius, int nb, float scale, int rec, VolumeDots& vd);
void putbumps(Vector gpos, float gradius, int nb, float scale, int rec, VolumeBlobs& vd);



int setChildParticle_Cauliflower(  Vector guideposition, float childradius, int nbcauliflowerclumps, float cauliflowerscale, int nbcauliflowerrecursions, VolumeDots& vd );
int setChildParticle_Cauliflower(  Vector guideposition, float childradius, int nbcauliflowerclumps, float cauliflowerscale, int nbcauliflowerrecursions, VolumeBlobs& vd );
//Several new things are going on in this algorith. First, because of the recursion, all of the emit() calls take place inside setChildParticle_Cauliflower(), and its return value is zero so that no additional emit() calls are made. The second new technique is recursion. With each recursive call to putbumps(), the next level of spheres are placed, with a sphere radius that is larger or smaller than the previous level by a factor of cauliflowerscale. By setting this parameter less than one, the bumps get smaller and smaller. In the example images above, the values used were:

/* cauliflower placement  */
//cauliflowerscale = 0.5;
//nbcauliflowerclumps = 20;
//nbcauliflowerrecursions = 3;


void Zero( VolumeGrid<float>& g );
void Zero( VolumeGrid<Color>& g );

}

#endif
