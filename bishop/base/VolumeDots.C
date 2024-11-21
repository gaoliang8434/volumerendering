
#include "VolumeDots.h"
#include "Color.h"
#include "ImplicitVolumeShapes.h"

#include <iostream>
#include <cstdlib>

using namespace std;
using namespace lux;


const Vector RandomUnitVector()
{
   Vector value( drand48()-0.5, drand48()-0.5, drand48()-0.5 );
   value.normalize();
   return value;
}

void lux::putbumps(Vector gpos, float gradius, int nb, float scale, int rec, VolumeDots& vd)
{
int p, trackrec;
float radius;
Vector pos;

/*  track the number of recursions  */
trackrec = rec-1;

/* loop over particles to clump onto surface of parent  */
p = 0;
while(p < nb)
{
/*  place the particle on the surface of its parent  */
pos = RandomUnitVector();
pos = gpos + gradius * pos;
vd.setAtt("position", pos);

/*  give the particle a reduced size  */
radius = gradius * scale;
vd.setAtt("radiusPP", radius);

vd.emit();

/*  now recurse to put particles on the surface of this particle  */
if(trackrec > 0){putbumps(pos, radius, nb, scale, trackrec, vd);}
p = p + 1;
}
return;
}

void lux::putbumps(Vector gpos, float gradius, int nb, float scale, int rec, VolumeBlobs& vd)
{
int p, trackrec;
float radius;
Vector pos;

/*  track the number of recursions  */
trackrec = rec-1;

/* loop over particles to clump onto surface of parent  */
p = 0;
while(p < nb)
{
/*  place the particle on the surface of its parent  */
pos = RandomUnitVector();
pos = gpos + gradius * pos;
vd.setAtt("position", pos);

/*  give the particle a reduced size  */
radius = gradius * scale;
vd.setAtt("radiusPP", radius);

vd.emit();

/*  now recurse to put particles on the surface of this particle  */
if(trackrec > 0){putbumps(pos, radius, nb, scale, trackrec, vd);}
p = p + 1;
}
return;
}



int lux::setChildParticle_Cauliflower(  Vector guideposition, float childradius, int nbcauliflowerclumps, float cauliflowerscale, int nbcauliflowerrecursions, VolumeDots& vd )
{
vd.emit(); /* emit guide particle */
/*  recursively places and emits particles  */
putbumps(guideposition, childradius, nbcauliflowerclumps,
      cauliflowerscale, nbcauliflowerrecursions, vd);
return 0;  /* dont emit because emission occured during recursion  */
}


int lux::setChildParticle_Cauliflower(  Vector guideposition, float childradius, int nbcauliflowerclumps, float cauliflowerscale, int nbcauliflowerrecursions, VolumeBlobs& vd )
{
vd.emit(); /* emit guide particle */
/*  recursively places and emits particles  */
putbumps(guideposition, childradius, nbcauliflowerclumps,
      cauliflowerscale, nbcauliflowerrecursions, vd);
return 0;  /* dont emit because emission occured during recursion  */
}




void VolumeDots::emit()
{
   Color white;
   white[0] = dotAttributes["color"][0];
   white[1] = dotAttributes["color"][1];
   white[2] = dotAttributes["color"][2];
   white[3] = 0.0;

   //cout << "Emit color " << white[0] << " " << white[1] << " " << white[2] << endl;


   Vector center = dotAttributes["position"];
   Vector radius = dotAttributes["radiusPP"];

   Vector LL = center - 2.0*radius;
   Vector UL = center + 2.0*radius;

   float skinthickness = 2.0* volume->dx();

   VolumeFloatPtr emitted = VolumeFloatPtr( new ClampVolume( new MultiplyVolume( new SphereVolume( center, radius[0] ), 1.0/skinthickness ), 0.0, 1.0 ) );

   int ixmin, ixmax, iymin, iymax, izmin, izmax;

   if( volume->getBox( LL, UL, ixmin, ixmax, iymin, iymax, izmin, izmax ) )
   {

   for( int iz=izmin;iz<=izmax;iz++ )
   {
      for( int iy=iymin;iy<=iymax;iy++ )
      {
         for( int ix=ixmin;ix<=ixmax;ix++ )
	 {
            float vvalue = emitted->eval(  volume->evalP(ix,iy,iz)  );
	    if( vvalue > 0 )
	    {
	       float olddensity  = volume->value(ix,iy,iz);
	       if( vvalue > olddensity )
	       {
	          volume->value(ix,iy,iz) = vvalue;
	          color->value(ix,iy,iz)  = white;
	       }
	       //float compdensity = olddensity * ( 1-vvalue ) + vvalue;
	       //Color compcolor   = ( color->value(ix,iy,iz) * olddensity  * (1-vvalue) + white * vvalue ) / compdensity;
	       //volume->value(ix,iy,iz) = compdensity;
	       //color->value(ix,iy,iz)  = compcolor;
	    }
	 }
      }
   }
   }

}



void VolumeBlobs::emit()
{
   Color white;
   white[0] = dotAttributes["color"][0];
   white[1] = dotAttributes["color"][1];
   white[2] = dotAttributes["color"][2];
   white[3] = 0.0;


   cout << "Emit color " << white[0] << " " << white[1] << " " << white[2] << endl;

   Vector center = dotAttributes["position"];
   Vector radius = dotAttributes["radiusPP"];


   Volume<float>* blob = new ScaleVolume( new TranslateVolume( pattern, center) , radius );
   volume = new BlinnBlendVolume( volume, blob );
   //volume = new UnionVolume( volume, blob );

}








void lux::Zero( VolumeGrid<float>& g )
{
   size_t n = g.nx() * g.ny() * g.nz();
   float* data = g.rawPtr();
   float defaultvalue = 0.0;
   for( size_t i=0;i<n;i++ ){ data[i] = defaultvalue; }
}
void lux::Zero( VolumeGrid<Color>& g )
{
   size_t n = g.nx() * g.ny() * g.nz();
   Color* data = g.rawPtr();
   Color defaultvalue(0,0,0,0);
   for( size_t i=0;i<n;i++ ){ data[i] = defaultvalue; }
}
