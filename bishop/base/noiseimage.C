
#include "Noise.h"
#include "PerlinNoise.h"
#include "FFTNoise.h"
#include <iostream>
#include <cmath>
using namespace std;
#include "CmdLineFind.h"
#include "ImageFiles.h"
using namespace lux;


int main( int argc, char** argv )
{
   CmdLineFind clf( argc, argv );


   Noise_t parms;
   parms.seed = clf.find("-seed", 485758, "Seed for the PRN");
   parms.octaves = clf.find("-octaves", 1.0f );
   parms.wavelength = clf.find("-wavelength", 1.0f );
   parms.fjump = clf.find("-fjump", 2.0f );
   parms.roughness = clf.find("-roughness", 0.5f );
   parms.time = clf.find("-time", 0.0f );
   clf.usage("-h");

   FractalSum<PerlinNoiseGustavson> perlin;
   perlin.setParameters( parms );
   FFTNoise fft;
   fft.setParameters( parms );




   Noise* noise = &perlin;
   if( clf.findFlag("-fft") ){ noise = &fft; }


    int imageWidth  = clf.find( "-NX", 800, "Image width");
    int imageHeight = clf.find( "-NY", 800, "Image height");

    float imageRangeX = clf.find("-rangeX", 10.0f );
    float imageRangeY = clf.find("-rangeY", 10.0f );


    Image image;
    image.reset( imageWidth, imageHeight );

   float scale = pow( 1.0+parms.roughness, parms.octaves-1.0);
   for( int j=0;j<imageHeight;j++ )
   {
      float y = j*(float)imageRangeY/(float)(imageHeight-1);
      for( int i=0;i<imageWidth;i++ )
      {
         float x = i*(float)imageRangeX/(float)(imageWidth-1);
	 Vector P(x,y,0.1);
	 float val = (noise->eval(P) + 0.5*scale )/scale; 
         image.value(i,j,0) = val;
         image.value(i,j,1) = val;
         image.value(i,j,2) = val;
         image.value(i,j,3) = 1.0;
      }
   }

   writeMagickImage( clf, image );

}
