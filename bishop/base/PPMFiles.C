

#include "ImageFiles.h"

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
using namespace lux;

float imagePlanePPMValue( float v, float dG, float dB )
{
   if( v == 0.0 ){ return 0.0; }
   return pow(v, dG) * dB;
}

void lux::writePPMImage(CmdLineFind& clf, Image& img )
{
   string filename = clf.find( "-ppmname", "", "Name of ppm file to write");
   if( filename == "" )
   {
      cout << "You need to give -ppmname the name of a file to write" << endl << flush;
      return;
   }
   int blackpoint = clf.find("-ppmblackpoint", 0, "Blackpoint");

   float displayBrightness = clf.find( "-ppmbrightness", 1.0f, "Scale brightness of image" );
   float displayGamma = clf.find( "-ppmgamma", 1.0f, "Gamma of image" );
   int ppmScale = 65535;
   float rgba[4];

   ofstream imagefile( filename.c_str() );
   imagefile << "P3 " << img.Width() << " " << img.Height() << " " << ppmScale << "\n";
   for( int j=0;j<img.Height();j++ )
   {
      for( int i=0;i<img.Width();i++ )
      {
         vector<float> pixel = img.pixel(i,img.Height()-j-1);
	 for( size_t c=0;c<3;c++ )
         {
	    rgba[c] = imagePlanePPMValue( pixel[c], displayGamma, displayBrightness ) * ppmScale;
	    if( rgba[c] > (float)ppmScale ){ rgba[c] = (float)ppmScale; }
	    if( rgba[c] < (float)blackpoint ){ rgba[c] = (float)blackpoint; }
         }
         imagefile << (int)rgba[0] << " " << (int)rgba[1] << " " << (int)rgba[2] << endl;
      }
   }




   imagefile.close();
   cout << "File " << filename << " written.\n";
}
