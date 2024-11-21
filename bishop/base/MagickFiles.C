


#include "ImageFiles.h"

#include <iostream>
#include <cmath>
#include <Magick++.h>
#include <list>
#include <ctime>
#include "LABLogo.h"

using namespace std;
using namespace lux;


float imagePlaneMagickValue( float v, float dG, float dB )
{
   if( v == 0.0 ){ return 0.0; }
   return pow(v, dG) * dB;
}

template <typename T> 
string tostr(const T& t) { std::ostringstream os; os<<t; return os.str(); }


void lux::writeMagickImage(CmdLineFind& clf, lux::Image& img)
{
   string filenameBase = clf.find( "-magickname", "", "Name of file to write via ImageMagick");
   if( filenameBase == "" )
   {
      cout << "You need to give -magickname the name of a file to write" << endl << flush;
      return;
   }
   vector<float> displayBrightnessSet = clf.findMultiple( "-magickbrightness", 1.0f, "Scale brightness of image" );
   vector<float> displayGammaSet = clf.findMultiple( "-magickgamma", 1.0f, "Gamma of image" );
   int blackpoint = clf.find("-magickblackpoint", 0, "Blackpoint");
   vector<string> tags = clf.findMultiple("-magicklabel", string(""), "Labels to embed in the image" );
   bool showAllOptions = clf.findFlag("-magickshowoptions", "Print all of the option values on the image");
   bool noOwnerAndFile = clf.findFlag("-magicknoowner", "suppress printing the owner and filename");
   bool noText = clf.findFlag("-magicknotext");

   for( size_t b=0;b<displayBrightnessSet.size();b++ )
   {
      float displayBrightness = displayBrightnessSet[b];
      for( size_t g=0;g<displayGammaSet.size();g++ )
      {
         float displayGamma = displayGammaSet[g];
         string filename =  filenameBase;
	 if( displayGammaSet.size() > 1 )
	 {
            filename =  "G" + tostr(displayGamma) + filename;
	 }
	 if( displayBrightnessSet.size() > 1 )
	 {
            filename = "B" + tostr(displayBrightness) + filename;
	 }

   	 Magick::Image image( Magick::Geometry(img.Width(), img.Height()), Magick::Color("white"));


   // fill image with the contents of img

   for( int j=0;j<img.Height();j++ )
   {
      for( int i=0;i<img.Width();i++ )
      {
         vector<float> pix = img.pixel(i,img.Height()-j-1);
	 for( size_t c=0;c<3;c++ )
         {
            pix[c] = imagePlaneMagickValue( pix[c], displayGamma, displayBrightness ) * 65535;
	    if( pix[c] > 65535.0 ){ pix[c] = 65535.0; }
	    if( pix[c] < (float)blackpoint ){ pix[c] = (float)blackpoint; }
         }

         image.pixelColor(i,j,Magick::Color( (int)pix[0],(int)pix[1],(int)pix[2],(int)pix[3] ) );
      }
   }

   string user = "tessendorf";
   time_t seconds = time(NULL);
   string timestamp = ctime( &seconds );
   string stampid = user + "  " + timestamp;

   if( !noText )
   {
   // embed label if one given
   list<Magick::Drawable> text_draw_list;
   int x = 50;
   int y = 25;
   int textspace = img.Height() - 150;
   if( !tags.empty() )
   {
      cout << "Embedding text in image:\n";
      // set the text to be drawn at specified position:
      for( size_t i=0;i<tags.size();i++ )
      {
         cout << "\t" << tags[i] << endl;
         text_draw_list.push_back( Magick::DrawableText(x, y, tags[i].c_str()));
	 y += 20;
	 if( y > textspace )
	 {
	    y =  25;
	    x += img.Width()/3;
	 }
      }
   }

   if( showAllOptions )
   {
      vector<string> finds = clf.listFinds();
      for( size_t i=0;i<finds.size();i++ )
      {
         text_draw_list.push_back( Magick::DrawableText(x, y, finds[i].c_str()));
	 y += 20;
	 if( y > textspace )
	 {
	    y = 25;
	    x += img.Width()/3;
	 }
      }
   }

   
   if( !noOwnerAndFile )
   {
      text_draw_list.push_back( Magick::DrawableText(50, img.Height()-15, stampid.c_str()));
      text_draw_list.push_back( Magick::DrawableText(img.Width() - 6*filename.size() - 50, img.Height()-15, filename.c_str()));
      text_draw_list.push_back( Magick::DrawableText(img.Width()/2 - 3, img.Height()-15, LABSHORTLOGO));
      text_draw_list.push_back( Magick::DrawableStrokeColor(Magick::Color("grey")));
      text_draw_list.push_back( Magick::DrawableFillColor(Magick::Color("white") ) );
   }
   if( !text_draw_list.empty() ){ image.draw( text_draw_list ); }
   }

   image.write(filename.c_str());
   cout << "File " << filename <<  " written by user " << user <<  " at " << timestamp <<endl;
   }
   }
}





void lux::writeMagickImage(CmdLineFind& clf, lux::Frame& img)
{
   string filename = clf.find( "-magickname", "", "Name of file to write via ImageMagick");
   if( filename == "" )
   {
      cout << "You need to give -magickname the name of a file to write" << endl << flush;
      return;
   }
   int blackpoint = clf.find("-magickblackpoint", 0, "Blackpoint");
   vector<string> tags = clf.findMultiple("-magicklabel", string(""), "Labels to embed in the image" );
   bool showAllOptions = clf.findFlag("-magickshowoptions", "Print all of the option values on the image");
   bool noOwnerAndFile = clf.findFlag("-magicknoowner", "suppress printing the owner and filename");
   bool noText = clf.findFlag("-magicknotext");

   Magick::Image image( Magick::Geometry(img.width, img.height), Magick::Color("white"));


   // fill image with the contents of img

   float pix[4];
   pix[3] = 1.0;
   long index = 0;
   for( int j=0;j<img.height;j++ )
   {
      for( int i=0;i<img.width;i++ )
      {
	 for( size_t c=0;c<3;c++ )
         {
            pix[c] = imagePlaneMagickValue( img.data[index], img.gamma, img.contrast ) * 65535;
	    if( pix[c] > 65535.0 ){ pix[c] = 65535.0; }
	    if( pix[c] < (float)blackpoint ){ pix[c] = (float)blackpoint; }
	    index++;
         }

         image.pixelColor(i,img.height - j - 1,Magick::Color( (int)pix[0],(int)pix[1],(int)pix[2],(int)pix[3] ) );
      }
   }

   string user = "tessendorf";
   time_t seconds = time(NULL);
   string timestamp = ctime( &seconds );
   string stampid = user + "  " + timestamp;

   if( !noText )
   {
   // embed label if one given
   list<Magick::Drawable> text_draw_list;
   int x = 50;
   int y = 25;
   int textspace = img.height - 150;
   if( !tags.empty() )
   {
      cout << "Embedding text in image:\n";
      // set the text to be drawn at specified position:
      for( size_t i=0;i<tags.size();i++ )
      {
         cout << "\t" << tags[i] << endl;
         text_draw_list.push_back( Magick::DrawableText(x, y, tags[i].c_str()));
	 y += 20;
	 if( y > textspace )
	 {
	    y =  25;
	    x += img.width/3;
	 }
      }
   }

   if( showAllOptions )
   {
      vector<string> finds = clf.listFinds();
      for( size_t i=0;i<finds.size();i++ )
      {
         text_draw_list.push_back( Magick::DrawableText(x, y, finds[i].c_str()));
	 y += 20;
	 if( y > textspace )
	 {
	    y = 25;
	    x += img.width/3;
	 }
      }
   }

   
   if( !noOwnerAndFile )
   {
      text_draw_list.push_back( Magick::DrawableText(50, img.height-15, stampid.c_str()));
      text_draw_list.push_back( Magick::DrawableText(img.width - 6*filename.size() - 50, img.height-15, filename.c_str()));
      text_draw_list.push_back( Magick::DrawableText(img.width/2 - 3, img.height-15, LABSHORTLOGO));
      text_draw_list.push_back( Magick::DrawableStrokeColor(Magick::Color("grey")));
      text_draw_list.push_back( Magick::DrawableFillColor(Magick::Color("white") ) );
   }
   if( !text_draw_list.empty() ){ image.draw( text_draw_list ); }
   }

   image.write(filename.c_str());
   cout << "File " << filename <<  " written by user " << user <<  " at " << timestamp <<endl;
}
