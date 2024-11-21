
#include "OIIOFiles.h"

#include <iostream>
#include <cmath>
#include <OpenImageIO/imageio.h> 
#include <list>
#include <ctime>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace std;
using namespace lux;
OIIO_NAMESPACE_USING


float imagePlaneValue( float v, float dG, float dB )
{
   if( v == 0.0 ){ return 0.0; }
   return pow(v, dG) * dB;
}


void lux::writeOIIOImage(CmdLineFind& clf, Image& img)
{
   string filename = clf.find( "-oiioname", "", "Name of file to write via oiio");
   if( filename == "" )
   {
      cout << "You need to give -oiioname the name of a file to write" << endl << flush;
      return;
   }
   float displayBrightness = clf.find( "-oiiobrightness", 1.0f, "Scale brightness of image" );
   float displayGamma = clf.find( "-oiiogamma", 1.0f, "Gamma of image" );

   float* imagedata = new float[ img.Width()* img.Height() * 3 ];

   // fill image with the contents of img

   long index = 0;
   for( int j=0;j<img.Height();j++ )
   {
      for( int i=0;i<img.Width();i++ )
      {
         vector<float> pix = img.pixel(i,img.Height() - j - 1);
	 for( size_t c=0;c<3;c++ )
         {
            pix[c] = imagePlaneValue( pix[c], displayGamma, displayBrightness );
	    imagedata[index++] = pix[c];
         }
      }
   }

   ImageOutput *out = ImageOutput::create (filename.c_str()); 
   if( !out )
   {
      cout << "Not able to write an image to file " << filename << endl;
   }
   else
   {
      ImageSpec spec (img.Width(), img.Height(), 3, TypeDesc::HALF); 
      spec.attribute("user", "lux");
      spec.attribute("writer", "OIIOFiles" );
      out->open (filename.c_str(), spec);
      out->write_image (TypeDesc::FLOAT, imagedata); 
      out->close (); 
      delete out;
   }
   delete[] imagedata;
}



void lux::writeOIIOImage( const char* fname, Image& img, float displayBrightness, float displayGamma  )
{
   float* imagedata = new float[ img.Width()* img.Height() * img.Depth() ];

   // fill image with the contents of img

   long index = 0;
   for( int j=0;j<img.Height();j++ )
   {
      for( int i=0;i<img.Width();i++ )
      {
         vector<float> pix = img.pixel(i,img.Height() - j - 1);
	 for( size_t c=0;c<(size_t)img.Depth();c++ )
         {
            pix[c] = imagePlaneValue( pix[c], displayGamma, displayBrightness );
	    imagedata[index++] = pix[c];
         }
      }
   }

   ImageOutput *out = ImageOutput::create (fname); 
   if( !out )
   {
      cout << "Not able to write an image to file " << fname << endl;
   }
   else
   {
      ImageSpec spec (img.Width(), img.Height(), img.Depth(), TypeDesc::FLOAT); 
      spec.attribute("user", "lux");
      spec.attribute("writer", "OIIOFiles" );
      out->open (fname, spec);
      out->write_image (TypeDesc::FLOAT, imagedata);
      out->close (); 
      delete out;
   }
   delete[] imagedata;
}



void lux::writeDeepImage(map<string, Image *> &images, const char *filename)
{
	ImageOutput *out = ImageOutput::create(filename);
	if (!out)
	{
		cout << "Not able to write images to file ";
		cout << filename << endl;
		return;
	}

	if (strcmp(out->format_name(), "openexr") != 0)
	{
		cout << "DeepImage format should be OpenEXR" << endl;
		cout << "Format provided: ";
		cout << out->format_name() << endl;
		return;
	}

	int num_channels = images.size();
	map<string, Image *>::iterator it;
	int width = images.begin()->second->Width();
	int height = images.begin()->second->Height();
	ImageSpec spec(width, height, num_channels, TypeDesc::HALF);
	float *imagedata = new float[width * height * num_channels];
	long long index = 0;
	spec.channelnames.clear();
	int cur_channel = 0;
	for (it = images.begin(); it != images.end(); it++)
	{
		index = 0;
		string fname = it->first;
		Image *img = it->second;
		spec.channelnames.push_back(fname);
		for (int j = 0; j < height; ++j)
		{
			for (int i = 0; i < width; ++i)
			{
				vector<float> pix = img->pixel(i, height - j - 1);
				/* Last two arguments should be displayGamma
				   and displayBrightness.  Ask about adding
				   in optional map of filename to gamma and
				   filename to brightness for unique values.
				*/
				pix[0] = imagePlaneValue(pix[0], 1.0, 1.0);
				long long ind = num_channels * index;
				ind += cur_channel;
				imagedata[ind] = pix[0];
				index++;
			}
		}
		cur_channel += 1;
	} // loop over images
	spec.attribute("user", "lux");
	spec.attribute("writer", "OIIOFiles");
	out->open(filename, spec, ImageOutput::Create);
	out->write_image(TypeDesc::FLOAT, imagedata);
	out->close();
	
	delete imagedata;
	delete out;
}




void lux::writeOIIOImage( const char* fname, Image& img, const vector<string>& keys, const vector<string>& tags,  float displayBrightness, float displayGamma )
{
   map<string,string> labels;
   size_t nblabels = keys.size();
   nblabels = ( nblabels < tags.size() )? nblabels : tags.size();
   for( size_t i=0;i<nblabels;i++ )
   {
      labels[ keys[i] ] = tags[i];
   }
   writeOIIOImage( fname, img, labels, displayBrightness, displayGamma );

}



void lux::writeOIIOImage( const char* fname, Image& img, const map<string,string>& labels, float displayBrightness, float displayGamma )
{
   float* imagedata = new float[ img.Width()* img.Height() * img.Depth() ];

   // fill image with the contents of img

   long index = 0;
   for( int j=0;j<img.Height();j++ )
   {
      for( int i=0;i<img.Width();i++ )
      {
         vector<float> pix = img.pixel(i,img.Height() - j - 1);
	 for( size_t c=0;c<(size_t)img.Depth();c++ )
         {
            pix[c] = imagePlaneValue( pix[c], displayGamma, displayBrightness );
	    imagedata[index++] = pix[c];
         }
      }
   }

   ImageOutput *out = ImageOutput::create (fname); 
   if( !out )
   {
      cout << "Not able to write an image to file " << fname << endl;
   }
   else
   {
      ImageSpec spec (img.Width(), img.Height(), img.Depth(), TypeDesc::FLOAT); 
      spec.attribute("user", "imageTools");
      spec.attribute("writer", "OIIOFiles" );
      if( labels.size() > 0 )
      {
         map<string,string>::const_iterator lab = labels.begin();
	 while( lab != labels.end() )
	 {
	    const string& name = lab->first;
	    const string& value = lab->second;
	    spec.attribute( name, value );
	    lab++;
	 }
      }
      out->open (fname, spec);
      out->write_image (TypeDesc::FLOAT, imagedata); 
      out->close (); 
      delete out;
   }
   delete[] imagedata;
}




void lux::readOIIOImage( const char* fname, Image& img  )
{
   int xres, yres, channels;
   ImageInput *in = ImageInput::create (fname);
   if (! in) {return;}
   ImageSpec spec;
   in->open (fname, spec);
   xres = spec.width;
   yres = spec.height;
   channels = spec.nchannels;
   float* pixels = new float[xres*yres*channels];
   in->read_image (TypeDesc::FLOAT, pixels);

   img.reset( xres, yres, channels );
   long index = 0;
   for( int j=0;j<yres;j++)
   {
      for( int i=0;i<xres;i++ )
      {
         for( int c=0;c<channels;c++ )
         {
	    img.value(i,img.Height() - j - 1,c) = pixels[index++];
         }
      }
   }
   in->close ();
   delete in;
   delete[] pixels;
}

void lux::readOIIOImage( const char* fname, Image& img, map<string,string>& labels )
{
   int xres, yres, channels;
   ImageInput *in = ImageInput::create (fname);
   if (! in) {return;}
   ImageSpec spec;
   in->open (fname, spec);
   xres = spec.width;
   yres = spec.height;
   channels = spec.nchannels;
   float* pixels = new float[xres*yres*channels];
   in->read_image (TypeDesc::FLOAT, pixels);

   img.reset( xres, yres, channels );
   long index = 0;
   for( int j=0;j<yres;j++)
   {
      for( int i=0;i<xres;i++ )
      {
         for( int c=0;c<channels;c++ )
         {
	    img.value(i,img.Height() - j - 1,c) = pixels[index++];
         }
      }
   }

   for( size_t i=0;i<spec.extra_attribs.size();i++)
   {
      const ParamValue& p = spec.extra_attribs[i];
      string name = p.name().c_str();
      string value = spec.metadata_val ( p, true);
      labels[name] = value;
   }

   in->close ();
   delete in;
   delete[] pixels;
}



void lux::readOIIOImage( const char* fname, Image& img, vector<string>& keys, vector<string>& tags )
{
   map<string,string> labels;
   readOIIOImage( fname, img, labels );
   map<string,string>::iterator p = labels.begin();
   keys.clear();
   tags.clear();
   while( p != labels.end() )
   {
      keys.push_back( p->first );
      tags.push_back( p->second );
      p++;
   }
}




void lux::printMetadata( const map<string,string>& meta )
{
   if( meta.empty() ){ return; }
   cout << "\n\nMetadata Labels\n==============================\n";
   map<string,string>::const_iterator p = meta.begin();
   size_t maxlength = 0;
   while( p != meta.end() )
   {
      maxlength = ( maxlength < p->first.size() ) ? p->first.size() : maxlength;
      p++;
   }
   maxlength += 1;
   p = meta.begin();
   while( p != meta.end() )
   {
      size_t extraspace = maxlength - p->first.size();
      cout << p->first;
      for( size_t i=0;i<extraspace;i++ ){ cout << " "; }
      cout << "---------> " << p->second << endl;
      p++;
   }
   cout << "\n==============================\n\n";

}


float* lux::convert( const Image& img, int channel )
{
    float* data = new float[img.Width() * img.Height() ];
    long index = 0;
    for( int y=0;y<img.Height();y++)
    {
       for( int x=0;x<img.Width();x++ )
       {
          data[index++] = img.value(x,y,channel);
       }
    }
    return data;
}

