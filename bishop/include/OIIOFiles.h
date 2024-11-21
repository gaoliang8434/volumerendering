

#ifndef __OIIOFILES_H__
#define __OIIOFILES_H__

#include "Image.h"
#include "CmdLineFind.h"


#include <map>
#include <string>

namespace lux{


void writeOIIOImage(lux::CmdLineFind& clf, Image& img);
void writeOIIOImage( const char* fname, Image& img, float brightness = 1.0, float gamma = 1.0 );



void writeOIIOImage( const char* fname, Image& img, const std::map<std::string,std::string>& labels, float brightness = 1.0, float gamma = 1.0 );
void writeOIIOImage( const char* fname, Image& img, const std::vector<std::string>& keys, const std::vector<std::string>& tags,  float brightness = 1.0, float gamma = 1.0 );
void readOIIOImage( const char* fname, Image& img );
void readOIIOImage( const char* fname, Image& img, std::map<std::string,std::string>& labels );
void readOIIOImage( const char* fname, Image& img, std::vector<std::string>& keys, std::vector<std::string>& tags );
void printMetadata( const std::map<std::string,std::string>& meta );




//! Stores image channels provided in images into an OpenEXR file.  Each
//! image should contain channel data in channel 0.  Each image name
//! should be the name of the corresponding channel.  For example, use
//! "R" for the red channel.  EXR filename should be stored in filename.
void writeDeepImage(std::map<std::string, Image *> &images, const char *filename);




//! Converts a channel of an image into a float array, for use with some other tools in python
float* convert( const Image& img, int channel );


}

#endif
