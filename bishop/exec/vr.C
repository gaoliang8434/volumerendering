//===============================================================================
//
//  Ray marcher for rendering volumes
//
//===============================================================================





// Standard utilities and system includes, plus project specific items
//===============================================================================
#include <algorithm>
using namespace std;

#include "CmdLineFind.h"
#include "RenderPatchAllocator.h"
#include "RayMarcher.h"
#include "Image.h"
#include "Volume.h"
#include "LABLogo.h"
#include "ImplicitColors.h"
#include "ImageFiles.h"
#include "ProgressMeter.h"
#include "Logger.h"
#include "version.h"
#include "OIIOFiles.h"

using namespace std;
using namespace lux;




void SetUpCamera( CmdLineFind& clf, Camera& cam );
void SetUpData( CmdLineFind& clf, RenderData& data );
void SetUpVolume( CmdLineFind& clf, RenderData& d );
void SetUpDSM( CmdLineFind& clf, RenderData& d );
RenderPatchAllocator* SetUpAllocator( CmdLineFind& clf, int width, int height, const Camera& cam );


void RenderLoop(CmdLineFind& clf, RenderPatchAllocator* rp, Image& img, RenderData& d );
void RenderBackground(CmdLineFind& clf, RenderPatchAllocator* rp, Image& img, RenderData& d );
//void writePPMImage( CmdLineFind& clf, Image& img );

// Main program
//===============================================================================
int main(int argc, char** argv)
{
    cout << versionString() << endl;
    Logger logger( "vr", argc, argv );
    cout << LABLOGOSTRING << flush;

    ProgressMeter vrtimer(1, "TOTAL EXECUTION TIME");

    CmdLineFind clf( argc, argv );

    Camera camera;
    SetUpCamera( clf, camera );

    int imageWidth  = clf.find( "-NX", 1920/2, "Image width");
    int imageHeight = clf.find( "-NY", 1080/2, "Image height");

    Image volumeImage;
    volumeImage.reset( imageWidth, imageHeight );
  
    RenderPatchAllocator* rpa = SetUpAllocator( clf, imageWidth, imageHeight, camera );

    RenderData rd;
    SetUpData( clf, rd );

    RenderLoop(clf, rpa, volumeImage, rd);
    RenderBackground( clf, rpa, volumeImage, rd );
    delete rpa;

#ifndef VRPPM
    //writeMagickImage( clf, volumeImage );
    writeOIIOImage( clf, volumeImage );
#endif
#ifdef VRPPM
    writePPMImage( clf, volumeImage );
#endif
    clf.usage("-help");
    clf.printFinds();
}




// Volume Renderer functions

void SetUpData( CmdLineFind& clf, RenderData& data)
{
    data.scatterCoefficient = clf.find( "-scatter", 1.0f, "Scatter coefficient" );
    data.maxPathlength= clf.find( "-maxpathlength", 30.0f, "Maximum length of a ray march path" );
    data.ds = clf.find( "-ds", 0.1f, "Size of ray march step (negative means use cell size)" );

    vector<float> inputColor;
    inputColor.push_back(1);
    inputColor.push_back(1);
    inputColor.push_back(1);
    vector<float> fieldcolor = clf.findArray( "-color", inputColor, "Color of the density under scattering");
    data.colorField = new ConstantColor( Color( fieldcolor[0], fieldcolor[1], fieldcolor[2], 1) );

    inputColor[0] = 0;
    inputColor[1] = 0;
    inputColor[2] = 0;
    fieldcolor = clf.findArray( "-emission", inputColor, "Color of the density under emission");
    data.ambientColorField = new ConstantColor( Color( fieldcolor[0], fieldcolor[1], fieldcolor[2], 1) );
    data.sparseGrid = 0;

    SetUpVolume(clf, data);
    SetUpDSM(clf, data );
}

void SetUpCamera( CmdLineFind& clf, Camera& cam )
{
   Vector eye = cam.eye();
   Vector view = cam.view();
   Vector up = cam.up();

   vector<float> vvalue;
   vvalue.push_back( 0.0 );
   vvalue.push_back( 0.0 );
   vvalue.push_back( 15.0 );

   vvalue = clf.findArray( "-eye", vvalue, "Position of the camera");
   eye = Vector( vvalue[0], vvalue[1], vvalue[2] );

   int turntable = clf.find("-turntable", -1, "Using a 120 frame turntable animation and set the frame number");
   if( turntable > 0 && turntable <=120 )
   {
      float angle = (turntable-1)*M_PI*2.0/120.0;
      float radius = (eye-Vector(0,eye[1],0)).magnitude();
      eye[0] = cos(angle)*radius;
      eye[2] = sin(angle)*radius;
   }


   vvalue.clear();
   view = Vector(0,0,0);
   vvalue.push_back( 0 );
   vvalue.push_back( 0 );
   vvalue.push_back( 0 );

   vvalue = clf.findArray( "-view", vvalue, "Position of view spot");
   view = Vector( vvalue[0], vvalue[1], vvalue[2] );

   view -= eye;

   vvalue.clear();
   vvalue.push_back( up[0] );
   vvalue.push_back( up[1] );
   vvalue.push_back( up[2] );

   vvalue = clf.findArray( "-up", vvalue, "Up direction vector of the camera");
   up = Vector( vvalue[0], vvalue[1], vvalue[2] );
  
   cam.setEyeViewUp( eye, view, up );
   cam.setFov( clf.find( "-fov", 60.0f, "Field of view" ) );
   cam.setAspectRatio( clf.find( "-aspect", (float)cam.aspectRatio(), "Aspect ratio" ) );
   cam.setNearPlane( clf.find("-near", 0.0f, "Camera new plane") );
   cam.setFarPlane( clf.find("-far", 1.0e6f, "Camera far plane") );
}






RenderPatchAllocator* SetUpAllocator( CmdLineFind& clf, int width, int height, const Camera& cam )
{
   int nbpatches = 100;
   nbpatches = clf.find("-nbpatches", nbpatches, "Number of patches to divide up the image plane into");
   int nbRaysPerPixel = clf.find("-raysperpixel", 1, "Number of rays per pixel (antialiasing)" );
   int nbpixels = (width*height);
   if( nbpixels%nbpatches != 0 ){ --nbpatches; }
   int psize = nbpixels/nbpatches;
   RenderPatchAllocator* rpa = new RenderPatchAllocator( width, height, cam );
   rpa->setPatches( psize, "base", 0, 0, nbRaysPerPixel );
   return rpa;
}



void RenderLoop(CmdLineFind& clf, RenderPatchAllocator* rpa, Image& img, RenderData& data)
{
   if( clf.findFlag("-dontrender") > 0 ){ return; }
   RenderPatchAllocator::PixelPatch pix;
   vector<Color> output;

   vector<int> selectPatches = clf.findMultiple( "-patch", -1, "Select a patch to render");
   int count =  rpa->nbUncompletedPatches();
   if( !selectPatches.empty() ){ count = selectPatches.size(); }
   ProgressMeter metr( count , "Render" );

   while( rpa->nbUncompletedPatches() != 0 )
   {
      bool renderThisOne = true;
      if( !selectPatches.empty() )
      {
         vector<int>::iterator check = find( selectPatches.begin(), selectPatches.end(), (int)( rpa->nbUncompletedPatches() ) );
	 renderThisOne = (check != selectPatches.end() );
      }
      data.startPosition.clear();
      data.startDirection.clear();
      rpa->popPatch( data.startPosition, data.startDirection, pix );
      if( output.size() != data.startPosition.size() )
      {
         output.clear();
	 output.resize( data.startPosition.size() );
      }

      if( renderThisOne )
      { 
         ssRayMarchAccumulation( data, output ); 

	 int nb = rpa->getRaysPerPixel();
         for( size_t i=0;i<pix.size();i++ )
         {
            vector<float>& pixel = img.pixel( pix[i] );
            for( size_t c=0;c<pixel.size();c++ ) { pixel[c] += output[i][c]/nb; }
         }
      }
      metr.update();
   }
}


void RenderBackground(CmdLineFind& clf, RenderPatchAllocator* rp, Image& img, RenderData& d )
{
/*
   bool renderbg = clf.findFlag("-renderbackground");
   if( !renderbg ) { return; }

   Vector backgroundTopColor = clf.find("-topcolor", Vector( 0.1, 0.1, 1.0 ) );
   Vector backgroundBottomColor = clf.find("-bottomcolor", backgroundTopColor*0.2 );
   d.backgroundTopColor = Color( backgroundTopColor[0], backgroundTopColor[1], backgroundTopColor[2],1.0) ;
   d.backgroundBottomColor = Color( backgroundBottomColor[0], backgroundBottomColor[1], backgroundBottomColor[2],1.0) ;

   RenderPatchAllocator::PixelPatch pix;
   vector<Color> output;

   vector<int> selectPatches = clf.findMultiple( "-patch", -1, "Select a patch to render");
   int count =  rpa->nbUncompletedPatches();
   if( !selectPatches.empty() ){ count = selectPatches.size(); }
   ProgressMeter metr( count , "RenderBackground" );

   while( rpa->nbUncompletedPatches() != 0 )
   {
      bool renderThisOne = true;
      if( !selectPatches.empty() )
      {
         vector<int>::iterator check = find( selectPatches.begin(), selectPatches.end(), (int)( rpa->nbUncompletedPatches() ) );
	 renderThisOne = (check != selectPatches.end() );
      }
      rpa->popPatch( data.startPosition, data.startDirection, pix );

      if( renderThisOne )
      { 
         ssBackgroundRender( data, output ); 

         for( size_t i=0;i<pix.size();i++ )
         {
            vector<float>& pixel = img.pixel( pix[i] );
            for( size_t c=0;c<pixel.size();c++ ) 
	    { 
	       pixel[c] += (1.0 - pixel[c][3])*output[i][c]; 
	    }
         }
      }
      metr.update();
   }

*/
}

