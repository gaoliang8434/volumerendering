//-------------------------------------------------------
//
//  TextureViewer.h
//
//
//--------------------------------------------------------

#include <iostream>
#include <string>
#include <cstdio>

#ifdef __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h> // GLUT support library.
#endif
#include "TextureViewer.h"
#include "ImageFiles.h"

using namespace std;

namespace lux{


TextureViewer::TextureViewer( CmdLineFind& clf) :
   pos_x  (20),
   pos_y  (40),
   translate_x  (0),
   translate_y  (0),
   pClf         (&clf)
{
   DepthBuffer();
}


TextureViewer::~TextureViewer(){}



void TextureViewer::Init( int * argc, char **argv )
{
   FrameViewer::Init( argc, argv );
   InitTextureList();
}


void TextureViewer::Display()
{
   if( CurrentFrame() < 0 || CurrentFrame() >= size() ){ return; }
   glClear(GL_DEPTH_BUFFER_BIT );
   glEnable(GL_DEPTH);

   // Insert texture data here
   SetTexture();
}

	

void TextureViewer::SetTexture()
{
   Frame& frame = GetFrame( CurrentFrame() );

   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();

   glBindTexture( GL_TEXTURE_2D, texturelist[CurrentFrame()] );
   glTranslatef( translate_x, translate_y, 0 ); 
  glBegin(GL_QUADS);
    glTexCoord2f(0.0,0.0);
    glVertex2f(0, 0);
    glTexCoord2f(1.0,0.0);
    glVertex2f(frame.width*GetZoom(), 0);
    glTexCoord2f(1.0,1.0);
    glVertex2f(frame.width*GetZoom(), frame.height*GetZoom());
    glTexCoord2f(0.0,1.0);
    glVertex2f(0, frame.height*GetZoom());
  glEnd();
   
}



void TextureViewer::InitTextureList()
{
   if( texturelist.size() > 0 ) { return; }

   texturelist.resize( size() );
   glEnable(GL_TEXTURE_2D);

   glGenTextures(size(), &(texturelist[0]));
   for( int i=0;i<size();i++ )
   {
      Frame& frame = GetFrame( i );
      glBindTexture(GL_TEXTURE_2D, texturelist[i] );
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexImage2D(GL_TEXTURE_2D, 0, frame.data_type, frame.width, frame.height, 0, frame.data_type, GL_FLOAT, &(frame.data[0]) );
   }
}

void TextureViewer::AddTextureToList()
{
   texturelist.push_back( 0 );
   glGenTextures( 1, &(texturelist[texturelist.size()-1]) );
   Frame& frame = GetFrame( size()-1 );
   glBindTexture(GL_TEXTURE_2D, texturelist[texturelist.size()-1] );
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
   glTexImage2D(GL_TEXTURE_2D, 0, frame.data_type, frame.width, frame.height, 0, frame.data_type, GL_FLOAT, &(frame.data[0]) );  
}


void TextureViewer::ResetTexture()
{
   Frame& frame = GetFrame( CurrentFrame() );
   int texture = texturelist[ CurrentFrame() ];
   glBindTexture(GL_TEXTURE_2D, texture );
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
   long tsize = frame.width * frame.height * 3;
   float *texturedata = new float[ tsize ];
   for( long i=0;i<tsize;i++ )
   {
      texturedata[i] = pow( (double)frame.data[i], (double)frame.gamma ) * frame.contrast; 
   }
   glTexImage2D(GL_TEXTURE_2D, 0, frame.data_type, frame.width, frame.height, 0, frame.data_type, GL_FLOAT, &(texturedata[0]) );
   delete[] texturedata;
}


void TextureViewer::Motion( int x, int y )
{
   FrameViewer::Motion(x,y);
   if( (_keystate == GLUT_ACTIVE_SHIFT ) ||
       (_keystate == GLUT_ACTIVE_CTRL  ) ||
       (_keystate == GLUT_ACTIVE_ALT   ) ) { return; }
   float dx = x - _mouse_x;
   float dy = y - _mouse_y;
   translate_x = pos_x + dx;
   translate_y = pos_y - dy;
}


void TextureViewer::Mouse( int button, int state, int x, int y )
{
   FrameViewer::Mouse( button, state, x, y );
   if( (_keystate == GLUT_ACTIVE_SHIFT ) ||
       (_keystate == GLUT_ACTIVE_CTRL  ) ||
       (_keystate == GLUT_ACTIVE_ALT   ) ) { return; }

   if( _mouse_state == GLUT_UP )
   {
      pos_x = translate_x;
      pos_y = translate_y;
   }
}



void TextureViewer::Keyboard( unsigned char key, int x, int y )
{
  Frame& frame = GetFrame( CurrentFrame() );
  switch( key )
   {
      case 'g':
        frame.gamma /= 1.3;
	cout << "Gamma " << frame.gamma << endl << flush;
	ResetTexture();
	glutPostRedisplay();
	break;
      case 'G':
        frame.gamma *= 1.3;
	cout << "Gamma " << frame.gamma << endl << flush;
	ResetTexture();
	glutPostRedisplay();
	break;
      case 'b':
        frame.contrast /= 1.3;
	cout << "Contrast " << frame.contrast << endl << flush;
	ResetTexture();
	glutPostRedisplay();
	break;
      case 'B':
        frame.contrast *= 1.3;
	cout << "Contrast " << frame.contrast << endl << flush;
	ResetTexture();
	glutPostRedisplay();
	break;
      case 'w':
        cout << "Writing file to disk....." << flush;
        writeMagickImage( *pClf, frame );
        cout << "DONE" << endl << flush;
	break;
      	
   }
   FrameViewer::Keyboard( key, x, y );
}






}

