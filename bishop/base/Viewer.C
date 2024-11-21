//-------------------------------------------------------
//
//  Viewer.C
//
//  This viewer is a wrapper of the Glut calls needed to
//  display one or more images.  Options for zooming, 
//  labeling the window frame, etc are available for 
//  derived classes to use.
//  
//
//  Copyright (c) 2003 Finelight Visual Technology, Inc.
//
//
//--------------------------------------------------------

//--------------------------------------------------------
//
//
//--------------------------------------------------------

#ifdef __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif

#include <iostream>
#include "Viewer.h"

using namespace std;
namespace lux{

void cbDisplayFunc()
{
   glClear( GL_COLOR_BUFFER_BIT );	
   Viewer::pViewer -> Display();
   glutSwapBuffers();
   glutPostRedisplay();
}


void cbKeyboardFunc( unsigned char key, int x, int y )
{
   Viewer::pViewer -> Keyboard( key, x, y );
}

void cbMotionFunc( int x, int y )
{
   
   Viewer::pViewer -> Motion( x, y );
}

void cbMouseFunc( int button, int state, int x, int y )
{
   Viewer::pViewer -> Mouse( button, state, x, y );
}

void cbReshapeFunc( int w, int h )
{
   Viewer::pViewer -> Reshape( w, h );
}


Viewer* Viewer::pViewer = 0;
	
Viewer::Viewer() : 
   _initialized    ( false ),
   _width          ( 512 ), 
   _height         ( 512 ),
   _title          ( string("Viewer") ),
   _display_mode   ( GLUT_DOUBLE | GLUT_RGBA ),
   _zoom_factor    ( 1.0 ),
   _zoom_increment ( 1.2 ),
   _mouse_x        ( 0 ),
   _mouse_y        ( 0 )
{}

Viewer::~Viewer(){}

void Viewer::Init( int * argc, char **argv )
{
   glutInit( argc, argv );
   glutInitDisplayMode( _display_mode );
   glutInitWindowSize( _width, _height );
   glutCreateWindow( _title.c_str() );
   glClearColor(0.3,0.3,0.3,0.0);
   pViewer = this;
   _initialized = true;

   _window_width = glutGet( GLUT_WINDOW_WIDTH );
   _window_height = glutGet( GLUT_WINDOW_HEIGHT );

   cameraAspect = (float)_window_width/(float)_window_height;
   cameraFov = 90.0;
   cameraEye = Vector( 0,0,1 );
   cameraView = Vector(0,0,0);
   cameraUp = Vector(0,1,0);
   cameraNear = 0.01;
   cameraFar = 100.0;

   glutDisplayFunc( &cbDisplayFunc );
   glutKeyboardFunc( &cbKeyboardFunc );
   glutMotionFunc( &cbMotionFunc );
   glutMouseFunc( &cbMouseFunc );
   glutReshapeFunc( &cbReshapeFunc );
}

void Viewer::MainLoop()
{
   glutMainLoop();
}


void Viewer::Display()
{
    glLoadIdentity();

    gluPerspective( cameraFov, cameraAspect, cameraNear, cameraFar );
    gluLookAt( cameraEye[0],  cameraEye[1],  cameraEye[2],
               cameraView[0], cameraView[1], cameraView[2],
               cameraUp[0],   cameraUp[1],   cameraUp[2]   );
}

void Viewer::DoubleBuffer() { _display_mode = _display_mode | GLUT_DOUBLE; }
void Viewer::SingleBuffer() { _display_mode = _display_mode | GLUT_SINGLE; }
void Viewer::Luminance() { _display_mode = _display_mode | GLUT_LUMINANCE; }
void Viewer::RGB() { _display_mode = _display_mode | GLUT_RGB; }
void Viewer::RGBA() { _display_mode = _display_mode | GLUT_RGBA; }
void Viewer::DepthBuffer()  { _display_mode = _display_mode | GLUT_DEPTH; }

void Viewer::Reshape( int w, int h )
{
   _window_width = w;
   _window_height = h;

   //gluPerspective( (float)h/(float)w, 1.0, 0.01, 10000.0 );
   glViewport( 0, 0, (GLsizei) _window_width, (GLsizei) _window_height );
   glMatrixMode( GL_PROJECTION );
   glLoadIdentity();
   //gluOrtho2D( 0.0, (GLfloat) _window_width, 0.0, (GLfloat) _window_height );
   glutPostRedisplay();
}

void Viewer::Keyboard( unsigned char key, int x, int y )
{
   switch (key)
   {
      case '+': case '=':
         _zoom_factor *= _zoom_increment;
	 glPixelZoom( _zoom_factor, _zoom_factor );
	 glutPostRedisplay();
	 break;

      case '-': case '_':
         _zoom_factor /= _zoom_increment;
	 glPixelZoom( _zoom_factor, _zoom_factor );
	 glutPostRedisplay();
	 break;	
	
      case 'c': case 'C':
	// this doesnt work because CTRL-C does not
	// come in as the 'c' character with modifier.
	if( glutGetModifiers() == GLUT_ACTIVE_CTRL )
	{
	   exit(0);
	}
        break;	
   }
}


void Viewer::Motion( int x, int y )
{
   if( (_keystate == GLUT_ACTIVE_SHIFT ) ||
       (_keystate == GLUT_ACTIVE_CTRL  ) ||
       (_keystate == GLUT_ACTIVE_ALT   ) ) { return; }
   float dx = x - _mouse_x;
   float dy = y - _mouse_y;
   float pos_x = _current_raster_pos[0] + dx;
   float pos_y = _current_raster_pos[1] - dy;
   glRasterPos2f( pos_x, pos_y ); 
   glutPostRedisplay();
}


void Viewer::Mouse( int button, int state, int x, int y )
{
   _mouse_x = x;
   _mouse_y = y;
   _keystate = glutGetModifiers();
   _button = button;
   _mouse_state = state;
   glGetFloatv( GL_CURRENT_RASTER_POSITION, _current_raster_pos );
}


void Viewer::usage()
{
   cout << "+,=                          Zoom in on image.\n";
   cout << "-,_                          Zoom out on image.\n";
   cout << "drag mouse                   Move image around in display window\n";
}

}


