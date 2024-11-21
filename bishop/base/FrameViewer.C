//-------------------------------------------------------
//
//  FrameViewer.h
//
//  Holds multiple frames and views them
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

#include "FrameViewer.h"

using namespace std;

namespace lux{

void cbTimedFunc( int v )
{
   Viewer::pViewer -> Timed(v);
}

void cbIdleFunc()
{
   Viewer::pViewer -> Idle();
}

void cbSpecialFunc( int key, int x, int y )
{
   Viewer::pViewer -> Special( key, x, y );
}
	
FrameViewer::FrameViewer() :
  _anim_direction(  1 ),
  _animate       ( true ),
  _averaging_number ( 50 ),
  _frame_rate    ( 0 ),
  _target_frame_rate (30),
  _current_frame ( -1 )
{}

FrameViewer::~FrameViewer()
{
   _frames.clear();
}

void FrameViewer::Display()
{
   if( _current_frame < 0 || _current_frame >= (int)_frames.size() ){ return; }
   Frame& frame = _frames[_current_frame];
   glDrawPixels( frame.width , frame.height, 
		 frame.data_type, GL_FLOAT, &(frame.data[0]) );
}

// this idle function advances to the next frame
void FrameViewer::Idle()
{
   // stop using idle to advance frames, and use a timed frame rate instead
   glutIdleFunc( NULL );
   Timed( 0 );
}
   
void FrameViewer::Timed( int v )
{
   if( _animate )
   {
      AdvanceFrame();
      SetFrameTitle();
      glutPostRedisplay();

      /*************************
      _clock_frame_counter++;
      if( _clock_frame_counter >= _averaging_number )
      {
         _clock_now = clock();
         _frame_rate = (float)CLOCKS_PER_SEC/(float)(_clock_now - _clock_last);
	 cout << "now, last " << _clock_now << ",  " << _clock_last;
         _clock_last = _clock_now;
	 _frame_rate *= (float)_averaging_number;
	 cout << "    rate " << _frame_rate << endl << flush;
	 _clock_frame_counter = 0;
      }
      ***************************/

      unsigned int milliseconds = 1000/_target_frame_rate;
      _frame_rate = (float)_target_frame_rate;
      if( size() > 1 ){ glutTimerFunc( milliseconds, &cbTimedFunc, 0 ); }
   }
   
}


void FrameViewer::push_back( const int w, const int h, const unsigned int type )
{
   Frame temp;
   temp.data_type = type;
   temp.width = w;
   temp.height = h;
   _frames.push_back( temp );
}

void FrameViewer::Init( int * argc, char **argv )
{
   Viewer::Init( argc, argv );
   glutIdleFunc( &cbIdleFunc );
   glutSpecialFunc( &cbSpecialFunc );
   _clock_last = clock();
   _clock_frame_counter = 0;
}


void FrameViewer::AdvanceFrame()
{
   if( _anim_direction > 0 )
   {   
      _current_frame++;
      if( _current_frame >= (int)_frames.size() ){ _current_frame = 0; }
   }
   else
   {
      _current_frame--;
      if( _current_frame < 0 ){ _current_frame = (int)_frames.size()-1; }
   }
}

void FrameViewer::Keyboard( unsigned char key, int x, int y )
{
   switch( key )
   {
      case ' ':
        _animate = !_animate;
	SetFrameTitle();
	glutPostRedisplay();
	if( _animate )
	{
           glutIdleFunc( &cbIdleFunc );
	}
	else
	{
	   glutIdleFunc( NULL );
	}
	break;
   }
   Viewer::Keyboard( key, x, y );
   SetFrameTitle();
}

void FrameViewer::Special( int key, int x, int y )
{
   switch( key )
   {
      case GLUT_KEY_RIGHT:
         _anim_direction = 1;
	 if( !_animate )
	 {
            AdvanceFrame();
            SetFrameTitle();
            glutPostRedisplay();
	 }
	 break;

      case GLUT_KEY_LEFT:
	 _anim_direction = -1;
	 if( !_animate )
	 {
            AdvanceFrame();
            SetFrameTitle();
            glutPostRedisplay();
	 }
	 break;

      case GLUT_KEY_UP:
	 _target_frame_rate++;
	 break;

      case GLUT_KEY_DOWN:
	 if( _target_frame_rate > 1 ){ _target_frame_rate--; }
	 break;
   }
   Viewer::Special( key, x, y );
}


void FrameViewer::Motion( int x, int y )
{
   Viewer::Motion( x, y );
   if( _keystate != GLUT_ACTIVE_SHIFT ) { return; }
   if( _animate ) { return; }

   float fraction = (float)(x - _scrub_x);
   fraction /= (float)GetWidth();
   int frame_move = (int)( fraction * _frames.size() );

   _current_frame = _scrub_frame + frame_move;
   if( _current_frame >= (int)_frames.size() ){ _current_frame = _frames.size()-1; }
   if( _current_frame < 0 ) { _current_frame = 0; }

  
   SetFrameTitle();
   glutPostRedisplay();   	   
}

void FrameViewer::SetFrameTitle()
{
   char title[200];
   if( _animate && (size() > 1) )
   {
      sprintf( title, "%s  Zoom %4.2f  Frame %4.4i  Rate %4.0f /sec", GetTitle().c_str(), GetZoom(), _current_frame+1, _frame_rate );
   }
   else
   {
      if( size() == 1 )
      {
         sprintf( title, "%s  Zoom %4.2f  %s", GetTitle().c_str(), GetZoom(), GetFrame( 0 ).label.c_str() );
      }
      else
      {
         sprintf( title, "%s  Zoom %4.2f  Frame %4.4i", GetTitle().c_str(), GetZoom(), _current_frame+1 );
      }
   }
   glutSetWindowTitle( title );
}

void FrameViewer::Mouse( int button, int state, int x, int y )
{
   Viewer::Mouse( button, state, x, y );
   _scrub_x = x;
   _scrub_y = y;
   _scrub_frame = _current_frame;
}

void FrameViewer::usage()
{
   cout << "FRAME VIEWER CONTROLS:\n";
   cout << "Spacebar                     stop/start animation.\n";
   cout << "Right/Left arrows            direction of animation.\n";
   cout << "Up/Down arrows               animation frame rate.\n";
   cout << "shift-drag mouse             scrub frames.\n";
   Viewer::usage();
}

void FrameViewer::Reshape( int w, int h )
{
   _window_width = w;
   _window_height = h;

   glViewport( 0, 0, (GLsizei) _window_width, (GLsizei) _window_height );
   glMatrixMode( GL_PROJECTION );
   glLoadIdentity();
   gluOrtho2D( 0.0, (GLfloat) _window_width, 0.0, (GLfloat) _window_height );
   glutPostRedisplay();
}



}


