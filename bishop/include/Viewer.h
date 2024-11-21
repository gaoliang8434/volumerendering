//-------------------------------------------------------
//
//  Viewer.h
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


#ifndef ____LUX_VIEWER_H____
#define ____LUX_VIEWER_H____

#include <string>
#include "Vector.h"

using namespace std;

namespace lux{

class Viewer
{
  public:

    Viewer();
    virtual ~Viewer();
    virtual void Init( int * argc, char **argv );
    virtual void MainLoop();
    
    void SetWidth( const int width ) { _width = width; }
    void SetHeight( const int height ) { _height = height; }

    const int& GetWidth() { return _width;  }
    const int& GetHeight() { return _height; }

    void SetTitle( const string& t ){ _title = t; }
    void SetTitle( const char * t ) { _title = t; }
    const string& GetTitle() { return _title; }

    void Luminance();
    void RGB();
    void RGBA();
    void SingleBuffer();
    void DoubleBuffer();
    void DepthBuffer();

    // Callback virtual functions
    virtual void Display();
    virtual void Keyboard( unsigned char key, int x, int y );
    virtual void Mouse( int button, int state, int x, int y );
    virtual void Motion( int x, int y );
    virtual void Special( int key, int x, int y ){}
    virtual void Idle(){}
    virtual void Reshape( int w, int h );
    virtual void Timed( int v ){}
    virtual void Menu( int item ){}

    const float GetZoom() const { return _zoom_factor; } 
    void SetZoom( const float z ){ _zoom_factor = z; }

    static Viewer* pViewer;

    void usage(); 
    
  private:

    bool _initialized;
    int _width, _height;
    string _title;
    unsigned int _display_mode;
    float _zoom_factor, _zoom_increment;  // for zooming the image


    
    // dont allow any of these
    Viewer( const Viewer& );
    Viewer& operator= (const Viewer&);

  protected:

    int _mouse_x, _mouse_y;
    int _keystate, _button;
    int _mouse_state;
    int _window_width, _window_height;
    float _current_raster_pos[4];
   
    Vector cameraEye;
    Vector cameraView;
    Vector cameraUp;
    float cameraFov;
    float cameraAspect;
    float cameraNear;
    float cameraFar;

 
};





}





#endif
