//-------------------------------------------------------
//
//  FrameViewer.h
//
//  Holds multiple frames and views them
//
//
//--------------------------------------------------------


#ifndef ____FVT_FRAMEVIEWER_H____
#define ____FVT_FRAMEVIEWER_H____

#include <vector>
#include <ctime>

#include "Frame.h"
#include "Viewer.h"

using namespace std;
namespace lux{

class FrameViewer : public Viewer
{
  public:

    FrameViewer();
    virtual ~FrameViewer();
    
    void push_back( const int w, const int h, const unsigned int type );
    Frame& GetFrame( const int i ) { return _frames[i]; }
    const int size() const { return _frames.size(); }

    virtual void Init( int * argc, char **argv );
    virtual void MainLoop(){ Viewer::MainLoop(); }

    // Callback virtual functions
    virtual void Display();
    virtual void Idle();
    virtual void Keyboard( unsigned char key, int x, int y );
    virtual void Special( int key, int x, int y );
    virtual void Motion( int x, int y );
    virtual void Mouse( int button, int state, int x, int y );
    virtual void Timed( int v );
    virtual void Reshape( int w, int h );

    void AdvanceFrame();
    const int CurrentFrame() const { return _current_frame; }
    virtual void SetFrameTitle();

    void SetFrameRate( unsigned int rate ) { _target_frame_rate = rate; }
    const unsigned int GetFrameRate() const { return _target_frame_rate; }

    void StopAnimation(){ _animate = false; }
    void StartAnimation(){ _animate = true; }
    
    
    void usage();
    
  private:

    vector<Frame> _frames;
    int _anim_direction;
    bool _animate;
    int _scrub_x, _scrub_y, _scrub_frame;

// Need to add a frame rate display
    int _clock_frame_counter, _averaging_number;
    clock_t _clock_now, _clock_last;
    float _frame_rate;
    unsigned int _target_frame_rate;

    
  protected:

    int _current_frame;

};


}

#endif
