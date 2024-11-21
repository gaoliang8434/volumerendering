//-------------------------------------------------------
//
//  Frame.h
//
//  Defines a container for a single frame that can be 
//  viewed in the Viewer.  
//
//
//--------------------------------------------------------


#ifndef ____FVT_FRAME_H____
#define ____FVT_FRAME_H____

#include <vector>
#include <string>


using namespace std;

namespace lux{

class Frame
{
  public:

    Frame(){}
   ~Frame(){ data.clear(); index.clear(); }
    
    int width, height;
    vector<float> data;
    vector<unsigned int> index;
    unsigned int data_type;
    float scale_low, scale_high;
    float contrast, gamma;
    string label;
    float aspect;
};




}

#endif
