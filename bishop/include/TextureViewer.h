//-------------------------------------------------------
//
//  TextureViewer.h
//
//  Uses Frame data as a mesh of 3D vertices
//
//  Copyright (c) 2004 Finelight Visual Technology, Inc.
//
//
//--------------------------------------------------------


#ifndef ____FVT_MESHVIEWER_H____
#define ____FVT_MESHVIEWER_H____

#include <vector>
using namespace std;
#include "FrameViewer.h"
#include "CmdLineFind.h"

namespace lux{

class TextureViewer : public FrameViewer
{

  public:

    TextureViewer( CmdLineFind& clf );
    virtual ~TextureViewer();

    virtual void Init( int * argc, char **argv );

    virtual void Display();
    virtual void Motion( int x, int y );
    virtual void Mouse( int button, int state, int x, int y );
    virtual void Keyboard( unsigned char key, int x, int y );

    virtual void AddTextureToList();

  protected:

    vector<GLuint> texturelist;
    float pos_x, pos_y;
    float translate_x, translate_y;
    
  private:

    void SetTexture();
    void ResetTexture();
    void InitTextureList();
    
    CmdLineFind* pClf;
};


}

#endif

