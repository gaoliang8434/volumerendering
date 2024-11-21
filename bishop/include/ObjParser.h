//-----------------------------------------------------------
//
//  ObjParser.h
//
//  Parses obj files and extracts geometry
//
//------------------------------------------------------------


#ifndef ____LUX_PLYPARSER_H____
#define ____LUX_PLYPARSER_H____

#include <string>
using namespace std;
#include "AsciiParser.h"
#include "Noise.h"

namespace lux
{

class TriangleGeometry;
	
class ObjParser
{
  public:

    ObjParser(){}
   ~ObjParser(){}

    const bool ParseFile( const string& filename );
    const bool Fill( TriangleGeometry& g );
    const bool Fill( AnchorChain& g );
    const bool List();

    const bool IsObjLine() const;
    const bool IsComment() const;
    const bool IsTexture() const;
    const bool IsVertex() const;
    const bool IsNormal() const;
    const bool IsFace() const;
    const bool IsGroup() const;
    
    const bool AdvanceToNextLine();
    const bool GetVertex( double& x, double& y, double& z );
    const bool GetNormal( double& x, double& y, double& z );
    const bool GetTextureCoordinate( double& x, double& y, double& z );
    const bool GetFace( int& x, int& y, int& z );
    const bool GetTexturedFace( int& x, int& y, int& z, int& xt, int& yt, int& zt );
    const bool GetTexturedNormaledFace( int& x, int& y, int& z, int& xt, int& yt, int& zt, int& xn, int& yn, int& zn );

    const bool GetObjLine();
    const int  GetVertexCount();
    const int  GetFaceCount();
    const bool GetHeaderEnd(){ return false; }

  private:

    AsciiParser parser;
    int nb_vertices;
    int nb_faces; 
};

}

#endif
