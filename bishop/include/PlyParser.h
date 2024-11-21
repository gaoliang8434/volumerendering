//-----------------------------------------------------------
//
//  PlyParser.h
//
//  Parses ply files and extracts geometry
//
//  Copyright (c) 2005, Finelight Visual Technology, Inc.
//
//------------------------------------------------------------


#ifndef ____LUX_PLYPARSER_H____
#define ____LUX_PLYPARSER_H____

#include <string>
using namespace std;
#include "AsciiParser.h"

namespace lux
{

class Geometry;
	
class PlyParser
{
  public:

    PlyParser(){}
   ~PlyParser(){}

    const bool ParseFile( const string& filename );
    const bool Fill( Geometry& g );
    const bool List();

    const bool IsPlyLine() const;
    const bool IsFormatLine() const;
    const bool IsComment() const;
    const bool IsElement() const;
    const bool IsVertex() const;
    const bool IsFace() const;
    const bool IsProperty() const;
    const bool IsHeaderEnd() const;
    
    const bool AdvanceToNextLine();
    const bool GetVertex( double& x, double& y, double& z );
    const bool GetFace( int& x, int& y, int& z );

    const bool GetPlyLine();
    const bool GetFormatLine();
    const int  GetVertexCount();
    const int  GetFaceCount();
    const bool GetHeaderEnd();

  private:

    AsciiParser parser;
    int nb_vertices;
    int nb_faces; 
};

}

#endif
