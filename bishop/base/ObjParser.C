//-----------------------------------------------------------
//
//  ObjParser.C
//
//  Parses obj files and extracts geometry
//
//------------------------------------------------------------

#include <iostream>
using namespace std;

#include "ObjParser.h"
#include "VolumeGeometry.h"

namespace lux
{

const bool ObjParser::ParseFile( const string& filename )
{
   if( ! parser.ParseFile( filename ) ) { return false; }
   nb_vertices = GetVertexCount();
   nb_faces    = GetFaceCount();
   return ( (nb_vertices > 0) && (nb_faces>0) );
}

	
const bool ObjParser::IsObjLine() const
{
   if( ! parser.IsToken() ) { return false; }
   if( ! parser.IsText() ){ return false; }
   return ( parser.TextValue() == "obj" );
}

const bool ObjParser::IsComment() const
{
   if( ! parser.IsText() ){ return false; }
   return ( parser.TextValue() == "#" );
}

const bool ObjParser::IsVertex() const
{
   if( ! parser.IsText() ){ return false; }
   return ( parser.TextValue() == "v" );
}

const bool ObjParser::IsNormal() const
{
   if( ! parser.IsText() ){ return false; }
   return ( parser.TextValue() == "vn" );
}


const bool ObjParser::IsTexture() const
{
   if( ! parser.IsText() ){ return false; }
   return ( parser.TextValue() == "vt" );
}


const bool ObjParser::IsFace() const
{
   if( ! parser.IsText() ){ return false; }
   return ( parser.TextValue() == "f" );
}

const bool ObjParser::IsGroup() const
{
   if( ! parser.IsText() ){ return false; }
   return ( parser.TextValue() == "g" );
}


const bool ObjParser::AdvanceToNextLine()
{
   while( !parser.IsEOL() )
   {
       if( !parser.GetToken() ) { return false; }
   }
   // move one token past EOL 
   if( !parser.GetToken() ) { return false; }
   return true;
}

const bool ObjParser::GetVertex( double& x, double& y, double& z )
{
   parser.GetToken();
   if( ! parser.IsFloat() && ! parser.IsInteger() ){ return false; }
   x = parser.FloatValue();
   parser.GetToken();
   if( ! parser.IsFloat() && ! parser.IsInteger() ){ return false; }
   y = parser.FloatValue();
   parser.GetToken();
   if( ! parser.IsFloat() && ! parser.IsInteger() ){ return false; }
   z = parser.FloatValue();
   return true;
}


const bool ObjParser::GetNormal( double& x, double& y, double& z )
{
   parser.GetToken();
   if( ! parser.IsFloat() && ! parser.IsInteger() ){ return false; }
   x = parser.FloatValue();
   parser.GetToken();
   if( ! parser.IsFloat() && ! parser.IsInteger() ){ return false; }
   y = parser.FloatValue();
   parser.GetToken();
   if( ! parser.IsFloat() && ! parser.IsInteger() ){ return false; }
   z = parser.FloatValue();
   return true;
}




const bool ObjParser::GetTextureCoordinate( double& x, double& y, double& z )
{
   parser.GetToken();
   if( ! parser.IsFloat() && ! parser.IsInteger() ){ return false; }
   x = parser.FloatValue();
   parser.GetToken();
   if( ! parser.IsFloat() && ! parser.IsInteger() ){ return false; }
   y = parser.FloatValue();
   parser.GetToken();
   if( ! parser.IsFloat() && ! parser.IsInteger() ){ z = 0; return true; }
   z = parser.FloatValue();
   return true;
}




const bool ObjParser::GetFace( int& x, int& y, int& z )
{
   parser.GetToken();
   if( ! parser.IsInteger() ){ return false; }
   x = parser.IntegerValue();
   parser.GetToken();
   // skip texture and normal info
   if( parser.IsSeparator() )
   {
      parser.GetToken();
      if( parser.IsInteger() ){ parser.GetToken(); }
      if( parser.IsSeparator() ){ parser.GetToken(); parser.GetToken(); }
   }

   if( ! parser.IsInteger() ){ return false; }
   y = parser.IntegerValue();
   parser.GetToken();
   // skip texture and normal info
   if( parser.IsSeparator() )
   { 
      parser.GetToken();
      if( parser.IsInteger() ){ parser.GetToken(); }
      if( parser.IsSeparator() ){ parser.GetToken(); parser.GetToken(); }
   }
   if( ! parser.IsInteger() ){  return false; }
   z = parser.IntegerValue();

   return true;
}



const bool ObjParser::GetTexturedFace( int& x, int& y, int& z, int& xt, int& yt, int& zt )
{
   parser.GetToken();
   if( ! parser.IsInteger() ){ return false; }
   x = parser.IntegerValue();
   parser.GetToken();
   if( parser.IsSeparator() )
   {
      parser.GetToken();
      if( parser.IsInteger() )
      {
         xt = parser.IntegerValue();
         parser.GetToken(); 
      }
      // skip normal
      if( parser.IsSeparator() ){ parser.GetToken(); parser.GetToken(); }
   }

   if( ! parser.IsInteger() ){ return false; }
   y = parser.IntegerValue();
   parser.GetToken();
   if( parser.IsSeparator() )
   { 
      parser.GetToken();
      if( parser.IsInteger() )
      { 
         yt = parser.IntegerValue();
         parser.GetToken(); 
      }
      // skip normal
      if( parser.IsSeparator() ){ parser.GetToken(); parser.GetToken(); }
   }
   if( ! parser.IsInteger() ){  return false; }
   z = parser.IntegerValue();
   parser.GetToken();
   if( parser.IsSeparator() )
   { 
      parser.GetToken();
      if( parser.IsInteger() )
      { 
         zt = parser.IntegerValue();
         parser.GetToken(); 
      }
      // skip normal
      //if( parser.IsSeparator() ){ parser.GetToken(); parser.GetToken(); }
   }

   return true;
}



const bool ObjParser::GetTexturedNormaledFace( int& x, int& y, int& z, int& xt, int& yt, int& zt, int& xn, int& yn, int& zn )
{
   parser.GetToken();
   if( ! parser.IsInteger() ){ return false; }
   x = parser.IntegerValue();
   parser.GetToken();
   if( parser.IsSeparator() )
   {
      parser.GetToken();
      if( parser.IsInteger() )
      {
         xt = parser.IntegerValue();
         parser.GetToken(); 
      }
      // normal
      if( parser.IsSeparator() )
      { 
         parser.GetToken(); 
        if( parser.IsInteger() )
        {
           xn = parser.IntegerValue();
           parser.GetToken(); 
        }
      }
   }

   if( ! parser.IsInteger() ){ return false; }
   y = parser.IntegerValue();
   parser.GetToken();
   if( parser.IsSeparator() )
   { 
      parser.GetToken();
      if( parser.IsInteger() )
      { 
         yt = parser.IntegerValue();
         parser.GetToken(); 
      }
      //  normal
      if( parser.IsSeparator() )
      { 
         parser.GetToken(); 
        if( parser.IsInteger() )
        {
           yn = parser.IntegerValue();
           parser.GetToken(); 
        }
      }
   }
   if( ! parser.IsInteger() ){  return false; }
   z = parser.IntegerValue();
   parser.GetToken();
   if( parser.IsSeparator() )
   { 
      parser.GetToken();
      if( parser.IsInteger() )
      { 
         zt = parser.IntegerValue();
         parser.GetToken(); 
      }
      // normal
      if( parser.IsSeparator() )
      { 
         parser.GetToken(); 
        if( parser.IsInteger() )
        {
           zn = parser.IntegerValue();
           parser.GetToken(); 
        }
      }
   }

   return true;
}





const bool ObjParser::GetObjLine()
{
   parser.Rewind();
   parser.GetToken();
   return IsObjLine();
}


const int  ObjParser::GetVertexCount()
{
   parser.Rewind();
   parser.GetToken();
   int count = 0;
   while( AdvanceToNextLine() )
   {
      if( IsVertex() ){ ++count; }
   }
   return count;
}

const int  ObjParser::GetFaceCount()
{
   parser.Rewind();
   parser.GetToken();
   int count = 0;
   while( AdvanceToNextLine() )
   {
      if( IsFace() ) { ++count; }
   }
   return count;
}


const bool ObjParser::List()
{
   parser.Rewind();
   parser.GetToken();
   
   while( AdvanceToNextLine() )
   {
      if( IsFace() )
      {
         int i,j,k;
	 GetFace( i,j,k );
      }
      else if ( IsVertex() )
      {
         double x, y, z;
	 GetVertex( x,y,z );
      }
    }
   return true;
}


const bool ObjParser::Fill( TriangleGeometry& g )
{
   parser.Rewind();
   parser.GetToken();
   bool hasTextureCoordinates = false; 
   bool hasNormals = false; 
   while( AdvanceToNextLine() )
   {
      if( IsFace() )
      {
         int i,j,k,it,jt,kt,in,jn,kn;
	 if( hasTextureCoordinates || hasNormals )
	 {
        if( hasNormals )
        {
	       if( !GetTexturedNormaledFace( i,j,k, it,jt,kt, in, jn, kn ) ) { cout << "Could not read face " << g.nbFaces() << endl; }
	       g.addTexturedNormaledFace( i-1,j-1,k-1,it-1,jt-1,kt-1, in-1, jn-1, kn-1 );
        }
        else
        {
	       if( !GetTexturedFace( i,j,k, it,jt,kt ) ) { cout << "Could not read face " << g.nbFaces() << endl; }
	       g.addTexturedFace( i-1,j-1,k-1,it-1,jt-1,kt-1 );
        }
	 }
	 else
	 {
	    if( !GetFace( i,j,k ) ) { cout << "Could not read face " << g.nbFaces() << endl; }
	    g.addFace( i-1,j-1,k-1 );
	 }
      }
      else if ( IsVertex() )
      {
         double x, y, z;
	 if( !GetVertex( x,y,z ) ) { cout << "Could not read vertex " << g.nbVertices() << endl; }
	 g.addVertex( Vector(x,y,z) );
      }
      else if ( IsTexture() )
      {
         hasTextureCoordinates = true;
         double x, y, z;
	     if( !GetTextureCoordinate( x,y,z ) ) { cout << "Could not read texture coordinate " << g.nbVertices() << endl; }
	     g.addTextureCoordinate( Vector(x,y,z) );
      }
      else if ( IsNormal() )
      {
         hasNormals = true;
         double x, y, z;
	     if( !GetNormal( x,y,z ) ) { cout << "Could not read normal " << g.nbVertices() << endl; }
	     g.addNormal( Vector(x,y,z) );
      }
      else
      {
         cout << "Token not recognized\n";
      }
   }
   return ( g.nbFaces() > 0 && g.nbVertices() > 0 );
}



const bool ObjParser::Fill( AnchorChain& g )
{
   parser.Rewind();
   parser.GetToken();
   while( AdvanceToNextLine() )
   {
      if ( IsVertex() )
      {
         double x, y, z;
	     if( !GetVertex( x,y,z ) ) { cout << "Could not read vertex " << g.size() << endl; }
         Noise_t point;
         point.P =  Vector(x,y,z);
         point.v = Vector(0,0,0);
         point.A = point.v;
         point.pscale = 0.1;
         point.translate = point.P;
	     g.push_back(point);
      }
   }
   return ( g.size() > 0 );
}





}

