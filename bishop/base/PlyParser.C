//-----------------------------------------------------------
//
//  PlyParser.C
//
//  Parses ply files and extracts geometry
//
//  Copyright (c) 2005, Finelight Visual Technology, Inc.
//
//------------------------------------------------------------

#include <iostream>
using namespace std;

#include "PlyParser.h"
#include "VolumeGeometry.h"


namespace lux
{

const bool PlyParser::ParseFile( const string& filename )
{
   if( ! parser.ParseFile( filename ) ) { return false; }
   if( ! GetPlyLine() ) { return false; }
   if( ! GetFormatLine() ){ return false; }
   nb_vertices = GetVertexCount();
   nb_faces    = GetFaceCount();
   return ( (nb_vertices > 0) && (nb_faces>0) );
}

	
const bool PlyParser::IsPlyLine() const
{
   if( ! parser.IsToken() ) { return false; }
   if( ! parser.IsText() ){ return false; }
   return ( parser.TextValue() == "ply" );
}

const bool PlyParser::IsFormatLine() const
{
   if( ! parser.IsText() ){ return false; }
   return ( parser.TextValue() == "format" );
}

const bool PlyParser::IsComment() const
{
   if( ! parser.IsText() ){ return false; }
   return ( parser.TextValue() == "comment" );
}

const bool PlyParser::IsElement() const
{
   if( ! parser.IsText() ){ return false; }
   return ( parser.TextValue() == "element" );
}

const bool PlyParser::IsVertex() const
{
   if( ! parser.IsText() ){ return false; }
   return ( parser.TextValue() == "vertex" );
}

const bool PlyParser::IsFace() const
{
   if( ! parser.IsText() ){ return false; }
   return ( parser.TextValue() == "face" );
}

const bool PlyParser::IsProperty() const
{
   if( ! parser.IsText() ){ return false; }
   return ( parser.TextValue() == "property" );
}

const bool PlyParser::IsHeaderEnd() const
{
   if( ! parser.IsText() ){ return false; }
   return ( parser.TextValue() == "end_header" );
}
    
const bool PlyParser::AdvanceToNextLine()
{
   while( !parser.IsEOL() )
   {
       if( !parser.GetToken() ) { return false; }
   }
   // move one token past EOL 
   if( !parser.GetToken() ) { return false; }
   return true;
}

const bool PlyParser::GetVertex( double& x, double& y, double& z )
{
   if( ! parser.IsFloat() ){ return false; }
   x = parser.FloatValue();
   parser.GetToken();
   if( ! parser.IsFloat() ){ return false; }
   y = parser.FloatValue();
   parser.GetToken();
   if( ! parser.IsFloat() ){ return false; }
   z = parser.FloatValue();
   return true;
}

const bool PlyParser::GetFace( int& x, int& y, int& z )
{
   if( ! parser.IsInteger() ){ return false; }
   parser.GetToken();
   if( ! parser.IsInteger() ){ return false; }
   x = parser.IntegerValue();
   parser.GetToken();
   if( ! parser.IsInteger() ){ return false; }
   y = parser.IntegerValue();
   parser.GetToken();
   if( ! parser.IsInteger() ){ return false; }
   z = parser.IntegerValue();
   return true;
}

const bool PlyParser::GetPlyLine()
{
   parser.Rewind();
   parser.GetToken();
   return IsPlyLine();
}

const bool PlyParser::GetFormatLine()
{
   parser.Rewind();
   parser.GetToken();
   while( AdvanceToNextLine() )
   {
      if( IsFormatLine() ) { return true; }
   }
   return false;
}

const int  PlyParser::GetVertexCount()
{
   parser.Rewind();
   parser.GetToken();
   while( AdvanceToNextLine() )
   {
      if( IsElement() )
      {
         if( parser.GetToken() )
	 {
            if( IsVertex() )
	    {
	       if( parser.GetToken() )
	       {
	          if( parser.IsInteger() )
		  {
		     return parser.IntegerValue();
		  }
		  return 0;
	       }
	       return 0;
	    }
	 }
      }
   }
   return 0;
}

const int  PlyParser::GetFaceCount()
{
   parser.Rewind();
   parser.GetToken();
   while( AdvanceToNextLine() )
   {
      if( IsElement() )
      {
         if( parser.GetToken() )
	 {
            if( IsFace() )
	    {
	       if( parser.GetToken() )
	       {
	          if( parser.IsInteger() )
		  {
		     return parser.IntegerValue();
		  }
		  return 0;
	       }
	       return 0;
	    }
	 }
      }
   }
   return 0;
}

const bool PlyParser::GetHeaderEnd()
{
   parser.Rewind();
   parser.GetToken();
   while( AdvanceToNextLine() )
   {
      if( IsHeaderEnd() ){ return true; }
   }
   return false;
}


const bool PlyParser::List()
{
   if( GetHeaderEnd() )
   {
      if( !AdvanceToNextLine() ) { return false; }
      double x,y,z;
      for( int i=0;i<nb_vertices;i++ )
      {
         if( GetVertex(x,y,z) )
	 {
               cout << "Vertex " << i << "   " << x << " " << y << " " << z << endl << flush;
	 }
	 else
	 {
               cout << "terminated prematurely\n" << flush;
	       return false;
	 }
	 if( !AdvanceToNextLine() ) { return false; }
      }

      int ix,iy,iz;
      for( int i=0;i<nb_faces;i++ )
      {
         if( GetFace(ix,iy,iz) )
	 {
               cout << "Face " << i << "   " << ix << " " << iy << " " << iz << endl << flush;
	 }
	 else
	 {
               cout << "terminated prematurely\n" << flush;
	       return false;
	 }
	 if( !AdvanceToNextLine() ) { return false; }
      }
   }
   else
   {
      return false;
   }
   return true;
}


const bool PlyParser::Fill( VolumeGeometry& g )
{

   // build a collection of vertices
   vector<Vector> vertex_list;
   
   if( GetHeaderEnd() )
   {
      if( !AdvanceToNextLine() ) { return false; }
      double x,y,z;
      for( int i=0;i<nb_vertices;i++ )
      {
         if( GetVertex(x,y,z) )
	 {
	    vertex_list.push_back( Vector(x,y,z)  );
	 }
	 else
	 {
	       return false;
	 }
	 if( !AdvanceToNextLine() ) { return false; }
      }
   }
   else
   {
      return false;
   }

   // now read the faces and build the GeoTriangles

   int ix,iy,iz;
   for( int i=0;i<nb_faces;i++ )
   {
      if( GetFace(ix,iy,iz) )
      {
         Triangle verts( vertex_list[ix], vertex_list[iy], vertex_list[iz] );
	 Triangle normals( verts.Normal(), verts.Normal(), verts.Normal() );
	 Triangle2 tex( Vector2(0,0), Vector2(0,1), Vector2(1,0) );
	 double theta = M_PI*0.5*random()/RAND_MAX;
	 double phi = M_PI*0.5*random()/RAND_MAX;
	 double cth = cos( theta );
	 double sth = sin( theta );
	 double sphi = sin( phi );
	 double cphi = cos( phi );
	 Color color( sth*cphi, cth, sth*sphi );
	 color *= ( 1 + verts.Normal() * Vector(0,1,0) )*0.6 + 0.1;
	 GeoTriangle gt( verts, normals, tex, color, &g );
	 g.triangles.push_back( gt );
      }
      else
      {
         return false;
      }
      if( !AdvanceToNextLine() ){ return false; }
   }



   
   return true;

	
}



}

