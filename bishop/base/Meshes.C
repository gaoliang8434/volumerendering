#include "Meshes.h"

#include "ObjParser.h"
#include "Grids.h"
#include "Fields.h"

using namespace lux;

Mesh lux::getMesh( const string fname )
{
   ObjParser parser;
   parser.ParseFile( fname );
   TriangleGeometry* tri = new TriangleGeometry();
   parser.Fill(*tri);
   cout << "Number of Vertices: " << tri->nbVertices() << endl;
   cout << "Number of Faces: " << tri->nbFaces() << endl;
   return Mesh( tri );
}


ScalarField lux::mesh2ls( Mesh& mesh, GridBox& gb )
{
   float maxdist = gb->Lx() + gb->Ly() + gb->Lz();
   ScalarGrid lsgrid = makeGrid( gb, -maxdist );
   RayMarchLevelSet( *mesh, lsgrid );  
   return gridded(lsgrid);
}

Mesh lux::ls2mesh( ScalarField& sf, GridBox& gb ) { return Mesh(); }

Vector lux::llc( Mesh& m ) { return m->LLC(); }
Vector lux::urc( Mesh& m ) { return m->URC(); }

