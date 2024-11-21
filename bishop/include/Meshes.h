
#ifndef __MESHES_H__
#define __MESHES_H__

#include "VolumeGeometry.h"



namespace lux
{

class ScalarField;
Mesh getMesh( const string fname );
ScalarField mesh2ls( Mesh& mesh, GridBox& gb );
Mesh ls2mesh( ScalarField& sf, GridBox& gb );

Vector llc( Mesh& m );
Vector urc( Mesh& m );

}

#endif
