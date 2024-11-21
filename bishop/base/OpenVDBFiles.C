#include "OpenVDBFiles.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetSphere.h>
#include "OpenVDBTypes.h"

using namespace std;
using namespace lux;



void lux::writeOpenVDB( const char* fname, const ScalarGrid& g )
{
   // register openvdb types with the system
   openvdb::initialize();
 
   /*** File Output ***/
   FloatGrid::Ptr grid = FloatGrid::create();
 
   // create a coord and accessor
   FloatGrid::Accessor accessor = grid->getAccessor();
   ProgressMeter pm( g->nx()*g->ny()*g->nz(), "writeOpenVDB" );
   for(int k=0;k<g->nz();k++)
   {
      for(int j=0;j<g->ny();j++)
      {
         for(int i=0;i<g->nx();i++)
         {
            Coord ijk(i, j, k);
            accessor.setValue(ijk, g->get(i,j,k));
            pm.update();
         }
      }
   }

 
   // set some optional meta data with the grid for identification.
   // There are 4 types, GRID_UNKNOWN, GRID_LEVEL_SET, GRID_FOG_VOLUME, and GRID_STAGGERED.
   //grid->setGridClass(openvdb::GRID_LEVEL_SET);
   grid->setGridClass(openvdb::GRID_LEVEL_SET);
 
   // set a name for your grid.
   grid->setName(fname);
 
   // create a VDB file object.
   openvdb::io::File file(fname);
 
   // add the grid pointer to a container (since we can write multiple grids out).
   openvdb::GridPtrVec grids;
   grids.push_back(grid);
 
   // write out the container and close the file.
   file.write(grids);
   file.close();
}



void lux::writeOpenVDB( const char* fname, const ScalarField& f, const GridBox& g )
{
   // register openvdb types with the system
   openvdb::initialize();
 
   /*** File Output ***/
   FloatGrid::Ptr grid = FloatGrid::create();
 
   // create a coord and accessor
   FloatGrid::Accessor accessor = grid->getAccessor();
   ProgressMeter pm( g->nx()*g->ny()*g->nz(), "writeOpenVDB" );
   for(int k=0;k<g->nz();k++)
   {
      for(int j=0;j<g->ny();j++)
      {
         for(int i=0;i<g->nx();i++)
         {
            Coord ijk(i, j, k);
            Vector P = g->evalP(i,j,k);
            accessor.setValue(ijk, f->eval(P));
            pm.update();
         }
      }
   }

 
   // set some optional meta data with the grid for identification.
   // There are 4 types, GRID_UNKNOWN, GRID_LEVEL_SET, GRID_FOG_VOLUME, and GRID_STAGGERED.
   //grid->setGridClass(openvdb::GRID_UNKNOWN);
   grid->setGridClass(openvdb::GRID_LEVEL_SET);
 
   // set a name for your grid.
   grid->setName(fname);
 
   // create a VDB file object.
   openvdb::io::File file(fname);
 
   // add the grid pointer to a container (since we can write multiple grids out).
   openvdb::GridPtrVec grids;
   grids.push_back(grid);
 
   // write out the container and close the file.
   file.write(grids);
   file.close();
}
