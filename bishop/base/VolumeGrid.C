
#include "Volume.h"
#include "VolumeGrid.h"

using namespace lux;




void lux::WriteFloatVolumeGrid( lux::VolumeGrid<float>& grid, const string& fname )
{
   ofstream floatout( fname.c_str() );
   lux::WriteVolumeGrid( grid, floatout );
}



void lux::ReadFloatVolumeGrid( lux::VolumeGrid<float>& grid, const string& fname )
{
   ifstream floatin( fname.c_str() );
   lux::ReadVolumeGrid( grid, floatin );
}


void lux::WriteVectorVolumeGrid( lux::VolumeGrid<Vector>& grid, const string& fname )
{
   ofstream floatout( fname.c_str() );
   lux::WriteVolumeGrid( grid, floatout );
}



void lux::ReadVectorVolumeGrid( lux::VolumeGrid<Vector>& grid, const string& fname )
{
   ifstream floatin( fname.c_str() );
   lux::ReadVolumeGrid( grid, floatin );
}


void lux::WriteColorVolumeGrid( lux::VolumeGrid<Color>& grid, const string& fname )
{
   ofstream floatout( fname.c_str() );
   lux::WriteVolumeGrid( grid, floatout );
}



void lux::ReadColorVolumeGrid( lux::VolumeGrid<Color>& grid, const string& fname )
{
   ifstream floatin( fname.c_str() );
   lux::ReadVolumeGrid( grid, floatin );
}



