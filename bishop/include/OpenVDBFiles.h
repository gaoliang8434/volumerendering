#ifndef __OPENVDBFILES_H__
#define __OPENVDBFILES_H__

#include "SparseGrid.h"
#include "Volume.h"

namespace lux{



void writeOpenVDB( const char* fname, const ScalarGrid& g );
void writeOpenVDB( const char* fname, const ScalarField& f, const GridBox& g );



}



#endif
