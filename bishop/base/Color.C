//*******************************************************************
//
//   Color.C
//
// 4D color class in the namespace lux
//
//
//
//*******************************************************************

#include "Color.h"

using namespace lux;


Color lux::exp(const Color& c ){ return Color( std::exp(c[0]), std::exp(c[1]), std::exp(c[2]), std::exp(c[3]) ); }



