
#include "version.h"
#include <sstream>
#include <string>


template <typename T> 
std::string tostr(const T& t) { std::stringstream os; os<<t; return os.str(); }

const std::string lux::versionString() 
{
   std::string v = "Bishop Revision: Compiled: " __DATE__ " " __TIME__ ;
   //std::string v = "Bishop Revision: " + tostr(__REVISION__) + "   Compiled: " __DATE__ " " __TIME__ ;
   return v;
}

