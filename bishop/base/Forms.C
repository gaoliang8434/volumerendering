
#include "Forms.h"

using namespace lux;

const Form lux::wedge( const Form& f1, const Form& f2 ) { return f1^f2; }

const Form lux::star( const Form& f ) { return f.star(); }
