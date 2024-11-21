//  Basic Viewer test program
//
//  Jerry Tessendorf
//
//  11 November, 2003
//

#include "Viewer.h"

using namespace std;
using namespace lux;


int main( int argc, char **argv )
{

    Viewer viewer;

    viewer.SetWidth( 512 );
    viewer.SetHeight( 512 );
    viewer.Init( &argc, argv );

    viewer.MainLoop();

	
};
