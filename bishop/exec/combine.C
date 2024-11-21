//===============================================================================
//
//  combining images of subpatches into a single frame
//
//===============================================================================





// Standard utilities and system includes, plus project specific items
//===============================================================================
using namespace std;

#include "LABLogo.h"
#include "ProgressMeter.h"
#include "Magick++.h"

using namespace std;
using namespace Magick;
using namespace lux;





// Main program
//===============================================================================
int main(int argc, char** argv)
{
    cout << LABLOGOSTRING << flush;


    Image combined;
    bool combinedInitialized = false;
    int width=0, height=0;

    for( int i=1;i<argc-1;i++ )
    {
       cout << "Processing file " << argv[i] << endl;
       Image input;
       input.read( argv[i] );
       if( !combinedInitialized )
       {
           width = input.baseColumns();
	   height = input.baseRows();
	   combined = input;
          combinedInitialized = true;
       }
       else
       {
          for( int j=0;j<height;j++ )
	  {
	     for( int i=0;i<width;i++ )
	     {
		ColorRGB incolor = input.pixelColor(i,j);
		ColorRGB combcolor = combined.pixelColor(i,j);
		if( i == width/2 )
		{
		   cout << incolor.red() << " " << incolor.green() << " " << incolor.blue() << "     " << combcolor.red() << " " << combcolor.green() << " " << combcolor.blue() << endl;
		}
		double red = incolor.red() + combcolor.red();
		double green = incolor.green() + combcolor.green();
		double blue = incolor.blue() + combcolor.blue();
		combcolor.red( red );
		combcolor.green( green );
		combcolor.blue( blue );
		combined.pixelColor(i,j,combcolor);
	     }
	  }
	  combined.syncPixels();
       }
    }

    cout << "Writing combined image to file " << argv[argc-1] << endl;
    combined.write( argv[argc-1] );

}

