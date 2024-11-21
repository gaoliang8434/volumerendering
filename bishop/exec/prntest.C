
#include "UniformPRN.h"
#include "GaussianPRN.h"
#include "LognormalPRN.h"
#include <iostream>
using namespace std;
#include "CmdLineFind.h"
using namespace lux;


int main( int argc, char** argv )
{
   CmdLineFind clf( argc, argv );
   UniformPRN uniform;
   GaussianPRN gauss;
   LognormalPRN lognormal;
   Noise_t parms;
   parms.seed = clf.find("-seed", 485758, "Seed for the PRN");

   clf.usage("-h");
   uniform.setParameters( parms );
   gauss.setParameters( parms );
   lognormal.setParameters( parms );

   PRN* prn = &uniform;
   if( clf.findFlag("-gauss" ) )
   { 
      prn = &gauss; 
   }
   else if( clf.findFlag("-lognormal" ) )
   { 
      prn = &lognormal; 
   }

   for( int i=0;i<10000;i++ )
   {
      cout << i << " " << prn->eval() << endl;
   }

}
