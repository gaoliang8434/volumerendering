

#include "Tracking.h"
#include <algorithm>

using namespace lux;
using namespace std;


const long lux::getBoundingBoxes( const ScalarGrid& g, vector<AARectangle>& boxes )
{
   long nbboxes = g->Size();

   if( nbboxes == 0 ){ return nbboxes; }

   boxes.clear();

   int i0,j0,k0, i1,j1,k1;
   for( int box=0;box<nbboxes;box++ )
   {
      g->blockBounds( box, i0, j0, k0, i1, j1, k1 );
      if( i0 >= 0 && i1 >= 0 && j0 >= 0 && j1 >= 0 && k0 >= 0 && k1 >= 0 )
      {
         Vector llc = g->evalP(i0,j0,k0);
         Vector urc = g->evalP(i1,j1,k1);
         boxes.push_back( AARectangle( llc, urc ) );
      }
   }
   return nbboxes;
}


void lux::findRayMarchBoxes( const ScalarGrid& g, const Vector X0, const Vector D, std::vector<TraceSegment>& sortedSegments, double tolerance )
{
   std::vector<AARectangle> goodBoxes;
   double dx = g->dx();
   dx = ( dx <= g->dy() ) ? dx : g->dy();
   dx = ( dx <= g->dz() ) ? dx : g->dz();
   double ds = g->blockSize() * dx;
   Vector offset( g->dx(), g->dy(), g->dz() );


   Vector P = X0;
   AARectangle global( g->llc(), g->urc() );
   std::vector<AARectangle> globe;
   globe.push_back(global);
   std::vector<double> startstop;
   findRayMarchBoxes( globe, X0, D, dx, startstop );
   if( startstop.size() != 2 ) { return; }

   findRayMarchBoxes( globe, X0, D, sortedSegments, tolerance );
   return;

   // march along ray and verify which boxes are hit
   P = X0 + startstop[0]*D;
   double s = startstop[0];
   int i0,j0,k0,i1,j1,k1;
   long lastindex = -1;
   long currentindex = -1;
   while( s < startstop[1] )
   {
      if( g->getGridIndex( P, i0, j0, k0 ) )
      {
         currentindex = g->index( i0,j0,k0 );
         if ( lastindex != currentindex )
	 {
            if( g->goodBlock( i0,j0,k0 ) )
	    {
	       g->blockBounds( currentindex, i0, j0, k0, i1, j1, k1 );
               if( i0 >= 0 && i1 >= 0 && j0 >= 0 && j1 >= 0 && k0 >= 0 && k1 >= 0 )
               {
                  Vector llc = g->evalP(i0,j0,k0) - offset/2.0;
                  Vector urc = g->evalP(i1,j1,k1) + offset/2.0;
                  goodBoxes.push_back( AARectangle( llc, urc ) );
	       }
	    }
	 }
         lastindex = currentindex;
      }
      P += D*ds;
      s += ds;
   }
   findRayMarchBoxes( goodBoxes, X0, D, sortedSegments, tolerance );
}






void lux::findRayMarchBoxes( const std::vector<AARectangle>& boxes, const Vector X0, const Vector D, std::vector<TraceSegment>& sortedSegments, double tolerance )
{
   // find the subset along the ray

   std::vector<double> hits;
   findRayMarchBoxes( boxes, X0, D, tolerance, hits );
   for( size_t i=0;i<hits.size(); i+=2 )
   {
      sortedSegments.push_back( TraceSegment( X0 + hits[i]*D, hits[i+1]-hits[i], hits[i] ) );
   }
}



void lux::findRayMarchBoxes( const std::vector<AARectangle>& boxes, const Vector X0, const Vector D, const double tolerance, std::vector<double>& sortedSegments )
{
   double near, far;
   std::vector<double> hitpoints;
   for( size_t b=0;b<boxes.size();b++ )
   {
      const AARectangle& box = boxes[b];
      if(  box.intersection( X0, D, near, far ) )
      {
	 hitpoints.push_back( near );
	 hitpoints.push_back( far );
      }
   }
   sortMergeHitPoints( hitpoints, sortedSegments, tolerance );
}

void lux::sortMergeHitPoints( std::vector<double>& hitpoints, std::vector<double>& sortedMerged, double tolerance )
{
   // first sort 
   std::sort( hitpoints.begin(), hitpoints.end() );
   sortedMerged.clear();

   // now merge close segments

   size_t i = 1;
   while( i < hitpoints.size() )
   {
      sortedMerged.push_back( hitpoints[i-1] );
      size_t icheck = i+1;
      while( (hitpoints[icheck] - hitpoints[icheck-1] <= tolerance) && ( icheck < hitpoints.size() ) )
      {
         icheck += 2;
      }

      if( icheck >= hitpoints.size() )
      { 
         icheck = hitpoints.size() - 1; 
      }
      sortedMerged.push_back( hitpoints[icheck] );
      i = icheck + 1;
   }
}
