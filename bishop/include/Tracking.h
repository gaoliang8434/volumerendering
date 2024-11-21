#ifndef __TRACKING_H__
#define __TRACKING_H__

#include "AARectangle.h"
#include "SparseGrid.h"
#include <vector>

namespace lux
{

class TraceSegment
{
   public:

   TraceSegment( const Vector& X, const double thickness, const double range ) :
      P  (X),
      thick (thickness),
      dist  (range)
   {}

   ~TraceSegment(){}

    const bool operator<( const TraceSegment& t ) const { return dist < t.dist; }

    const Vector& pos() const { return P; }
    const double& thickness() const { return thick; }
    const double& range() const { return dist; }

  private:

    Vector P;
    double thick, dist;
};

const long getBoundingBoxes( const ScalarGrid& g, std::vector<AARectangle>& boxes ); 
void findRayMarchBoxes( const ScalarGrid& g, const Vector X0, const Vector D, std::vector<TraceSegment>& sortedSegments, double tolerance );
void findRayMarchBoxes( const std::vector<AARectangle>& boxes, const Vector X0, const Vector D, std::vector<TraceSegment>& sortedSegments, double tolerance );
void findRayMarchBoxes( const std::vector<AARectangle>& boxes, const Vector X0, const Vector D, double tolerance, std::vector<double>& sortedSegments );
void sortMergeHitPoints( std::vector<double>& hitpoints, std::vector<double>& sortedMerged, double tolerance );

}

#endif
