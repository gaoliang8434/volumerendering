
#ifndef ___INTERVALTREE_H____
#define ___INTERVALTREE_H____

#include <vector>
#include "AARectangle.h"
#include <memory>

namespace lux {

struct IntervalData
{
   double tmin;
   double tmax;
   bool status;
};


class IntervalTree
{
  public:

    IntervalTree( const Vector& llc, const Vector& urc, const int lvl, const int maxlvl, const int minobj );
   ~IntervalTree();

    void Divide();

    void addObject( const AARectangle* o );
    void addObject( AABB& o );

    const IntervalData interval_intersect( const Vector& start, const Vector& direction ) const;

    const size_t nbObjects() const;

    const std::vector<AABB>& objects() const { return object_list; };

    const AABB& object( size_t i ) const { return object_list[i]; };

  private:

    AARectangle aabb;
    IntervalTree* node1;
    IntervalTree* node2;

    int level;
    int max_levels;
    int min_objects;

    std::vector<AABB> object_list;
};

typedef std::shared_ptr<IntervalTree> IntervalSetBase;

class IntervalSet : public IntervalSetBase
{
  public:

    IntervalSet( const Vector& llc, const Vector& urc, const int lvl, const int maxlvl, const int minobj );
   ~IntervalSet();

}; 

IntervalSet makeIntervalSet( const Vector& llc, const Vector& urc, const int lvl, const int maxlvl, const int minobj );
void AddIntervalSetAABB( IntervalSet& is, AABB& aabb );
void DivideIntervalSet( IntervalSet& is );

}

#endif
