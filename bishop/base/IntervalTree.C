
#include "IntervalTree.h"
#include <iostream>
using namespace std;

using namespace lux;

IntervalTree::IntervalTree( const Vector& llc, const Vector& urc, const int lvl, const int maxlvl, const int minobj ) :
   aabb        (AARectangle(llc,urc)),
   node1       (0),
   node2       (0),
   level       (lvl),
   max_levels  (maxlvl),
   min_objects (minobj)
   {}


IntervalTree::~IntervalTree()
{
   if( node1 ){ delete node1; }
   if( node2 ){ delete node2; }
}

void IntervalTree::addObject( const AARectangle* o ) 
{
   Vector llc = o->llc();
   Vector urc = o->urc(); 
   object_list.push_back(makeAABB(llc,urc)); 
}
void IntervalTree::addObject( AABB& o ) { object_list.push_back(o); }

const size_t IntervalTree::nbObjects() const
{
   size_t total = object_list.size();
   if( node1 ){ total += node1->nbObjects(); }
   if( node2 ){ total += node2->nbObjects(); }
   return total;
}



const IntervalData IntervalTree::interval_intersect( const Vector& start, const Vector& direction ) const
{
   IntervalData hit;
   hit.status = false;
   double tmin0, tmax0;
   if( aabb.isInside(start) || aabb.intersection( start, direction, tmin0, tmax0 ) )
   {
      if( node1 == 0 && node2 == 0 )
      {
         // leaf level.  Intesection against object_list and return results
         for( size_t i=0;i<object_list.size();i++ )
         {
            double tmin, tmax;
            bool result = object_list[i]->intersection( start, direction, tmin, tmax );
            if( result )
            {
               if( !hit.status )
	       {
	          hit.tmin = tmin;
                  hit.tmax = tmax;
	          hit.status = true;
                  //cout << "First object hit tmin, tmax " << tmin << "   " << tmax << endl;
	       }
               else 
               {
                  //cout << "Additional hit tmin, tmax " << tmin << "   " << tmax << "      current tmin tmax " << hit.tmin << "   " << hit.tmax << endl;
                  if( tmin < hit.tmin )
                  {
                     hit.tmin = tmin;
                  }
                  if( tmin <= hit.tmax && tmax > hit.tmax )
                  {
                     hit.tmax = tmax;
                  }
               }
            }
         }
         return hit;
      }
      else
      {
         if( node1 != 0 && node2 != 0 )
         {
            const IntervalData hit1 = node1->interval_intersect( start, direction );
            const IntervalData hit2 = node2->interval_intersect( start, direction );
            if( !hit1.status && !hit2.status ){ return hit; }
            if( hit1.status && !hit2.status ){ return hit1; }
            if( !hit1.status && hit2.status ){ return hit2; }
            if( hit1.tmin < hit2.tmin )
            {
               return hit1; 
            }
            else
            {
               return hit2; 
            }
         }
         else if( node1 == 0 )
         {
            return node2->interval_intersect( start, direction );
         }
         else if( node2 == 0 )
         {
            return node1->interval_intersect( start, direction );
         }
      }
   }
   return hit;
}



void IntervalTree::Divide()
{
   if( level >= max_levels || object_list.size() <= (size_t)min_objects ){ return; }

   AARectangle aabb1( aabb.llc(), aabb.urc() ); 
   AARectangle aabb2( aabb.llc(), aabb.urc() ); 
   aabb.split( level%3, aabb1, aabb2 );


   if( node1 != 0 ){ delete node1; }
   node1 = new IntervalTree( aabb1.llc(), aabb1.urc(), level+1, max_levels, min_objects );
   // move objects into this node:
   for( size_t i=0;i<object_list.size();i++ )
   {
      if( aabb1.intersects( *object_list[i] ) )
      {
         node1->addObject( object_list[i] );
      }
   }
   if( node1->nbObjects() == 0 ){ delete node1; node1 = 0; }

   if( node2 != 0 ){ delete node2; }
   node2 = new IntervalTree( aabb2.llc(), aabb2.urc(), level+1, max_levels, min_objects );
   // move objects into this node:
   for( size_t i=0;i<object_list.size();i++ )
   {
      if( aabb2.intersects( *object_list[i] ) )
      {
         node2->addObject( object_list[i] );
      }
   }
   if( node2->nbObjects() == 0 ){ delete node2; node2 = 0; }

   //object_list.clear();  // dont need this list anymore, so dont take up the space.
   if( node1 ) { node1->Divide(); }
   if( node2 ) { node2->Divide(); }
}


IntervalSet::IntervalSet(  const Vector& llc, const Vector& urc, const int lvl, const int maxlvl, const int minobj ) :
   IntervalSetBase( new IntervalTree( llc, urc, lvl, maxlvl, minobj ) )
{}

IntervalSet::~IntervalSet(){}


IntervalSet lux::makeIntervalSet( const Vector& llc, const Vector& urc, const int lvl, const int maxlvl, const int minobj )
{
   return IntervalSet( llc, urc, lvl, maxlvl, minobj );
}


void lux::AddIntervalSetAABB( IntervalSet& is, AABB& aabb )
{
   is->addObject(aabb);
}

void lux::DivideIntervalSet( IntervalSet& is ){ is->Divide(); }
