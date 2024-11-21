

#include "OVDBGrid.h"
#include "Volume.h"



using namespace lux;

ScalarOGrid::ScalarOGrid( OVDBGrid<OVDB_Grid_Style<float> >* f )  : ScalarOGridBase(f) {}
ScalarOGrid::~ScalarOGrid(){}


ScalarOGrid ScalarOGrid::operator+=( const ScalarField& e2 )
{
   for( OVDBGrid<OVDB_Grid_Style<float> >::Iterator iter = (*this)->getIterator(); iter.test(); ++iter)
   {
      openvdb::Coord ijk = iter.getCoord();
	  iter.setValue( *iter + e2->eval( (*this)->evalP(ijk[0],ijk[1],ijk[2]) ) );
   }
   return *this;
}

ScalarOGrid ScalarOGrid::operator-=( const ScalarField& e2 )
{
   for( OVDBGrid<OVDB_Grid_Style<float> >::Iterator iter = (*this)->getIterator(); iter.test(); ++iter)
   {
      openvdb::Coord ijk = iter.getCoord();
	  iter.setValue( *iter - e2->eval( (*this)->evalP(ijk[0],ijk[1],ijk[2]) ) );
   }
   return *this;
}






ScalarOGrid ScalarOGrid::operator-()
{
   for( OVDBGrid<OVDB_Grid_Style<float> >::Iterator iter = (*this)->getIterator(); iter.test(); ++iter)
   {
	  iter.setValue( -(*iter) );
   }
   return *this;
}

ScalarOGrid ScalarOGrid::operator*=( const ScalarField& e2 )
{
   for( OVDBGrid<OVDB_Grid_Style<float> >::Iterator iter = (*this)->getIterator(); iter.test(); ++iter)
   {
      openvdb::Coord ijk = iter.getCoord();
	  iter.setValue( *iter * e2->eval( (*this)->evalP(ijk[0],ijk[1],ijk[2]) ) );
   }
   return *this;
}

ScalarOGrid ScalarOGrid::operator/=( const ScalarField& e2 )
{
   for( OVDBGrid<OVDB_Grid_Style<float> >::Iterator iter = (*this)->getIterator(); iter.test(); ++iter)
   {
      openvdb::Coord ijk = iter.getCoord();
	  iter.setValue( *iter / e2->eval( (*this)->evalP(ijk[0],ijk[1],ijk[2]) ) );
   }
   return *this;
}





















VectorOGrid::VectorOGrid( OVDBGrid<OVDB_Grid_Style<Vector> >* f )  : VectorOGridBase(f) {}

VectorOGrid::~VectorOGrid(){}


VectorOGrid VectorOGrid::operator+=( const VectorField& e2 )
{
   for( OVDBGrid<OVDB_Grid_Style<Vector> >::Iterator iter = (*this)->getIterator(); iter.test(); ++iter)
   {
      openvdb::Coord ijk = iter.getCoord();
	  iter.setValue( equals( equals(*iter) + e2->eval( (*this)->evalP(ijk[0],ijk[1],ijk[2]) ) ) );
   }
   return *this;
}

VectorOGrid VectorOGrid::operator-=( const VectorField& e2 )
{
   for( OVDBGrid<OVDB_Grid_Style<Vector> >::Iterator iter = (*this)->getIterator(); iter.test(); ++iter)
   {
      openvdb::Coord ijk = iter.getCoord();
	  iter.setValue( equals( equals(*iter) - e2->eval( (*this)->evalP(ijk[0],ijk[1],ijk[2]) ) ) );
   }
   return *this;
}

VectorOGrid VectorOGrid::operator-()
{
   for( OVDBGrid<OVDB_Grid_Style<Vector> >::Iterator iter = (*this)->getIterator(); iter.test(); ++iter)
   {
	  iter.setValue( equals(-equals(*iter) ) );
   }
   return *this;
}

VectorOGrid VectorOGrid::operator*=( const ScalarField& e2 )
{
   for( OVDBGrid<OVDB_Grid_Style<Vector> >::Iterator iter = (*this)->getIterator(); iter.test(); ++iter)
   {
      openvdb::Coord ijk = iter.getCoord();
	  iter.setValue( equals( equals(*iter) * e2->eval( (*this)->evalP(ijk[0],ijk[1],ijk[2]) ) ) );
   }
   return *this;
}

VectorOGrid VectorOGrid::operator/=( const ScalarField& e2 )
{
   for( OVDBGrid<OVDB_Grid_Style<Vector> >::Iterator iter = (*this)->getIterator(); iter.test(); ++iter)
   {
      openvdb::Coord ijk = iter.getCoord();
	  iter.setValue( equals( equals(*iter) / e2->eval( (*this)->evalP(ijk[0],ijk[1],ijk[2]) ) ) );
   }
   return *this;
}







ColorOGrid::ColorOGrid( OVDBGrid<OVDB_Grid_Style<Color> >* f ): ColorOGridBase(f) {}
ColorOGrid::~ColorOGrid(){}


ColorOGrid ColorOGrid::operator+=( const ColorField& e2 )
{
   for( OVDBGrid<OVDB_Grid_Style<Color> >::Iterator iter = (*this)->getIterator(); iter.test(); ++iter)
   {
      openvdb::Coord ijk = iter.getCoord();
	  iter.setValue( equals( equals(*iter) + e2->eval( (*this)->evalP(ijk[0],ijk[1],ijk[2]) ) ) );
   }
   return *this;
}

ColorOGrid ColorOGrid::operator*=( const ColorField& e2 )
{
   for( OVDBGrid<OVDB_Grid_Style<Color> >::Iterator iter = (*this)->getIterator(); iter.test(); ++iter)
   {
      openvdb::Coord ijk = iter.getCoord();
	  iter.setValue( equals( equals(*iter) * e2->eval( (*this)->evalP(ijk[0],ijk[1],ijk[2]) ) ) );
   }
   return *this;
}


ColorOGrid ColorOGrid::operator*=( const ScalarField& e2 )
{
   for( OVDBGrid<OVDB_Grid_Style<Color> >::Iterator iter = (*this)->getIterator(); iter.test(); ++iter)
   {
      openvdb::Coord ijk = iter.getCoord();
	  iter.setValue( equals( equals(*iter) * e2->eval( (*this)->evalP(ijk[0],ijk[1],ijk[2]) ) ) );
   }
   return *this;
}

ColorOGrid ColorOGrid::operator/=( const ScalarField& e2 )
{
   for( OVDBGrid<OVDB_Grid_Style<Color> >::Iterator iter = (*this)->getIterator(); iter.test(); ++iter)
   {
      openvdb::Coord ijk = iter.getCoord();
	  iter.setValue( equals( equals(*iter) / e2->eval( (*this)->evalP(ijk[0],ijk[1],ijk[2]) ) ) );
   }
   return *this;
}

