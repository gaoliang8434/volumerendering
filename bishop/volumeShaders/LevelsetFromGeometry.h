
#include "Vector.h"
#include "Volume.h"
#include "ImplicitVolumeShapes.h"
#include "GridVolumes.h"
#include "VolumeGeometry.h"

namespace lux
{

class LevelsetFromGeometry : public Volume<float>
{

  public:

    LevelsetFromGeometry(  const lux::TriangleGeometry& geom, int nx, int ny, int nz, float lx, float ly, float lz )
    {
       lux::SignedDistance* sdf = new lux::SignedDistance(geom );
       lux::Vector origin = geom.LLC();
       lux::Vector urc = geom.URC();
       lux::Vector LL = (urc-origin);
       origin += LL/2.0;
       LL = lux::Vector( lx, ly, lz );
       origin -= LL/2.0;

       ls.init( nx, ny, nz, lx, ly, lz, origin );
       ls.setOutsideValue( -LL.magnitude() );
       lux::Sample( &ls, sdf );
       delete sdf;
    };

    ~LevelsetFromGeometry(){}

    const lux::VolumeGrid<float>* getGrid() const { return &ls; }

    const float eval( const Vector& P ) const
    {
       return ls.eval(P);
       //return sdf->eval(P);
    };

  private:

    lux::VolumeGrid<float> ls;
};

ScalarField geom2ls(  const lux::TriangleGeometry& geom, int nx, int ny, int nz, float lx, float ly, float lz )
   { return ScalarField(new LevelsetFromGeometry( geom, nx, ny, nz, lx, ly, lz )  ); }

}
