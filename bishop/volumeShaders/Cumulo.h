
#include "Vector.h"
#include "Volume.h"
#include "ImplicitVolumeShapes.h"
#include "GridVolumes.h"
#include "VolumeGeometry.h"
#include "VolumeGrid.h"

namespace lux
{

typedef Volume<float> FloatVolume;
typedef vector<FloatVolume* > FloatVolumeArray;

class Cumulo: public Volume<float>
{

  public:

    Cumulo( const Volume<float>* base_sdf, const Volume<float>* displacement, float dx, int iterations ) 
    {
       addElement( base_sdf );
       addElement( displacement );

       X = new ImplicitPointVectorVolume( base_sdf, dx, iterations );
       N = new WarpVolume( displacement, X );
       cumulo_sdf = new AddVolume( cumulo_sdf, N );
    };

    ~Cumulo(){}

    const float eval( const Vector& P ) const
    {
       return cumulo_sdf->eval(P); 
    };

  private:

    const Volume<Vector>* X;
    const Volume<float>* N;
    const Volume<float>* cumulo_sdf;
};


class GriddedCumulo: public Volume<float>
{

  public:

    GriddedCumulo( const Volume<float>* base_sdf, const FloatVolumeArray& displacements, float dx, int iterations,
            int nx, int ny, int nz, const Vector llc, const Vector dims ) 
    {
       addElement( base_sdf );
       for( size_t i=0;i<displacements.size();i++ ){ addElement( displacements[i] ); }

       cumulo_sdf = base_sdf;
       for( size_t i=0;i<displacements.size();i++ )
       {
          Volume<Vector>* X = new ImplicitPointVectorVolume( cumulo_sdf, dx, iterations );
	  VolumeGrid<Vector>* gX = new VolumeGrid<Vector>();
	  gX->init( nx, ny, nz, dims[0], dims[1], dims[2], llc );
	  Sample( gX, X );
	  delete X;
          Volume<Vector>* XX = new GriddedVectorVolume( gX );
	  Volume<float>* N = new WarpVolume( displacements[i], XX );
	  cumulo_sdf = new AddVolume( cumulo_sdf, N );
       }
    };

    ~GriddedCumulo(){}

    const float eval( const Vector& P ) const
    {
       return cumulo_sdf->eval(P); 
    };

  private:

    const Volume<float>* cumulo_sdf;
};




}
