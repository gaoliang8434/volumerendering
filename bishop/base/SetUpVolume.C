


#include <algorithm>
using namespace std;

#include "CmdLineFind.h"
#include "Volume.h"
#include "VolumeGrid.h"
#include "ImplicitVolumeShapes.h"
#include "MoreImplicitVolumes.h"
#include "Noise.h"
#include "PerlinNoise.h"
#include "ImplicitColors.h"
#include "RayMarcher.h"
#include "VolumeGeometry.h"
#include "ObjParser.h"
#include "FFTNoise.h"
#include "VolumeDots.h"
#include "UniformPRN.h"
#include "LognormalPRN.h"
#include "Stamp.h"
#include "GridVolumes.h"
#include "Ballistics.h"
#include "PoissonSolvers.h"
#include "cfd.h"
#include "SparseGrid.h"
#include "Tracking.h"


using namespace std;
using namespace lux;


Volume<float>* SetUpPyroSphere( CmdLineFind& clf )
{
   vector<float> tvalue;
   tvalue.push_back( 0.0 );
   tvalue.push_back( 0.0 );
   tvalue.push_back( 0.0 );

   tvalue = clf.findArray( "-translate", tvalue, "Noise translation");
   Vector trans( tvalue[0], tvalue[1], tvalue[2] );
   float freq = clf.find("-freq", 1.432234f, "Noise frequency");
   float amp = clf.find("-amp", 1.5f, "Noise amplitude" );
   float rough = clf.find("-rough", 0.5f, "Noise roughness" );
   float octaves = clf.find("-octaves", 3.0f, "Noise octaves" );
   float pyrogamma = clf.find("-pyrogamma", 0.33333f, "Noise octaves" );
   float pyrotime = clf.find("-time", 0.0f, "Noise time" );
   float timescale = clf.find("-timescale", 1.0f, "Scale factor on time");
   pyrotime *= timescale;


   Volume<float>* v1 = new PyroclasticVolume( Vector(0,0,0), 2.5, amp, octaves, freq, rough, trans, pyrotime, pyrogamma ); 
   return v1;
}




Volume<float>* SetUpSphere( CmdLineFind& clf )
{
   float rad = clf.find( "-radius", 1.0f, "Radius of sphere" );
   Vector center = clf.find("-center", Vector(0,0,0), "Center of sphere");
   Volume<float>* v1 = new SphereVolume( center, rad ); 
   return v1;
}


Volume<float>* SetUpBretzel2( CmdLineFind& clf )
{
   float c = clf.find("-bretzel2", 1.0f );
   Volume<float>* v1 = new Bretzel2Volume( c ); 
   return v1;
}


Volume<float>* SetUpMobiusStrip( CmdLineFind& clf )
{
   float rad = clf.find( "-radius", 1.0f, "Radius of strip" );
   float thick = clf.find( "-thick", 0.0f, "Thickness of strip" );
   Vector center = clf.find("-center", Vector(0,0,0), "Center of strip");
   Vector axis = clf.find("-axis", Vector(0,1,0), "Axis of strip");
   Volume<float>* v1 = new MobiusStripVolume( center, axis, rad, thick ); 
   return v1;
}


Volume<float>* SetUpCylinder( CmdLineFind& clf )
{
   float rad = clf.find( "-radius", 1.0f, "Radius of cylinder" );
   float length = clf.find( "-length", 5.0f, "length of cylinder" );
   Vector center = clf.find("-center", Vector(0,0,0), "Center of cylinder");
   Vector axis = clf.find("-axis", Vector(0,1,0), "Axis of cylinder");
   Volume<float>* v1 = new ImplicitCylinder( center, axis, length, rad ); 
   return v1;
}



Volume<float>* SetUpTube( CmdLineFind& clf )
{
   float rad = clf.find( "-outerradius", 1.0f, "Outer radius of tube" );
   float innerrad = clf.find( "-innerradius", 0.3f, "Inner radius of tube" );
   float length = clf.find( "-length", 5.0f, "length of cylinder" );
   Vector center = clf.find("-center", Vector(0,0,0), "Center of cylinder");
   Vector axis = clf.find("-axis", Vector(0,1,0), "Axis of cylinder");
   Volume<float>* v1 = new ImplicitCylinder( center, axis, length, rad ); 
   Volume<float>* v2 = new ImplicitCylinder( center, axis, length, innerrad );
   Volume<float>* v3 = new CutoutVolume( v1, v2 );
   return v3;
}



Volume<float>* SetUpEllipse( CmdLineFind& clf )
{
   float majorrad = clf.find( "-majorradius", 3.0f, "Major radius of ellipse" );
   float minorrad = clf.find( "-minorradius", 1.0f, "Minor radius of ellipse" );
   float gamma = clf.find( "-gamma", 0.0f, "Gamma on the ellipse" );

   vector<float> tvalue;
   tvalue.push_back( 1.0 );
   tvalue.push_back( 0.0 );
   tvalue.push_back( 0.0 );

   tvalue = clf.findArray( "-axis", tvalue, "Ellipse axis");
   Vector axis( tvalue[0], tvalue[1], tvalue[2] );

   Volume<float>* v1 = new EllipseVolume( Vector(0,0,0), axis, majorrad, minorrad ); 

   if( gamma > 0.0 )
   {
      v1 = new MultiplyVolume( v1, new GammaVolume( v1, gamma-1 ) );
   }

   return v1;
}

Volume<float>* SetUpJack( CmdLineFind& clf )
{
   float majorrad = clf.find( "-majorradius", 3.0f, "Major radius of ellipse" );
   float minorrad = clf.find( "-minorradius", 1.0f, "Minor radius of ellipse" );
   Vector center(0,0,0);


   Vector axis0( 1,1,1 );
   Vector axis1( 1,-1,0 );
   Vector axis2 = axis0^axis1;
   axis0.normalize();
   axis1.normalize();
   axis2.normalize();
   Volume<float>* v1 = new EllipseVolume( center, axis0, majorrad*1.3, minorrad );
   v1 = new UnionVolume( v1, new EllipseVolume( center, axis1, majorrad, minorrad ) );
   v1 = new UnionVolume( v1, new EllipseVolume( center, axis2, majorrad, minorrad ) );

   v1 = new UnionVolume( v1, new SphereVolume( center+majorrad*axis1, majorrad/3.0 ) );
   v1 = new UnionVolume( v1, new SphereVolume( center-majorrad*axis1, majorrad/3.0 ) );

   v1 = new UnionVolume( v1, new SphereVolume( center+majorrad*axis2, majorrad/3.0 ) );
   v1 = new UnionVolume( v1, new SphereVolume( center-majorrad*axis2, majorrad/3.0 ) );


   return v1;
}

Volume<float>* SetUpBlendJack( CmdLineFind& clf )
{
   float majorrad = clf.find( "-majorradius", 3.0f, "Major radius of ellipse" );
   float minorrad = clf.find( "-minorradius", 1.0f, "Minor radius of ellipse" );
   float blend = clf.find( "-blendfactor", 1.0f, "Blending scale factor" );
   float ballradius = majorrad/3.0;
   ballradius = clf.find( "-ballradius", ballradius );
   Vector center(0,0,0);


   Vector axis0( 1,1,1 );
   Vector axis1( 1,-1,0 );
   Vector axis2 = axis0^axis1;
   axis0.normalize();
   axis1.normalize();
   axis2.normalize();

   std::vector<Volume<float>*> vs;

   vs.push_back( new MultiplyVolume( new EllipseVolume( center, axis0, majorrad*1.3, minorrad ), blend ) );
   vs.push_back( new MultiplyVolume( new EllipseVolume( center, axis1, majorrad, minorrad ), blend ) );
   vs.push_back( new MultiplyVolume( new EllipseVolume( center, axis2, majorrad, minorrad ), blend ) );

   vs.push_back( new MultiplyVolume( new SphereVolume( center+majorrad*axis1, ballradius ), blend ) );
   vs.push_back( new MultiplyVolume( new SphereVolume( center-majorrad*axis1, ballradius ), blend ) );

   vs.push_back( new MultiplyVolume( new SphereVolume( center+majorrad*axis2, ballradius ), blend ) );
   vs.push_back( new MultiplyVolume( new SphereVolume( center-majorrad*axis2, ballradius ), blend ) );

   Volume<float>* v1 = new MultiBlendVolume( vs );

   return v1;
}




Volume<float>* SetUpBox( CmdLineFind& clf )
{
   float length = clf.find("-length", 1.0f, "Size of box");
   float power = clf.find("-smoothing", 3.0f, "Control how smooth the corners are");
   
   Volume<float>* v1 = new CsgBoxVolume( Vector(0,0,0), length, power); 
   return v1;
}


Volume<float>* SetUpCone( CmdLineFind& clf )
{
   float length = clf.find("-length", 1.0f, "Size of cone");
   float angle = clf.find("-angle", 30.0f, "Angle of cone in degrees");

   vector<float> tvalue;
   tvalue.push_back( 0.0 );
   tvalue.push_back( 1.0 );
   tvalue.push_back( 0.0 );

   tvalue = clf.findArray( "-axis", tvalue, "axis of cone");
   Vector axis( tvalue[0], tvalue[1], tvalue[2] );

   tvalue.clear();
   tvalue.push_back( 0.0 );
   tvalue.push_back( -length*0.5 );
   tvalue.push_back( 0.0 );

   tvalue = clf.findArray( "-tip", tvalue, "position of tip of cone");
   Vector tip( tvalue[0], tvalue[1], tvalue[2] );
   
   Volume<float>* v1 = new ConeVolume( tip, axis, length, angle); 
   return v1;
}





Volume<float>* SetUpSteinerPatch( CmdLineFind& clf )
{
   Volume<float>* v1 = new SteinerPatchVolume(); 
   return v1;
}


Volume<float>* SetUpIcosahedron( CmdLineFind& clf )
{
   Volume<float>* v1 = new IcosahedronVolume(); 
   return v1;
}




Volume<float>* SetUpTorus( CmdLineFind& clf )
{
   float majorrad = clf.find( "-majorradius", 3.0f, "Major radius of torus" );
   float minorrad = clf.find( "-minorradius", 0.5f, "Minor radius of torus" );


   vector<float> tvalue;
   tvalue.push_back( 0.0 );
   tvalue.push_back( 0.0 );
   tvalue.push_back( 0.0 );

   tvalue = clf.findArray( "-toruscenter", tvalue, "Location of the torus");
   Vector center( tvalue[0], tvalue[1], tvalue[2] );

   tvalue.clear();
   tvalue.push_back( 1.0 );
   tvalue.push_back( 1.0 );
   tvalue.push_back( 1.0 );

   tvalue = clf.findArray( "-torusaxis", tvalue, "Axis of the torus");
   Vector axis( tvalue[0], tvalue[1], tvalue[2] );
   axis.normalize();

   Volume<float>* v1 = new TorusVolume( center, axis, majorrad, minorrad ); 
   return v1;
}


Volume<float>* SetUpTorii( CmdLineFind& clf )
{
   int nbtorii = clf.find( "-nbtorii", 2, "Number of torii" );
   float majorrad = clf.find( "-majorradius", 3.0f, "Major radius of torus" );
   float minorrad = clf.find( "-minorradius", 0.5f, "Minor radius of torus" );
   float blendscale = clf.find( "-blendscale", 1.0f, "Blending scale" );
   bool doBlend = clf.findFlag("-blendtorii" );

   if( doBlend ){ cout << "Blending torii\n" << flush; }

   Volume<float>* tor = 0;
   for( int i=0;i<nbtorii;i++ )
   {
      Vector center( drand48()-0.5, drand48()-0.5, drand48()-0.5 );
      center *= 4.0;
      Vector axis( drand48()-0.5, drand48()-0.5, drand48()-0.5 );
      axis.normalize();
      float minrad = minorrad * drand48();
      float majrad = majorrad * drand48() + minrad;
      if( i==0 )
      { 
         tor = new TorusVolume( center, axis, majrad, minrad ); 
      }
      else
      {
         if( doBlend )
	 {
            tor = new BlinnBlendVolume( tor, new MultiplyVolume( new TorusVolume( center, axis, majrad, minrad ), blendscale ) );
	 }
	 else
	 {
            tor = new UnionVolume( tor, new TorusVolume( center, axis, majrad, minrad ) );
	 }
      }
   }
   return tor;
}

Volume<float>* SetUpPyramid( CmdLineFind& clf )
{
   
   vector<float> tvalue;
   tvalue.push_back( 0.0 );
   tvalue.push_back( 0.0 );
   tvalue.push_back( 0.0 );

   tvalue = clf.findArray( "-pyramidcenter", tvalue, "Location of the pyramid");
   Vector center( tvalue[0], tvalue[1], tvalue[2] );

   float length = clf.find( "-pyramidlength", 1.0f, "Length of the sides of the pyramid" );

   Volume<float>* v = ProcessPyramid( length, center );

   return v;
}


Volume<float>* SetUpObj( CmdLineFind& clf )
{
    string filename = clf.find( "-objfilename", "", "Name of obj file with triangular geometry" );
    ObjParser p;
    if( !p.ParseFile( filename ) )
    {
       cout << "Could not read file " << filename << endl << flush;
       return SetUpSphere( clf );
    }

    TriangleGeometry g;
    float scaling = clf.find("-objscale", 1.0f, "Scale vertices of obj file" );
    g.setScaling( scaling );
    if( !p.Fill( g ) )
    {
       cout << "Could not read geometry from file " << filename << endl << flush;
       return SetUpSphere( clf );
    }

    /*
    cout << "VERTICES:\n";
    for( int i=0;i<g.nbVertices();i++ )
    {
       const Vector& v = g.getVertex(i);
       cout << i << "  " << v[0] << " " << v[1] << " " << v[2] << endl;
    }

    cout << "TRIANGLES:\n";
    int iif, jjf, kkf;
    for( int i=0;i<g.nbFaces();i++ )
    {
       g.getFace(i, iif,jjf,kkf);
       cout << i << "  " << iif << " " << jjf << " " << kkf << endl;
    }
    */





    Vector dims = g.URC() - g.LLC();
    Vector center = (g.URC() + g.LLC())*0.5;
    dims *= 2.0;
    Vector llc = center - dims*0.5;
    Vector urc = llc + dims;
    int nx = 100;
    int ny = nx * ( dims[1]/dims[0] );
    int nz = nx * ( dims[2]/dims[0] );

    cout << "Obj faces: " << g.nbFaces() << endl;
    cout << "Obj vertices: " << g.nbVertices() << endl;
    cout << "Obj BB:\n";
    cout << llc[0] << " " <<  llc[1] << " " <<  llc[2] << "  X  " <<  urc[0] << " " <<  urc[1] << " " <<  urc[2] << endl << flush;

    nx  = clf.find("-objnx", nx, "Grid poins for obj levelset");
    ny  = clf.find("-objny", ny, "Grid poins for obj levelset");
    nz  = clf.find("-objnz", nz, "Grid poins for obj levelset");
    
    cout << "Level set dimensions:\n";
    cout << nx << " " <<  ny << " " <<  nz << endl << flush;

    VolumeGrid<float>* levelset = new VolumeGrid<float>();
    levelset->init( nx, ny, nz, dims[0], dims[1], dims[2], llc );
    levelset->setClearValue( -dims.magnitude() );

    float lsthresh = clf.find("-lsthresh", 0.1f, "Threshold for filtering level set intersections");
    RayMarchLevelSet( g, *levelset, lsthresh );
    Volume<float>* v = new GriddedVolume(levelset);

   return v;
}



Volume<float>* SetUpBlinnSpheres( CmdLineFind& clf )
{
   float rad = clf.find( "-radius", 1.0f, "Radius of Blinn spheres" );
   float separation = rad*1.5;
   separation = clf.find( "-separation", separation, "Separation distance of Blinn spheres" );
   float scale = clf.find( "-blendscale", 1.0f, "Blending scale of spheres" );
   Volume<float>* v1 = new SphereVolume( Vector(-separation*0.5,0,0), rad );
   v1 = new MultiplyVolume( v1, 1.0/scale );
   Volume<float>* v2 = new SphereVolume( Vector( separation*0.5,0,0), rad ); 
   v2 = new MultiplyVolume( v2, 1.0/scale );
   Volume<float>* v3 = new BlinnBlendVolume( v1, v2 );
   return v3;
}



Volume<float>* SetUpBlinnPyroSpheres( CmdLineFind& clf )
{

   vector<float> tvalue;
   tvalue.push_back( 0.0 );
   tvalue.push_back( 0.0 );
   tvalue.push_back( 0.0 );

   tvalue = clf.findArray( "-translate", tvalue, "Noise translation");
   Vector trans( tvalue[0], tvalue[1], tvalue[2] );
   float freq = clf.find("-freq", 1.432234f, "Noise frequency");
   float amp = clf.find("-amp", 1.5f, "Noise amplitude" );
   float rough = clf.find("-rough", 0.5f, "Noise roughness" );
   int octaves = clf.find("-octaves", 3, "Number of pyro octaves");

   float rad = clf.find( "-radius", 1.0f, "Radius of Blinn spheres" );
   float separation = rad*1.5;
   separation = clf.find( "-separation", separation, "Separation distance of Blinn spheres" );
   float scale = clf.find( "-blendscale", 1.0f, "Blending scale of spheres" );
   Volume<float>* v1 = new PyroclasticVolume( Vector(-separation*0.5,0,0), rad, amp, octaves, freq, rough, trans, 0.0 ); 
   v1 = new MultiplyVolume( v1, 1.0/scale );
   Volume<float>* v2 = new PyroclasticVolume( Vector(separation*0.5,0,0), rad, amp, octaves, freq, rough, trans + Vector(separation,0,0), 0.0 ); 
   v2 = new MultiplyVolume( v2, 1.0/scale );
   Volume<float>* v3 = new BlinnBlendVolume( v1, v2 );
   return v3;
}


Volume<float>* SetUpCutoutSpheres( CmdLineFind& clf )
{

   float rad = clf.find( "-radius", 1.0f, "Radius of cutout spheres" );
   float separation = 0.75*rad;
   separation = clf.find( "-separation", separation, "Separation distance of cutout spheres" );
   Volume<float>* v1 = new SphereVolume( Vector(-separation*0.5,0,0), rad ); 
   Volume<float>* v2 = new SphereVolume( Vector(separation*0.5,0,0), rad ); 
   Volume<float>* v3 = new CutoutVolume( v1, v2 );
   return v3;
}


Volume<float>* SetUpCSG( CmdLineFind& clf )
{

   float rad = clf.find( "-radius", 2.0f, "Radius of csg spheres" );
   float separation = 1.5*rad;
   separation = clf.find( "-separation", separation, "Separation distance of csg spheres" );
   Volume<float>* v1 = new SphereVolume( Vector(-separation*0.5,0,0), rad ); 
   Volume<float>* v2 = new SphereVolume( Vector(separation*0.5,0,0), rad ); 

   Volume<float>* v3 = 0;
   vector<string> csgMenu;
   csgMenu.push_back( "union" );
   csgMenu.push_back( "intersect" );
   csgMenu.push_back( "blend" );
   csgMenu.push_back( "cutout" );

   string csgtype = clf.findMenu( "-csg", csgMenu, "Select a predefined csg operation" );

   if( csgtype == "intersect" )
   {
      v3 = new IntersectionVolume( v1, v2 );
   }
   else if( csgtype == "blend" )
   {
      v3 = new BlinnBlendVolume( v1, v2 );
   }
   else if( csgtype == "cutout" )
   {
      v3 = new CutoutVolume( v1, v2 );
   }
   else
   {
      v3 = new UnionVolume( v1, v2 );
   }

   return v3;
}



Volume<float>* SetUpCSGSphereTorus( CmdLineFind& clf )
{

   Volume<float>* v2 = SetUpSphere( clf ); 
   Volume<float>* v1 = SetUpTorus( clf ); 

   Volume<float>* v3 = 0;
   vector<string> csgMenu;
   csgMenu.push_back( "union" );
   csgMenu.push_back( "intersect" );
   csgMenu.push_back( "blend" );
   csgMenu.push_back( "cutout" );
   csgMenu.push_back( "othercutout" );

   string csgtype = clf.findMenu( "-csg", csgMenu, "Select a predefined csg operation" );

   if( csgtype == "intersect" )
   {
      v3 = new IntersectionVolume( v1, v2 );
   }
   else if( csgtype == "blend" )
   {
      v3 = new BlinnBlendVolume( v1, v2 );
   }
   else if( csgtype == "cutout" )
   {
      v3 = new CutoutVolume( v1, v2 );
   }
   else if( csgtype == "othercutout" )
   {
      v3 = new CutoutVolume( v2, v1 );
   }
   else
   {
      v3 = new UnionVolume( v1, v2 );
   }

   return v3;
}






/*
Volume<float>* SetUpFFT( CmdLineFind& clf )
{
   Noise_t parms;
   parms.fftNbGridPoints = 256;
   parms.fftPower = clf.find("-fftpower", parms.fftPower );
   return new MultiplyVolume( new MaskVolume(SetUpSphere( clf )), new FFTNoiseVolume( parms.fftPower, parms.fftLowCutoff, parms.fftHighCutoff, parms.fftLength, parms.fftNbGridPoints) );
}
*/

/*
void SetUpStampNoise( CmdLineFind& clf, RenderData& data )
{
   int nx = clf.find( "-stampnx", 100, "Volume size in X direction" );
   int ny = clf.find( "-stampny", 100, "Volume size in Y direction" );
   int nz = clf.find( "-stampnz", 100, "Volume size in Z direction" );

   float dx = clf.find( "-stampdx", 0.1f, "Volume cell size in X direction");
   float dy = clf.find( "-stampdy", 0.1f, "Volume cell size in Y direction");
   float dz = clf.find( "-stampdz", 0.1f, "Volume cell size in Z direction");

   float x0 = clf.find( "-stampllcx", -5.0f, "Volume llc in X direction" );
   float y0 = clf.find( "-stampllcy", -5.0f, "Volume llc in Y direction" );
   float z0 = clf.find( "-stampllcz", -5.0f, "Volume llc in Z direction" );

   VolumeGrid<float>* densityGrid = new VolumeGrid<float>;
   densityGrid->init( nx, ny, nz, dx*nx, dy*ny, dz*nz, Vector(x0,y0,z0) );

   FractalSum<PerlinNoise> pn;
   float radius = clf.find("-stampradius", 2.0f );

   Noise_t parms;
   parms.wavelength= 3.0/radius;
   pn.setParameters(parms);

   StampNoise( *densityGrid, &pn, Vector(0,0,0), radius );
 
   data.densityField = new GriddedVolume( densityGrid );
}
*/

void SetUpCauliflower( CmdLineFind& clf, RenderData& d )
{
   int nx = clf.find( "-dotnx", 100, "Volume size in X direction" );
   int ny = clf.find( "-dotny", 100, "Volume size in Y direction" );
   int nz = clf.find( "-dotnz", 100, "Volume size in Z direction" );

   float dx = clf.find( "-dotdx", 0.1f, "Volume cell size in X direction");
   float dy = clf.find( "-dotdy", 0.1f, "Volume cell size in Y direction");
   float dz = clf.find( "-dotdz", 0.1f, "Volume cell size in Z direction");

   float x0 = clf.find( "-dotllcx", -5.0f, "Volume llc in X direction" );
   float y0 = clf.find( "-dotllcy", -5.0f, "Volume llc in Y direction" );
   float z0 = clf.find( "-dotllcz", -5.0f, "Volume llc in Z direction" );

   VolumeGrid<float>* densityGrid = new VolumeGrid<float>;
   VolumeGrid<Color>* colorGrid = new VolumeGrid<Color>;

   densityGrid->init( nx, ny, nz, dx*nx, dy*ny, dz*nz, Vector(x0,y0,z0) );
   colorGrid->init( nx, ny, nz, dx*nx, dy*ny, dz*nz, Vector(x0,y0,z0) );




   VolumeDots dots( densityGrid, colorGrid  ); 
   Vector guideposition;
   float childradius = clf.find( "-childradius", 0.5f, "Size of cauliflower children");;
   int nbcauliflowerclumps = clf.find("-nbcauliflowerclumps", 5, "Number of cauliflower clumps");
   float cauliflowerscale = clf.find("-cauliflowerscale", 0.65f );
   int nbcauliflowerrecursions = clf.find("-nbcauliflowerrecursions", 4 );
   int nbSpheres = clf.find("-nbcauliflowers", 10, "Number of cauliflower clusters");
   float sphereRadius = clf.find("-dotradius", 1.0f, "Maximum radius of any of the spheres");

   UniformPRN unoise;
   LognormalPRN lnoise;
   Noise_t parms;
   parms.seed = clf.find("-seed", 48473, "Cauliflower seed" );
   unoise.setParameters( parms );

   Vector dims = Vector( densityGrid->Lx(), densityGrid->Ly(), densityGrid->Lz() );
   dims *= 0.7;

   for( int i=0;i<nbSpheres;i++ )
   {
      float cx = (unoise.eval()-0.5) * dims[0];
      float cy = (unoise.eval()-0.5) * dims[1];
      float cz = (unoise.eval()-0.5) * dims[2];

      guideposition = Vector( cx, cy, cz );
      dots.setAtt( "position", guideposition );

      float rad = sphereRadius;
      childradius = rad;
      dots.setAtt("radiusPP", rad );

      float r = unoise.eval();
      float g = unoise.eval();
      float b = unoise.eval();

      dots.setAtt("color", Vector(r,g,b) );

      setChildParticle_Cauliflower( guideposition, childradius, nbcauliflowerclumps, cauliflowerscale, nbcauliflowerrecursions, dots );
   }

   d.densityField = new GriddedVolume( densityGrid );
   d.colorField = new GriddedColor( colorGrid );
}

Volume<float>* SetUpBoxTorus( CmdLineFind& clf, RenderData& data )
{
   Volume<float>* box = new CsgBoxVolume( Vector(0,0,0), 4.0, 3.0 );
   Volume<float>* torus = new TorusVolume( Vector( 2.0,0,0) , Vector(0.3,0,1), 2.5, 0.5 );

   Volume<float>* toruswithboxcutout = new CutoutVolume( torus, box );

   Volume<Color>* boxcolor = new ConstantColor( Color(0,0.7,0,1) );  
   Volume<Color>* toruscolor = new ConstantColor( Color(0,0,0.8,1) ); 

   Volume<Color>* boxcolorclipped = new VolumeMultiplyColor( boxcolor, new MaskVolume(box) );
   Volume<Color>* toruscolorclipped = new VolumeMultiplyColor( toruscolor, new MaskVolume(toruswithboxcutout) );

   float thresh = clf.find("-clamp", 0.1f, "clamp level");

   Vector rot = clf.find("-rotate", Vector(-M_PI/6.0,0,0) );
   data.colorField = new RotateColor( new AddColor( boxcolorclipped, toruscolorclipped ), rot ) ;
   Volume<float>* v = new RotateVolume( new UnionVolume( box, torus ), rot );
   data.densityField = new ClampVolume( new MultiplyVolume( v, 1.0/thresh), 0.0, 1.0 );
   return v;
}







void SetUpFireBall( CmdLineFind& clf, RenderData& data )
{
   vector<float> tvalue;
   tvalue.push_back( 0.0 );
   tvalue.push_back( 0.0 );
   tvalue.push_back( 0.0 );

   tvalue = clf.findArray( "-translate", tvalue, "Noise translation");
   Vector trans( tvalue[0], tvalue[1], tvalue[2] );
   float freq = clf.find("-freq", 1.432234f, "Noise frequency");
   float amp = clf.find("-amp", 1.5f, "Noise amplitude" );
   float rough = clf.find("-rough", 0.5f, "Noise roughness" );
   float octaves = clf.find("-octaves", 3.0f, "Noise octaves" );
   float pyrogamma = clf.find("-pyrogamma", 0.33333f, "Noise octaves" );
   float pyrotime = clf.find("-time", 0.0f, "Noise time");
   pyrotime *= clf.find("-timescale", 1.0f );

   Volume<float>* v1 = new PyroclasticVolume( Vector(0,0,0), 2.5, amp, octaves, freq, rough, trans, pyrotime, pyrogamma ); 
   float thresh = clf.find("-fireballclamp", 0.1f, "Clamp threshold for fireball");
   data.densityField = new ClampVolume( new MultiplyVolume( v1, 1.0/thresh ), 0.0, 1.0 );

   float burnthickness = clf.find("-fireballthickness", 0.3f );
   Volume<float>* burnfield = new ClampVolume( new MultiplyVolume( v1, 1.0/burnthickness ), 0.0, 1.0 );

   vector<int> NNN;
   NNN.push_back( 300 );
   NNN.push_back( 300 );
   NNN.push_back( 300 );
   NNN = clf.findArray( "-fireballNxNyNz", NNN, "Number of grid points in fireball volume" );

   Vector lX = clf.find( "-fireballize", Vector(10,10,10), "Volume size of fireball volume");
   Vector X0 = clf.find( "-fireballorigin", Vector(-5,-5,-5), "Origin of fireball volume");

   VolumeGrid<float>* Kgrid = new VolumeGrid<float>;

   Kgrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );

   float gritlevel = clf.find("-grit", 0.1f );
   Volume<float>* K = new MultiplyVolume( new DivergenceVectorVolume( new UnitVectorVolume( new GradientVectorVolume( burnfield ) ) ), gritlevel );


   Volume<float>* g1 = new ClampVolume( new ExpVolume( K ), 0.01, 100.0);
   lux::Sample( Kgrid, g1 );
   Blur( *Kgrid );
   Blur( *Kgrid );
   Blur( *Kgrid );
   g1 = new GriddedVolume( Kgrid );


   Noise_t aparms;
   aparms.translate = clf.find( "-vtranslate", Vector(0,0,0), "Noise translation");
   aparms.wavelength = clf.find("-vfreq", 1.432234f, "Noise frequency");
   aparms.amplitude = clf.find("-vamp", 1.5f, "Noise amplitude" );
   aparms.roughness = clf.find("-vrough", 0.5f, "Noise roughness" );
   aparms.octaves = clf.find("-voctaves", 3.0f, "Noise octaves" );

   Noise * anoise = new FractalSum<PerlinNoiseGustavson>();
   anoise->setParameters(aparms);

   Volume<Vector>* velocity = new NoiseVectorVolume( anoise );

   float dt = clf.find("-dt", 0.1f,"Advection time step");
   int nbSteps = clf.find("-nbsteps", 0, "Number of advection steps");

   float ddt = dt/nbSteps;

   for( int i=0;i<nbSteps;i++ )
   {
      g1 = new AdvectVolume( g1, velocity, ddt );
      v1 = new AdvectVolume( v1, velocity, ddt );
   }




   data.ambientColorField = new VolumeMultiplyColor( new AddColor( new VolumeMultiplyColor( new ConstantColor( Color(0.2,0.16,0,1) ), g1 ),  new ConstantColor( Color( 0.64, 0.2, 0, 1) )  ), new MaskVolume( v1 ) ) ;
   data.colorField = new ConstantColor( Color(0,0,0,1) );
}



void SetUpFireBallAndSmoke( CmdLineFind& clf, RenderData& data )
{
   vector<float> tvalue;
   tvalue.push_back( 0.0 );
   tvalue.push_back( 0.0 );
   tvalue.push_back( 0.0 );

   tvalue = clf.findArray( "-translate", tvalue, "Noise translation");
   Vector trans( tvalue[0], tvalue[1], tvalue[2] );
   float freq = clf.find("-freq", 1.432234f, "Noise frequency");
   float amp = clf.find("-amp", 1.5f, "Noise amplitude" );
   float rough = clf.find("-rough", 0.5f, "Noise roughness" );
   float octaves = clf.find("-octaves", 3.0f, "Noise octaves" );
   float pyrogamma = clf.find("-pyrogamma", 0.33333f, "Noise octaves" );
   float pyroradius = clf.find("-fireradius", 2.5f );
   float pyrotime = clf.find("-pyrotime", 0.0f, "Noise time");

   Volume<float>* v1 = new PyroclasticVolume( Vector(0,0,0), pyroradius, amp, octaves, freq, rough, trans, pyrotime, pyrogamma ); 
   float thresh = clf.find("-fireballclamp", 0.1f, "Clamp threshold for fireball");
   data.densityField = new ClampVolume( new MultiplyVolume( v1, 1.0/thresh ), 0.0, 1.0 );

   float burnthickness = clf.find("-fireballthickness", 0.3f );
   Volume<float>* burnfield = new ClampVolume( new MultiplyVolume( v1, 1.0/burnthickness ), 0.0, 1.0 );

   vector<int> NNN;
   NNN.push_back( 300 );
   NNN.push_back( 300 );
   NNN.push_back( 300 );
   NNN = clf.findArray( "-fireballNxNyNz", NNN, "Number of grid points in fireball volume" );

   Vector lX = clf.find( "-fireballize", Vector(10,10,10), "Volume size of fireball volume");
   Vector X0 = clf.find( "-fireballorigin", Vector(-5,-5,-5), "Origin of fireball volume");

   VolumeGrid<float>* Kgrid = new VolumeGrid<float>;

   Kgrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );

   float gritlevel = clf.find("-grit", 0.1f );
   Volume<float>* K = new MultiplyVolume( new DivergenceVectorVolume( new UnitVectorVolume( new GradientVectorVolume( burnfield ) ) ), gritlevel );


   Volume<float>* g1 = new ClampVolume( new ExpVolume( K ), 0.01, 100.0);
   lux::Sample( Kgrid, g1 );
   Blur( *Kgrid );
   Blur( *Kgrid );
   Blur( *Kgrid );
   g1 = new GriddedVolume( Kgrid );

   data.ambientColorField = new VolumeMultiplyColor( new AddColor( new VolumeMultiplyColor( new ConstantColor( Color(0.1,0.08,0,1) ), g1 ),  new ConstantColor( Color( 0.32, 0.06, 0, 1) )  ), new MaskVolume( v1 ) ) ;
   data.colorField = new ConstantColor( Color(0,0,0,1) );



   float smokegamma = clf.find( "-smokegamma", (float)(pyrogamma*2.0),"Smoke gamma" );
   float smokeradius = clf.find("-smokeradius", pyroradius);
   float smokeamp = clf.find("-smokeamp", (float)(amp*1.3));
   float smokefreq = clf.find("-smokefreq", freq);
   float smokerough = clf.find("-smokerough", rough);
   float smokeoctaves = clf.find("-smokeoctaves", octaves );
   float smoketime = clf.find("-smoketime", 0.0f, "Noise time");
   Volume<float>* smokefield = new PyroclasticVolume( Vector(0,0,0), smokeradius, smokeamp, smokeoctaves, smokefreq, smokerough, trans, smoketime, smokegamma ); 

   float smokeclamp = clf.find( "-smokeclamp", 0.1f );
   float densityScale = clf.find("-smokedensity", 10.0f );
   Volume<float>* smokedensity = new MultiplyVolume( new ClampVolume( new MultiplyVolume( smokefield, 1.0/smokeclamp ), 0.0, 1.0 ), densityScale );



   Volume<float>* smokemask = new MaskVolume( smokefield );
   Volume<float>* smokemaskcomplement = new MaskVolume( new NegateVolume(smokefield) );

   data.densityField =  new AddVolume( new MultiplyVolume( data.densityField, smokemaskcomplement ), smokedensity );
   data.ambientColorField = new VolumeMultiplyColor( data.ambientColorField , smokemaskcomplement ) ;
   data.colorField = new VolumeMultiplyColor( new ConstantColor( Color(1,1,1,1) ), smokemask );

}





void SetUpFireBallIntoSmoke( CmdLineFind& clf, RenderData& data )
{
   vector<float> tvalue;
   tvalue.push_back( 0.0 );
   tvalue.push_back( 0.0 );
   tvalue.push_back( 0.0 );

   float trans = clf.find( "-radialtranslate", 0.0f, "Noise translation");
   float freq = clf.find("-freq", 1.432234f, "Noise frequency");
   float amp = clf.find("-amp", 1.5f, "Noise amplitude" );
   float rough = clf.find("-rough", 0.5f, "Noise roughness" );
   float octaves = clf.find("-octaves", 3.0f, "Noise octaves" );
   float pyrogamma = clf.find("-pyrogamma", 0.33333f, "Noise octaves" );
   float pyroradius = clf.find("-fireradius", 2.5f );
   float pyrotime = clf.find("-pyrotime", 0.0f, "Noise time");

   Volume<float>* v1 = new RadialPyroclasticVolume( Vector(0,0,0), pyroradius, amp, octaves, freq, rough, trans, pyrotime, pyrogamma ); 
   float thresh = clf.find("-fireballclamp", 0.1f, "Clamp threshold for fireball");


   float smokegamma = clf.find( "-smokegamma", (float)(pyrogamma*2.0),"Smoke gamma" );
   float smokeradius = clf.find("-smokeradius", (float)(5.0*pyroradius) );
   float smokeamp = clf.find("-smokeamp", (float)(amp*1.3));
   float smokefreq = clf.find("-smokefreq", freq);
   float smokerough = clf.find("-smokerough", rough);
   float smokeoctaves = clf.find("-smokeoctaves", octaves );
   float smoketime = clf.find("-smoketime", 0.0f, "Noise time");
   Volume<float>* smokefield = new CutoutVolume( v1, new SphereVolume( Vector(0,0,0), smokeradius )  );
   v1 = new CutoutVolume(  v1, smokefield );
   float smokeclamp = clf.find( "-smokeclamp", 0.1f );
   float densityScale = clf.find("-smokedensity", 10.0f );
   Volume<float>* smokedensity = new MultiplyVolume( new ClampVolume( new MultiplyVolume( smokefield, 1.0/smokeclamp ), 0.0, 1.0 ), densityScale );


   float burnthickness = clf.find("-fireballthickness", 0.3f );
   Volume<float>* burnfield = new ClampVolume( new MultiplyVolume( v1, 1.0/burnthickness ), 0.0, 1.0 );

   vector<int> NNN;
   NNN.push_back( 300 );
   NNN.push_back( 300 );
   NNN.push_back( 300 );
   NNN = clf.findArray( "-fireballNxNyNz", NNN, "Number of grid points in fireball volume" );

   Vector lX = clf.find( "-fireballize", Vector(10,10,10), "Volume size of fireball volume");
   Vector X0 = clf.find( "-fireballorigin", Vector(-5,-5,-5), "Origin of fireball volume");

   VolumeGrid<float>* Kgrid = new VolumeGrid<float>;

   Kgrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );

   float gritlevel = clf.find("-grit", 0.1f );
   Volume<float>* K = new MultiplyVolume( new DivergenceVectorVolume( new UnitVectorVolume( new GradientVectorVolume( burnfield ) ) ), gritlevel );


   Volume<float>* g1 = new ClampVolume( new ExpVolume( K ), 0.01, 100.0);
   lux::Sample( Kgrid, g1 );
   Blur( *Kgrid );
   Blur( *Kgrid );
   Blur( *Kgrid );
   g1 = new GriddedVolume( Kgrid );



   data.densityField = new AddVolume( new ClampVolume( new MultiplyVolume( v1, 1.0/thresh ), 0.0, 1.0 ), smokedensity );
   data.ambientColorField = new AddColor( new VolumeMultiplyColor( new ConstantColor( Color(0.1,0.08,0,1) ), g1 ),  new ConstantColor( Color( 0.32, 0.06, 0, 1) ));
   data.ambientColorField = new VolumeMultiplyColor( data.ambientColorField, new MaskVolume( new MultiplyVolume( smokefield, -1.0) ) );
   data.colorField = new VolumeMultiplyColor( new ConstantColor( Color(1,1,1,1) ), new MaskVolume( smokefield ) );

}









void SetUpStampedParticleSpheres( CmdLineFind& clf, RenderData& data )
{
   vector<int> NNN;
   NNN.push_back( 100 );
   NNN.push_back( 100 );
   NNN.push_back( 100 );
   NNN = clf.findArray( "-stampNxNyNz", NNN, "Number of grid points in stamped volume" );

   Vector lX = clf.find( "-stampsize", Vector(10,10,10), "Volume size of stamped volume");
   Vector X0 = clf.find( "-stamporigin", Vector(-5,-5,-5), "Origin of stamped volume");

   VolumeGrid<float>* densityGrid = new VolumeGrid<float>;
   VolumeGrid<Color>* colorGrid = new VolumeGrid<Color>;

   densityGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
   colorGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );

   ParticleGroupA particles;
   int nbParticles = clf.find("-nbparticles", 10, "Number of particles to stamp");
   float pscale = clf.find("-pscale", 3.0f, "Max particle pscale");
   float fade = clf.find("-fade", 0.1f, "Relative fade distance");
   float offset = clf.find("-offset", 1.0f, "Relative offset");
   float opacity = clf.find("-opacity", 1.0f, "Opacity");

   UniformPRN prn;
   for( int i=0;i<nbParticles;i++ )
   {
      Particle p;
      p.pscale() = prn.eval() * pscale;
      p.P() = Vector( i*lX[0]/nbParticles, 0.5*lX[1]*(1.0-prn.eval()*0.1)/0.95, 0.5*lX[2]*(1.0-prn.eval()*0.1)/0.95);
      if (p.P()[0] < p.pscale() ){ p.P()[0] = p.pscale(); }
      if (p.P()[1] < p.pscale() ){ p.P()[1] = p.pscale(); }
      if (p.P()[2] < p.pscale() ){ p.P()[2] = p.pscale(); }
      if (p.P()[0] > lX[0] - p.pscale() ){ p.P()[0] = lX[0] - p.pscale(); }
      if (p.P()[1] > lX[1] - p.pscale() ){ p.P()[1] = lX[1] - p.pscale(); }
      if (p.P()[2] > lX[2] - p.pscale() ){ p.P()[2] = lX[2] - p.pscale(); }
      p.P() += X0;
      //p.Cd() = Color( prn.eval(), prn.eval(), prn.eval(), 1.0 );
      p.Cd() = Color( 1.0, 1.0, 1.0, 1.0 );
      p.id() = i;
      p.octaves() = 6.0 * prn.eval();
      p.roughness() = prn.eval();
      p.freq() = 1.0 + prn.eval();
      p.translate() = Vector( prn.eval(), prn.eval(), prn.eval() ) * 100.0;
      p.offset() = -prn.eval() * offset;
      p.fade() = fade * prn.eval();
      p.opacity() = opacity; 
      particles.push_back( p );
   }





   FractalSum<PerlinNoise> noise;
   StampNoiseAndColor( *densityGrid, *colorGrid, &noise, particles );

   data.densityField = new GriddedVolume( densityGrid );
   data.colorField = new GriddedColor( colorGrid );

}



void SetUpStampSparseFileWisps( CmdLineFind& clf, RenderData& data )
{
   vector<int> NNNN;
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   vector<int> NNN = clf.findArray( "-stampNxNyNz", NNNN, "Number of grid points in stamped volume" );

   int sparsesize = clf.find("-sparsesize", 16 );
   Vector lX = clf.find( "-stampsize", Vector(10,10,10), "Volume size of stamped volume");
   Vector X0 = clf.find( "-stamporigin", Vector(-5,-5,-5), "Origin of stamped volume");

   float blurScale = clf.find("-blurscale", 1.0f );

   vector<string> files = clf.findMultiple("-particlefile", "" );

   vector<Vector> typeColVectors = clf.findMultiple("-particlecolor", Vector(1,1,1) );
   vector<Color> typeColors;
   for( size_t i=0;i<typeColVectors.size();i++ )
   {
      Color col( typeColVectors[i][0],  typeColVectors[i][1],  typeColVectors[i][2], 1 );
      typeColors.push_back( col );
   }
   //typeColors.push_back( Color(1,1,1,1) ); // temporary hack to test ideas
   //typeColors.push_back( Color(0.3412, 0.4745, 0.5843, 1) );
   //typeColors.push_back( Color(0.5294, 0.3882, 0.3333,1)*2.0 );

   //vector<float> colorScale;
   //for( size_t c=0;c<typeColors.size();c++){ colorScale.push_back( 1.0 ); }
   //colorScale = clf.findArray( "-colorscale", colorScale, "Scale factors for the color of each particle type");
   //for( size_t c=0;c<typeColors.size();c++){ typeColors[c] *= colorScale[c]; }


   vector<float> opacity;
   opacity = clf.findMultiple( "-opacityscale", 1.0f, "Scale factors for the opacity of each particle type");

   SparseGrid* densityGrid = new SparseGrid;
   densityGrid->setPartitionSize( sparsesize );

   string densityFileName = clf.find("-readdensity", "" );
   if( densityFileName != "" )
   {
      ifstream input( densityFileName.c_str() );
      ReadVolumeGrid( *densityGrid, input );
      input.close();
      cout << "Density read from file " << densityFileName << endl;
   }
   else
   {
      densityGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
   }

   SparseColorGrid* colorGrid = new SparseColorGrid();
   colorGrid->setPartitionSize( sparsesize );

   string colorFileName = clf.find("-readcolor", "" );
   if( colorFileName != "" )
   {
      ifstream input( colorFileName.c_str() );
      ReadVolumeGrid( *colorGrid, input );
      input.close();
      cout << "Color read from file " << colorFileName << endl;
   }
   else
   {
      colorGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
   }

   if( densityFileName == "" || colorFileName == "" )
   {
      StampPointWisps( *densityGrid, *colorGrid, files, typeColors, opacity, blurScale );
      int blur = clf.find("-blurgrid", 0 );
      if( blur > 0 )
      {
             for( int b=0;b<blur;b++ ){Blur( *densityGrid ); Blur( *colorGrid ); }
      }
   }


   densityFileName = clf.find("-writedensity", "" );
   if( densityFileName != "" )
   {
      ofstream output( densityFileName.c_str() );
      WriteVolumeGrid( *densityGrid, output );
      output.close();
      cout << "Density written to file " << densityFileName << endl;
   }

   colorFileName = clf.find("-writecolor", "" );
   if( colorFileName != "" )
   {
      ofstream output( colorFileName.c_str() );
      WriteVolumeGrid( *colorGrid, output );
      output.close();
      cout << "Color written to file " << densityFileName << endl;
   }

   if( clf.findFlag( "-usegridlines" ) )
   {
      int gridSpacing = clf.find("-gridlinespacing", 100 );
      Vector gridCol = clf.find("-gridlinecolor", Vector(1,1,1) );
      Color gridColor = Color( gridCol[0], gridCol[1], gridCol[2], 1.0 );
      StampGridPattern( *densityGrid, *colorGrid, gridSpacing, gridColor );
   }

   data.densityField = new GriddedVolume( densityGrid );
   data.ambientColorField = new GriddedColor( colorGrid );
   data.colorField = new ConstantColor( Color(0,0,0,1) );

   float thresh = clf.find( "-clamp", -1.0f, "Value to clamp density at. Negative value suppresses clamp." );
   if( thresh > 0.0 ){ data.densityField = new ClampVolume( data.densityField, 0.0, thresh ); }




   if( clf.findFlag("-useboundingbox" ) )
   {
      cout << "Getting bounding boxes\n";
      //long bbs = getBoundingBoxes( *densityGrid, data.boundingBoxes );
      data.sparseGrid = densityGrid;
      //cout << "Found " << bbs << " bounding boxes" << endl;
   }

   data.nbDensitySamples = clf.find("-nbdensitysamples", 1 );

}





void SetUpStampSparsePointWisps( CmdLineFind& clf, RenderData& data )
{
   vector<int> NNNN;
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   vector<int> NNN = clf.findArray( "-stampNxNyNz", NNNN, "Number of grid points in stamped volume" );

   Vector lX = clf.find( "-stampsize", Vector(10,10,10), "Volume size of stamped volume");
   Vector X0 = clf.find( "-stamporigin", Vector(-5,-5,-5), "Origin of stamped volume");

   SparseGrid* densityGrid = new SparseGrid;
   densityGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );

   ParticleGroupA particles;
   int nbParticles = clf.find("-nbparticles", 1, "Number of particles to stamp");
   float pscale = clf.find("-pscale", 3.0f, "Max particle pscale");
   int nbwisps = clf.find("-nbwisps", 1000000, "Number of wisps per particle");
   float dispscale = clf.find("-wispdisplacement", 1.0f, "Wisp displacement scale factor" );
   float correlation = clf.find("-wispcorrelation", 0.0f, "Wisp walk correlation" );
   float wispoctaves = clf.find("-wispoctaves", 1.0f, "Number of octaves for wisps");
   float wisproughness = clf.find("-wisproughness", 0.5f, "Number of octaves for wisps");
   float wispfreq = clf.find("-wispfreq", 1.0f, "Number of octaves for wisps");

   float octaves = clf.find("-octaves", 3.0f,"Octaves for wisp displacement");
   float roughness = clf.find("-roughness", 0.5f, "Roughness for wisp displacement");
   float freq = clf.find("-freq", 1.5f, "Freq for wisp displacement");
   float opacity = clf.find("-opacity", 1.0f );
   float shutter = clf.find("-shutter", 0.0f );

   Vector vel = clf.find("-velocity", Vector(0,0,0), "Point wisp velocity");
   Vector accel = clf.find("-accel", Vector(0,0,0), "Point wisp acceleration");

   if( nbParticles == 1 )
   {
      Particle p;
      p.pscale() = pscale;
      p.P() = Vector( 0,0,0 );
      p.Cd() = Color( 1.0, 1.0, 1.0, 1.0 );
      p.id() = 1;
      p.accel() = accel;
      p.v() = vel;
      p.octaves() = octaves;
      p.roughness() = roughness;
      p.freq() = freq;
      p.opacity() = opacity;
      p.nbWisps() = nbwisps;
      p.wispDisplacementScale() = dispscale;
      p.wispCorrelation() = correlation;
      p.wispOctaves() = wispoctaves;
      p.wispRoughness() = wisproughness;
      p.wispFreq() = wispfreq;
      p.shutter() = shutter;
      particles.push_back( p );
   }
   else
   {




      Noise_t prnparms;
      prnparms.gaussianstandarddeviation = clf.find("-noisespread", 1.5f );
      prnparms.lognormalmean = 1.0;
      prnparms.seed = 48575;
      LognormalPRN lnprn;
      lnprn.setParameters( prnparms );

      float time = clf.find("-time", 1.0f, "Time for simple missile trail" );
      float timeScale = clf.find("-timescale", 1.0f, "Time scaling for perlin noise input" );
      float translatescale = clf.find("-transscale", 1.0f );
      Vector V0( -8.3, 10.0, 0);
      Vector accel( 0,-9.8,0);
      Vector P0 (5,-5,0);

      BallisticString( P0, V0, accel, 0.0, time, 1.2, nbParticles, particles );

      for( size_t ip=0;ip<particles.size();ip++ )
      {
         Particle& p = particles[ip];
	 const Vector& V = p.v();
	 float f = p.age()/p.lifetime();
         float af = 3.0 * f;
	 if( f > 1.0 ){ f = 1.0; }
	 if( af > 1.0 ){ af = 1.0; }
	 float gscale = lnprn.eval();
	 p.pscale() = pscale * ( 0.2 + f)/1.2; //(gscale + 0.2) * pscale;
         p.octaves() = octaves;
         p.roughness() = roughness;
	 p.opacity() = opacity * tanh(f*70.0);
	 Vector u = V.unitvector();
         p.translate() = -u * translatescale * time + ((p.P()-P0)*u)*u ; 
         p.Cd() = Color( 1.0, 1.0, 1.0, 1.0 );
	 p.age() *= timeScale;

         p.nbWisps() = nbwisps;
         p.wispDisplacementScale() = dispscale * p.pscale();
         p.wispCorrelation() = correlation;
         p.freq() = freq;
         p.wispOctaves() = wispoctaves;
         p.wispRoughness() = wisproughness;
         p.wispFreq() = wispfreq;
         p.shutter() = shutter;
	 p.normal() = V.unitvector();
	 p.right() = Vector(0.0,1.0,0.0);
	 p.right() -= (p.right() * p.normal() ) * p.normal();
	 p.right().normalize();
	 p.up() = p.normal() ^ p.right();
      }
   }

   StampPointWisps( *densityGrid, particles );

   Volume<float>* dfield = new MaskVolume( new GriddedVolume( densityGrid ) );
   Volume<Color>* whitefield = new VolumeMultiplyColor( new ConstantColor( Color(1,1,1,1) ), dfield );
   data.densityField = new GriddedVolume( densityGrid );
   data.colorField = whitefield;

}







void SetUpFlameWisps( CmdLineFind& clf, RenderData& data )
{
   vector<int> NNNN;
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   vector<int> NNN = clf.findArray( "-stampNxNyNz", NNNN, "Number of grid points in stamped volume" );

   Vector lX = clf.find( "-stampsize", Vector(10,10,10), "Volume size of stamped volume");
   Vector X0 = clf.find( "-stamporigin", Vector(-5,-5,-5), "Origin of stamped volume");

   int sparsesize = clf.find("-sparsesize", 16 );


   SparseGrid* densityGrid = new SparseGrid( sparsesize );
   densityGrid->setPartitionSize( sparsesize );

   SparseColorGrid* colorGrid = new SparseColorGrid();
   colorGrid->setPartitionSize( sparsesize );


   string densityFileName = clf.find("-readdensity", "" );
   if( densityFileName != "" )
   {
      ifstream input( densityFileName.c_str() );
      ReadVolumeGrid( *densityGrid, input );
      input.close();
      cout << "Density read from file " << densityFileName << endl;
   }
   else
   {
      densityGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
   }

   string colorFileName = clf.find("-readcolor", "" );
   if( colorFileName != "" )
   {
      ifstream input( colorFileName.c_str() );
      ReadVolumeGrid( *colorGrid, input );
      input.close();
      cout << "Color read from file " << colorFileName << endl;
   }
   else
   {
      colorGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
   }


  
   Noise_t rparms, gparms, bparms;


   rparms.octaves = clf.find("-flameoctaves", 2.0f ); 
   rparms.roughness = clf.find("-flameroughness", 0.5f ); 
   rparms.wavelength= clf.find("-flamefreq", 2.0f ); 
   rparms.translate = clf.find("-flameshift", Vector(0.4332,83.444,-0.12847) ); 
   float flamegamma = clf.find("-flamegamma", 1.0f );
   float flamebrightness = clf.find("-flamebrightness", 1.0f );

   gparms = rparms;
   gparms.translate = Vector( rparms.translate[2], -rparms.translate[0],  rparms.translate[0]  );

   bparms = rparms;
   bparms.translate = -rparms.translate;


   Noise * nr = new FractalSum<PerlinNoiseGustavson>();
   nr->setParameters(rparms);
   Noise * ng = new FractalSum<PerlinNoiseGustavson>();
   ng->setParameters(gparms);
   Noise * nb = new FractalSum<PerlinNoiseGustavson>();
   nb->setParameters(bparms);

   Volume<Color>* baseColor = new FloatMultiplyColor( new FloatGammaColor( new ComponentColor( new AbsoluteVolume(new NoiseVolume( nr )), 
                                                                                               new AbsoluteVolume(new NoiseVolume( ng )), 
						                                               new AbsoluteVolume(new NoiseVolume( nb )) ), flamegamma ), flamebrightness ); 






   ParticleGroupA particles;
   int nbParticles = clf.find("-nbparticles", 1, "Number of particles to stamp");
   float pscale = clf.find("-pscale", 3.0f, "Max particle pscale");
   int nbwisps = clf.find("-nbwisps", 1000000, "Number of wisps per particle");
   float dispscale = clf.find("-wispdisplacement", 1.0f, "Wisp displacement scale factor" );
   float correlation = clf.find("-wispcorrelation", 0.0f, "Wisp walk correlation" );
   float wispoctaves = clf.find("-wispoctaves", 1.0f, "Number of octaves for wisps");
   float wisproughness = clf.find("-wisproughness", 0.5f, "Number of octaves for wisps");
   float wispfreq = clf.find("-wispfreq", 1.0f, "Number of octaves for wisps");

   float octaves = clf.find("-octaves", 3.0f,"Octaves for wisp displacement");
   float roughness = clf.find("-roughness", 0.5f, "Roughness for wisp displacement");
   float freq = clf.find("-freq", 1.5f, "Freq for wisp displacement");
   float opacity = clf.find("-opacity", 1.0f );
   float shutter = clf.find("-shutter", 0.5f );
   Vector translate = clf.find("-translate", Vector(0,0,0) );

   Vector vel = clf.find("-velocity", Vector(0,0,0), "Point wisp velocity");
   Vector accel = clf.find("-accel", Vector(0,0,0), "Point wisp acceleration");

   if( nbParticles == 1 )
   {
      Particle p;
      p.pscale() = pscale;
      p.P() = Vector( 0,0,0 );
      p.Cd() = Color( 1.0, 1.0, 1.0, 1.0 );
      p.id() = 1;
      p.accel() = accel;
      p.v() = vel;
      p.octaves() = octaves;
      p.roughness() = roughness;
      p.freq() = freq;
      p.opacity() = opacity;
      p.translate() = translate;
      p.nbWisps() = nbwisps;
      p.wispDisplacementScale() = dispscale;
      p.wispCorrelation() = correlation;
      p.wispOctaves() = wispoctaves;
      p.wispRoughness() = wisproughness;
      p.wispFreq() = wispfreq;
      p.shutter() = shutter;
      particles.push_back( p );
   }
   else
   {

      Noise_t prnparms;
      prnparms.gaussianstandarddeviation = clf.find("-noisespread", 1.5f );
      prnparms.lognormalmean = 1.0;
      prnparms.seed = 48575;
      LognormalPRN lnprn;
      lnprn.setParameters( prnparms );

      float time = clf.find("-time", 1.0f, "Time for simple missile trail" );
      float timeScale = clf.find("-timescale", 1.0f, "Time scaling for perlin noise input" );
      float translatescale = clf.find("-transscale", 1.0f );
      Vector V0( -8.3, 10.0, 0);
      Vector accel( 0,-9.8,0);
      Vector P0 (5,-5,0);

      BallisticString( P0, V0, accel, 0.0, time, 1.2, nbParticles, particles );

      for( size_t ip=0;ip<particles.size();ip++ )
      {
         Particle& p = particles[ip];
	 const Vector& V = p.v();
	 float f = p.age()/p.lifetime();
         float af = 3.0 * f;
	 if( f > 1.0 ){ f = 1.0; }
	 if( af > 1.0 ){ af = 1.0; }
	 float gscale = lnprn.eval();
	 p.pscale() = pscale * ( 0.2 + f)/1.2; //(gscale + 0.2) * pscale;
         p.octaves() = octaves;
         p.roughness() = roughness;
	 p.opacity() = opacity * tanh(f*70.0);
	 Vector u = V.unitvector();
         p.translate() = -u * translatescale * time + ((p.P()-P0)*u)*u ; 
         p.Cd() = Color( 1.0, 1.0, 1.0, 1.0 );
	 p.age() *= timeScale;

         p.nbWisps() = nbwisps;
         p.wispDisplacementScale() = dispscale * p.pscale();
         p.wispCorrelation() = correlation;
         p.freq() = freq;
         p.wispOctaves() = wispoctaves;
         p.wispRoughness() = wisproughness;
         p.wispFreq() = wispfreq;
         p.shutter() = shutter;
	 p.normal() = V.unitvector();
	 p.right() = Vector(0.0,1.0,0.0);
	 p.right() -= (p.right() * p.normal() ) * p.normal();
	 p.right().normalize();
	 p.up() = p.normal() ^ p.right();
      }
   }

   if( densityFileName == "" || colorFileName == "" )
   {
      StampFlameWisps( *densityGrid, *colorGrid, particles, *baseColor );
   }




   densityFileName = clf.find("-writedensity", "" );
   if( densityFileName != "" )
   {
      cout << "Writing density to file " << densityFileName << endl;
      ofstream output( densityFileName.c_str() );
      WriteVolumeGrid( *densityGrid, output );
      output.close();
      cout << "Density written to file " << densityFileName << endl;
   }

   colorFileName = clf.find("-writecolor", "" );
   if( colorFileName != "" )
   {
      ofstream output( colorFileName.c_str() );
      WriteVolumeGrid( *colorGrid, output );
      output.close();
      cout << "Color written to file " << densityFileName << endl;
   }






   data.densityField = new ConstantVolume( 0.0 );
   data.ambientDensityField = new GriddedVolume( densityGrid );
   data.ambientColorField = new GriddedColor( colorGrid );
   data.colorField = new ConstantColor( Color(0,0,0,1) );

}










void SetUpStampPointWisps( CmdLineFind& clf, RenderData& data )
{
   vector<int> NNNN;
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   vector<int> NNN = clf.findArray( "-stampNxNyNz", NNNN, "Number of grid points in stamped volume" );

   Vector lX = clf.find( "-stampsize", Vector(10,10,10), "Volume size of stamped volume");
   Vector X0 = clf.find( "-stamporigin", Vector(-5,-5,-5), "Origin of stamped volume");

   VolumeGrid<float>* densityGrid = new VolumeGrid<float>;
   densityGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );

   ParticleGroupA particles;
   int nbParticles = clf.find("-nbparticles", 1, "Number of particles to stamp");
   float pscale = clf.find("-pscale", 3.0f, "Max particle pscale");
   int nbwisps = clf.find("-nbwisps", 1000000, "Number of wisps per particle");
   float dispscale = clf.find("-wispdisplacement", 1.0f, "Wisp displacement scale factor" );
   float correlation = clf.find("-wispcorrelation", 0.0f, "Wisp walk correlation" );
   float wispoctaves = clf.find("-wispoctaves", 1.0f, "Number of octaves for wisps");
   float wisproughness = clf.find("-wisproughness", 0.5f, "Number of octaves for wisps");
   float wispfreq = clf.find("-wispfreq", 1.0f, "Number of octaves for wisps");

   float octaves = clf.find("-octaves", 3.0f,"Octaves for wisp displacement");
   float roughness = clf.find("-roughness", 0.5f, "Roughness for wisp displacement");;
   float freq = clf.find("-freq", 1.5f, "Freq for wisp displacement");
   float opacity = clf.find("-opacity", 1.0f );
   float shutter = clf.find("-shutter", 0.0f );

   Vector vel = clf.find("-velocity", Vector(0,0,0), "Point wisp velocity");
   Vector accel = clf.find("-accel", Vector(0,0,0), "Point wisp acceleration");

   if( nbParticles == 1 )
   {
      Particle p;
      p.pscale() = pscale;
      p.P() = Vector( 0,0,0 );
      p.Cd() = Color( 1.0, 1.0, 1.0, 1.0 );
      p.id() = 1;
      p.accel() = accel;
      p.v() = vel;
      p.octaves() = octaves;
      p.roughness() = roughness;
      p.freq() = freq;
      p.opacity() = opacity;
      p.nbWisps() = nbwisps;
      p.wispDisplacementScale() = dispscale;
      p.wispCorrelation() = correlation;
      p.wispOctaves() = wispoctaves;
      p.wispRoughness() = wisproughness;
      p.wispFreq() = wispfreq;
      p.shutter() = shutter;
      particles.push_back( p );
   }
   else
   {




      Noise_t prnparms;
      prnparms.gaussianstandarddeviation = clf.find("-noisespread", 1.5f );
      prnparms.lognormalmean = 1.0;
      prnparms.seed = 48575;
      LognormalPRN lnprn;
      lnprn.setParameters( prnparms );

      float time = clf.find("-time", 1.0f, "Time for simple missile trail" );
      float timeScale = clf.find("-timescale", 1.0f, "Time scaling for perlin noise input" );
      float translatescale = clf.find("-transscale", 1.0f );
      Vector V0( -8.3, 10.0, 0);
      Vector accel( 0,-9.8,0);
      Vector P0 (5,-5,0);

      BallisticString( P0, V0, accel, 0.0, time, 1.2, nbParticles, particles );

      for( size_t ip=0;ip<particles.size();ip++ )
      {
         Particle& p = particles[ip];
	 const Vector& V = p.v();
	 float f = p.age()/p.lifetime();
         float af = 3.0 * f;
	 if( f > 1.0 ){ f = 1.0; }
	 if( af > 1.0 ){ af = 1.0; }
	 float gscale = lnprn.eval();
	 p.pscale() = pscale * ( 0.2 + f)/1.2; //(gscale + 0.2) * pscale;
         p.octaves() = octaves;
         p.roughness() = roughness;
	 p.opacity() = opacity * tanh(f*70.0);
	 Vector u = V.unitvector();
         p.translate() = -u * translatescale * time + ((p.P()-P0)*u)*u ; 
         p.Cd() = Color( 1.0, 1.0, 1.0, 1.0 );
	 p.age() *= timeScale;

         p.nbWisps() = nbwisps;
         p.wispDisplacementScale() = dispscale * p.pscale();
         p.wispCorrelation() = correlation;
         p.freq() = freq;
         p.wispOctaves() = wispoctaves;
         p.wispRoughness() = wisproughness;
         p.wispFreq() = wispfreq;
         p.shutter() = shutter;
	 p.normal() = V.unitvector();
	 p.right() = Vector(0.0,1.0,0.0);
	 p.right() -= (p.right() * p.normal() ) * p.normal();
	 p.right().normalize();
	 p.up() = p.normal() ^ p.right();
      }




/*
      float time = clf.find("-time", 1.0f, "Time for simple missile trail" );
      Vector V0( -8.3, 10.0, 0);
      Vector accel( 0,-9.8,0);
      Vector P0 (5,-5,0);
      float T = 1.0;
      float dT = T/nbParticles;
      UniformPRN prn;
      Noise_t prnparms;
      prnparms.gaussianstandarddeviation = clf.find("-noisespread", 1.5f );
      prnparms.lognormalmean = 1.0;
      LognormalPRN lnprn;
      lnprn.setParameters( prnparms );

      float t = 0;
      int ii = 0;
      while( t <=time )
      {
         Particle p;
	 float jitteredtime = t + dT*(prn.eval()-0.5);
	 p.P() = P0 + jitteredtime*V0 + 0.5*jitteredtime*jitteredtime*accel;
	 Vector V = V0 + jitteredtime*accel;
	 float f = (time-t)/T;
	 if( f > 1.0 ){f = 1.0; }

	 float gscale = lnprn.eval();
	 //p.pscale() = gscale * pscale;
         p.octaves() = (2.5 + 1.0*f); // *gscale;
	 p.pscale() = (0.1 + sqrt( f ))*pscale / 1.1;
         p.roughness() = 0.7 - 0.4*f;
         p.freq() = 1.4;

         p.translate() = p.P() + Vector( 0, time, time )* (V.magnitude())*0.1;
         //p.translate() = p.P()*2.0 + V * time;
         p.nbWisps() = nbwisps * (f + 0.01)/1.01;
         p.wispDisplacementScale() = sqrt(f) * dispscale * (1.0 + prn.eval());
         p.wispCorrelation() = correlation;
         p.freq() = freq;
         p.opacity() = opacity * tanh( 100.0*f );
         p.wispCorrelation() = correlation;
         p.wispOctaves() = wispoctaves;
         p.wispRoughness() = wisproughness*(1.0 + 0.5*f);
         p.wispFreq() = wispfreq;
         p.shutter() = shutter;
         p.Cd() = Color( 1.0, 1.0, 1.0, 1.0 );
         p.id() = ii;

	 p.normal() = V.unitvector();
	 p.right() = Vector(prn.eval(), prn.eval(), prn.eval() );
	 p.right() -= (p.right() * p.normal() ) * p.normal();
	 p.right().normalize();
	 p.up() = p.normal() ^ p.right();

         Vector blur = Vector( 0, dT, dT )* (V.magnitude())*0.1;
	 p.v() = blur[0] * p.normal() + blur[1] * p.right() + blur[2] * p.up();


         particles.push_back( p );
         t += dT;
	 ++ii;
      }
*/
   }



   string densityReadName = clf.find("-readdensity", "" );
   if( densityReadName != "" )
   {
      ifstream output( densityReadName.c_str() );
      ReadVolumeGrid( *densityGrid, output );
      output.close();
      cout << "Density file " << densityReadName << " read\n";
   }
   else
   {
      StampPointWisps( *densityGrid, particles );
      //StampPointWisps( *densityGrid, *colorGrid, particles );
      string densityWriteName = clf.find("-writedensity", "" );
      if( densityWriteName != "" )
      {
         ofstream output( densityWriteName.c_str() );
         WriteVolumeGrid( *densityGrid, output );
         output.close();
         cout << "Density file " << densityWriteName << " written\n";
      }
   }



   Volume<float>* dfield = new MaskVolume( new GriddedVolume( densityGrid ) );
   Volume<Color>* whitefield = new VolumeMultiplyColor( new ConstantColor( Color(1,1,1,1) ), dfield );
   data.densityField = new GriddedVolume( densityGrid );
   data.colorField = whitefield;

}






void SetUpStampNoise( CmdLineFind& clf, RenderData& data )
{
   vector<int> NNNN;
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   vector<int> NNN = clf.findArray( "-stampNxNyNz", NNNN, "Number of grid points in stamped volume" );

   Vector lX = clf.find( "-stampsize", Vector(10,10,10), "Volume size of stamped volume");
   Vector X0 = clf.find( "-stamporigin", Vector(-5,-5,-5), "Origin of stamped volume");

   VolumeGrid<float>* densityGrid = new VolumeGrid<float>;
   VolumeGrid<Color>* colorGrid = new VolumeGrid<Color>;

   densityGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
   colorGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );

   ParticleGroupA particles;
   int nbParticles = clf.find("-nbparticles", 1, "Number of particles to stamp");
   float pscale = clf.find("-pscale", 3.0f, "Max particle pscale");
   float octaves = clf.find("-octaves", 3.0f);
   float roughness = clf.find("-roughness", 0.5f);
   float freq = clf.find("-freq", 1.5f);
   float offset = clf.find("-offset", 0.0f, "offset for biasing noise");
   float fade = clf.find("-fade", 0.0f, "offset for biasing noise");
   float opacity = clf.find("-opacity", 1.0f, "opacity for noise");
   float shutter = clf.find("-shutter", 0.5f, "opacity for noise");

   Vector vel = clf.find("-velocity", Vector(0,0,0), "Point wisp velocity");
   Vector accel = clf.find("-accel", Vector(0,0,0), "Point wisp acceleration");

   float time = clf.find("-time", 1.0f, "Time for simple missile trail" );
   if( nbParticles == 1 )
   {
      Particle p;
      p.pscale() = pscale;
      p.P() = Vector( 0,0,0 );
      p.Cd() = Color( 1.0, 1.0, 1.0, 1.0 );
      p.id() = 1;
      p.octaves() = octaves;
      p.roughness() = roughness;
      p.freq() = freq;
      p.fade() = fade;
      p.offset() = offset;
      p.opacity() = opacity;
      p.shutter() = shutter;
      p.v() = vel;
      p.accel() = accel;
      particles.push_back( p );
   }
   else
   {

      Vector V0( -10.0, 10.0, 0);
      Vector accel( 0,-9.8,0);
      Vector P0 (5,-5,0);
      float T = 1.0;
      float dT = T/nbParticles;
      UniformPRN prn;
      LognormalPRN lnprn;

      Noise_t prnparms;
      prnparms.gaussianstandarddeviation = clf.find("-noisespread", 1.5f );
      prnparms.lognormalmean = 1.0;
      lnprn.setParameters( prnparms );

      float t = 0;
      int ii = 0;
      while( t <=time )
      {
         Particle p;
	 float jitteredtime = t + dT*(prn.eval()-0.5);
	 p.P() = P0 + jitteredtime*V0 + 0.5*jitteredtime*jitteredtime*accel;
	 Vector V = V0 + jitteredtime*accel;
	 float f = (time-t)/T;
	 if( f > 1.0 ){f = 1.0; }


         p.octaves() = 2.0 + 1.5*f;
	// p.pscale() = (0.1 + sqrt( f ))*pscale / 1.1;
	 p.pscale() = lnprn.eval() * pscale;
         p.roughness() = 0.7 - 0.4*f;
         p.freq() = freq;
         p.fade() = fade;
         //p.translate() = Vector( prn.eval(), prn.eval(), prn.eval() ) + V * time;
	 Vector u = V.unitvector();
	 p.v() = V;
	 p.accel() = accel;
         //p.translate() = (p.P()*u)*u + (((p.P() - (u*p.P())*u)).unitvector())*V.magnitude()*time ;
         p.translate() = Vector( prn.eval(), prn.eval(), prn.eval() ) + (((p.P() - (u*p.P())*u)).unitvector())*V.magnitude()*time*0.1 ;
         p.Cd() = Color( 1.0, 1.0, 1.0, 1.0 );
         p.id() = ii;
         p.offset() = offset;
         p.opacity() = opacity;


	 p.v() =  (((p.P() - (u*p.P())*u)).unitvector())*V.magnitude()*dT*0.1 ;
	 p.shutter() = shutter;

         particles.push_back( p );
         t += dT;
	 ++ii;
      }
   }


   FractalSum<PerlinNoiseGustavson> noise;


   StampNoiseAndColor( *densityGrid, *colorGrid, &noise, particles );

   string densityWriteName = clf.find("-writedensity", "" );
   if( densityWriteName != "" )
   {
      ofstream output( densityWriteName.c_str() );
      WriteVolumeGrid( *densityGrid, output );
      output.close();
      cout << "Density file " << densityWriteName << " written\n";
   }

   string colorWriteName = clf.find("-writecolor", "" );
   if( colorWriteName != "" )
   {
      ofstream output( colorWriteName.c_str() );
      WriteVolumeGrid( *colorGrid, output );
      output.close();
      cout << "Color file " << colorWriteName << " written\n";
   }




   data.densityField = new GriddedVolume( densityGrid );
   data.colorField = new GriddedColor( colorGrid );

}




void SetUpStampNoiseLine( CmdLineFind& clf, RenderData& data )
{
   vector<int> NNNN;
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   vector<int> NNN = clf.findArray( "-stampNxNyNz", NNNN, "Number of grid points in stamped volume" );

   Vector lX = clf.find( "-stampsize", Vector(10,10,10), "Volume size of stamped volume");
   Vector X0 = clf.find( "-stamporigin", Vector(-5,-5,-5), "Origin of stamped volume");

   VolumeGrid<float>* densityGrid = new VolumeGrid<float>;
   VolumeGrid<Color>* colorGrid = new VolumeGrid<Color>;

   densityGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
   colorGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );

   ParticleGroupA particles;
   int nbParticles = clf.find("-nbparticles", 100, "Number of particles to stamp");
   float pscale = clf.find("-pscale", 3.0f, "Max particle pscale");
   float octaves = clf.find("-octaves", 3.0f);
   float roughness = clf.find("-roughness", 0.5f);
   float freq = clf.find("-freq", 1.5f);
   float offset = clf.find("-offset", 0.0f, "offset for biasing noise");
   float fade = clf.find("-fade", 0.0f, "offset for biasing noise");
   float opacity = clf.find("-opacity", 1.0f, "opacity for noise");
   float shutter = clf.find("-shutter", 0.5f, "opacity for noise");

   Vector V0 = clf.find("-velocity", Vector(-8.3,10,0), "Point wisp velocity");
   Vector accel = clf.find("-accel", Vector(0,-9.8,0), "Point wisp acceleration");

   float time = clf.find("-time", 1.0f, "Time for simple missile trail" );
   float translatescale = clf.find("-translatescale", 1.0f );
   float correlation = clf.find("-correlation", 0.9f );





      Noise_t prnparms;
      prnparms.gaussianstandarddeviation = clf.find("-noisespread", 1.5f );
      prnparms.lognormalmean = 1.0;
      prnparms.seed = 84747;
      UniformPRN prn;
      prn.setParameters( prnparms );

      float timeScale = clf.find("-timescale", 1.0f, "Time scaling for perlin noise input" );
      Vector P0 (5,-5,0);

      BallisticString( P0, V0, accel, 0.0, 1.2, 1.2, nbParticles, particles );

      for( size_t ip=0;ip<particles.size();ip++ )
      {
	 Particle& p = particles[ip];
	 const Vector& V = p.v();
	 p.pscale() = pscale;
         p.octaves() = octaves;
         p.roughness() = roughness;
	 p.opacity() = opacity;
	 Vector u = V.unitvector();
         p.translate() = -u * translatescale * time + ((p.P()-P0)*u)*u ; 
         p.Cd() = Color( 1.0, 1.0, 1.0, 1.0 );
	 p.age() *= timeScale;
         p.freq() = freq;
         p.shutter() = shutter;
	 p.Cd() = Color( prn.eval(), prn.eval(), prn.eval(), 1.0 );
	 if( ip > 0 )
	 {
	    p.Cd() = particles[ip-1].Cd() * correlation + p.Cd() * (1.0 - correlation);
	 }
      }



   FractalSum<PerlinNoiseGustavson> noise;


   StampNoiseAndColor( *densityGrid, *colorGrid, &noise, particles );

   string densityWriteName = clf.find("-writedensity", "" );
   if( densityWriteName != "" )
   {
      ofstream output( densityWriteName.c_str() );
      WriteVolumeGrid( *densityGrid, output );
      output.close();
      cout << "Density file " << densityWriteName << " written\n";
   }

   string colorWriteName = clf.find("-writecolor", "" );
   if( colorWriteName != "" )
   {
      ofstream output( colorWriteName.c_str() );
      WriteVolumeGrid( *colorGrid, output );
      output.close();
      cout << "Color file " << colorWriteName << " written\n";
   }




   data.densityField = new GriddedVolume( densityGrid );
   data.colorField = new GriddedColor( colorGrid );

}





void SetUpStampObjParticles( CmdLineFind& clf, RenderData& data )
{
    string filename = clf.find( "-objfilename", "", "Name of obj file with triangular geometry" );
    ObjParser p;
    if( !p.ParseFile( filename ) )
    {
       cout << "Could not read file " << filename << endl << flush;
       return;
    }

    TriangleGeometry g;
    float scaling = clf.find("-objscale", 1.0f, "Scale vertices of obj file" );
    g.setScaling( scaling );
    if( !p.Fill( g ) )
    {
       cout << "Could not read geometry from file " << filename << endl << flush;
       return;
    }

    Vector llc = g.LLC();
    Vector urc = g.URC();
    cout << "Obj faces: " << g.nbFaces() << endl;
    cout << "Obj vertices: " << g.nbVertices() << endl;
    cout << "Obj BB:\n";
    cout << llc[0] << " " <<  llc[1] << " " <<  llc[2] << "  X  " <<  urc[0] << " " <<  urc[1] << " " <<  urc[2] << endl << flush;

    float octaves = clf.find( "-octaves", 2.5f, "Octaves for Obj particles");
    float freq = clf.find("-freq", 1.5f);
    float rough = clf.find( "-roughness", 0.65f, "Roughness for Obj particles");
    float offset = clf.find("-offset", 0.0f, "offset for biasing noise");
    float fade = clf.find("-fade", 0.0f, "offset for biasing noise");
    float opacity = clf.find("-opacity", 1.0f, "particle opacity");

    ParticleGroupA particles;
    cout << "Converting vertices to particles\n";
    Geometry2Particles( g, particles );
    Noise_t parms;
    UniformPRN noise;
    for( size_t i=0;i<particles.size();i++ )
    {
       parms.seed = particles[i].id();
       noise.setParameters( parms);
       particles[i].Cd() = Color( 1,1,1, 1 );
       //particles[i].Cd() = Color( noise.eval(), noise.eval(), noise.eval(), 1 );
       particles[i].octaves() = octaves;
       particles[i].roughness() = rough;
       particles[i].freq() = freq;
       particles[i].offset() = offset;
       particles[i].fade() = fade;
       particles[i].opacity() = opacity;
    }

   vector<int> NNNN;
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   vector<int> NNN = clf.findArray( "-stampNxNyNz", NNNN, "Number of grid points in stamped volume" );

   Vector lX = clf.find( "-stampsize", Vector(10,10,10), "Volume size of stamped volume");
   Vector X0 = clf.find( "-stamporigin", Vector(-5,-5,-5), "Origin of stamped volume");

   VolumeGrid<float>* densityGrid = new VolumeGrid<float>;
   VolumeGrid<Color>* colorGrid = new VolumeGrid<Color>;

   densityGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
   colorGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );

    FractalSum<PerlinNoise> fsp;
    StampNoiseAndColor( *densityGrid, *colorGrid, &fsp, particles );

   data.densityField = new GriddedVolume( densityGrid );
   data.colorField = new GriddedColor( colorGrid );
}






void SetUpStampPyro( CmdLineFind& clf, RenderData& data )
{
   vector<int> NNNN;
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   vector<int> NNN = clf.findArray( "-stampNxNyNz", NNNN, "Number of grid points in stamped volume" );

   Vector lX = clf.find( "-stampsize", Vector(10,10,10), "Volume size of stamped volume");
   Vector X0 = clf.find( "-stamporigin", Vector(-5,-5,-5), "Origin of stamped volume");

   VolumeGrid<float>* densityGrid = new VolumeGrid<float>;
   VolumeGrid<Color>* colorGrid = new VolumeGrid<Color>;

   densityGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
   colorGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );

   ParticleGroupA particles;
   int nbParticles = clf.find("-nbparticles", 1, "Number of particles to stamp");
   float pscale = clf.find("-pscale", 3.0f, "Max particle pscale");
   float pyroamplitude = clf.find("-pyroamplitude", 0.3f, "Pyro displacement scale factor" );
   float pyrogamma = clf.find("-pyrogamma", 1.0f, "Pyro gamma of displacements" );
   float octaves = clf.find("-octaves", 3.0f,"Octaves for wisp displacement");
   float roughness = clf.find("-roughness", 0.5f, "Roughness for wisp displacement");;
   float freq = clf.find("-freq", 1.5f, "Freq for wisp displacement");
   float pyrodensity = clf.find("-pyrodensity", 1.0f );
   float translatescale = clf.find("-transscale", 1.0f );

   if( nbParticles == 1 )
   {
      Particle p;
      p.pscale() = pscale;
      p.P() = Vector( 0,0,0 );
      p.Cd() = Color( 1.0, 1.0, 1.0, 1.0 );
      p.id() = 1;
      p.octaves() = octaves;
      p.roughness() = roughness;
      p.freq() = freq;
      p.pyroAmplitude() = pyroamplitude;
      p.pyroGamma() = pyrogamma;
      p.pyroDensity() = pyrodensity;
      particles.push_back( p );
   }
   else
   {



      Noise_t prnparms;
      prnparms.gaussianstandarddeviation = clf.find("-noisespread", 1.5f );
      prnparms.lognormalmean = 1.0;
      LognormalPRN lnprn;
      lnprn.setParameters( prnparms );

      float time = clf.find("-time", 1.0f, "Time for simple missile trail" );
      float timeScale = clf.find("-timescale", 1.0f, "Time scaling for perlin noise input" );
      Vector V0( -8.3, 10.0, 0);
      Vector accel( 0,-9.8,0);
      Vector P0 (5,-5,0);

      BallisticString( P0, V0, accel, 0.0, time, 1.2, nbParticles, particles );

      for( size_t ip=0;ip<particles.size();ip++ )
      {
         Particle& p = particles[ip];
	 const Vector& V = p.v();
	 float f = p.age()/p.lifetime();
	 if( f > 1.0 ){ f = 1.0; }
         p.pyroAmplitude() = pyroamplitude * (0.2 + f)/1.2;
	 p.pscale() = pscale *(0.2 + f)/1.2;
         p.octaves() = octaves;
         p.roughness() = roughness;
         p.freq() = freq;
	 Vector u = V.unitvector();
         p.translate() = -(p.P()*u)*u * translatescale ; 
         p.pyroDensity() = pyrodensity;
         p.Cd() = Color( 1.0, 1.0, 1.0, 1.0 );
	 p.age() *= timeScale;
      }
   }


   string densityWriteBaseName = clf.find("-writedensity", "", "Base name of files to write density into" );
   string densityReadBaseName =  clf.find("-readdensity", "", "Base name of files to read density from" );

       if( densityReadBaseName != "" )
       {
          string filen = densityReadBaseName;
	  ifstream input( filen.c_str() );
          ReadVolumeGrid( *densityGrid, input );
	  input.close();
	  cout << "Density file " << filen << " read\n";
	  Volume<Color>* whitefield = new VolumeMultiplyColor( new ConstantColor( Color(1,1,1,1) ), new MaskVolume( new GriddedVolume( densityGrid ) ) );
	  lux::Sample( colorGrid, whitefield );
       }
       else
       {
          StampPyro( *densityGrid, *colorGrid, particles );
          if( densityWriteBaseName != "" )
          {
             string filen = densityWriteBaseName;
	     ofstream output( filen.c_str() );
             WriteVolumeGrid( *densityGrid, output );
	     output.close();
	     cout << "Density file " << filen << " written\n";
          }
       }



   data.densityField = new GriddedVolume( densityGrid );
   data.colorField = new GriddedColor( colorGrid );

}








Volume<float>* SetUpPyroTorus( CmdLineFind& clf )
{
   Volume<float> * torus = SetUpObj( clf );
   FractalSum<PerlinNoise>* noise = new FractalSum<PerlinNoise>();

   Noise_t parms;
   parms.octaves = clf.find( "-octaves", 5.2f );
   parms.wavelength = clf.find("-freq", 3.0f );
   noise->setParameters( parms );
   float amplitude = clf.find("-pyroTorusAmplitude", 0.5f );
   Volume<float>* nz = new NoiseVolume( noise );
   Volume<float>* anz = new AbsoluteVolume(nz);
   Volume<float>* snz = new MultiplyVolume( anz, amplitude );
   int nbiter = clf.find("-nbiterations", 5 );
   float marchX = clf.find("-dx" , 0.5f );
   Volume<float>* displ = new ImplicitFunctionDisplacement( torus, snz, marchX, nbiter  );

   vector<int> NNNN;
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   vector<int> NNN = clf.findArray( "-levelsetNxNyNz", NNNN, "Number of grid points in levelset volume" );

   Vector lX = clf.find( "-levelsetsize", Vector(10,10,10), "Volume size of levelset volume");
   Vector X0 = clf.find( "-levelsetorigin", Vector(-5,-5,-5), "Origin of levelset volume");

   VolumeGrid<float>* levelsetGrid = new VolumeGrid<float>;
   levelsetGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );

   lux::Sample( levelsetGrid, displ );

   return new GriddedVolume( levelsetGrid );
}


void SetUpImplicitHumanoid( CmdLineFind& clf, RenderData& data )
{
   // box torso
   Vector torsoDimensions( 1,2,1);
   Vector torsoP( 0,0,0);
   float torsoPwr = clf.find( "-torsopower", 1.0f, "Power for torso box");
   Volume<float>* torso = new ScaleVolume( new CsgBoxVolume( torsoP, 1.0, torsoPwr ), torsoDimensions );

   float armlength = torsoDimensions[0]*1.5;
   float armthickness = armlength/5.0;
   Vector armP = torsoP + Vector(1,0,0)*torsoDimensions[0];
   Vector upperarmaxis(1,0,0);
   Vector lowerarmaxis(1,1,0);
   Volume<float>* leftupperarm = new EllipseVolume( armP, upperarmaxis , armlength, armthickness ); 
   Volume<float>* leftlowerarm = new EllipseVolume( armP + upperarmaxis*armlength*0.75+lowerarmaxis*armlength*0.5 , lowerarmaxis , armlength, armthickness ); 

   armP = torsoP - Vector(1,0,0)*torsoDimensions[0];
   upperarmaxis = Vector(-1,0,0);
   lowerarmaxis = Vector(-1,1,0);
   Volume<float>* rightupperarm = new EllipseVolume( armP, upperarmaxis , armlength, armthickness ); 
   Volume<float>* rightlowerarm = new EllipseVolume( armP + upperarmaxis*armlength*0.75+lowerarmaxis*armlength*0.5 , lowerarmaxis , armlength, armthickness ); 

   armP = torsoP - Vector(0,1,0)*torsoDimensions[1] - Vector(1,0,0)*torsoDimensions[0]*0.75;
   upperarmaxis = Vector(-1,-1,0);
   lowerarmaxis = Vector(-0.3,-1,0);
   Volume<float>* rightupperleg = new EllipseVolume( armP, upperarmaxis , armlength, armthickness ); 
   Volume<float>* rightlowerleg = new EllipseVolume( armP + upperarmaxis*armlength*0.5+lowerarmaxis*armlength*0.5 , lowerarmaxis , armlength, armthickness ); 

   armP = torsoP - Vector(0,1,0)*torsoDimensions[1] + Vector(1,0,0)*torsoDimensions[0]*0.75;
   upperarmaxis = Vector(1,-1,0);
   lowerarmaxis = Vector(0.3,-1,0);
   Volume<float>* leftupperleg = new EllipseVolume( armP, upperarmaxis , armlength, armthickness ); 
   Volume<float>* leftlowerleg = new EllipseVolume( armP + upperarmaxis*armlength*0.5+lowerarmaxis*armlength*0.5 , lowerarmaxis , armlength, armthickness ); 


    Vector neckP = torsoP + Vector(0,1,0)*torsoDimensions[1]*0.5;
    float neckLength = clf.find("-necklength", 2.0f, "Length of neck part");
    Vector neckAxis(0,1,0);
    float neckRadius = clf.find("-neckthickness", 0.4f, "Thickness of neck part");
    Volume<float>* neck = new ImplicitCylinder( neckP, neckAxis, neckLength, neckRadius );

    float headRadius = clf.find("-headsize", 1.5f, "size of head part");
    Vector headP = neckP + Vector(0,1,0)*(neckLength*0.25 + headRadius*0.75);
    Volume<float>* head = new SphereVolume( headP, headRadius );

    float eyeRadius = clf.find("-eyesize", 0.3f, "size of eye part");
    Vector lefteyeP = headP + (Vector(0.35,0.35,1).unitvector())*headRadius*0.75 ;
    Volume<float>* lefteye = new IntersectionVolume( new SphereVolume( lefteyeP-Vector(0,0,eyeRadius*0.25), eyeRadius ), new SphereVolume( lefteyeP+Vector(0,0,eyeRadius*0.25), eyeRadius )  );
    Vector righteyeP = headP + (Vector(-0.35,0.35,1).unitvector())*headRadius*0.75 ;
    //Volume<float>* righteye = new SphereVolume( righteyeP, eyeRadius );
    Volume<float>* righteye = new IntersectionVolume( new SphereVolume( righteyeP-Vector(0,0,eyeRadius*0.25), eyeRadius ), new SphereVolume( righteyeP+Vector(0,0,eyeRadius*0.25), eyeRadius )  );

    Vector noseP = headP + (Vector(0,0.0,1).unitvector())*headRadius + Vector(0,0,0.4);
    Vector noseAxis = Vector(0,0,-1);
    float noseLength = clf.find( "-noselength", 0.5f );
    float noseWidth  = clf.find( "-nosewidth", 0.3f );
    float tana = noseWidth/noseLength;
    float angle = atan(tana)*0.5*180.0/M_PI;
    Volume<float>* nose = new ConeVolume( noseP, noseAxis, noseLength, angle );

   float smileLength = clf.find("-smilelength", 0.75f );
   float smileThickness = clf.find("-smilethickness", 0.2f);
   Vector smileAxis = Vector( 0,-0.3,1);
   Vector smileP = headP + (smileAxis+Vector(0,0.15,0)).unitvector()*headRadius + Vector(0,0,-0.7);
   Volume<float>* smile = new TorusVolume( smileP, smileAxis, smileLength, smileThickness );
   smile = new CutoutVolume( smile, new PlaneVolume(smileP, Vector(0,1,0) )  );

   float blendstrength = clf.find( "-blendstrength", 5, "Strength of blending of implicit humanoid");
   vector<Volume<float>*> parts;
   parts.push_back( new MultiplyVolume(torso, blendstrength) );
   parts.push_back( new MultiplyVolume(leftupperarm, blendstrength) );
   parts.push_back( new MultiplyVolume(leftlowerarm, blendstrength ) );
   parts.push_back( new MultiplyVolume(rightupperarm, blendstrength) );
   parts.push_back( new MultiplyVolume(rightlowerarm, blendstrength ) );
   parts.push_back( new MultiplyVolume(rightupperleg, blendstrength) );
   parts.push_back( new MultiplyVolume(rightlowerleg, blendstrength ) );
   parts.push_back( new MultiplyVolume(leftupperleg, blendstrength) );
   parts.push_back( new MultiplyVolume(leftlowerleg, blendstrength ) );
   parts.push_back( new MultiplyVolume(neck, blendstrength*2.0 ) );
   parts.push_back( new MultiplyVolume(head, blendstrength ) );
   Volume<float>* body = new MultiBlendVolume( parts );

   float pantthickness = clf.find("-pantthickness", 11.2f );
   Volume<float>* clothes = new AddVolume( body, new ConstantVolume(pantthickness) );
   clothes = new CutoutVolume( clothes, body );
   Volume<float>* pants = new CutoutVolume( clothes,  new PlaneVolume( Vector( 0,-0.5,0), Vector(0,1,0) ) );
   pants = new CutoutVolume( pants,  new PlaneVolume( Vector( 0,-3.5,0), Vector(0,-1,0) ) );
   Volume<Color>* pantscolor = new VolumeMultiplyColor( new ConstantColor( Color(1,0,1,1) ), new MaskVolume(pants) );
   
   Volume<float>* shirt = new CutoutVolume( clothes,  new PlaneVolume( Vector( 0,-0.45,0), Vector(0,-1,0) ) );
   shirt = new CutoutVolume( shirt,  new PlaneVolume( Vector( 0,1.2,0), Vector(0,1,0) ) );
   shirt = new CutoutVolume( shirt,  new PlaneVolume( Vector( -1.4,0,0), Vector(-1,0,0) ) );
   shirt = new CutoutVolume( shirt,  new PlaneVolume( Vector( 1.4,0,0), Vector(1,0,0) ) );
   Volume<Color>* shirtcolor = new VolumeMultiplyColor( new ConstantColor( Color(0.4,1,0.4,1) ), new MaskVolume(shirt) );


   // now cut out eyes and a smile
   body = new CutoutVolume( body, righteye );
   body = new CutoutVolume( body, lefteye );
   body = new CutoutVolume( body, nose );
   body = new CutoutVolume( body, smile );


   Volume<Color>* bodycolor = new ConstantColor( Color( 1,0.6,0,1) );
   bodycolor = new VolumeMultiplyColor( bodycolor, new MaskVolume( body ) );

   float ambientlevel = clf.find("-ambientlevel", 0.1f );
   Volume<Color>* torsocolor = new FloatMultiplyColor( bodycolor, ambientlevel );

   body = new UnionVolume( body, righteye );
   body = new UnionVolume( body, lefteye );
   body = new UnionVolume(body, nose );
   body = new UnionVolume(body, smile );

   bodycolor = new AddColor( bodycolor , new VolumeMultiplyColor( new ConstantColor( Color(1,1,1,1) ), new MaskVolume(righteye) ) );
   bodycolor = new AddColor( bodycolor , new VolumeMultiplyColor( new ConstantColor( Color(1,1,1,1) ), new MaskVolume(lefteye) ) );
   bodycolor = new AddColor( bodycolor , new VolumeMultiplyColor( new ConstantColor( Color(0,0,1,1) ), new MaskVolume(nose) ) );
   bodycolor = new AddColor( bodycolor , new VolumeMultiplyColor( new ConstantColor( Color(1,0,0,1) ), new MaskVolume(smile) ) );


   float hairLength = clf.find("-hairlength", 1.5f );
   float hairThickness = clf.find("-hairthickness", 0.05f);
   Vector hairAxis = Vector( 0,0,1);
   Vector hairP = headP + (hairAxis^Vector(0,1,0)).unitvector()*hairLength;
   Volume<float>* hairstrand = new TorusVolume( hairP, hairAxis, hairLength, hairThickness );

   int nbhairs = clf.find("-nbhairs", 1 );
   UniformPRN rng;
   for( int i=1;i<nbhairs;i++ )
   {
       hairAxis = Vector( rng.eval()-0.5,0,rng.eval()-0.5).unitvector();
       hairP = headP + (hairAxis^Vector(0,1,0)).unitvector()*hairLength;
       hairstrand = new UnionVolume( hairstrand, new TorusVolume( hairP, hairAxis, hairLength, hairThickness ) );
   }
   Volume<float>* hairclipper = new PlaneVolume( hairP, Vector( 0,-1,0) );
   hairstrand = new CutoutVolume( hairstrand, hairclipper );
   hairstrand = new CutoutVolume( hairstrand, body );

   Volume<Color>* haircolor = new VolumeMultiplyColor( new ConstantColor( Color(0,1,0,1) ), new MaskVolume( hairstrand ) );
   body = new UnionVolume( body, hairstrand );
   bodycolor = new AddColor( bodycolor, haircolor );

   if( clf.findFlag("-putpantsonhim") )
   {
      cout << "\nPUTTING PANTS ON THE HUMANOID\n";
      body = new UnionVolume( body, pants );
      bodycolor = new AddColor( bodycolor, pantscolor );
      body = new UnionVolume( body, shirt );
      bodycolor = new AddColor( bodycolor, shirtcolor );
   }

   if( clf.findFlag("-crown") )
   {
      cout << "\nCROWN\n";
      Volume<float>* crown = new TranslateVolume( new ScaleVolume( new IcosahedronVolume( ), 0.3*0.9 ), Vector(0,3.75,0.2) );
      crown = new CutoutVolume( crown, body );
      Volume<Color>* crowncolor = new VolumeMultiplyColor( new ConstantColor( Color(1.0,1,0.1,1) ), new MaskVolume(crown) );

      body = new UnionVolume( body, crown );
      bodycolor = new AddColor( bodycolor, crowncolor );
   }

   // set up option for sampling in grid


   float clamp = clf.find("-humanoidclamp", 0.02f, "Clamp for implicit humanoid");
   data.densityField = new ClampVolume( new MultiplyVolume( body, 1.0/clamp ), 0.0, 1.0 );

   bool isLit = clf.findFlag("-islit" );
   if( isLit )
   {
      
      data.colorField = bodycolor;
      data.ambientColorField = torsocolor;
   }
   else
   {
      data.colorField = new ConstantColor( Color(0,0,0,1) );
      data.ambientColorField = bodycolor;
   }
}


void SetUpOlympicRings( CmdLineFind& clf, RenderData& data )
{
   float majorradius = 1.0;
   float minorradius = clf.find("-olympicringminorradius",0.10f);
   Vector axis = clf.find("-yellowring", Vector(0.3,-0.2,1) );
   Vector center( -(majorradius+minorradius)*1.1, -majorradius/2.0, 0 );
   Volume<float>* v0 = new TorusVolume( center, axis, majorradius, minorradius );
   Volume<Color>* cd0 = new VolumeMultiplyColor( new ConstantColor( Color( 1,1,0,1) ),  new MaskVolume(v0)  );

   Volume<float>* v = v0;
   Volume<Color>* Cd = cd0;

   axis = clf.find("-greenring", Vector(0,0.3,1) );
   center = Vector( (majorradius+minorradius)*1.1, -majorradius/2.0, 0 );
   Volume<float>* v1 = new TorusVolume( center, axis, majorradius, minorradius );
   Volume<Color>* cd1 = new VolumeMultiplyColor( new ConstantColor( Color( 0,1,0,1) ),  new MaskVolume(v1)  );

   v = new UnionVolume( v, v1 );
   Cd = new AddColor( Cd, cd1 );


   center = Vector( 0, majorradius/2.0, 0 );
   axis = clf.find("-blackring", Vector(0,-0.1,1) );
   Volume<float>* v2 = new TorusVolume( center, axis, majorradius, minorradius );
   Volume<Color>* cd2 = new VolumeMultiplyColor( new ConstantColor( Color( 0.5,0.5,0.5,1) ),  new MaskVolume(v2)  );

   v = new UnionVolume( v, v2 );
   Cd = new AddColor( Cd, cd2 );

   center = Vector( 0, majorradius/2.0, 0 ) + Vector( -(majorradius+minorradius)*2.2,0,0);
   axis = clf.find("-bluering", Vector(-0.4,0,1) );
   Volume<float>* v3 = new TorusVolume( center, axis, majorradius, minorradius );
   Volume<Color>* cd3 = new VolumeMultiplyColor( new ConstantColor( Color( 0,0,1,1) ),  new MaskVolume(v3)  );

   v = new UnionVolume( v, v3 );
   Cd = new AddColor( Cd, cd3 );

   center = Vector( 0, majorradius/2.0, 0 ) + Vector( (majorradius+minorradius)*2.2,0,0);
   axis = clf.find("-redring", Vector(-0.5,0.3,1) );
   Volume<float>* v4 = new TorusVolume( center, axis, majorradius, minorradius );
   Volume<Color>* cd4 = new VolumeMultiplyColor( new ConstantColor( Color( 1,0,0,1) ),  new MaskVolume(v4)  );

   v = new UnionVolume( v, v4 );
   Cd = new AddColor( Cd, cd4 );




   data.ambientColorField = Cd;
   float thresh = clf.find("-olympicringsclamp", 0.01f );
   data.densityField =  new ClampVolume( new MultiplyVolume( v, 1.0/thresh ), 0.0, 1.0);
}





Volume<float>* SetUpCutTeapot( CmdLineFind& clf )
{
    string filename = clf.find( "-objfilename", "~/projects/models/teapot.obj", "Name of obj file with triangular geometry" );
    ObjParser p;
    if( !p.ParseFile( filename ) )
    {
       cout << "Could not read file " << filename << endl << flush;
       return SetUpSphere( clf );
    }

    TriangleGeometry g;
    float scaling = clf.find("-objscale", 0.03f, "Scale vertices of obj file" );
    g.setScaling( scaling );
    if( !p.Fill( g ) )
    {
       cout << "Could not read geometry from file " << filename << endl << flush;
       return SetUpSphere( clf );
    }



    Vector dims = g.URC() - g.LLC();
    Vector center = (g.URC() + g.LLC())*0.5;
    dims *= 2.0;
    Vector llc = center - dims*0.5;
    Vector urc = llc + dims;
    int nx = 300;
    int ny = nx * ( dims[1]/dims[0] );
    int nz = nx * ( dims[2]/dims[0] );

    nx  = clf.find("-objnx", nx, "Grid poins for obj levelset");
    ny  = clf.find("-objny", ny, "Grid poins for obj levelset");
    nz  = clf.find("-objnz", nz, "Grid poins for obj levelset");
    
    cout << "Level set dimensions:\n";
    cout << nx << " " <<  ny << " " <<  nz << endl << flush;

    VolumeGrid<float>* levelset = new VolumeGrid<float>();
    levelset->init( nx, ny, nz, dims[0], dims[1], dims[2], llc );
    levelset->setClearValue( -dims.magnitude() );

    float lsthresh = clf.find("-lsthresh", 0.1f, "Threshold for filtering level set intersections");
    RayMarchLevelSet( g, *levelset, lsthresh );
    Volume<float>* v = new GriddedVolume(levelset);


    // insert knives (spheres) to cut things up

    Volume<float>* cutTeapot = v;
    vector<Vector> knives = clf.findMultiple("-knife", Vector(0,0,0), "Knives for cutting geometry");
    for( size_t i=0;i<knives.size();i++ )
    {
       Volume<float>* knife = new SphereVolume( knives[i], 0.75 );
       cutTeapot = new CutoutVolume( cutTeapot, knife );
    }

   return cutTeapot;
}



Volume<float>* SetUpCutAndHandledTeapot( CmdLineFind& clf )
{
   Volume<float>* cutTeapot = SetUpCutTeapot( clf );
   Vector handleP = clf.find("-handleP", Vector(0,0,0) );
   Vector handleD = clf.find("-handleD", Vector(0,1,0) );
   float handleLength = clf.find("-handlelength", 1.0f );

   Volume<float>* handle = new TorusVolume( handleP, handleD, handleLength, 0.2 );
   Volume<float>* teapot = new UnionVolume( cutTeapot, handle );


   Vector spoutP = clf.find("-spoutP", Vector(0,0,0) );
   Vector spoutD = clf.find("-spoutD", Vector(0,1,0) );
   float spoutlength = clf.find("-spoutlength", 1.0f );
   float spoutwidth = clf.find("-spoutwidth", 0.5f );

   Volume<float>* spout = new TranslateVolume( new ImplicitCylinder( Vector(0,0,0), spoutD, spoutlength, spoutwidth ), spoutP );
   Volume<float>* insidespout = new TranslateVolume( new ImplicitCylinder( Vector(0,0,0), spoutD, spoutlength, spoutwidth*0.5 ), spoutP );
   spout = new CutoutVolume( spout, insidespout );
   teapot = new UnionVolume( teapot, spout );


   return teapot;
}




Volume<float>* SetUpCutDisplacedTeapot( CmdLineFind& clf )
{
   Volume<float>* cutTeapot = SetUpCutAndHandledTeapot( clf );
   FractalSum<PerlinNoise>* noise = new FractalSum<PerlinNoise>();

   Noise_t parms;
   parms.octaves = clf.find( "-octaves", 5.2f );
   parms.wavelength = clf.find("-freq", 3.0f );
   parms.roughness = clf.find("-roughness", 0.5f );
   noise->setParameters( parms );
   float amplitude = clf.find("-displaceamplitude", 0.1f );
   Volume<float>* nz = new NoiseVolume( noise );
   Volume<float>* anz = new AbsoluteVolume(nz);
   Volume<float>* snz = new MultiplyVolume( anz, amplitude );

   // Use Cylinders to trim the top and bottom of the displacement

   
   float rad = 10.0;
   float length = 10.0;
   Vector center = Vector( 0,length*0.5-0.5,0 );
   Vector axis = Vector(0,1,0);
   Volume<float>* cyl1 = new ClampVolume( new MultiplyVolume( new ImplicitCylinder( center, axis, length, rad ), 10.0 ), 0.0, 1.0 ); 

   center = Vector( 0,-length*0.5+0.0,0 );
   axis = -axis;
   Volume<float>* cyl2 = new ClampVolume( new MultiplyVolume( new ImplicitCylinder( center, axis, length, rad ), 10.0 ), 0.0, 1.0 ); 

   snz = new MultiplyVolume( snz, cyl1 );
   snz = new MultiplyVolume( snz, cyl2 );


   int nbiter = clf.find("-nbiterations", 5 );
   float marchX = clf.find("-dx" , 0.5f );
   Volume<float>* displ = new ImplicitFunctionDisplacement( cutTeapot, snz, marchX, nbiter  );
   vector<int> NNNN;
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   vector<int> NNN = clf.findArray( "-levelsetNxNyNz", NNNN, "Number of grid points in levelset volume" );

   Vector lX = clf.find( "-levelsetsize", Vector(10,10,10), "Volume size of levelset volume");
   Vector X0 = clf.find( "-levelsetorigin", Vector(-5,-5,-5), "Origin of levelset volume");

   VolumeGrid<float>* levelsetGrid = new VolumeGrid<float>;
   levelsetGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );

   lux::Sample( levelsetGrid, displ );

   return new GriddedVolume( levelsetGrid );
}


Volume<float>* SetUpAdvectedSphere( CmdLineFind& clf )
{
   Volume<float>* volume = SetUpSphere( clf );

   Noise_t parms;
   parms.translate = clf.find( "-translate", Vector(0,0,0), "Noise translation");
   parms.wavelength = clf.find("-freq", 1.432234f, "Noise frequency");
   parms.amplitude = clf.find("-amp", 1.5f, "Noise amplitude" );
   parms.roughness = clf.find("-rough", 0.5f, "Noise roughness" );
   parms.octaves = clf.find("-octaves", 3.0f, "Noise octaves" );

   Noise * noise = new FractalSum<PerlinNoise>();
   noise->setParameters(parms);

   Volume<Vector>* velocity = new NoiseVectorVolume( noise );


   bool makeDivFree = clf.findFlag("-divfree" );
   if( makeDivFree )
   {
      vector<int> NNNN;
      NNNN.push_back( 100 );
      NNNN.push_back( 100 );
      NNNN.push_back( 100 );
      vector<int> NNN = clf.findArray( "-divfreeNxNyNz", NNNN, "Number of grid points in divfree velocity" );

      Vector lX = clf.find( "-divfreesize", Vector(10,10,10), "Volume size of divfree velocity");
      Vector X0 = clf.find( "-divfreeorigin", Vector(-5,-5,-5), "Origin of divfree velocity");

      VolumeGrid<Vector>* divfreevel = new VolumeGrid<Vector>;
      divfreevel->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
      GaussSeidelDivFree( *divfreevel, velocity, 1000, 0.00001 );
      velocity = new GriddedVectorVolume( divfreevel );
   }



   float dt = clf.find("-dt", 0.1f,"Advection time step");
   int nbSteps = clf.find("-nbsteps", 1, "Number of advection steps");

   float ddt = dt/nbSteps;

   for( int i=0;i<nbSteps;i++ )
   {
      volume = new AdvectVolume( volume, velocity, ddt );
   }
   return volume;
}


Volume<float>* SetUpOtherAdvectedSphere( CmdLineFind& clf )
{
   Volume<float>* volume = SetUpSphere( clf );

   Noise_t parms;
   parms.translate = clf.find( "-translate", Vector(0,0,0), "Noise translation");
   parms.wavelength = clf.find("-freq", 1.432234f, "Noise frequency");
   parms.amplitude = clf.find("-amp", 1.5f, "Noise amplitude" );
   parms.roughness = clf.find("-rough", 0.5f, "Noise roughness" );
   parms.octaves = clf.find("-octaves", 3.0f, "Noise octaves" );

   Noise * noise = new FractalSum<PerlinNoiseGustavson>();
   noise->setParameters(parms);

   Volume<Vector>* velocity = new NoiseSampleVectorVolume( noise );

   bool makeDivFree = clf.findFlag("-divfree" );
   if( makeDivFree )
   {
      vector<int> NNNN;
      NNNN.push_back( 100 );
      NNNN.push_back( 100 );
      NNNN.push_back( 100 );
      vector<int> NNN = clf.findArray( "-divfreeNxNyNz", NNNN, "Number of grid points in divfree velocity" );

      Vector lX = clf.find( "-divfreesize", Vector(10,10,10), "Volume size of divfree velocity");
      Vector X0 = clf.find( "-divfreeorigin", Vector(-5,-5,-5), "Origin of divfree velocity");

      VolumeGrid<Vector>* divfreevel = new VolumeGrid<Vector>;
      divfreevel->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
      GaussSeidelDivFree( *divfreevel, velocity, 1000, 0.00001 );
      velocity = new GriddedVectorVolume( divfreevel );
   }


   float dt = clf.find("-dt", 0.1f,"Advection time step");
   int nbSteps = clf.find("-nbsteps", 1, "Number of advection steps");

   float ddt = dt/nbSteps;

   for( int i=0;i<nbSteps;i++ )
   {
      volume = new AdvectVolume( volume, velocity, ddt );
   }
   return volume;
}


Volume<float>* SetUpBFECCAdvectedSphere( CmdLineFind& clf )
{
   Volume<float>* volume = SetUpSphere( clf );

   Noise_t parms;
   parms.translate = clf.find( "-translate", Vector(0,0,0), "Noise translation");
   parms.wavelength = clf.find("-freq", 1.432234f, "Noise frequency");
   parms.amplitude = clf.find("-amp", 1.5f, "Noise amplitude" );
   parms.roughness = clf.find("-rough", 0.5f, "Noise roughness" );
   parms.octaves = clf.find("-octaves", 3.0f, "Noise octaves" );

   Noise * noise = new FractalSum<PerlinNoiseGustavson>();
   noise->setParameters(parms);

   Volume<Vector>* velocity = new NoiseVectorVolume( noise );

   bool makeDivFree = clf.findFlag("-divfree" );
   if( makeDivFree )
   {
      vector<int> NNNN;
      NNNN.push_back( 100 );
      NNNN.push_back( 100 );
      NNNN.push_back( 100 );
      vector<int> NNN = clf.findArray( "-divfreeNxNyNz", NNNN, "Number of grid points in divfree velocity" );

      Vector lX = clf.find( "-divfreesize", Vector(10,10,10), "Volume size of divfree velocity");
      Vector X0 = clf.find( "-divfreeorigin", Vector(-5,-5,-5), "Origin of divfree velocity");

      VolumeGrid<Vector>* divfreevel = new VolumeGrid<Vector>;
      divfreevel->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
      GaussSeidelDivFree( *divfreevel, velocity, 1000, 0.00001 );
      velocity = new GriddedVectorVolume( divfreevel );
   }


   float dt = clf.find("-dt", 0.1f,"Advection time step");
   int nbSteps = clf.find("-nbsteps", 1, "Number of advection steps");
   int nbCorrections = clf.find("-nbcorrections", 1, "Number of advection steps");

   float ddt = dt/nbSteps;

   for( int i=0;i<nbSteps;i++ )
   {
      volume = new BFECCAdvectVolume( volume, velocity, ddt, nbCorrections );
   }
   return volume;
}




Volume<float>* SetUpSelmaSphere( CmdLineFind& clf )
{
   Volume<float>* volume = new ClampVolume( new MultiplyVolume(SetUpSphere( clf ), 3.0) , 0, 1 );

   Noise_t parms;
   parms.translate = clf.find( "-translate", Vector(0,0,0), "Noise translation");
   parms.wavelength = clf.find("-freq", 1.432234f, "Noise frequency");
   parms.amplitude = clf.find("-amp", 1.5f, "Noise amplitude" );
   parms.roughness = clf.find("-rough", 0.5f, "Noise roughness" );
   parms.octaves = clf.find("-octaves", 3.0f, "Noise octaves" );

   Noise * noise = new FractalSum<PerlinNoiseGustavson>();
   noise->setParameters(parms);

   Volume<Vector>* velocity = new NoiseVectorVolume( noise );

   bool makeDivFree = clf.findFlag("-divfree" );
   if( makeDivFree )
   {
      vector<int> NNNN;
      NNNN.push_back( 100 );
      NNNN.push_back( 100 );
      NNNN.push_back( 100 );
      vector<int> NNN = clf.findArray( "-divfreeNxNyNz", NNNN, "Number of grid points in divfree velocity" );

      Vector lX = clf.find( "-divfreesize", Vector(10,10,10), "Volume size of divfree velocity");
      Vector X0 = clf.find( "-divfreeorigin", Vector(-5,-5,-5), "Origin of divfree velocity");

      VolumeGrid<Vector>* divfreevel = new VolumeGrid<Vector>;
      divfreevel->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
      GaussSeidelDivFree( *divfreevel, velocity, 1000, 0.00001 );
      velocity = new GriddedVectorVolume( divfreevel );
   }


   float dt = clf.find("-dt", 0.1f,"Selma time step");
   int nbSteps = clf.find("-nbsteps", 1, "Number of selma steps");

   float ddt = dt/nbSteps;


      vector<int> NNNN;
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   vector<int> NNN = clf.findArray( "-selmaNxNyNz", NNNN, "Number of grid points in selma volume" );

   Vector lX = clf.find( "-selmasize", Vector(10,10,10), "Volume size of selma volume");
   Vector X0 = clf.find( "-selmaorigin", Vector(-5,-5,-5), "Origin of selma volume");

   Volume<Vector>* selmaMap = new IdentityVectorVolume();

   VolumeGrid<Vector>* xGrid0 = 0;
   VolumeGrid<Vector>* xGrid1 = 0;

   for( int i=0;i<nbSteps;i++ )
   {
      if( i%2 == 0 )
      {
         xGrid0 = new VolumeGrid<Vector>;
         xGrid0->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
         selmaMap = StampSELMA( *xGrid0, selmaMap, velocity, ddt );
	 if( xGrid1 ) { delete xGrid1; xGrid1 = 0; }
      }
      else
      {
         xGrid1 = new VolumeGrid<Vector>;
         xGrid1->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
         selmaMap = StampSELMA( *xGrid1, selmaMap, velocity, ddt );
	 if( xGrid0 ) { delete xGrid0; xGrid0 = 0;}
      }
   }
   volume = new WarpVolume( volume, selmaMap );
   return volume;
}


Volume<float>* SetUpBFECCSelmaSphere( CmdLineFind& clf )
{
   Volume<float>* volume = new ClampVolume( new MultiplyVolume(SetUpSphere( clf ), 10.0) , 0, 1 );

   Noise_t parms;
   parms.translate = clf.find( "-translate", Vector(0,0,0), "Noise translation");
   parms.wavelength = clf.find("-freq", 1.432234f, "Noise frequency");
   parms.amplitude = clf.find("-amp", 1.5f, "Noise amplitude" );
   parms.roughness = clf.find("-rough", 0.5f, "Noise roughness" );
   parms.octaves = clf.find("-octaves", 3.0f, "Noise octaves" );

   Noise * noise = new FractalSum<PerlinNoiseGustavson>();
   noise->setParameters(parms);

   Volume<Vector>* velocity = new NoiseVectorVolume( noise );

   bool makeDivFree = clf.findFlag("-divfree" );
   if( makeDivFree )
   {
      vector<int> NNNN;
      NNNN.push_back( 100 );
      NNNN.push_back( 100 );
      NNNN.push_back( 100 );
      vector<int> NNN = clf.findArray( "-divfreeNxNyNz", NNNN, "Number of grid points in divfree velocity" );

      Vector lX = clf.find( "-divfreesize", Vector(10,10,10), "Volume size of divfree velocity");
      Vector X0 = clf.find( "-divfreeorigin", Vector(-5,-5,-5), "Origin of divfree velocity");

      VolumeGrid<Vector>* divfreevel = new VolumeGrid<Vector>;
      divfreevel->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
      GaussSeidelDivFree( *divfreevel, velocity, 1000, 0.00001 );
      velocity = new GriddedVectorVolume( divfreevel );
   }


   float dt = clf.find("-dt", 0.1f,"Selma time step");
   int nbSteps = clf.find("-nbsteps", 1, "Number of selma steps");

   float ddt = dt/nbSteps;


      vector<int> NNNN;
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   vector<int> NNN = clf.findArray( "-selmaNxNyNz", NNNN, "Number of grid points in selma volume" );

   Vector lX = clf.find( "-selmasize", Vector(10,10,10), "Volume size of selma volume");
   Vector X0 = clf.find( "-selmaorigin", Vector(-5,-5,-5), "Origin of selma volume");

   Volume<Vector>* selmaMap = new IdentityVectorVolume();

   VolumeGrid<Vector>* xGrid0 = 0;
   VolumeGrid<Vector>* xGrid1 = 0;

   for( int i=0;i<nbSteps;i++ )
   {
      if( i%2 == 0 )
      {
         xGrid0 = new VolumeGrid<Vector>;
         xGrid0->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
         selmaMap = StampBFECCSELMA( *xGrid0, selmaMap, velocity, ddt );
	 if( xGrid1 ) { delete xGrid1; xGrid1 = 0; }
      }
      else
      {
         xGrid1 = new VolumeGrid<Vector>;
         xGrid1->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
         selmaMap = StampBFECCSELMA( *xGrid1, selmaMap, velocity, ddt );
	 if( xGrid0 ) { delete xGrid0; xGrid0 = 0;}
      }
   }
   volume = new WarpVolume( volume, selmaMap );
   return volume;
}




/*

Volume<float>* SetUpCFDSphere( CmdLineFind& clf )
{
   Volume<float>* volume = new ClampVolume( new MultiplyVolume(SetUpSphere( clf ), 0.1) , 0, 1 );

   float dt = clf.find("-dt", 0.1f,"Time step");
   int nbSteps = clf.find("-nbsteps", 1, "Number of time steps");
   int nbStepsPer = clf.find("-nbstepsperupdate", 1, "Number of steps per update");


   vector<int> NNNN;
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   vector<int> NNN = clf.findArray( "-cfdNxNyNz", NNNN, "Number of grid points in cfd volume" );

   Vector lX = clf.find( "-cfdsize", Vector(10,10,10), "Volume size of selma volume");
   Vector X0 = clf.find( "-cfdorigin", Vector(-5,-5,-5), "Origin of selma volume");

   Vector gravity = clf.find("-gravity", Vector( 0, -10, 0 ) );
   float couple = clf.find("-coupling", 1.0f );
   float refdensity = clf.find("-refdensity", 0.0f );
   float errortolerance = clf.find("-error", 0.0001f );
   int nbprojections = clf.find("-nbprojections", 1000000 );

   int noisestart = clf.find("-noisestart", 1 );
   int noiseend = clf.find("-noiseend", 1000 );
   int blendlength = clf.find("-blendlength", 7 );

   string advectionmethod = clf.find("-advectionmethod", "semilagrangian" );

   GasSystem cfdgas( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0, gravity, couple, refdensity ); 
   cfdgas.InitializeDensity( volume );
   cfdgas.setStepsPerUpdate( nbStepsPer );
   cfdgas.setProjectionTolerance( errortolerance );
   cfdgas.setNbDivFreeIterations( nbprojections );

   vector<string> mapMenu;
   mapMenu.push_back("none");
   mapMenu.push_back("ga");
   mapMenu.push_back("ga3");
   mapMenu.push_back("selma");
   string mapchoice = clf.findMenu("-mapping", mapMenu );

   vector<int> selmaNNN;
   Vector selmaX;
   Vector selmaX0;
   if( mapchoice == "selma" )
   {
      selmaNNN = clf.findArray( "-selmaNxNyNz", NNN, "Number of grid points in selma volume" );
      selmaX = clf.find( "-selmasize", lX, "Volume size of selma volume");
      selmaX0 = clf.find( "-selmaorigin", X0, "Origin of selma volume");
      cout << "using selma map\n";
   }

   VolumeGrid<Vector>* xGrid0 = 0;
   VolumeGrid<Vector>* xGrid1 = 0;
   VolumeGrid<Vector>* xNoiseGrid0 = 0;
   VolumeGrid<Vector>* xNoiseGrid1 = 0;
   Volume<Vector>* selmaMap = new IdentityVectorVolume();
   Volume<Vector>* selmaNoiseMap = new IdentityVectorVolume();



   bool injectVelocityNoise = clf.findFlag("-injectvelocitynoise" );
   if( injectVelocityNoise )
   {
      Noise_t parms;
      parms.translate = clf.find( "-velocitytranslate", Vector(0,0,0), "Noise translation");
      parms.wavelength = clf.find("-velocityfreq", 1.432234f, "Noise frequency");
      parms.amplitude = clf.find("-velocityamp", 1.5f, "Noise amplitude" );
      parms.roughness = clf.find("-velocityrough", 0.5f, "Noise roughness" );
      parms.octaves = clf.find("-velocityoctaves", 3.0f, "Noise octaves" );

      Noise * noise = new FractalSum<PerlinNoiseGustavson>();
      noise->setParameters(parms);

      cfdgas.InitializeVelocity( new MultiplyVectorVolume( new NoiseVectorVolume( noise ), new MaskVolume( volume ) ) );
      cout << "Noise injected into the velocity" << endl;
   }


   bool injectMapNoise = clf.findFlag("-injectmapnoise" );
   Noise* mapnoise = 0;
   Noise_t mapparms;
   if( injectMapNoise )
   {
      mapparms.translate = clf.find( "-maptranslate", Vector(0,0,0), "Noise translation");
      mapparms.wavelength = clf.find("-mapfreq", 1.432234f, "Noise frequency");
      mapparms.amplitude = clf.find("-mapamp", 1.5f, "Noise amplitude" );
      mapparms.roughness = clf.find("-maprough", 0.5f, "Noise roughness" );
      mapparms.octaves = clf.find("-mapoctaves", 3.0f, "Noise octaves" );

      mapnoise = new FractalSum<PerlinNoiseGustavson>();
      mapnoise->setParameters(mapparms);

      cout << "Noise injected into the map" << endl;
   }



   vector< VolumeGrid<Vector>* > ga3U;
   if( mapchoice == "ga3" )
   {
      ga3U.push_back( new VolumeGrid<Vector> );
      ga3U.push_back( new VolumeGrid<Vector> );
      ga3U.push_back( new VolumeGrid<Vector> );
      ga3U[0]->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
      ga3U[1]->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
      ga3U[2]->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
   }


   for( int i=0;i<nbSteps;i++ )
   {
      cout << "\n\n============= STEP " << i << " =====================\n\n" << flush;
      cfdgas.update(dt, advectionmethod);
      if( mapchoice == "selma" )
      {
         if( mapnoise )
	 {
	    mapparms.time = dt * i;
            mapnoise->setParameters(mapparms);
	    selmaNoiseMap = new AdvectVectorVolume( selmaNoiseMap, new MultiplyVectorVolume( new NoiseVectorVolume( mapnoise ),  mapparms.amplitude ) , dt );
	 }
         if( i%2 == 0 )
         {
            xGrid0 = new VolumeGrid<Vector>;
            xGrid0->init( selmaNNN[0], selmaNNN[1], selmaNNN[2], selmaX[0], selmaX[1], selmaX[2], selmaX0 );
	    xGrid0->setClearValue( Vector(0,0,0) );
            selmaMap = StampSELMA( *xGrid0, selmaMap, new GriddedVectorVolume( &cfdgas.getVelocityGrid() ), dt );
	    if( xGrid1 ) { delete xGrid1; xGrid1 = 0; }
            xNoiseGrid0 = new VolumeGrid<Vector>;
            xNoiseGrid0->init( selmaNNN[0], selmaNNN[1], selmaNNN[2], selmaX[0], selmaX[1], selmaX[2], selmaX0 );
	    xNoiseGrid0->setClearValue( Vector(0,0,0) );
            selmaNoiseMap = StampSELMA( *xNoiseGrid0, selmaNoiseMap, new GriddedVectorVolume( &cfdgas.getVelocityGrid() ), dt );
	    if( xNoiseGrid1 ) { delete xNoiseGrid1; xNoiseGrid1 = 0; }
         }
         else
         {
            xGrid1 = new VolumeGrid<Vector>;
            xGrid1->init( selmaNNN[0], selmaNNN[1], selmaNNN[2], selmaX[0], selmaX[1], selmaX[2], selmaX0 );
	    xGrid1->setClearValue( Vector(0,0,0) );
            selmaMap = StampSELMA( *xGrid1, selmaMap, new GriddedVectorVolume( &cfdgas.getVelocityGrid() ), dt );
	    if( xGrid0 ) { delete xGrid0; xGrid0 = 0;}
            xNoiseGrid1 = new VolumeGrid<Vector>;
            xNoiseGrid1->init( selmaNNN[0], selmaNNN[1], selmaNNN[2], selmaX[0], selmaX[1], selmaX[2], selmaX0 );
	    xNoiseGrid1->setClearValue( Vector(0,0,0) );
            selmaNoiseMap = StampSELMA( *xNoiseGrid1, selmaNoiseMap, new GriddedVectorVolume( &cfdgas.getVelocityGrid() ), dt );
	    if( xNoiseGrid0 ) { delete xNoiseGrid0; xNoiseGrid0 = 0; }
         }
      }
      else if( mapchoice == "ga" )
      {
         VolumeGrid<Vector>* xGrid = new VolumeGrid<Vector>;
	 const VolumeGrid<Vector>& gg = cfdgas.getVelocityGrid();
         xGrid->init( gg.nx(), gg.ny(), gg.nz(), gg.Lx(), gg.Ly(), gg.Lz(), gg.llc() );
	 Sample( xGrid, new GriddedVectorVolume( &cfdgas.getVelocityGrid() ) );
         selmaMap = new AdvectVectorVolume( selmaMap, new GriddedVectorVolume( xGrid ), dt );
      }
      else if( mapchoice == "ga3" )
      {
         int ii = i%3;
         Sample( ga3U[ii], new GriddedVectorVolume( &cfdgas.getVelocityGrid() ) );
      }
   }

   bool injectDensityNoise = clf.findFlag("-injectdensitynoise" );
   if( injectDensityNoise )
   {
      Noise_t parms;
      parms.translate = clf.find( "-densitytranslate", Vector(0,0,0), "Noise translation");
      parms.wavelength = clf.find("-densityfreq", 1.432234f, "Noise frequency");
      parms.amplitude = clf.find("-densityamp", 1.5f, "Noise amplitude" );
      parms.roughness = clf.find("-densityrough", 0.5f, "Noise roughness" );
      parms.octaves = clf.find("-densityoctaves", 3.0f, "Noise octaves" );

      Noise * noise = new FractalSum<PerlinNoiseGustavson>();
      noise->setParameters(parms);

      volume = new ClampVolume( new AddVolume( volume, new MultiplyVolume( new NoiseVolume( noise ),  new MaskVolume(volume) ) ), 0.0, 100.0 );
      cout << "Noise injected into the density field" << endl;
   }

   if( mapchoice == "ga" )
   {
      volume = new MultiplyVolume( new WarpVolume( volume, selmaMap ), new DetGradMapVolume( selmaMap ) );
      cout << "Constructed " + mapchoice + " advected volume\n";
   }
   else if( mapchoice == "selma" )
   {
      if( nbSteps >= noisestart && nbSteps <= noiseend )
      {
         float alpha = 1.0;
	 if( nbSteps-noisestart < blendlength ){ alpha = (float)(nbSteps-noisestart)/blendlength; }
	 if( noiseend-nbSteps < blendlength ){ alpha = (float)(noiseend-nbSteps)/blendlength; }
	 cout << "\bSelma Blend Alpha = " << alpha << endl << endl;
	 selmaMap = new AddVectorVolume( new MultiplyVectorVolume( selmaMap, 1.0-alpha ), new MultiplyVectorVolume( selmaNoiseMap, alpha ) );
      }
      volume = new MultiplyVolume( new WarpVolume( volume, selmaMap ), new DetGradMapVolume( selmaMap ) );
      cout << "Constructed " + mapchoice + " advected volume\n";
   }
   else if( mapchoice == "ga3" )
   {
      int ii0 = nbSteps%3;
      int ii1 = (nbSteps-1)%3;
      int ii2 = (nbSteps-2)%3;
      selmaMap = new IdentityVectorVolume;
      if( ii2 >= 0 ){ selmaMap = new AdvectVectorVolume( selmaMap, new GriddedVectorVolume( ga3U[ii2] ), dt ); }
      if( ii1 >= 0 ){ selmaMap = new AdvectVectorVolume( selmaMap, new GriddedVectorVolume( ga3U[ii1] ), dt ); }
      selmaMap = new AdvectVectorVolume( selmaMap, new GriddedVectorVolume( ga3U[ii0] ), dt );
      volume = new MultiplyVolume( new WarpVolume( volume, selmaMap ), new DetGradMapVolume( selmaMap ) );
      cout << "Constructed " + mapchoice + " advected volume\n";
   }
   else
   {
      // copy into a grid we will keep
      VolumeGrid<float>* dengrid = new VolumeGrid<float>();
      dengrid->init(  NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
      Volume<float>* cfddensity = new GriddedVolume( &cfdgas.getDensityGrid() );
      Sample( dengrid, cfddensity );
      volume = new GriddedVolume( dengrid );
   }

   int sweeten = clf.find("-sweeten", 0 );
   float sweetendt = clf.find("-sweetendt", dt );
   if( sweeten > 0 )
   {
      VolumeGrid<Vector>* velgrid = new VolumeGrid<Vector>();
      velgrid->init(  NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
      Volume<Vector>* velg = new GriddedVectorVolume( &cfdgas.getVelocityGrid() );
      Sample( velgrid, velg );
      Volume<Vector>* vel = new GriddedVectorVolume( velgrid );
      for( int i=0;i<sweeten;i++ )
      {
         volume = new AdvectVolume( volume, vel, sweetendt );
      }
   }


   return volume;
}

*/



void SetUpCFDObstacle( CmdLineFind& clf, RenderData& data )
{
   Volume<float>* volume = new ClampVolume( new MultiplyVolume(SetUpSphere( clf ), 0.1) , 0, 1 );


   Volume<float>* obstacle = new MaskVolume( SetUpObj( clf ) );
   Volume<float>* obstaclecomplement = new SubtractVolume( new ConstantVolume(1.0), obstacle );


   float dt = clf.find("-dt", 0.1f,"Time step");
   int nbSteps = clf.find("-nbsteps", 1, "Number of time steps");
   int nbStepsPer = clf.find("-nbstepsperupdate", 1, "Number of steps per update");


   vector<int> NNNN;
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   vector<int> NNN = clf.findArray( "-cfdNxNyNz", NNNN, "Number of grid points in cfd volume" );

   Vector lX = clf.find( "-cfdsize", Vector(10,10,10), "Volume size of selma volume");
   Vector X0 = clf.find( "-cfdorigin", Vector(-5,-5,-5), "Origin of selma volume");

   Vector gravity = clf.find("-gravity", Vector( 0, -10, 0 ) );
   float couple = clf.find("-coupling", 1.0f );
   float refdensity = clf.find("-refdensity", 0.0f );
   float errortolerance = clf.find("-error", 0.0001f );
   int nbprojections = clf.find("-nbprojections", 1000000 );

   int noisestart = clf.find("-noisestart", 1 );
   int noiseend = clf.find("-noiseend", 1000 );
   int blendlength = clf.find("-blendlength", 7 );

   string advectionmethod = clf.find("-advectionmethod", "semilagrangian" );

   GasSystem cfdgas( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0, gravity, couple, refdensity ); 
   cfdgas.InitializeDensity( volume );
   cfdgas.setStepsPerUpdate( nbStepsPer );
   cfdgas.setProjectionTolerance( errortolerance );
   cfdgas.setNbDivFreeIterations( nbprojections );

   vector<string> mapMenu;
   mapMenu.push_back("none");
   mapMenu.push_back("ga");
   mapMenu.push_back("ga3");
   mapMenu.push_back("selma");
   string mapchoice = clf.findMenu("-mapping", mapMenu );

   vector<int> selmaNNN;
   Vector selmaX;
   Vector selmaX0;
   if( mapchoice == "selma" )
   {
      selmaNNN = clf.findArray( "-selmaNxNyNz", NNN, "Number of grid points in selma volume" );
      selmaX = clf.find( "-selmasize", lX, "Volume size of selma volume");
      selmaX0 = clf.find( "-selmaorigin", X0, "Origin of selma volume");
      cout << "using selma map\n";
   }

   VolumeGrid<Vector>* xGrid0 = 0;
   VolumeGrid<Vector>* xGrid1 = 0;
   VolumeGrid<Vector>* xNoiseGrid0 = 0;
   VolumeGrid<Vector>* xNoiseGrid1 = 0;
   Volume<Vector>* selmaMap = new IdentityVectorVolume();
   Volume<Vector>* selmaNoiseMap = new IdentityVectorVolume();



   bool injectVelocityNoise = clf.findFlag("-injectvelocitynoise" );
   if( injectVelocityNoise )
   {
      Noise_t parms;
      parms.translate = clf.find( "-velocitytranslate", Vector(0,0,0), "Noise translation");
      parms.wavelength = clf.find("-velocityfreq", 1.432234f, "Noise frequency");
      parms.amplitude = clf.find("-velocityamp", 1.5f, "Noise amplitude" );
      parms.roughness = clf.find("-velocityrough", 0.5f, "Noise roughness" );
      parms.octaves = clf.find("-velocityoctaves", 3.0f, "Noise octaves" );

      Noise * noise = new FractalSum<PerlinNoiseGustavson>();
      noise->setParameters(parms);

      cfdgas.InitializeVelocity( new MultiplyVectorVolume( new NoiseVectorVolume( noise ), new MaskVolume( volume ) ) );
      cout << "Noise injected into the velocity" << endl;
   }


   bool injectMapNoise = clf.findFlag("-injectmapnoise" );
   Noise* mapnoise = 0;
   Noise_t mapparms;
   if( injectMapNoise )
   {
      mapparms.translate = clf.find( "-maptranslate", Vector(0,0,0), "Noise translation");
      mapparms.wavelength = clf.find("-mapfreq", 1.432234f, "Noise frequency");
      mapparms.amplitude = clf.find("-mapamp", 1.5f, "Noise amplitude" );
      mapparms.roughness = clf.find("-maprough", 0.5f, "Noise roughness" );
      mapparms.octaves = clf.find("-mapoctaves", 3.0f, "Noise octaves" );

      mapnoise = new FractalSum<PerlinNoiseGustavson>();
      mapnoise->setParameters(mapparms);

      cout << "Noise injected into the map" << endl;
   }



   vector< VolumeGrid<Vector>* > ga3U;
   if( mapchoice == "ga3" )
   {
      ga3U.push_back( new VolumeGrid<Vector> );
      ga3U.push_back( new VolumeGrid<Vector> );
      ga3U.push_back( new VolumeGrid<Vector> );
      ga3U[0]->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
      ga3U[1]->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
      ga3U[2]->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
   }


   for( int i=0;i<nbSteps;i++ )
   {
      cout << "\n\n============= STEP " << i << " =====================\n\n" << flush;
      cfdgas.update(dt, advectionmethod);
      if( mapchoice == "selma" )
      {
         if( mapnoise )
	 {
	    mapparms.time = dt * i;
            mapnoise->setParameters(mapparms);
	    selmaNoiseMap = new AdvectVectorVolume( selmaNoiseMap, new MultiplyVectorVolume( new NoiseVectorVolume( mapnoise ),  mapparms.amplitude ) , dt );
	 }
	 Volume<Vector>* obstU = new MultiplyVectorVolume( new GriddedVectorVolume( &cfdgas.getVelocityGrid() ), obstaclecomplement );
         if( i%2 == 0 )
         {
            xGrid0 = new VolumeGrid<Vector>;
            xGrid0->init( selmaNNN[0], selmaNNN[1], selmaNNN[2], selmaX[0], selmaX[1], selmaX[2], selmaX0 );
	    xGrid0->setClearValue( Vector(0,0,0) );
            selmaMap = StampSELMA( *xGrid0, selmaMap, obstU, dt );
	    if( xGrid1 ) { delete xGrid1; xGrid1 = 0; }
            xNoiseGrid0 = new VolumeGrid<Vector>;
            xNoiseGrid0->init( selmaNNN[0], selmaNNN[1], selmaNNN[2], selmaX[0], selmaX[1], selmaX[2], selmaX0 );
	    xNoiseGrid0->setClearValue( Vector(0,0,0) );
            selmaNoiseMap = StampSELMA( *xNoiseGrid0, selmaNoiseMap, obstU, dt );
	    if( xNoiseGrid1 ) { delete xNoiseGrid1; xNoiseGrid1 = 0; }
         }
         else
         {
            xGrid1 = new VolumeGrid<Vector>;
            xGrid1->init( selmaNNN[0], selmaNNN[1], selmaNNN[2], selmaX[0], selmaX[1], selmaX[2], selmaX0 );
	    xGrid1->setClearValue( Vector(0,0,0) );
            selmaMap = StampSELMA( *xGrid1, selmaMap, obstU, dt );
	    if( xGrid0 ) { delete xGrid0; xGrid0 = 0;}
            xNoiseGrid1 = new VolumeGrid<Vector>;
            xNoiseGrid1->init( selmaNNN[0], selmaNNN[1], selmaNNN[2], selmaX[0], selmaX[1], selmaX[2], selmaX0 );
	    xNoiseGrid1->setClearValue( Vector(0,0,0) );
            selmaNoiseMap = StampSELMA( *xNoiseGrid1, selmaNoiseMap, obstU, dt );
	    if( xNoiseGrid0 ) { delete xNoiseGrid0; xNoiseGrid0 = 0; }
         }
         selmaMap = new AddVectorVolume(  new MultiplyVectorVolume( new IdentityVectorVolume(), obstacle ), new MultiplyVectorVolume( selmaMap , obstaclecomplement) );
         selmaNoiseMap = new AddVectorVolume(  new MultiplyVectorVolume( new IdentityVectorVolume(), obstacle ), new MultiplyVectorVolume( selmaNoiseMap , obstaclecomplement) );

      }
      else if( mapchoice == "ga" )
      {
         VolumeGrid<Vector>* xGrid = new VolumeGrid<Vector>;
	 const VolumeGrid<Vector>& gg = cfdgas.getVelocityGrid();
         xGrid->init( gg.nx(), gg.ny(), gg.nz(), gg.Lx(), gg.Ly(), gg.Lz(), gg.llc() );
	 Sample( xGrid, new GriddedVectorVolume( &cfdgas.getVelocityGrid() ) );
         selmaMap = new AdvectVectorVolume( selmaMap, new GriddedVectorVolume( xGrid ), dt );
      }
      else if( mapchoice == "ga3" )
      {
         int ii = i%3;
         Sample( ga3U[ii], new GriddedVectorVolume( &cfdgas.getVelocityGrid() ) );
      }
   }

   bool injectDensityNoise = clf.findFlag("-injectdensitynoise" );
   if( injectDensityNoise )
   {
      Noise_t parms;
      parms.translate = clf.find( "-densitytranslate", Vector(0,0,0), "Noise translation");
      parms.wavelength = clf.find("-densityfreq", 1.432234f, "Noise frequency");
      parms.amplitude = clf.find("-densityamp", 1.5f, "Noise amplitude" );
      parms.roughness = clf.find("-densityrough", 0.5f, "Noise roughness" );
      parms.octaves = clf.find("-densityoctaves", 3.0f, "Noise octaves" );

      Noise * noise = new FractalSum<PerlinNoiseGustavson>();
      noise->setParameters(parms);

      volume = new ClampVolume( new AddVolume( volume, new MultiplyVolume( new NoiseVolume( noise ),  new MaskVolume(volume) ) ), 0.0, 100.0 );
      cout << "Noise injected into the density field" << endl;
   }

   if( mapchoice == "ga" )
   {
      volume = new WarpVolume( volume, selmaMap );
      cout << "Constructed " + mapchoice + " advected volume\n";
   }
   else if( mapchoice == "selma" )
   {
      if( nbSteps >= noisestart && nbSteps <= noiseend )
      {
         float alpha = 1.0;
	 if( nbSteps-noisestart < blendlength ){ alpha = (float)(nbSteps-noisestart)/blendlength; }
	 if( noiseend-nbSteps < blendlength ){ alpha = (float)(noiseend-nbSteps)/blendlength; }
	 cout << "\bSelma Blend Alpha = " << alpha << endl << endl;
	 selmaMap = new AddVectorVolume( new MultiplyVectorVolume( selmaMap, 1.0-alpha ), new MultiplyVectorVolume( selmaNoiseMap, alpha ) );
      }
      volume = new WarpVolume( volume, selmaMap );
      cout << "Constructed " + mapchoice + " advected volume\n";
   }
   else if( mapchoice == "ga3" )
   {
      int ii0 = nbSteps%3;
      int ii1 = (nbSteps-1)%3;
      int ii2 = (nbSteps-2)%3;
      selmaMap = new IdentityVectorVolume;
      if( ii2 >= 0 ){ selmaMap = new AdvectVectorVolume( selmaMap, new GriddedVectorVolume( ga3U[ii2] ), dt ); }
      if( ii1 >= 0 ){ selmaMap = new AdvectVectorVolume( selmaMap, new GriddedVectorVolume( ga3U[ii1] ), dt ); }
      selmaMap = new AdvectVectorVolume( selmaMap, new GriddedVectorVolume( ga3U[ii0] ), dt );
      volume = new WarpVolume( volume, selmaMap );
      cout << "Constructed " + mapchoice + " advected volume\n";
   }
   else
   {
      // copy into a grid we will keep
      VolumeGrid<float>* dengrid = new VolumeGrid<float>();
      dengrid->init(  NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
      Volume<float>* cfddensity = new GriddedVolume( &cfdgas.getDensityGrid() );
      Sample( dengrid, cfddensity );
      volume = new GriddedVolume( dengrid );
   }

   int sweeten = clf.find("-sweeten", 0 );
   float sweetendt = clf.find("-sweetendt", dt );
   if( sweeten > 0 )
   {
      VolumeGrid<Vector>* velgrid = new VolumeGrid<Vector>();
      velgrid->init(  NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
      Volume<Vector>* velg = new GriddedVectorVolume( &cfdgas.getVelocityGrid() );
      Sample( velgrid, velg );
      Volume<Vector>* vel = new GriddedVectorVolume( velgrid );
      for( int i=0;i<sweeten;i++ )
      {
         volume = new AdvectVolume( volume, vel, sweetendt );
      }
   }

  
   Vector objcolor = clf.find("-objcolor", Vector( 1.5, 0.75, 0.5 ) );
   float objdensity = clf.find("-objdensity", 3.0f );
   float thresh = clf.find("-clamp", 0.1f );
   volume = new ClampVolume( new MultiplyVolume( volume, 1.0/thresh ), 0.0, 1.0 );
   data.colorField = new AddColor( new VolumeMultiplyColor( new ConstantColor( Color(1,1,1,1) ), obstaclecomplement ),   new VolumeMultiplyColor( new ConstantColor( Color( objcolor[0], objcolor[1], objcolor[2],1) ), obstacle ) );
   data.densityField = new AddVolume( volume, new MultiplyVolume( obstacle, objdensity ) );
}










Volume<float>* SetUpGradCrossGradAdvectedSphere( CmdLineFind& clf )
{
   Volume<float>* volume = SetUpSphere( clf );

   

   Noise_t parms;
   parms.translate = clf.find( "-translate", Vector(0,0,0), "Noise translation");
   parms.wavelength = clf.find("-freq", 1.432234f, "Noise frequency");
   parms.amplitude = clf.find("-amp", 1.5f, "Noise amplitude" );
   parms.roughness = clf.find("-rough", 0.5f, "Noise roughness" );
   parms.octaves = clf.find("-octaves", 3.0f, "Noise octaves" );

   Noise * noise = new FractalSum<PerlinNoise>();
   noise->setParameters(parms);

   Volume<float>* noise1 = new NoiseVolume( noise );
   Volume<float>* noise2 = new TranslateVolume( new NoiseVolume( noise ), Vector( 34.374764, -1.84757, 100.8848) );

   Volume<Vector>* velocity = new CrossProductVectorVolume( new GradientVectorVolume( noise1 ), new GradientVectorVolume( noise2 ) );

   float dt = clf.find("-dt", 0.1f,"Advection time step");
   int nbSteps = clf.find("-nbsteps", 1, "Number of advection steps");

   float ddt = dt/nbSteps;

   for( int i=0;i<nbSteps;i++ )
   {
      volume = new AdvectVolume( volume, velocity, ddt );
   }
   return volume;
}




Volume<float>* SetUpCurlAdvectedSphere( CmdLineFind& clf )
{
   Volume<float>* volume = SetUpSphere( clf );

   Noise_t parms;
   parms.translate = clf.find( "-translate", Vector(0,0,0), "Noise translation");
   parms.wavelength = clf.find("-freq", 1.432234f, "Noise frequency");
   parms.amplitude = clf.find("-amp", 1.5f, "Noise amplitude" );
   parms.roughness = clf.find("-rough", 0.5f, "Noise roughness" );
   parms.octaves = clf.find("-octaves", 3.0f, "Noise octaves" );

   Noise * noise = new FractalSum<PerlinNoiseGustavson>();
   noise->setParameters(parms);

   Volume<Vector>* velocity = new NoiseVectorVolume( noise );
   velocity =  new CurlVectorVolume( velocity );

   float dt = clf.find("-dt", 0.1f,"Advection time step");
   int nbSteps = clf.find("-nbsteps", 1, "Number of advection steps");

   float ddt = dt/nbSteps;

   for( int i=0;i<nbSteps;i++ )
   {
      volume = new AdvectVolume( volume, velocity, ddt );
   }
   return volume;
}


Volume<float>* SetUpBamf( CmdLineFind& clf )
{
   vector<int> NNNN;
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   vector<int> NNN = clf.findArray( "-bamfNxNyNz", NNNN, "Number of grid points in bamf volume" );
   Vector lX = clf.find( "-bamfsize", Vector(10,10,10), "Volume size of bamf volume");
   Vector X0 = clf.find( "-bamforigin", Vector(-5,-5,-5), "Origin of bamf volume");

   VolumeGrid<float>* bamfGrid = new VolumeGrid<float>;
   bamfGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );



   Vector bamfEmissionPoint = clf.find("-bamfemissionpoint", Vector(0,0,0) );
   float bamfEmissionRate = clf.find("-bamfemissionrate", 400.0f );
   float time = clf.find("-time", 1.0f );
   float dt = clf.find("-dt", 0.04f );

  Noise_t parms;
   parms.translate = clf.find( "-translate", Vector(0,0,0), "Noise translation");
   parms.wavelength = clf.find("-freq", 1.432234f, "Noise frequency");
   parms.amplitude = clf.find("-amp", 1.5f, "Noise amplitude" );
   parms.roughness = clf.find("-rough", 0.5f, "Noise roughness" );
   parms.octaves = clf.find("-octaves", 3.0f, "Noise octaves" );

   Noise * noise = new FractalSum<PerlinNoiseGustavson>();
   noise->setParameters(parms);

   Volume<Vector>* velocity = new NoiseVectorVolume( noise );

   return new GriddedVolume( bamfGrid );
}



void SetUpGriddedVolume( CmdLineFind& clf, RenderData& data )
{
   vector<int> NNNN;
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   NNNN.push_back( 100 );
   vector<int> NNN = clf.findArray( "-stampNxNyNz", NNNN, "Number of grid points in stamped volume" );

   Vector lX = clf.find( "-stampsize", Vector(10,10,10), "Volume size of stamped volume");
   Vector X0 = clf.find( "-stamporigin", Vector(-5,-5,-5), "Origin of stamped volume");

   VolumeGrid<float>* densityGrid = new VolumeGrid<float>;
   VolumeGrid<Color>* colorGrid = new VolumeGrid<Color>;

   densityGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );
   colorGrid->init( NNN[0], NNN[1], NNN[2], lX[0], lX[1], lX[2], X0 );

   Sample( densityGrid, data.densityField );
   Sample( colorGrid, data.colorField );


   data.densityField = new GriddedVolume( densityGrid );
   data.colorField   = new GriddedColor( colorGrid );

}






//Volume<float>* SetUpVolume( CmdLineFind& clf )
void SetUpVolume( CmdLineFind& clf, RenderData& data )
{
   vector<string> volumeMenu;
   volumeMenu.push_back( "sphere" );
   volumeMenu.push_back( "pyrosphere" );
   volumeMenu.push_back( "blinnspheres" );
   volumeMenu.push_back( "blinnpyrospheres" );
   volumeMenu.push_back( "pyramid" );
   volumeMenu.push_back( "torus" );
   volumeMenu.push_back( "torii" );
   volumeMenu.push_back( "obj" );
   volumeMenu.push_back( "cutoutspheres" );
   volumeMenu.push_back( "csg" );
   volumeMenu.push_back( "csgspheretorus" );
   volumeMenu.push_back( "steinerpatch" );
   volumeMenu.push_back( "fft" );
   volumeMenu.push_back( "cauliflower" );
   volumeMenu.push_back( "stampnoise" );
   volumeMenu.push_back( "ellipse" );
   volumeMenu.push_back( "jack" );
   volumeMenu.push_back( "blendjack" );
   volumeMenu.push_back( "box" );
   volumeMenu.push_back( "cone" );
   volumeMenu.push_back( "icosahedron" );
   volumeMenu.push_back( "boxtorus" );
   volumeMenu.push_back( "stampedparticles" );
   volumeMenu.push_back( "pointwisps" );
   volumeMenu.push_back( "sparsepointwisps" );
   volumeMenu.push_back( "pyro" );
   volumeMenu.push_back( "mobiusstrip" );
   volumeMenu.push_back( "cylinder" );
   volumeMenu.push_back( "tube" );
   volumeMenu.push_back( "pyrotorus" );
   volumeMenu.push_back( "bretzel2" );
   volumeMenu.push_back( "olympicrings" );
   volumeMenu.push_back( "cutteapot" );
   volumeMenu.push_back( "stampobjparticles" );
   volumeMenu.push_back( "fireball" );
   volumeMenu.push_back( "fireballandsmoke" );
   volumeMenu.push_back( "fireballinsmoke" );
   volumeMenu.push_back( "implicithumanoid" );
   volumeMenu.push_back( "displacedteapot" );
   volumeMenu.push_back( "cutandhandledteapot" );
   volumeMenu.push_back( "advectedsphere" );
   volumeMenu.push_back( "advectedothersphere" );
   volumeMenu.push_back( "bfeccadvectedsphere" );
   volumeMenu.push_back( "selmasphere" );
   volumeMenu.push_back( "bfeccselmasphere" );
   volumeMenu.push_back( "advectedgradcrossgradsphere" );
   volumeMenu.push_back( "curladvectedsphere" );
   volumeMenu.push_back( "stampnoiseline" );
   volumeMenu.push_back( "cfdsphere" );
   volumeMenu.push_back( "cfdobstacle" );

   vector<string> volumetypes = clf.findMultiple( "-volume", volumeMenu[0], "Select a predefined volume" );

   vector<string> blendMenu;
   blendMenu.push_back("union");
   blendMenu.push_back("blinn");
   blendMenu.push_back("intersect");
   blendMenu.push_back("cutout");
   vector<string> blendType = clf.findMultiple( "-blendmethod", blendMenu[0], "Select method of blending multiple volumes");

   Volume<float>* assembledresult = 0;

   for( size_t i=0;i<volumetypes.size();i++  )
   {
      string volumetype = volumetypes[i];
      Volume<float>* result = 0;
   if( volumetype == "pyrosphere" )
   {
      result = SetUpPyroSphere( clf );
   }
   else if( volumetype == "blinnspheres" )
   {
      result = SetUpBlinnSpheres( clf );
   }
   else if( volumetype == "blinnpyrospheres" )
   {
      result = SetUpBlinnPyroSpheres( clf );
   }
   else if( volumetype == "pyramid" )
   {
      result = SetUpPyramid( clf );
   }
   else if( volumetype == "torus" )
   {
      result = SetUpTorus( clf );
   }
   else if( volumetype == "torii" )
   {
      result = SetUpTorii( clf );
   }
   else if( volumetype == "obj" )
   {
      result = SetUpObj( clf );
   }
   else if( volumetype == "cutoutspheres" )
   {
      result = SetUpCutoutSpheres( clf );
   }
   else if( volumetype == "csg" )
   {
      result = SetUpCSG( clf );
   }
   else if( volumetype == "csgspheretorus" )
   {
      result = SetUpCSGSphereTorus( clf );
   }
   else if( volumetype == "steinerpatch" )
   {
      result = SetUpSteinerPatch( clf );
   }
   //else if( volumetype == "fft" )
   //{
   //   result = SetUpFFT( clf );
   //}
   else if( volumetype == "ellipse" )
   {
      result = SetUpEllipse( clf );
   }
   else if( volumetype == "jack" )
   {
      result = SetUpJack( clf );
   }
   else if( volumetype == "blendjack" )
   {
      result = SetUpBlendJack( clf );
   }
   else if( volumetype == "box" )
   {
      result = SetUpBox( clf );
   }
   else if( volumetype == "cone" )
   {
      result = SetUpCone( clf );
   }
   else if( volumetype == "icosahedron" )
   {
      result = SetUpIcosahedron( clf );
   }
   else if( volumetype == "boxtorus" )
   {
      result = SetUpBoxTorus( clf, data );
   }
   else if( volumetype == "mobiusstrip" )
   {
      result = SetUpMobiusStrip( clf );
   }
   else if( volumetype == "cylinder" )
   {
      result = SetUpCylinder( clf );
   }
   else if( volumetype == "tube" )
   {
      result = SetUpTube( clf );
   }
   else if( volumetype == "pyrotorus" )
   {
      result = SetUpPyroTorus( clf );
   }
   else if( volumetype == "bretzel2" )
   {
      result = SetUpBretzel2( clf );
   }
   else if( volumetype == "cutteapot" )
   {
      result = SetUpCutTeapot( clf );
   }
   else if( volumetype == "displacedteapot" )
   {
      result = SetUpCutDisplacedTeapot( clf );
   }
   else if( volumetype == "cutandhandledteapot" )
   {
      result = SetUpCutAndHandledTeapot( clf );
   }
   else if( volumetype == "advectedsphere" )
   {
      result = SetUpAdvectedSphere( clf );
   }
   else if( volumetype == "advectedothersphere" )
   {
      result = SetUpOtherAdvectedSphere( clf );
   }
   else if( volumetype == "bfeccadvectedsphere" )
   {
      result = SetUpBFECCAdvectedSphere( clf );
   }
   else if( volumetype == "selmasphere" )
   {
      result = SetUpSelmaSphere( clf );
   }
   //else if( volumetype == "cfdsphere" )
   //{
   //   result = SetUpCFDSphere( clf );
   //}
   else if( volumetype == "bfeccselmasphere" )
   {
      result = SetUpBFECCSelmaSphere( clf );
   }
   else if( volumetype == "advectedgradcrossgradsphere" )
   {
      result = SetUpGradCrossGradAdvectedSphere( clf );
   }
   else if( volumetype == "curladvectedsphere" )
   {
      result = SetUpCurlAdvectedSphere( clf );
   }
   else if( volumetype == "sphere" )
   {
      // default choice
      result = SetUpSphere( clf );
   }



      if( assembledresult )
      {
         if( blendType[i-1] == "union" )
	 {
	    assembledresult = new UnionVolume( assembledresult, result );
	 }
	 else if( blendType[i-1] == "blinn" )
	 {
	    assembledresult = new BlinnBlendVolume( assembledresult, result );
	 }
	 else if( blendType[i-1] == "intersect" )
	 {
	    assembledresult = new IntersectionVolume( assembledresult, result );
	 }
	 else if( blendType[i-1] == "cutout" )
	 {
	    assembledresult = new CutoutVolume( assembledresult, result );
	 }
      }
      else
      {
         assembledresult = result;
      }

   }

   bool noclamp = clf.findFlag("-noclamp" );
   if( !noclamp )
   {
      float thresh = clf.find( "-clamp", 0.1f, "Value to clamp density at. Negative value suppresses clamp." );
      if( thresh > 0.0 ){ assembledresult = new ClampVolume( new MultiplyVolume( assembledresult, 1.0/thresh ), 0.0, 1.0 ); }
   }
   data.densityField = assembledresult;

   for( size_t i=0;i<volumetypes.size();i++  )
   {
      string volumetype = volumetypes[i];
      if( volumetype == "stampedparticles" )
      {
         SetUpStampedParticleSpheres( clf, data );
      }
      else if( volumetype == "pointwisps" )
      {
         SetUpStampPointWisps( clf, data );
      }
      else if( volumetype == "filewisps" )
      {
         SetUpStampSparseFileWisps( clf, data );
      }
      else if( volumetype == "flamewisps" )
      {
         SetUpFlameWisps( clf, data );
      }
      else if( volumetype == "sparsepointwisps" )
      {
         SetUpStampSparsePointWisps( clf, data );
      }
      else if( volumetype == "pyro" )
      {
         SetUpStampPyro( clf, data );
      }
      else if( volumetype == "cauliflower" )
      {
         SetUpCauliflower( clf, data );
      }
      else if( volumetype == "stampnoise" )
      {
         SetUpStampNoise( clf, data );
      }
      else if( volumetype == "stampnoiseline" )
      {
         SetUpStampNoiseLine( clf, data );
      }
      else if( volumetype == "olympicrings" )
      {
         SetUpOlympicRings( clf, data );
      }
      else if( volumetype == "fireball" )
      {
         SetUpFireBall( clf, data );
      }
      else if( volumetype == "fireballandsmoke" )
      {
         SetUpFireBallAndSmoke( clf, data );
      }
      else if( volumetype == "fireballinsmoke" )
      {
         SetUpFireBallIntoSmoke( clf, data );
      }
      else if( volumetype == "stampobjparticles" )
      {
         SetUpStampObjParticles( clf, data );
      }
      else if( volumetype == "implicithumanoid" )
      {
         SetUpImplicitHumanoid( clf, data );
      }
      else if( volumetype == "cfdobstacle" )
      {
         SetUpCFDObstacle( clf, data );
      }
   }
 
   vector<string>::iterator jitter = find( volumetypes.begin(), volumetypes.end(), string("jitter") );
   if( jitter != volumetypes.end() )
   {
      float jitterradius = clf.find("-jitterradius", 0.1f,"Radius of jitter sampling");
      int nbjittersamples = clf.find("-nbjitter", 1, "Number of jittered samples");
      data.densityField = new JitterSampleVolume( data.densityField, jitterradius, nbjittersamples );
      data.colorField = new JitterSampleColor( data.colorField, jitterradius, nbjittersamples );
      data.ambientColorField = new JitterSampleColor( data.ambientColorField, jitterradius, nbjittersamples );
      cout << "Jitter sampling applied.\n";
   }

   if( clf.findFlag("-gridit") )
   {
      SetUpGriddedVolume( clf, data );
   }
}


void SetUpDSM( CmdLineFind& clf, RenderData& data )
{
   data.dsmField.clear();
   bool noDsm = clf.findFlag( "-nodsm", "Dont calculate a DSM");
   if( noDsm ) { return; }

   vector<int> nxnynz;
   nxnynz.push_back(100);
   nxnynz.push_back(100);
   nxnynz.push_back(100);
   nxnynz = clf.findArray( "-dsmNxNyNz", nxnynz, "DSM size in each direction" );

   Vector L = clf.find( "-dsmlength", Vector(10,10,10), "DSM length in each direction");
   Vector X0 = clf.find( "-dsmllc", Vector(-5,-5,-5), "DSM lower left corner");

   vector<Vector> lightPs;
   vector<Vector> lightCds;
   vector<string> lightMenu;
   lightMenu.push_back("custom");
   lightMenu.push_back("keyrimfill");
   lightMenu.push_back("keyfill");
   lightMenu.push_back("whitekeyrimfill");
   lightMenu.push_back("whitekeyfill");
   string lightstyle = clf.findMenu("-lightstyle", lightMenu, "Use a pre-built lighting arrangement" );
   if( lightstyle == "keyrimfill" )
   {
      lightCds.push_back( Vector( 0,0,1 ) );
      lightCds.push_back( Vector( 0.3,0,0 ) );
      lightCds.push_back( Vector( 0,1,0 ) );

      lightPs.push_back( Vector( 0, 95, 20 ) );
      lightPs.push_back( Vector( 0, -95, 0 ) );
      lightPs.push_back( Vector( 0, 0, -95 ) );
   }
   else if( lightstyle == "keyfill"  )
   {
      lightCds.push_back( Vector( 0,0,1 ) );
      lightCds.push_back( Vector( 0.3,0,0 ) );

      lightPs.push_back( Vector( 0, 95, 20 ) );
      lightPs.push_back( Vector( 0, -95, 0 ) );
   }
   else if( lightstyle == "whitekeyrimfill"  )
   {
      lightCds.push_back( Vector( 1,1,1 ) );
      lightCds.push_back( Vector( 0.3,0.3,0.3 ) );
      lightCds.push_back( Vector( 1,1,1 ) );

      lightPs.push_back( Vector( 0, 95, 20 ) );
      lightPs.push_back( Vector( 0, -95, 0 ) );
      lightPs.push_back( Vector( 0, 0, -95 ) );
   }
   else if( lightstyle == "whitekeyfill"  )
   {
      lightCds.push_back( Vector( 1,1,1 ) );
      lightCds.push_back( Vector( 0.3,0.3,0.3 ) );

      lightPs.push_back( Vector( 0, 95, 20 ) );
      lightPs.push_back( Vector( 0, -95, 0 ) );
   }
   else
   {
      lightPs = clf.findMultiple( "-lightP", Vector( 0, 95, 0 ), "Light position");
      lightCds = clf.findMultiple( "-lightCd", Vector( 1,1,1 ), "Light color");
   }

   const int blur = clf.find("-blurdsm", 0 );

   size_t nblights = lightPs.size();
   if( nblights > lightCds.size() ){ nblights = lightCds.size(); }

   string dsmWriteBaseName = clf.find("-writedsm", "", "Base name of files to write dsms into" );
   string dsmReadBaseName =  clf.find("-readdsm", "", "Base name of files to read dsms from" );

   for( size_t i=0;i<nblights;i++ )
   {
       VolumeGrid<float>* dsmGrid = new VolumeGrid<float>;
       if( dsmReadBaseName != "" )
       {
          string filen = dsmReadBaseName + "." + tostr(i) + ".dsm";
	  ifstream input( filen.c_str() );
          ReadVolumeGrid( *dsmGrid, input );
	  input.close();
	  cout << "DSM file " << filen << " read\n";
       }
       else
       {
          dsmGrid->init( nxnynz[0], nxnynz[1], nxnynz[2], L[0], L[1], L[2], X0 );
          RayMarchDSMAccumulation( data.densityField, lightPs[i], data.ds, *dsmGrid );
          if( blur > 0 )
          {
             for( int b=0;b<blur;b++ ){Blur( *dsmGrid ); }
          }
          if( dsmWriteBaseName != "" )
          {
             string filen = dsmWriteBaseName + "." + tostr(i) + ".dsm";
	     ofstream output( filen.c_str() );
             WriteVolumeGrid( *dsmGrid, output );
	     output.close();
	     cout << "DSM file " << filen << " written\n";
          }
       }


       data.dsmField.push_back( new GriddedVolume(dsmGrid) );
       data.lightColor.push_back( Color( lightCds[i][0], lightCds[i][1], lightCds[i][2], 1 ) );
   }

   data.lightPosition = lightPs;
}
