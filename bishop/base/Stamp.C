

#include "Stamp.h"
#include "ProgressMeter.h"
#include "GridVolumes.h"
#include "ImplicitVectorShapes.h"

using namespace lux;


void lux::StampNoise( VolumeGrid<float>& grid, Noise* noise, const Vector& center, const float radius )
{
   int nx = grid.nx();
   int ny = grid.ny();
   int nz = grid.nz();

   ProgressMeter meter( nx*ny*nz, "StampNoise" );

   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
	    Vector P = grid.evalP(i,j,k);
	    float distance = (P-center).magnitude();
	    if( distance <= radius )
	    {
	       float nv = noise->eval( P ); 
	       if( nv > grid.value(i,j,k) ) { grid.value(i,j,k) = nv; } 
	    }
	    meter.update();
         }
      }
   }
}


// Motion blur version
// For motion blur, it has to assume that it is summing values
void lux::StampNoise( VolumeGrid<float>& grid, Noise* noise, const Vector& center, const float radius, const Vector& vel, const Vector& accel, const float timestep, const int seed )
{
   int nx = grid.nx();
   int ny = grid.ny();
   int nz = grid.nz();

   float dx = grid.dx();
   dx = ( dx < grid.dy() ) ? dx : grid.dy();
   dx = ( dx < grid.dz() ) ? dx : grid.dz();

   float approxDistance = timestep * vel.magnitude() + 0.5*timestep*timestep*accel.magnitude();
   int nbsamples = 1 + (int)(  approxDistance/dx  );



   UniformPRN prn;
   Noise_t parms;
   parms.seed = seed + 646383;
   prn.setParameters(parms);

   int ii,jj,kk;

   ProgressMeter meter( nx*ny*nz, "StampNoise" );
   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
	    Vector P = grid.evalP(i,j,k);
	    float distance = (P-center).magnitude();
	    if( distance <= radius )
	    {
	       float nv = noise->eval( P ) / nbsamples; 

	       // accumulate along the streak
	       if( nv > 0 )
	       {
               for( int q=0;q<nbsamples;q++ )
	       {
	          float local = prn.eval() * timestep;
		  Vector PP = P + local*vel + (0.5*local*local)*accel;
                  if( grid.getGridIndex( PP, ii,jj,kk ) ) { grid.value(ii,jj,kk) += nv; }
	       }
	       }
	    }
	    meter.update();
         }
      }
   }
}





void lux::StampNoiseAndColor( VolumeGrid<float>& grid, VolumeGrid<Color>& cgrid,  Noise* noise, const Vector& center, const float radius, const float fade, const Color& Cd )
{
   int ix, iix, iy,iiy, iz,iiz;
   Vector R(radius, radius,radius);
   if( grid.getBox( center-R, center+R, ix, iix, iy, iiy, iz, iiz ) )
   {

   for( int k=iz;k<=iiz;k++ )
   {
      for( int j=iy;j<=iiy;j++ )
      {
         for( int i=ix;i<=iix;i++ )
         {
	    Vector P = grid.evalP(i,j,k);
	    float distance = (P-center).magnitude();
	    if( distance <= radius )
	    {
	       float fadefactor = (1.0-distance/radius)/fade;
	       if( fade <= 0.0 ) { fadefactor = 1.0; }
	       if( fadefactor > 1.0 ){ fadefactor = 1.0; }
	       if( fadefactor < 0.0 ){ fadefactor = 0.0; }
	       float nv = noise->eval( P ) * fadefactor; 
	       if( nv > grid.value(i,j,k) ) 
	       {
	          grid.value(i,j,k) = nv;
		  cgrid.value(i,j,k) = Cd;
	       } 
	    }
         }
      }
   }
   }
}



void lux::StampNoiseAndColor( VolumeGrid<float>& grid, VolumeGrid<Color>& cgrid,  Noise* noise, const Vector& center, const float radius, const float fade, const Color& Cd, const Vector& vel, const Vector& accel, const float timestep, const int seed  )
{
   float dx = grid.dx();
   dx = ( dx < grid.dy() ) ? dx : grid.dy();
   dx = ( dx < grid.dz() ) ? dx : grid.dz();

   float approxDistance = timestep * vel.magnitude() + 0.5*timestep*timestep*accel.magnitude();
   int nbsamples = 1 + (int)(  approxDistance/dx  + 0.5 );

   float dt = timestep/nbsamples;

   UniformPRN prn;
   Noise_t parms;
   parms.seed = seed + 646383;
   prn.setParameters(parms);

   int ii,jj,kk;
   Color blurCd = Cd/nbsamples;


   int ix, iix, iy,iiy, iz,iiz;
   Vector R(radius, radius,radius);
   if( grid.getBox( center-R, center+R, ix, iix, iy, iiy, iz, iiz ) )
   {
   for( int k=iz;k<=iiz;k++ )
   {
      for( int j=iy;j<=iiy;j++ )
      {
         for( int i=ix;i<=iix;i++ )
         {
	    Vector P = grid.evalP(i,j,k);
            grid.getGridIndex( P, ii,jj,kk );
	    float distance = (P-center).magnitude();
	    if( distance <= radius )
	    {
	       float fadefactor = (1.0-distance/radius)/fade;
	       if( fade <= 0.0 ) { fadefactor = 1.0; }
	       if( fadefactor > 1.0 ){ fadefactor = 1.0; }
	       if( fadefactor < 0.0 ){ fadefactor = 0.0; }
	       float nv = noise->eval( P ) * fadefactor/nbsamples; 
	       
	       // accumulate along the streak
	       if( nv > 0 )
	       {
               for( int q=0;q<nbsamples;q++ )
	       {
	          float local = q*dt;
		  Vector PP = P + local*vel + (0.5*local*local)*accel;
                  if( grid.getGridIndex( PP, ii,jj,kk ) )
		  { 
		     grid.value(ii,jj,kk) += nv; 
		     cgrid.value(ii,jj,kk) += blurCd;
		  }
	       }
	       }
	    }
         }
      }
   }
   }
}







void lux::StampNoise( VolumeGrid<float>& grid, Noise* noise, const ParticleGroupA& particles )
{
   Noise_t parms;
   for( size_t p=0;p<particles.size();p++ )
   {
      const Vector& center = particles[p].P();
      const float& radius = particles[p].pscale();
      parms.octaves = particles[p].octaves();
      parms.roughness = particles[p].roughness();
      parms.frequency = particles[p].freq();
      parms.fjump = particles[p].fjump();
      parms.offset = particles[p].offset();
      parms.translate = particles[p].translate();
      parms.translate = particles[p].translate();
      noise->setParameters( parms );
      if( (particles[p].v().magnitude() == 0 && particles[p].accel().magnitude() == 0) || particles[p].shutter() == 0 )
      {
         StampNoise( grid, noise, center, radius );
      }
      else
      {
         float timestep = particles[p].shutter()/24.0;
         StampNoise( grid, noise, center, radius, particles[p].v(), particles[p].accel(), timestep, particles[p].id() );
      }
   }
}

void lux::StampNoiseAndColor( VolumeGrid<float>& grid, VolumeGrid<Color>& cgrid,  Noise* noise, const ParticleGroupA& particles )
{
   Noise_t parms;
   ProgressMeter meter( particles.size(), "StampNoiseAndColor" );
   for( size_t p=0;p<particles.size();p++ )
   {
      parms.octaves = particles[p].octaves();
      parms.roughness = particles[p].roughness();
      parms.frequency = particles[p].freq();
      parms.fjump = particles[p].fjump();
      parms.offset = particles[p].offset();
      parms.translate = particles[p].translate();
      parms.offset = particles[p].offset();
      noise->setParameters( parms );
      if( (particles[p].v().magnitude() == 0 && particles[p].accel().magnitude() == 0) || particles[p].shutter() == 0 )
      {
         StampNoiseAndColor( grid, cgrid, noise, particles[p].P(), particles[p].pscale() , particles[p].fade(), particles[p].Cd()*particles[p].opacity() );
      }
      else
      {
         float timestep = particles[p].shutter()/24.0;
         StampNoiseAndColor( grid, cgrid, noise, particles[p].P(), particles[p].pscale() , particles[p].fade(), particles[p].Cd()*particles[p].opacity(), 
	                     particles[p].v(), particles[p].accel(), timestep, particles[p].id()  );
      }
      meter.update();
   }
}

void lux::StampField( VolumeGrid<float>& grid, Volume<float>* field )
{
   int nx = grid.nx();
   int ny = grid.ny();
   int nz = grid.nz();

   ProgressMeter meter( nx*ny*nz, "StampField" );

   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
#pragma omp parallel for
         for( int i=0;i<nx;i++ )
         {
	    float nv = field->eval( grid.evalP(i,j,k) ); 
	    if( nv > grid.value(i,j,k) ) { grid.value(i,j,k) = nv; } 
	    meter.update();
         }
      }
   }
}


void lux::StampPointWisps( VolumeGrid<float>& grid, VolumeGrid<Color>& cgrid, const ParticleGroupA& particles )
{

   PointWispWanderer pww;
   int metertotal = 0;

   float dx = grid.dx();
   if( dx > grid.dy() ){ dx = grid.dy(); }
   if( dx > grid.dz() ){ dx = grid.dz(); }

   for( size_t i=0;i<particles.size();i++ )
   {
      metertotal += particles[i].nbWisps();
   }

   ProgressMeter meter( metertotal, "StampPointWisps");
   for( size_t p=0;p<particles.size();p++ )
   {
      pww.reset(particles[p]);
      int nb = particles[p].nbWisps();
      Vector V = particles[p].v();
      Vector A = particles[p].accel();
      float shutter = particles[p].shutter();
      float timestep = shutter/24.0;

         for( int i=0;i<nb;i++ )
         {
            const Vector P = pww.step();
            StampBlurredWisps( grid, cgrid, P,  timestep, V, A, particles[p].opacity(), particles[p].Cd(), particles[p].id() );
            meter.update();
         }
   }
}



void lux::StampPointWisps( VolumeGrid<float>& grid, const ParticleGroupA& particles )
{

   PointWispWanderer pww;
   int metertotal = 0;

   float dx = grid.dx();
   if( dx > grid.dy() ){ dx = grid.dy(); }
   if( dx > grid.dz() ){ dx = grid.dz(); }

   for( size_t i=0;i<particles.size();i++ )
   {
      metertotal += particles[i].nbWisps();
   }

   ProgressMeter meter( metertotal, "StampPointWisps");
   for( size_t p=0;p<particles.size();p++ )
   {
      pww.reset(particles[p]);
      int nb = particles[p].nbWisps();
      Vector V = particles[p].v();
      Vector A = particles[p].accel();
      float opacity = particles[p].opacity();
      float shutter = particles[p].shutter();
      int id = particles[p].id();
      const float timestep = shutter/24.0;
      for( int i=0;i<nb;i++ )
      {
            const Vector P = pww.step();
            // now stamp a dot into the grid
	    StampBlurredWisps( grid, P, timestep, V, A, opacity, id );

            meter.update();
      }
   }
}


void lux::StampPointWisps( SparseGrid& grid, const ParticleGroupA& particles )
{

   PointWispWanderer pww;
   int metertotal = 0;

   float dx = grid.dx();
   if( dx > grid.dy() ){ dx = grid.dy(); }
   if( dx > grid.dz() ){ dx = grid.dz(); }

   for( size_t i=0;i<particles.size();i++ )
   {
      metertotal += particles[i].nbWisps();
   }

   ProgressMeter meter( metertotal, "StampSparsePointWisps");
   for( size_t p=0;p<particles.size();p++ )
   {
      pww.reset(particles[p]);
      int nb = particles[p].nbWisps();
      Vector V = particles[p].v();
      Vector A = particles[p].accel();
      float opacity = particles[p].opacity();
      float shutter = particles[p].shutter();
      int id = particles[p].id();
      const float timestep = shutter/24.0;
      for( int i=0;i<nb;i++ )
      {
            const Vector P = pww.step();
            // now stamp a dot into the grid
	    StampBlurredWisps( grid, P, timestep, V, A, opacity, id );

            meter.update();
      }
   }
}



void lux::StampPointWisps( SparseGrid& grid, SparseColorGrid& cgrid, const ParticleGroupA& particles )
{

   PointWispWanderer pww;
   int metertotal = 0;

   float dx = grid.dx();
   if( dx > grid.dy() ){ dx = grid.dy(); }
   if( dx > grid.dz() ){ dx = grid.dz(); }

   for( size_t i=0;i<particles.size();i++ )
   {
      metertotal += particles[i].nbWisps();
   }

   ProgressMeter meter( metertotal, "StampSparsePointWisps");
   for( size_t p=0;p<particles.size();p++ )
   {
      pww.reset(particles[p]);
      int nb = particles[p].nbWisps();
      Vector V = particles[p].v();
      Vector A = particles[p].accel();
      float opacity = particles[p].opacity();
      float shutter = particles[p].shutter();
      int id = particles[p].id();
      Color cd = particles[p].Cd();
      const float timestep = shutter/24.0;
      for( int i=0;i<nb;i++ )
      {
            const Vector P = pww.step();
            // now stamp a dot into the grid
	    StampBlurredWisps( grid, cgrid, P, timestep, V, A, opacity, cd, id );

            meter.update();
      }
   }
}






void lux::StampPointWisps( SparseGrid& grid, SparseColorGrid& cgrid, const vector<string> files, const vector<Color> typeColors, vector<float>& opacity, float blurScale )
{

   int metertotal = 0;

   float dx = grid.dx();
   if( dx > grid.dy() ){ dx = grid.dy(); }
   if( dx > grid.dz() ){ dx = grid.dz(); }

   //for( size_t i=0;i<particles.size();i++ )
   //{
   //   metertotal += particles[i].nbWisps();
   //}

   ProgressMeter meter( 1, "Stamp Sparse Point Wisp Files");

   Vector V(0,0,0);
   Vector A(0,0,0);
   float timestep = 1.0/24.0;
   Vector P;
   int particleType = 0;
   int id = 0;
   float vx,vy,vz;
   float px,py,pz;
   float mass;


   Vector centerOfMass;
   float llcx=100000000000000, llcy=10000000000, llcz=1000000000, urcx=-100000000000, urcy=-100000000, urcz=-10000000;
   int particleTypePopulation[3] = {0,0,0};

   int pupdate = 100;

   for( size_t ff=0;ff<files.size();ff++ )
   {
      ifstream particles;
      particles.open( files[ff].c_str() );
      if ( particles )
      {
         while( !particles.eof() )
         {
	    particles >> id >> px >> py >> pz >> vx >> vy >> vz >> mass >> particleType;
	    P = Vector(px,py,pz);
	    V = Vector(vx,vy,vz) * blurScale;
	    float opac = opacity[ff] * mass;
            // now stamp a dot into the grid
	    StampBlurredWisps( grid, cgrid, P, timestep, V, A, opac, typeColors[ff], id );
            ++metertotal;

	    centerOfMass += P;
	    llcx = ( llcx > px ) ? px : llcx;
	    llcy = ( llcy > py ) ? py : llcy;
	    llcz = ( llcz > pz ) ? pz : llcz;
	    urcx = ( urcx < px ) ? px : urcx;
	    urcy = ( urcy < py ) ? py : urcy;
	    urcz = ( urcz < pz ) ? pz : urcz;
	    particleTypePopulation[particleType-1]++;

	    if( ( metertotal % pupdate ) == 0 )
	    {
	       cout << "\rParticle " << metertotal << flush;
	    }
         }
      }
   }
   cout << "\n\nTotal particles: " << metertotal << endl << flush;
   centerOfMass /= metertotal;
   cout << "Center of mass: " << centerOfMass[0] << " " << centerOfMass[1] << " " << centerOfMass[2] << endl;
   cout << "Bounding box:  [ " <<  llcx << " " << llcy << " " << llcz << " ]  X  [ " << urcx << " " << urcy << " " << urcz << " ]\n";
   cout << "Particle distribution: " << particleTypePopulation[0] << " (" << (int)(100.0*particleTypePopulation[0]/metertotal) << "%)    "   << particleTypePopulation[1] << " (" << (int)(100.0*particleTypePopulation[1]/metertotal) << "%)    "    << particleTypePopulation[2] << " (" << (int)(100.0*particleTypePopulation[2]/metertotal) << "%)    " <<endl;

   cgrid.normalize( grid );
}











void lux::StampFlameWisps( SparseGrid& grid, SparseColorGrid& cgrid, const ParticleGroupA& particles, const Volume<Color>& baseColor )
{

   PointWispWanderer pww;
   int metertotal = 0;

   float dx = grid.dx();
   if( dx > grid.dy() ){ dx = grid.dy(); }
   if( dx > grid.dz() ){ dx = grid.dz(); }

   for( size_t i=0;i<particles.size();i++ )
   {
      metertotal += particles[i].nbWisps();
   }

   ProgressMeter meter( metertotal, "StampFlameWisps");
   for( size_t p=0;p<particles.size();p++ )
   {
      pww.reset(particles[p]);
      int nb = particles[p].nbWisps();
      Vector V = particles[p].v();
      Vector A = particles[p].accel();
      float opacity = particles[p].opacity();
      float shutter = particles[p].shutter();
      int id = particles[p].id();
      const float timestep = shutter/24.0;
      for( int i=0;i<nb;i++ )
      {
            const Vector P = pww.step();
      	    const Color cd = baseColor.eval( pww.pos() );
            // now stamp a dot into the grid
	    StampBlurredWisps( grid, cgrid, P, timestep, V, A, opacity, cd, id );

            meter.update();
      }
   }
}






void PointWispWanderer::reset( const Particle& p )
{
   part = p;
   Noise_t parms;
   parms.seed = part.id();
   walkNoise.setParameters(parms);

   parms.octaves = part.octaves();
   parms.frequency = part.freq();
   parms.translate = part.translate();
   parms.offset = part.offset();
   parms.fjump = part.fjump();
   parms.roughness = part.roughness();
   displacementNoise.setParameters( parms );

   parms.octaves = part.wispOctaves();
   parms.frequency = part.wispFreq();
   parms.translate = part.wispTranslate();
   parms.offset = part.wispOffset();
   parms.fjump = part.wispFjump();
   parms.roughness = part.wispRoughness();
   transformNoise.setParameters( parms );

   d1 = Vector( walkNoise.eval()-0.5, walkNoise.eval()-0.5, walkNoise.eval()-0.5 );
   d2 = Vector( walkNoise.eval()-0.5, walkNoise.eval()-0.5, walkNoise.eval()-0.5 );
   walkPosition = Vector( walkNoise.eval()-0.5, walkNoise.eval()-0.5, walkNoise.eval()-0.5 );
}




void PointWispWanderer::reset( const Anchor& p )
{
   part.P() = p.P;
   part.Cd() = p.Cd;  
   part.id() = p.seed;
   part.pscale() = p.pscale;
   part.octaves() = p.octaves;
   part.roughness() = p.roughness;
   part.freq() = p.frequency;
   part.fjump() = p.fjump;
   part.offset() = p.offset;
   part.translate() = p.translate;
   part.fade() = p.falloff;
   part.opacity() = p.amplitude;
   part.v() = p.v; 
   part.accel() = p.A;  
   part.normal() = p.tangent; 
   part.right() = p.normal; 
   part.up() = p.binormal;   
   part.nbWisps() = p.nbWisps; 
   part.wispOctaves() = p.wispOctaves;
   part.wispRoughness() = p.wispRoughness; 
   part.wispFreq() = p.wispFreq;  
   part.wispFjump() = p.wispFjump; 
   part.wispOffset() = p.wispOffset; 
   part.wispTranslate() = p.wispTranslate;
   part.shutter() = p.shutter;    
   part.framerate() = p.frameRate; 
   part.lifetime() = p.lifeTime; 
   part.age() = p.age;     
   part.wispCorrelation() = p.wispCorrelation; 
   part.wispRadialGroup() = p.wispRadialGroup;
   part.wispDisplacementScale() = p.wispDisplacementScale;   
   


   walkNoise.setParameters(p);
   displacementNoise.setParameters( p );
   transformNoise.setParameters( p );

   d1 = Vector( walkNoise.eval()-0.5, walkNoise.eval()-0.5, walkNoise.eval()-0.5 );
   d2 = Vector( walkNoise.eval()-0.5, walkNoise.eval()-0.5, walkNoise.eval()-0.5 );
   walkPosition = Vector( walkNoise.eval()-0.5, walkNoise.eval()-0.5, walkNoise.eval()-0.5 );
}


const Vector PointWispWanderer::step()
{
   // do a correlated step
   walkPosition = walkPosition * part.wispCorrelation()
                + Vector( walkNoise.eval()-0.5, walkNoise.eval()-0.5, walkNoise.eval()-0.5 ) * (1.0-part.wispCorrelation());

   // find the unit circle
   Vector wP = walkPosition.unitvector();
   float radius = pow( fabs(transformNoise.eval( walkPosition ) ), part.wispRadialGroup()  ) ;
   wP *= radius;
   wP = wP[0] * part.normal() + wP[1] * part.right() + wP[2] * part.up();

   Vector displacement( displacementNoise.eval(wP), 
                        displacementNoise.eval(wP+d1),  
			displacementNoise.eval(wP+d2) );
   displacement *= part.wispDisplacementScale();
   wP += displacement;
   wP *= part.pscale();
   wP += part.P();

   return wP;
}






void lux::StampPyro( VolumeGrid<float>& grid, VolumeGrid<Color>& cgrid, const ParticleGroupA& particles )
{

   ProgressMeter meter( particles.size(), "StampPyro");
   for( size_t p=0;p<particles.size();p++ )
   {
      StampPyro( grid, cgrid, particles[p] );
      meter.update();
   }
}


void lux::StampPyro(VolumeGrid<float>& grid, VolumeGrid<Color>& cgrid, const Particle& part )
{
   // no motion blur yet
   int ixmin,iymin,izmin;
   int ixmax,iymax,izmax;

   float cellsize = grid.dx();
   cellsize = ( cellsize < grid.dy() ) ? cellsize : grid.dy();
   cellsize = ( cellsize < grid.dz() ) ? cellsize : grid.dz();

   const float outtest = part.pscale() + part.pyroAmplitude()*pow( (double)(1.0+part.roughness()), (double)part.octaves() );
   Vector llc = part.P() - Vector( outtest, outtest, outtest );
   Vector urc = part.P() + Vector( outtest, outtest, outtest );

   const float intest = part.pscale();

   if( !grid.getBox( llc, urc, ixmin, ixmax, iymin, iymax, izmin, izmax ) ){ return; }

   FractalSum<PerlinNoiseGustavson> noise;
   Noise_t parms;
   
   parms.octaves = part.octaves();
   parms.roughness = part.roughness();
   parms.frequency = part.freq();
   parms.translate = part.translate();
   parms.seed = part.id();
   parms.time = part.age();

   noise.setParameters( parms );

   for( int k=izmin;k<=izmax;k++ )
   {
      for( int j=iymin;j<=iymax;j++ )
      {
         for( int i=ixmin;i<=ixmax;i++ )
         {
	     Vector P = grid.evalP(i,j,k) - part.P();
	     float d = P.magnitude();
	     float val = grid.value(i,j,k);
	     Color cval = cgrid.value(i,j,k);
             if( d < intest )
	     {
		 if( part.pyroDensity()*part.opacity() > val )
		 { 
		    val = part.pyroDensity() * part.opacity(); 
		    cval = part.Cd();
		 }
	     }
	     else
	     {
		float displace = fabs( noise.eval( P.unitvector() ) ) * part.pyroAmplitude();
		if( d < part.pscale() + displace )
		{
		   float trim = 1.0;
		   if(  part.pscale() + displace - d < 4.0 *cellsize )
		   {
		      trim = (part.pscale() + displace - d)/(4.0*cellsize);
		   }
		   if( part.pyroDensity()*trim*part.opacity() > val )
		   { 
		      val = part.pyroDensity()*trim*part.opacity(); 
		      cval = part.Cd();
		   }
		}
	     }
	     grid.value(i,j,k) = val;
	     cgrid.value(i,j,k) = cval;
         }
      }
   }
}


UniformPRN  rn;
void lux::StampBlurredWisps( VolumeGrid<float>& grid, VolumeGrid<Color>& cgrid, const Vector& P, const float timestep, const Vector& velocity, const Vector& acceleration, const float opacity, const Color& cd, int seed )
{
   float dx = grid.dx();
   dx = ( dx < grid.dy() ) ? dx : grid.dy();
   dx = ( dx < grid.dz() ) ? dx : grid.dz();

   float approxDistance = timestep * velocity.magnitude() + 0.5*timestep*timestep*acceleration.magnitude();
   int nbsamples = 1 + (int)(  approxDistance/dx  );

   if( nbsamples > 1 )
   {
      Noise_t parms;
      parms.seed = seed;
      rn.setParameters(parms);
   }

   int ix,iy,iz;

   float gvalue = opacity/nbsamples;
   Color cdvalue = cd*opacity/nbsamples;


   for( int q=0;q<nbsamples;q++ )
   {
      float local = rn.eval() * timestep;
      if( q==0 ){ local = 0; }

      Vector position = P + local * velocity + (0.5*local*local)*acceleration;

      if( grid.getGridIndex( position, ix,iy,iz ) )
      { 
         Vector P000 = grid.evalP(ix,iy,iz);
         position -= P000;

         float wx = position[0]/grid.dx();
         float wy = position[1]/grid.dy();
         float wz = position[2]/grid.dz();

         float vx = 1.0-wx;
         float vy = 1.0-wy;
         float vz = 1.0-wz;

         float w000 = vx*vy*vz;
         float w100 = wx*vy*vz;
         float w010 = vx*wy*vz;
         float w001 = vx*vy*wz;
         float w110 = wx*wy*vz;
         float w011 = vx*wy*wz;
         float w101 = wx*vy*wz;
         float w111 = wx*wy*wz;

         bool tx = ( ix < grid.nx()-1 );
         bool ty = ( iy < grid.ny()-1 );
         bool tz = ( iz < grid.nz()-1 );

         grid.value(ix,iy,iz) += gvalue*w000; cgrid.value(ix,iy,iz) += cdvalue*w000;
         if( tx ){ grid.value(ix+1,iy,iz) += gvalue*w100; cgrid.value(ix+1,iy,iz) += cdvalue*w100; }
         if( ty ){ grid.value(ix,iy+1,iz) += gvalue*w010; cgrid.value(ix,iy+1,iz) += cdvalue*w010; }
         if( tz ){ grid.value(ix,iy,iz+1) += gvalue*w001; cgrid.value(ix,iy,iz+1) += cdvalue*w001; }
         if( tx && ty ){ grid.value(ix+1,iy+1,iz) += gvalue*w110; cgrid.value(ix+1,iy+1,iz) += cdvalue*w110; }
         if( tx && tz ){ grid.value(ix+1,iy,iz+1) += gvalue*w101; cgrid.value(ix+1,iy,iz+1) += cdvalue*w101; }
         if( ty && tz ){ grid.value(ix,iy+1,iz+1) += gvalue*w011; cgrid.value(ix,iy+1,iz+1) += cdvalue*w011; }
         if( tx && ty && tz ){ grid.value(ix+1,iy+1,iz+1) += gvalue*w111; cgrid.value(ix+1,iy+1,iz+1) += cdvalue*w111; }
      }
   }
}

void lux::StampBlurredWisps( VolumeGrid<float>& grid, const Vector& P, const float timestep, const Vector& velocity, const Vector& acceleration, const float opacity, int seed )
{
   float dx = grid.dx();
   dx = ( dx < grid.dy() ) ? dx : grid.dy();
   dx = ( dx < grid.dz() ) ? dx : grid.dz();

   float approxDistance = timestep * velocity.magnitude() + 0.5*timestep*timestep*acceleration.magnitude();
   int nbsamples = 1 + (int)(  approxDistance/dx  );

   if( nbsamples > 1 )
   {
      Noise_t parms;
      parms.seed = seed;
      rn.setParameters(parms);
   }

   int ix,iy,iz;

   float gvalue = opacity/nbsamples;
   for( int q=0;q<nbsamples;q++ )
   {
      float local = rn.eval() * timestep;
      if( q==0 ){ local = 0; }

      Vector position = P + local * velocity + (0.5*local*local)*acceleration;

      if( grid.getGridIndex( position, ix,iy,iz ) )
      { 
         Vector P000 = grid.evalP(ix,iy,iz);
         position -= P000;

         float wx = position[0]/grid.dx();
         float wy = position[1]/grid.dy();
         float wz = position[2]/grid.dz();

         float vx = 1.0-wx;
         float vy = 1.0-wy;
         float vz = 1.0-wz;

         float w000 = vx*vy*vz;
         float w100 = wx*vy*vz;
         float w010 = vx*wy*vz;
         float w001 = vx*vy*wz;
         float w110 = wx*wy*vz;
         float w011 = vx*wy*wz;
         float w101 = wx*vy*wz;
         float w111 = wx*wy*wz;

         bool tx = ( ix < grid.nx()-1 );
         bool ty = ( iy < grid.ny()-1 );
         bool tz = ( iz < grid.nz()-1 );

         grid.value(ix,iy,iz) += gvalue*w000;
         if( tx ){ grid.value(ix+1,iy,iz) += gvalue*w100; }
         if( ty ){ grid.value(ix,iy+1,iz) += gvalue*w010; }
         if( tz ){ grid.value(ix,iy,iz+1) += gvalue*w001; }
         if( tx && ty ){ grid.value(ix+1,iy+1,iz) += gvalue*w110; }
         if( tx && tz ){ grid.value(ix+1,iy,iz+1) += gvalue*w101; }
         if( ty && tz ){ grid.value(ix,iy+1,iz+1) += gvalue*w011; }
         if( tx && ty && tz ){ grid.value(ix+1,iy+1,iz+1) += gvalue*w111; }
      }
   }
}









void lux::StampBlurredWisps( SparseGrid& grid, const Vector& P, const float timestep, const Vector& velocity, const Vector& acceleration, const float opacity, int seed )
{
   float dx = grid.dx();
   dx = ( dx < grid.dy() ) ? dx : grid.dy();
   dx = ( dx < grid.dz() ) ? dx : grid.dz();

   float approxDistance = timestep * velocity.magnitude() + 0.5*timestep*timestep*acceleration.magnitude();
   int nbsamples = 1 + (int)(  approxDistance/dx  );

   if( nbsamples > 1 )
   {
      Noise_t parms;
      parms.seed = seed;
      rn.setParameters(parms);
   }

   int ix,iy,iz;

   float gvalue = opacity/nbsamples;
   for( int q=0;q<nbsamples;q++ )
   {
      float local = rn.eval() * timestep;
      if( q==0 ){ local = 0; }

      Vector position = P + local * velocity + (0.5*local*local)*acceleration;

      if( grid.getGridIndex( position, ix,iy,iz ) )
      { 
         Vector P000 = grid.evalP(ix,iy,iz);
         position -= P000;

         float wx = position[0]/grid.dx();
         float wy = position[1]/grid.dy();
         float wz = position[2]/grid.dz();

         float vx = 1.0-wx;
         float vy = 1.0-wy;
         float vz = 1.0-wz;

         float w000 = vx*vy*vz;
         float w100 = wx*vy*vz;
         float w010 = vx*wy*vz;
         float w001 = vx*vy*wz;
         float w110 = wx*wy*vz;
         float w011 = vx*wy*wz;
         float w101 = wx*vy*wz;
         float w111 = wx*wy*wz;

         bool tx = ( ix < grid.nx()-1 );
         bool ty = ( iy < grid.ny()-1 );
         bool tz = ( iz < grid.nz()-1 );

         grid.set(grid.get(ix,iy,iz) + gvalue*w000, ix, iy, iz );
         if( tx ){ grid.set( grid.get(ix+1,iy,iz) + gvalue*w100, ix+1, iy, iz ); }
         if( ty ){ grid.set( grid.get(ix,iy+1,iz) + gvalue*w010, ix, iy+1, iz ); }
         if( tz ){ grid.set( grid.get(ix,iy,iz+1) + gvalue*w001, ix, iy, iz+1 ); }
         if( tx && ty ){ grid.set( grid.get(ix+1,iy+1,iz) + gvalue*w110, ix+1, iy+1, iz ); }
         if( tx && tz ){ grid.set( grid.get(ix+1,iy,iz+1) + gvalue*w101, ix+1, iy, iz+1 ); }
         if( ty && tz ){ grid.set( grid.get(ix,iy+1,iz+1) + gvalue*w011, ix, iy+1, iz+1 ); }
         if( tx && ty && tz ){ grid.set( grid.get(ix+1,iy+1,iz+1) + gvalue*w111, ix+1, iy+1, iz+1 ); }
      }
   }
}



void lux::StampBlurredWisps( SparseGrid& grid, SparseColorGrid& cgrid, const Vector& P, const float timestep, const Vector& velocity, const Vector& acceleration, const float opacity, const Color& cd, int seed )
{
   float dx = grid.dx();
   dx = ( dx < grid.dy() ) ? dx : grid.dy();
   dx = ( dx < grid.dz() ) ? dx : grid.dz();

   float streaklength  = (velocity + timestep*acceleration).magnitude();
   float nbsamples = streaklength/dx;
   if( nbsamples < 1.0 ){ nbsamples = 1.0; }
   int ix,iy,iz;

   float gvalue = opacity/nbsamples;
   Color cdvalue = cd*opacity/nbsamples;

   float streakTime = 0.0;
   float dt = timestep/nbsamples;
   int samplesTaken = 0;
   while( streakTime <= timestep && samplesTaken < nbsamples )
   {
      float vmag = ( velocity + streakTime*acceleration).magnitude();
      if( vmag != 0 )
      {
         dt = dx / vmag;
      }
      dt = (dt < (timestep-streakTime) )? dt : timestep-streakTime;
      streakTime += dt;
      ++samplesTaken;

      Vector position = P + streakTime * velocity + (0.5*streakTime*streakTime)*acceleration;

      if( grid.getGridIndex( position, ix,iy,iz ) )
      { 
         Vector P000 = grid.evalP(ix,iy,iz);
         position -= P000;

         float wx = position[0]/grid.dx();
         float wy = position[1]/grid.dy();
         float wz = position[2]/grid.dz();

         float vx = 1.0-wx;
         float vy = 1.0-wy;
         float vz = 1.0-wz;

         float w000 = vx*vy*vz;
         float w100 = wx*vy*vz;
         float w010 = vx*wy*vz;
         float w001 = vx*vy*wz;
         float w110 = wx*wy*vz;
         float w011 = vx*wy*wz;
         float w101 = wx*vy*wz;
         float w111 = wx*wy*wz;

         bool tx = ( ix < grid.nx()-1 );
         bool ty = ( iy < grid.ny()-1 );
         bool tz = ( iz < grid.nz()-1 );

         grid.set(grid.get(ix,iy,iz) + gvalue*w000, ix, iy, iz );
         if( tx ){ grid.set( grid.get(ix+1,iy,iz) + gvalue*w100, ix+1, iy, iz ); }
         if( ty ){ grid.set( grid.get(ix,iy+1,iz) + gvalue*w010, ix, iy+1, iz ); }
         if( tz ){ grid.set( grid.get(ix,iy,iz+1) + gvalue*w001, ix, iy, iz+1 ); }
         if( tx && ty ){ grid.set( grid.get(ix+1,iy+1,iz) + gvalue*w110, ix+1, iy+1, iz ); }
         if( tx && tz ){ grid.set( grid.get(ix+1,iy,iz+1) + gvalue*w101, ix+1, iy, iz+1 ); }
         if( ty && tz ){ grid.set( grid.get(ix,iy+1,iz+1) + gvalue*w011, ix, iy+1, iz+1 ); }
         if( tx && ty && tz ){ grid.set( grid.get(ix+1,iy+1,iz+1) + gvalue*w111, ix+1, iy+1, iz+1 ); }


         cgrid.set(cgrid.get(ix,iy,iz) + cdvalue*w000, ix, iy, iz );
         if( tx ){ cgrid.set( cgrid.get(ix+1,iy,iz) + cdvalue*w100, ix+1, iy, iz ); }
         if( ty ){ cgrid.set( cgrid.get(ix,iy+1,iz) + cdvalue*w010, ix, iy+1, iz ); }
         if( tz ){ cgrid.set( cgrid.get(ix,iy,iz+1) + cdvalue*w001, ix, iy, iz+1 ); }
         if( tx && ty ){ cgrid.set( cgrid.get(ix+1,iy+1,iz) + cdvalue*w110, ix+1, iy+1, iz ); }
         if( tx && tz ){ cgrid.set( cgrid.get(ix+1,iy,iz+1) + cdvalue*w101, ix+1, iy, iz+1 ); }
         if( ty && tz ){ cgrid.set( cgrid.get(ix,iy+1,iz+1) + cdvalue*w011, ix, iy+1, iz+1 ); }
         if( tx && ty && tz ){ cgrid.set( cgrid.get(ix+1,iy+1,iz+1) + cdvalue*w111, ix+1, iy+1, iz+1 ); }
      }
   }
}



void lux::StampParticles( VolumeGrid<float>& grid, const ParticleGroupA& particles, const float timestep )
{
   float dx = grid.dx();
   dx = ( dx < grid.dy() ) ? dx : grid.dy();
   dx = ( dx < grid.dz() ) ? dx : grid.dz();


   for( size_t p=0;p<particles.size();p++ )
   {
      const Vector& V = particles[p].v();
      const Vector& A = particles[p].accel();
      const Vector& P = particles[p].P();
      float approxDistance = timestep * V.magnitude() + 0.5*timestep*timestep*A.magnitude();
      int nbsamples = 1 + (int)(  approxDistance/dx  );

      if( nbsamples > 1 )
      {
         Noise_t parms;
         parms.seed = particles[p].id();
         rn.setParameters(parms);
      }

      float dt = timestep/nbsamples;

      int ix,iy,iz;

      float gvalue = particles[p].opacity()/nbsamples;
      for( int q=0;q<nbsamples;q++ )
      {
         float local = q * dt;
         Vector position = P + local*V + (0.5*local*local)*A;
         if( grid.getGridIndex( position, ix,iy,iz ) )
         { 
         Vector P000 = grid.evalP(ix,iy,iz);
         position -= P000;

         float wx = position[0]/grid.dx();
         float wy = position[1]/grid.dy();
         float wz = position[2]/grid.dz();

         float vx = 1.0-wx;
         float vy = 1.0-wy;
         float vz = 1.0-wz;

         float w000 = vx*vy*vz;
         float w100 = wx*vy*vz;
         float w010 = vx*wy*vz;
         float w001 = vx*vy*wz;
         float w110 = wx*wy*vz;
         float w011 = vx*wy*wz;
         float w101 = wx*vy*wz;
         float w111 = wx*wy*wz;

         bool tx = ( ix < grid.nx()-1 );
         bool ty = ( iy < grid.ny()-1 );
         bool tz = ( iz < grid.nz()-1 );

         grid.value(ix,iy,iz) += gvalue*w000;
         if( tx ){ grid.value(ix+1,iy,iz) += gvalue*w100; }
         if( ty ){ grid.value(ix,iy+1,iz) += gvalue*w010; }
         if( tz ){ grid.value(ix,iy,iz+1) += gvalue*w001; }
         if( tx && ty ){ grid.value(ix+1,iy+1,iz) += gvalue*w110; }
         if( tx && tz ){ grid.value(ix+1,iy,iz+1) += gvalue*w101; }
         if( ty && tz ){ grid.value(ix,iy+1,iz+1) += gvalue*w011; }
         if( tx && ty && tz ){ grid.value(ix+1,iy+1,iz+1) += gvalue*w111; }
         }
      }
   }
}



void lux::StampParticles( VolumeGrid<float>& grid, VolumeGrid<Color>& cgrid, const ParticleGroupA& particles, const float timestep )
{
   float dx = grid.dx();
   dx = ( dx < grid.dy() ) ? dx : grid.dy();
   dx = ( dx < grid.dz() ) ? dx : grid.dz();


   for( size_t p=0;p<particles.size();p++ )
   {
      const Vector& V = particles[p].v();
      const Vector& A = particles[p].accel();
      const Vector& P = particles[p].P();
      float approxDistance = timestep * V.magnitude() + 0.5*timestep*timestep*A.magnitude();
      int nbsamples = 1 + (int)(  approxDistance/dx  );

      float dt = timestep/nbsamples;

      if( nbsamples > 1 )
      {
         Noise_t parms;
         parms.seed = particles[p].id();
         rn.setParameters(parms);
      }

      int ix,iy,iz;

      float gvalue = particles[p].opacity()/nbsamples;
      Color cdvalue = particles[p].Cd() * particles[p].opacity()/nbsamples;
      for( int q=0;q<nbsamples;q++ )
      {
         float local = q * dt;
         Vector position = P + local*V + (0.5*local*local)*A;
         if( grid.getGridIndex( position, ix,iy,iz ) )
         { 
         Vector P000 = grid.evalP(ix,iy,iz);
         position -= P000;

         float wx = position[0]/grid.dx();
         float wy = position[1]/grid.dy();
         float wz = position[2]/grid.dz();

         float vx = 1.0-wx;
         float vy = 1.0-wy;
         float vz = 1.0-wz;

         float w000 = vx*vy*vz;
         float w100 = wx*vy*vz;
         float w010 = vx*wy*vz;
         float w001 = vx*vy*wz;
         float w110 = wx*wy*vz;
         float w011 = vx*wy*wz;
         float w101 = wx*vy*wz;
         float w111 = wx*wy*wz;

         bool tx = ( ix < grid.nx()-1 );
         bool ty = ( iy < grid.ny()-1 );
         bool tz = ( iz < grid.nz()-1 );

         grid.value(ix,iy,iz) += gvalue*w000; cgrid.value(ix,iy,iz) += cdvalue*w000;
         if( tx ){ grid.value(ix+1,iy,iz) += gvalue*w100; cgrid.value(ix+1,iy,iz) += cdvalue*w100; }
         if( ty ){ grid.value(ix,iy+1,iz) += gvalue*w010; cgrid.value(ix,iy+1,iz) += cdvalue*w010; }
         if( tz ){ grid.value(ix,iy,iz+1) += gvalue*w001; cgrid.value(ix,iy,iz+1) += cdvalue*w001; }
         if( tx && ty ){ grid.value(ix+1,iy+1,iz) += gvalue*w110; cgrid.value(ix+1,iy+1,iz) += cdvalue*w110; }
         if( tx && tz ){ grid.value(ix+1,iy,iz+1) += gvalue*w101; cgrid.value(ix+1,iy,iz+1) += cdvalue*w101; }
         if( ty && tz ){ grid.value(ix,iy+1,iz+1) += gvalue*w011; cgrid.value(ix,iy+1,iz+1) += cdvalue*w011; }
         if( tx && ty && tz ){ grid.value(ix+1,iy+1,iz+1) += gvalue*w111; cgrid.value(ix+1,iy+1,iz+1) += cdvalue*w111; }
         }
      }
   }
}





Volume<Vector>* lux::StampSELMA( VolumeGrid<Vector>& Xgrid, Volume<Vector>* X, Volume<Vector>* velocity, float dt )
{
   Volume<Vector>* Xadvected = new AdvectVectorVolume( X, velocity, dt );
   Xadvected = new SubtractVectorVolume( Xadvected, new IdentityVectorVolume() );
   lux::Sample( &Xgrid, Xadvected );
   return new AddVectorVolume( new IdentityVectorVolume(), new GriddedVectorVolume( &Xgrid ) );
}

/*
Volume<Vector>* lux::StampBFECCSELMA( VolumeGrid<Vector>& Xgrid, Volume<Vector>* X, Volume<Vector>* velocity, float dt )
{
   Volume<Vector>* Xadvected = new BFECCAdvectVectorVolume( X, velocity, dt );
   Xadvected = new SubtractVectorVolume( Xadvected, new IdentityVectorVolume() );
   lux::Sample( &Xgrid, Xadvected );
   return new AddVectorVolume( new IdentityVectorVolume(), new GriddedVectorVolume( &Xgrid ) );
}
*/

void lux::StampGridPattern( SparseGrid& grid, SparseColorGrid& cgrid, int spacing, const Color& col )
{
   int nx = grid.nx();
   int ny = grid.ny();
   int nz = grid.nz();
   long nb =  nx*ny*nz*3/(spacing*spacing) + 2*nx*ny/spacing + 2*nx*nz/spacing + 2*ny*nz/spacing;
   ProgressMeter meter( nb, "Stamping Gridlines" );

   // Laydown x-axis grid lines
   for( int k=0;k<nz;k+=spacing)
   {
      for( int j=0;j<ny;j+=spacing)
      {
         for( int i=0;i<nx;i++)
	 {
	    grid.set( 1.0, i,j,k );
	    cgrid.set( col, i,j,k );
	    meter.update();
	 }
      }
   }

   // Laydown y-axis grid lines
   for( int k=0;k<nz;k+=spacing)
   {
      for( int i=0;i<nx;i+=spacing)
      {
         for( int j=0;j<ny;j++)
	 {
	    grid.set( 1.0, i,j,k );
	    cgrid.set( col, i,j,k );
	    meter.update();
	 }
      }
   }

   // Laydown z-axis grid lines
   for( int i=0;i<nx;i+=spacing)
   {
      for( int j=0;j<ny;j+=spacing)
      {
         for( int k=0;k<nz;k++)
	 {
	    grid.set( 1.0, i,j,k );
	    cgrid.set( col, i,j,k );
	    meter.update();
	 }
      }
   }

   //for( int k=0;k<nz;k+=spacing)
   //{
      for( int j=0;j<ny;j+=spacing)
      {
         for( int i=0;i<nx;i++)
	 {
	    grid.set( 1.0, i,j,nz-1 );
	    cgrid.set( col, i,j,nz-1 );
	    meter.update();
	 }
      }
   //}


   for( int k=0;k<nz;k+=spacing)
   {
      //for( int j=0;j<ny;j+=spacing)
      //{
         for( int i=0;i<nx;i++)
	 {
	    grid.set( 1.0, i,ny-1,k );
	    cgrid.set( col, i,ny-1,k );
	    meter.update();
	 }
      //}
   }

   //for( int k=0;k<nz;k+=spacing)
   //{
      for( int i=0;i<nx;i+=spacing)
      {
         for( int j=0;j<ny;j++)
	 {
	    grid.set( 1.0, i,j,nz-1 );
	    cgrid.set( col, i,j,nz-1 );
	    meter.update();
	 }
      }
   //}

   for( int k=0;k<nz;k+=spacing)
   {
      //for( int i=0;i<nx;i+=spacing)
      //{
         for( int j=0;j<ny;j++)
	 {
	    grid.set( 1.0, nx-1,j,k );
	    cgrid.set( col, nx-1,j,k );
	    meter.update();
	 }
      //}
   }

   //for( int i=0;i<nx;i+=spacing)
   //{
      for( int j=0;j<ny;j+=spacing)
      {
         for( int k=0;k<nz;k++)
	 {
	    grid.set( 1.0, nx-1,j,k );
	    cgrid.set( col, nx-1,j,k );
	    meter.update();
	 }
      }
   //}

   for( int i=0;i<nx;i+=spacing)
   {
      //for( int j=0;j<ny;j+=spacing)
      //{
         for( int k=0;k<nz;k++)
	 {
	    grid.set( 1.0, i,ny-1,k );
	    cgrid.set( col, i,ny-1,k );
	    meter.update();
	 }
      //}
   }

}

void lux::StampGridPattern( SparseGrid& grid, int spacing )
{
   int nx = grid.nx();
   int ny = grid.ny();
   int nz = grid.nz();

   ProgressMeter meter( nx*ny*nz*3/(spacing*spacing), "Stamping Gridlines" );

   // Laydown x-axis grid lines
   for( int k=0;k<nz;k+=spacing)
   {
      for( int j=0;j<ny;j+=spacing)
      {
         for( int i=0;i<nx;i++)
	 {
	    grid.set( 1.0, i,j,k );
	 }
      }
   }

   // Laydown y-axis grid lines
   for( int k=0;k<nz;k+=spacing)
   {
      for( int i=0;i<nx;i+=spacing)
      {
         for( int j=0;j<ny;j++)
	 {
	    grid.set( 1.0, i,j,k );
	 }
      }
   }

   // Laydown z-axis grid lines
   for( int i=0;i<nx;i+=spacing)
   {
      for( int j=0;j<ny;j+=spacing)
      {
         for( int k=0;k<nz;k++)
	 {
	    grid.set( 1.0, i,j,k );
	 }
      }
   }


   //for( int k=0;k<nz;k+=spacing)
   //{
      for( int j=0;j<ny;j+=spacing)
      {
         for( int i=0;i<nx;i++)
	 {
	    grid.set( 1.0, i,j,nz-1 );
	 }
      }
   //}


   for( int k=0;k<nz;k+=spacing)
   {
      //for( int j=0;j<ny;j+=spacing)
      //{
         for( int i=0;i<nx;i++)
	 {
	    grid.set( 1.0, i,ny-1,k );
	 }
      //}
   }

   //for( int k=0;k<nz;k+=spacing)
   //{
      for( int i=0;i<nx;i+=spacing)
      {
         for( int j=0;j<ny;j++)
	 {
	    grid.set( 1.0, i,j,nz-1 );
	 }
      }
   //}

   for( int k=0;k<nz;k+=spacing)
   {
      //for( int i=0;i<nx;i+=spacing)
      //{
         for( int j=0;j<ny;j++)
	 {
	    grid.set( 1.0, nx-1,j,k );
	 }
      //}
   }

   //for( int i=0;i<nx;i+=spacing)
   //{
      for( int j=0;j<ny;j+=spacing)
      {
         for( int k=0;k<nz;k++)
	 {
	    grid.set( 1.0, nx-1,j,k );
	 }
      }
   //}

   for( int i=0;i<nx;i+=spacing)
   {
      //for( int j=0;j<ny;j+=spacing)
      //{
         for( int k=0;k<nz;k++)
	 {
	    grid.set( 1.0, i,ny-1,k );
	 }
      //}
   }

}





void SplineWispWanderer::reset( const Particle& p0, const Particle& p1 )
{
   part0 = p0;
   part1 = p1;
   parms.seed = part0.id();
   walkNoise.setParameters(parms);
   d1 = Vector( walkNoise.eval()-0.5, walkNoise.eval()-0.5, walkNoise.eval()-0.5 );
   d2 = Vector( walkNoise.eval()-0.5, walkNoise.eval()-0.5, walkNoise.eval()-0.5 );
   walkPosition = Vector( 2.0*walkNoise.eval()-1.0, 2.0*walkNoise.eval()-1.0, walkNoise.eval() );
}

const Vector SplineWispWanderer::step()
{
   float corr = part0.wispCorrelation() + walkPosition[2] * ( part1.wispCorrelation() - part0.wispCorrelation() );
   // do a correlated step
   walkPosition = walkPosition * corr
                + Vector( 2.0*walkNoise.eval()-1.0, 2.0*walkNoise.eval()-1.0, 2.0*walkNoise.eval()-1.0 ) * (1.0-corr);
   // reflection boundary conditions
   if( walkPosition[2] < 0 ){ walkPosition[2] = -walkPosition[2]; }
   if( walkPosition[2] > 1 ){ walkPosition[2] = walkPosition[2]-int(walkPosition[2]); }

   // Interpolate parameters.
   Interpolate( part0, part1, walkPosition[2], interpolated );
   // set parm values
   
   parms.octaves = interpolated.octaves();
   parms.frequency = interpolated.freq();
   parms.translate = interpolated.translate();
   parms.offset = interpolated.offset();
   parms.fjump = interpolated.fjump();
   parms.roughness = interpolated.roughness();
   displacementNoise.setParameters( parms );

   parms.octaves = interpolated.wispOctaves();
   parms.frequency = interpolated.wispFreq();
   parms.translate = interpolated.wispTranslate();
   parms.offset = interpolated.wispOffset();
   parms.fjump = interpolated.wispFjump();
   parms.roughness = interpolated.wispRoughness();
   transformNoise.setParameters( parms );

   // find the unit cylinder
   float radius = pow( fabs(transformNoise.eval( walkPosition ) ), interpolated.wispRadialGroup()  ) ;
   float factor = radius/sqrt( walkPosition[0]*walkPosition[0] + walkPosition[1]*walkPosition[1]  );
   Vector wP = walkPosition;
   wP[0] *= factor;
   wP[1] *= factor;
   wP = wP[0] * interpolated.right() + wP[1] * interpolated.up();
   wP *= interpolated.pscale();
   //wP = wP[2] * interpolated.normal() * (part0.P()-part1.P()).magnitude() + (wP[0] * interpolated.right() + wP[1] * interpolated.up())*interpolated.pscale();

   Vector displacement( displacementNoise.eval(wP), 
                        displacementNoise.eval(wP+d1),  
			displacementNoise.eval(wP+d2) );
   displacement *= interpolated.wispDisplacementScale();
   wP += interpolated.P();
   wP += displacement;

   return wP;
}

void lux::Interpolate( const Particle& p0, const Particle& p1, const float z, Particle& interp )
{
    interp.P() = p0.P() + z*( p1.P() - p0.P() );
    interp.Cd() = p0.Cd() + z*( p1.Cd() - p0.Cd() );
    interp.id() = p0.id();
    interp.pscale() = p0.pscale() + z*( p1.pscale() - p0.pscale() );
    interp.octaves() = p0.octaves() + z*( p1.octaves() - p0.octaves() );
    interp.roughness() = p0.roughness() + z*( p1.roughness() - p0.roughness() );
    interp.freq() = p0.freq() + z*( p1.freq() - p0.freq() );
    interp.fjump() = p0.fjump() + z*( p1.fjump() - p0.fjump() );
    interp.offset() = p0.offset() + z*( p1.offset() - p0.offset() );
    interp.translate() = p0.translate() + z*( p1.translate() - p0.translate() );
    interp.fade() = p0.fade() + z*( p1.fade() - p0.fade() );
    interp.opacity() = p0.opacity() + z*( p1.opacity() - p0.opacity() );

    interp.v() = p0.v() + z*( p1.v() - p0.v() );
    interp.accel() = p0.accel() + z*( p1.accel() - p0.accel() );
    interp.shutter() = p0.shutter() + z*( p1.shutter() - p0.shutter() );
    interp.framerate() = p0.framerate() + z*( p1.framerate() - p0.framerate() );
    interp.lifetime() = p0.lifetime() + z*( p1.lifetime() - p0.lifetime() );
    interp.age() = p0.age() + z*( p1.age() - p0.age() );

    // additional parameters for point wisps
    interp.nbWisps() = p0.nbWisps();
    interp.wispOctaves() = p0.wispOctaves() + z*( p1.wispOctaves() - p0.wispOctaves() );
    interp.wispRoughness() = p0.wispRoughness() + z*( p1.wispRoughness() - p0.wispRoughness() );
    interp.wispFreq() = p0.wispFreq() + z*( p1.wispFreq() - p0.wispFreq() );
    interp.wispFjump() = p0.wispFjump() + z*( p1.wispFjump() - p0.wispFjump() );
    interp.wispOffset() = p0.wispOffset() + z*( p1.wispOffset() - p0.wispOffset() );
    interp.wispTranslate() = p0.wispTranslate() + z*( p1.wispTranslate() - p0.wispTranslate() );
    interp.wispCorrelation() = p0.wispCorrelation() + z*( p1.wispCorrelation() - p0.wispCorrelation() );
    interp.wispRadialGroup() = p0.wispRadialGroup() + z*( p1.wispRadialGroup() - p0.wispRadialGroup() );
    interp.wispDisplacementScale() = p0.wispDisplacementScale() + z*( p1.wispDisplacementScale() - p0.wispDisplacementScale() );

    //additional parameters for pyroclastics
    interp.pyroAmplitude() = p0.pyroAmplitude() + z*( p1.pyroAmplitude() - p0.pyroAmplitude() );
    interp.pyroGamma() = p0.pyroGamma() + z*( p1.pyroGamma() - p0.pyroGamma() );
    interp.pyroDensity() = p0.pyroDensity() + z*( p1.pyroDensity() - p0.pyroDensity() );

    // rotate axes
    Vector axis = (p0.normal()^p1.normal());
    if( axis.magnitude() > 0 ){ axis.normalize(); }
    float cosangle = (p0.normal().unitvector() * p1.normal().unitvector());
    float angle = acos( cosangle ) * z;
    cosangle = cos(angle);
    float sinangle = sin(angle);
    interp.normal() = p0.normal() * cosangle + axis * (axis*p0.normal())*(1.0-cosangle)  + (axis^p0.normal()) * sinangle;

    axis = (p0.right()^p1.right());
    if( axis.magnitude() > 0 ){ axis.normalize(); }
    cosangle = (p0.right().unitvector() * p1.right().unitvector());
    angle = acos( cosangle ) * z;
    cosangle = cos(angle);
    sinangle = sin(angle);
    interp.right() = p0.right() * cosangle + axis * (axis*p0.right())*(1.0-cosangle)  + (axis^p0.right()) * sinangle;


    axis = (p0.up()^p1.up());
    if( axis.magnitude() > 0 ){ axis.normalize(); }
    cosangle = (p0.up().unitvector() * p1.up().unitvector());
    angle = acos( cosangle ) * z;
    cosangle = cos(angle);
    sinangle = sin(angle);
    interp.up() = p0.up() * cosangle + axis * (axis*p0.up())*(1.0-cosangle)  + (axis^p0.up()) * sinangle;


}



void lux::StampSplineWisps( ScalarGrid& grid, ColorGrid& cgrid, const ParticleGroupA& particles )
{
   SplineWispWanderer pww;
   int metertotal = 0;

   for( size_t i=0;i<particles.size()-1;i++ )
   {
      metertotal += particles[i].nbWisps();
   }

   ProgressMeter meter( metertotal, "StampSplineWisps");
   long count = 0;
   for( size_t p=0;p<particles.size()-1;p++ )
   {
      pww.reset(particles[p], particles[p+1]);
      int nb = particles[p].nbWisps();

         for( int i=0;i<nb;i++ )
         {
            const Vector P = pww.step();
	    const Particle& pinterp = pww.interpolatedParticle();
            StampBlurredWisps( grid, cgrid, P,  pinterp.shutter()/24.0, pinterp.v(), pinterp.accel(), pinterp.opacity(), pinterp.Cd(), pinterp.id() );
            meter.update();
	    ++count;
         }
   }
}





void lux::StampNoise( ScalarGrid& grid, Noise* noise, const Vector& center, const float radius, const float fade )
{
   Vector llc = center - Vector(radius,radius,radius);
   Vector urc = center + Vector(radius,radius,radius);
   int il,iu,jl,ju,kl,ku;
   if( grid->getBox( llc, urc, il,iu,jl,ju,kl,ku) )
   {


   for( int k=kl;k<=ku;k++ )
   {
      for( int j=jl;j<=ju;j++ )
      {
         for( int i=il;i<=iu;i++ )
         {
	        Vector P = grid->evalP(i,j,k);
	        float distance = (P-center).magnitude();
	        if( distance <= radius )
	        {
               float fadefactor = std::pow((1.0-distance/radius), fade);
	           if( fade <= 0.0 ) { fadefactor = 1.0; }
	           if( fadefactor > 1.0 ){ fadefactor = 1.0; }
	           if( fadefactor < 0.0 ){ fadefactor = 0.0; }
	           float nv = noise->eval( P ) * fadefactor; 
	           if( nv > grid->get(i,j,k) ) { grid->set(i,j,k, nv); } 
	        }
         }
      }
   }
   }
}


// Motion blur version
// For motion blur, it has to assume that it is summing values
void lux::StampNoise( ScalarGrid& grid, Noise* noise, const Vector& center, const float radius, const Vector& vel, const Vector& accel, const float timestep, const float fade, const int seed )
{
   int nx = grid->nx();
   int ny = grid->ny();
   int nz = grid->nz();

   float dx = grid->dx();
   dx = ( dx < grid->dy() ) ? dx : grid->dy();
   dx = ( dx < grid->dz() ) ? dx : grid->dz();

   float approxDistance = timestep * vel.magnitude() + 0.5*timestep*timestep*accel.magnitude();
   int nbsamples = 1 + (int)(  approxDistance/dx  );



   UniformPRN prn;
   Noise_t parms;
   parms.seed = seed + 646383;
   prn.setParameters(parms);

   int ii,jj,kk;

   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
	    Vector P = grid->evalP(i,j,k);
	    float distance = (P-center).magnitude();
	    if( distance <= radius )
	    {
           float fadefactor = (1.0-distance/radius)/fade;
	       if( fade <= 0.0 ) { fadefactor = 1.0; }
	       if( fadefactor > 1.0 ){ fadefactor = 1.0; }
	       if( fadefactor < 0.0 ){ fadefactor = 0.0; }
	       float nv = noise->eval( P ) * fadefactor / nbsamples; 

	       // accumulate along the streak
	       if( nv > 0 )
	       {
               for( int q=0;q<nbsamples;q++ )
	       {
	          float local = prn.eval() * timestep;
		      Vector PP = P + local*vel + (0.5*local*local)*accel;
                  if( grid->getGridIndex( PP, ii,jj,kk ) ) { grid->set(ii,jj,kk, grid->get(ii,jj,kk) + nv); }
	       }
	       }
	    }
         }
      }
   }
}



void lux::StampNoise( ScalarGrid& grid, const AnchorChain& particles )
{
   Noise* noise = new FractalSum<PerlinNoiseGustavson>();
   ProgressMeter meter( particles.size(), "StampNoise" );
   for( size_t p=0;p<particles.size();p++ )
   {
      noise->setParameters( particles[p] );
      const Vector& center = particles[p].P;
      const float radius = particles[p].pscale;
      const float fade = particles[p].falloff;
      if( (particles[p].v.magnitude() == 0 && particles[p].A.magnitude() == 0) || particles[p].shutter == 0 )
      {
         StampNoise( grid, noise, center, radius, fade );
      }
      else
      {
         float timestep = particles[p].shutter/24.0;
         StampNoise( grid, noise, center, radius, particles[p].v, particles[p].A, timestep, fade, p);
      }
	  meter.update();
   }
   delete noise;
}




void lux::StampPointWisps( ScalarGrid& grid, const AnchorChain& particles )
{

   PointWispWanderer pww;
   int metertotal = 0;

   for( size_t i=0;i<particles.size();i++ )
   {
      metertotal += particles[i].nbWisps;
   }

   ProgressMeter meter( metertotal, "StampSparsePointWisps");
   for( size_t p=0;p<particles.size();p++ )
   {
      pww.reset(particles[p]);
      int nb = particles[p].nbWisps;
      const Vector& V = particles[p].v;
      const Vector& A = particles[p].A;
      float opacity = particles[p].amplitude;
      float shutter = particles[p].shutter;
      int id = particles[p].seed;
      const float timestep = shutter*particles[p].frameRate;
      for( int i=0;i<nb;i++ )
      {
            const Vector P = pww.step();
            // now stamp a dot into the grid
        if( (particles[p].v.magnitude() == 0 && particles[p].A.magnitude() == 0) || particles[p].shutter == 0 )
        {
	       StampWisps( grid, P, opacity );
        }
        else
        {
	       StampBlurredWisps( grid, P, timestep, V, A, opacity, id );
        }
        meter.update();
      }
   }
}

void lux::StampWisps( ScalarGrid& grid, const Vector& P, const float opacity )
{
   int ix,iy,iz;
   if( grid->getGridIndex( P, ix,iy,iz ) )
   { 
         Vector position = P - grid->evalP(ix,iy,iz);
         float wx = position[0]/grid->dx();
         float wy = position[1]/grid->dy();
         float wz = position[2]/grid->dz();

         float vx = 1.0-wx;
         float vy = 1.0-wy;
         float vz = 1.0-wz;

         float w000 = vx*vy*vz;
         float w100 = wx*vy*vz;
         float w010 = vx*wy*vz;
         float w001 = vx*vy*wz;
         float w110 = wx*wy*vz;
         float w011 = vx*wy*wz;
         float w101 = wx*vy*wz;
         float w111 = wx*wy*wz;

         bool tx = ( ix < grid->nx()-1 );
         bool ty = ( iy < grid->ny()-1 );
         bool tz = ( iz < grid->nz()-1 );

         grid->set(ix,iy,iz, grid->get(ix,iy,iz) + opacity*w000 );
         if( tx ){ grid->set( ix+1,iy,iz, grid->get(ix+1,iy,iz) + opacity*w100 ); }
         if( ty ){ grid->set( ix,iy+1,iz, grid->get(ix,iy+1,iz) + opacity*w010 ); }
         if( tz ){ grid->set( ix,iy,iz+1, grid->get(ix,iy,iz+1) + opacity*w001 ); }
         if( tx && ty ){ grid->set( ix+1,iy+1,iz, grid->get(ix+1,iy+1,iz) + opacity*w110 ); }
         if( tx && tz ){ grid->set( ix+1,iy,iz+1, grid->get(ix+1,iy,iz+1) + opacity*w101 ); }
         if( ty && tz ){ grid->set( ix,iy+1,iz+1, grid->get(ix,iy+1,iz+1) + opacity*w011 ); }
         if( tx && ty && tz ){ grid->set( ix+1,iy+1,iz+1, grid->get(ix+1,iy+1,iz+1) + opacity*w111 ); }
   }
}


void lux::StampWisps( ScalarGrid& grid, ColorGrid& cgrid, const Vector& P, const float opacity, const Color& Cd )
{
   int ix,iy,iz;
   if( grid->getGridIndex( P, ix,iy,iz ) )
   { 
         Vector position = P - grid->evalP(ix,iy,iz);
         float wx = position[0]/grid->dx();
         float wy = position[1]/grid->dy();
         float wz = position[2]/grid->dz();

         float vx = 1.0-wx;
         float vy = 1.0-wy;
         float vz = 1.0-wz;

         float w000 = vx*vy*vz;
         float w100 = wx*vy*vz;
         float w010 = vx*wy*vz;
         float w001 = vx*vy*wz;
         float w110 = wx*wy*vz;
         float w011 = vx*wy*wz;
         float w101 = wx*vy*wz;
         float w111 = wx*wy*wz;

         bool tx = ( ix < grid->nx()-1 );
         bool ty = ( iy < grid->ny()-1 );
         bool tz = ( iz < grid->nz()-1 );

         grid->set(ix,iy,iz, grid->get(ix,iy,iz) + opacity*w000 );
         if( tx ){ grid->set( ix+1,iy,iz, grid->get(ix+1,iy,iz) + opacity*w100 ); }
         if( ty ){ grid->set( ix,iy+1,iz, grid->get(ix,iy+1,iz) + opacity*w010 ); }
         if( tz ){ grid->set( ix,iy,iz+1, grid->get(ix,iy,iz+1) + opacity*w001 ); }
         if( tx && ty ){ grid->set( ix+1,iy+1,iz, grid->get(ix+1,iy+1,iz) + opacity*w110 ); }
         if( tx && tz ){ grid->set( ix+1,iy,iz+1, grid->get(ix+1,iy,iz+1) + opacity*w101 ); }
         if( ty && tz ){ grid->set( ix,iy+1,iz+1, grid->get(ix,iy+1,iz+1) + opacity*w011 ); }
         if( tx && ty && tz ){ grid->set( ix+1,iy+1,iz+1, grid->get(ix+1,iy+1,iz+1) + opacity*w111 ); }


         cgrid->set(ix,iy,iz, cgrid->get(ix,iy,iz) + Cd*w000 );
         if( tx ){ cgrid->set( ix+1,iy,iz, cgrid->get(ix+1,iy,iz) + Cd*w100 ); }
         if( ty ){ cgrid->set( ix,iy+1,iz, cgrid->get(ix,iy+1,iz) + Cd*w010 ); }
         if( tz ){ cgrid->set( ix,iy,iz+1, cgrid->get(ix,iy,iz+1) + Cd*w001 ); }
         if( tx && ty ){ cgrid->set( ix+1,iy+1,iz, cgrid->get(ix+1,iy+1,iz) + Cd*w110 ); }
         if( tx && tz ){ cgrid->set( ix+1,iy,iz+1, cgrid->get(ix+1,iy,iz+1) + Cd*w101 ); }
         if( ty && tz ){ cgrid->set( ix,iy+1,iz+1, cgrid->get(ix,iy+1,iz+1) + Cd*w011 ); }
         if( tx && ty && tz ){ cgrid->set( ix+1,iy+1,iz+1, cgrid->get(ix+1,iy+1,iz+1) + Cd*w111 ); }
   }
}



void lux::StampBlurredWisps( ScalarGrid& grid, const Vector& P, const float timestep, const Vector& velocity, const Vector& acceleration, const float opacity, int seed )
{
   float dx = grid->dx();
   dx = ( dx < grid->dy() ) ? dx : grid->dy();
   dx = ( dx < grid->dz() ) ? dx : grid->dz();

   float approxDistance = timestep * velocity.magnitude() + 0.5*timestep*timestep*acceleration.magnitude();
   int nbsamples = 1 + (int)(  approxDistance/dx  );

   if( nbsamples > 1 )
   {
      Noise_t parms;
      parms.seed = seed;
      rn.setParameters(parms);
   }

   int ix,iy,iz;

   float gvalue = opacity/nbsamples;
   for( int q=0;q<nbsamples;q++ )
   {
      float local = rn.eval() * timestep;
      if( q==0 ){ local = 0; }

      Vector position = P + local * velocity + (0.5*local*local)*acceleration;

      if( grid->getGridIndex( position, ix,iy,iz ) )
      { 
         Vector P000 = grid->evalP(ix,iy,iz);
         position -= P000;

         float wx = position[0]/grid->dx();
         float wy = position[1]/grid->dy();
         float wz = position[2]/grid->dz();

         float vx = 1.0-wx;
         float vy = 1.0-wy;
         float vz = 1.0-wz;

         float w000 = vx*vy*vz;
         float w100 = wx*vy*vz;
         float w010 = vx*wy*vz;
         float w001 = vx*vy*wz;
         float w110 = wx*wy*vz;
         float w011 = vx*wy*wz;
         float w101 = wx*vy*wz;
         float w111 = wx*wy*wz;

         bool tx = ( ix < grid->nx()-1 );
         bool ty = ( iy < grid->ny()-1 );
         bool tz = ( iz < grid->nz()-1 );

         grid->set(ix,iy,iz, grid->get(ix,iy,iz) + gvalue*w000 );
         if( tx ){ grid->set( ix+1,iy,iz, grid->get(ix+1,iy,iz) + gvalue*w100 ); }
         if( ty ){ grid->set( ix,iy+1,iz, grid->get(ix,iy+1,iz) + gvalue*w010 ); }
         if( tz ){ grid->set( ix,iy,iz+1, grid->get(ix,iy,iz+1) + gvalue*w001 ); }
         if( tx && ty ){ grid->set( ix+1,iy+1,iz, grid->get(ix+1,iy+1,iz) + gvalue*w110 ); }
         if( tx && tz ){ grid->set( ix+1,iy,iz+1, grid->get(ix+1,iy,iz+1) + gvalue*w101 ); }
         if( ty && tz ){ grid->set( ix,iy+1,iz+1, grid->get(ix,iy+1,iz+1) + gvalue*w011 ); }
         if( tx && ty && tz ){ grid->set( ix+1,iy+1,iz+1, grid->get(ix+1,iy+1,iz+1) + gvalue*w111 ); }
      }
   }
}




void lux::StampBlurredWisps( ScalarGrid& grid, ColorGrid& cgrid, const Vector& P, const float timestep, const Vector& velocity, const Vector& acceleration, const float opacity, const Color& Cd, int seed )
{
   float dx = grid->dx();
   dx = ( dx < grid->dy() ) ? dx : grid->dy();
   dx = ( dx < grid->dz() ) ? dx : grid->dz();

   float approxDistance = timestep * velocity.magnitude() + 0.5*timestep*timestep*acceleration.magnitude();
   int nbsamples = 1 + (int)(  approxDistance/dx  );

   if( nbsamples > 1 )
   {
      Noise_t parms;
      parms.seed = seed;
      rn.setParameters(parms);
   }

   int ix,iy,iz;

   float gvalue = opacity/nbsamples;
   Color Cdvalue = Cd/nbsamples;
   for( int q=0;q<nbsamples;q++ )
   {
      float local = rn.eval() * timestep;
      if( q==0 ){ local = 0; }

      Vector position = P + local * velocity + (0.5*local*local)*acceleration;

      if( grid->getGridIndex( position, ix,iy,iz ) )
      { 
         Vector P000 = grid->evalP(ix,iy,iz);
         position -= P000;

         float wx = position[0]/grid->dx();
         float wy = position[1]/grid->dy();
         float wz = position[2]/grid->dz();

         float vx = 1.0-wx;
         float vy = 1.0-wy;
         float vz = 1.0-wz;

         float w000 = vx*vy*vz;
         float w100 = wx*vy*vz;
         float w010 = vx*wy*vz;
         float w001 = vx*vy*wz;
         float w110 = wx*wy*vz;
         float w011 = vx*wy*wz;
         float w101 = wx*vy*wz;
         float w111 = wx*wy*wz;

         bool tx = ( ix < grid->nx()-1 );
         bool ty = ( iy < grid->ny()-1 );
         bool tz = ( iz < grid->nz()-1 );

         grid->set(ix,iy,iz, grid->get(ix,iy,iz) + gvalue*w000 );
         if( tx ){ grid->set( ix+1,iy,iz, grid->get(ix+1,iy,iz) + gvalue*w100 ); }
         if( ty ){ grid->set( ix,iy+1,iz, grid->get(ix,iy+1,iz) + gvalue*w010 ); }
         if( tz ){ grid->set( ix,iy,iz+1, grid->get(ix,iy,iz+1) + gvalue*w001 ); }
         if( tx && ty ){ grid->set( ix+1,iy+1,iz, grid->get(ix+1,iy+1,iz) + gvalue*w110 ); }
         if( tx && tz ){ grid->set( ix+1,iy,iz+1, grid->get(ix+1,iy,iz+1) + gvalue*w101 ); }
         if( ty && tz ){ grid->set( ix,iy+1,iz+1, grid->get(ix,iy+1,iz+1) + gvalue*w011 ); }
         if( tx && ty && tz ){ grid->set( ix+1,iy+1,iz+1, grid->get(ix+1,iy+1,iz+1) + gvalue*w111 ); }


         cgrid->set(ix,iy,iz, cgrid->get(ix,iy,iz) + Cdvalue*w000 );
         if( tx ){ cgrid->set( ix+1,iy,iz, cgrid->get(ix+1,iy,iz) + Cdvalue*w100 ); }
         if( ty ){ cgrid->set( ix,iy+1,iz, cgrid->get(ix,iy+1,iz) + Cdvalue*w010 ); }
         if( tz ){ cgrid->set( ix,iy,iz+1, cgrid->get(ix,iy,iz+1) + Cdvalue*w001 ); }
         if( tx && ty ){ cgrid->set( ix+1,iy+1,iz, cgrid->get(ix+1,iy+1,iz) + Cdvalue*w110 ); }
         if( tx && tz ){ cgrid->set( ix+1,iy,iz+1, cgrid->get(ix+1,iy,iz+1) + Cdvalue*w101 ); }
         if( ty && tz ){ cgrid->set( ix,iy+1,iz+1, cgrid->get(ix,iy+1,iz+1) + Cdvalue*w011 ); }
         if( tx && ty && tz ){ cgrid->set( ix+1,iy+1,iz+1, cgrid->get(ix+1,iy+1,iz+1) + Cdvalue*w111 ); }
      }
   }
}






void lux::StampPointWisps( ScalarGrid& grid, ColorGrid& cgrid,  const AnchorChain& particles )
{

   PointWispWanderer pww;
   int metertotal = 0;

   for( size_t i=0;i<particles.size();i++ )
   {
      metertotal += particles[i].nbWisps;
   }

   ProgressMeter meter( metertotal, "StampSparsePointWisps");
   for( size_t p=0;p<particles.size();p++ )
   {
      pww.reset(particles[p]);
      int nb = particles[p].nbWisps;
      const Vector& V = particles[p].v;
      const Vector& A = particles[p].A;
      float opacity = particles[p].amplitude;
      float shutter = particles[p].shutter;
      int id = particles[p].seed;
      Color Cd = particles[p].Cd;
      const float timestep = shutter*particles[p].frameRate;
      for( int i=0;i<nb;i++ )
      {
            const Vector P = pww.step();
            // now stamp a dot into the grid
        if( (particles[p].v.magnitude() == 0 && particles[p].A.magnitude() == 0) || particles[p].shutter == 0 )
        {
	       StampWisps( grid, cgrid, P, opacity, Cd );
        }
        else
        {
	       StampBlurredWisps( grid, cgrid, P, timestep, V, A, opacity, Cd, id );
        }
        meter.update();
      }
   }
}


