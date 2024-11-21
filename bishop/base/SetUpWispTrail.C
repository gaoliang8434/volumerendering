
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
         p.translate() = -(p.P()*u)*u * translatescale ; 
         p.Cd() = Color( 1.0, 1.0, 1.0, 1.0 );
	 p.age() *= timeScale;

         p.nbWisps() = nbwisps;
         p.wispDisplacementScale() = dispscale * (0.2 + f )/1.2;
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


