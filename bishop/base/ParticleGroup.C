
#include "ParticleGroup.h"

using namespace lux;


ParticleGroup::ParticleGroup()
{
   // create standard set of attributes:
   //   id, pos, vel, ci
   int_attributes["id"] = PGAttribute<int>( "id", -1 );
   vector_attributes["pos"] = PGAttribute<Vector>( "pos", Vector(0,0,0) );
   vector_attributes["vel"] = PGAttribute<Vector>( "vel", Vector(0,0,0) );
   vector_attributes["translate"] = PGAttribute<Vector>( "translate", Vector(0,0,0) );
   color_attributes["ci"] = PGAttribute<Color>( "ci", Color(0,0,0,0) );
   float_attributes["pscale"] = PGAttribute<float>( "pscale", 1.0 );
   float_attributes["density"] = PGAttribute<float>( "density", 0.0 );
   float_attributes["roughness"] = PGAttribute<float>( "roughness", 0.5 );
   float_attributes["octaves"] = PGAttribute<float>( "octaves", 1.0 );
   float_attributes["fjump"] = PGAttribute<float>( "fjump", 2.2 );
   float_attributes["frequency"] = PGAttribute<float>( "frequency", 1.0 );
}

ParticleGroup::~ParticleGroup(){}


void ParticleGroup::create_attr( const std::string& nam, const int& def )
{
   if( int_attributes.find(nam) != int_attributes.end() ){ return; }
   int_attributes[nam] = PGAttribute<int>( nam, def );
   int_attributes[nam].expand_to( nb_particles() );
}

void ParticleGroup::create_attr( const std::string& nam, const float& def )
{
   if( float_attributes.find(nam) != float_attributes.end() ){ return; }
   float_attributes[nam] = PGAttribute<float>( nam, def );
   float_attributes[nam].expand_to( nb_particles() );
}

void ParticleGroup::create_attr( const std::string& nam, const Vector& def )
{
   if( vector_attributes.find(nam) != vector_attributes.end() ){ return; }
   vector_attributes[nam] = PGAttribute<Vector>( nam, def );
   vector_attributes[nam].expand_to( nb_particles() );
}

void ParticleGroup::create_attr( const std::string& nam, const Color& def )
{
   if( color_attributes.find(nam) != color_attributes.end() ){ return; }
   color_attributes[nam] = PGAttribute<Color>( nam, def );
   color_attributes[nam].expand_to( nb_particles() );
}


const size_t ParticleGroup::add_particle()
{
   return add_particles(1);
}

const size_t ParticleGroup::add_particles( const size_t nb )
{
   size_t current_size = nb_particles();
   size_t add_size = current_size + nb;

   for( std::map<std::string,PGAttribute<int> >::iterator a = int_attributes.begin(); a != int_attributes.end(); a++ )
   {
      a->second.expand_to(add_size);
   }
   for( std::map<std::string,PGAttribute<float> >::iterator a = float_attributes.begin(); a != float_attributes.end(); a++ )
   {
      a->second.expand_to(add_size);
   }
   for( std::map<std::string,PGAttribute<Vector> >::iterator a = vector_attributes.begin(); a != vector_attributes.end(); a++ )
   {
      a->second.expand_to(add_size);
   }
   for( std::map<std::string,PGAttribute<Color> >::iterator a = color_attributes.begin(); a != color_attributes.end(); a++ )
   {
      a->second.expand_to(add_size);
   }
   return add_size-1; // return the index of the last particle
}

size_t ParticleGroup::nb_particles() const 
{ 
   std::map<std::string,PGAttribute<int> >::const_iterator a = int_attributes.find("id");
   return a->second.size();
}

const int& ParticleGroup::get_int_attr( const std::string& nam, const size_t p ) const
{
   std::map<std::string,PGAttribute<int> >::const_iterator a = int_attributes.find(nam);
   return a->second.get(p);
}

const float& ParticleGroup::get_float_attr( const std::string& nam, const size_t p ) const
{
   std::map<std::string,PGAttribute<float> >::const_iterator a = float_attributes.find(nam);
   return a->second.get(p);
}

const Vector& ParticleGroup::get_vector_attr( const std::string& nam, const size_t p ) const
{
   std::map<std::string,PGAttribute<Vector> >::const_iterator a = vector_attributes.find(nam);
   return a->second.get(p);
}

const Color& ParticleGroup::get_color_attr( const std::string& nam, const size_t p ) const
{
   std::map<std::string,PGAttribute<Color> >::const_iterator a = color_attributes.find(nam);
   return a->second.get(p);
}

    
const int& ParticleGroup::id( const size_t p ) const
{
   std::map<std::string,PGAttribute<int> >::const_iterator a = int_attributes.find("id");
   return a->second.get(p);
}

const Vector& ParticleGroup::pos( const size_t p ) const
{
   std::map<std::string,PGAttribute<Vector> >::const_iterator a = vector_attributes.find("pos");
   return a->second.get(p);
}

const Vector& ParticleGroup::vel( const size_t p ) const
{
   std::map<std::string,PGAttribute<Vector> >::const_iterator a = vector_attributes.find("vel");
   return a->second.get(p);
}

const Color& ParticleGroup::ci( const size_t p ) const
{
   std::map<std::string,PGAttribute<Color> >::const_iterator a = color_attributes.find("ci");
   return a->second.get(p);
}

const float& ParticleGroup::pscale( const size_t p ) const
{
   std::map<std::string,PGAttribute<float> >::const_iterator a = float_attributes.find("pscale");
   return a->second.get(p);
}

const float& ParticleGroup::density( const size_t p ) const
{
   std::map<std::string,PGAttribute<float> >::const_iterator a = float_attributes.find("density");
   return a->second.get(p);
}

const float& ParticleGroup::roughness( const size_t p ) const
{
   std::map<std::string,PGAttribute<float> >::const_iterator a = float_attributes.find("roughness");
   return a->second.get(p);
}

const float& ParticleGroup::octaves( const size_t p ) const
{
   std::map<std::string,PGAttribute<float> >::const_iterator a = float_attributes.find("octaves");
   return a->second.get(p);
}

const float& ParticleGroup::fjump( const size_t p ) const
{
   std::map<std::string,PGAttribute<float> >::const_iterator a = float_attributes.find("fjump");
   return a->second.get(p);
}

const float& ParticleGroup::frequency( const size_t p ) const
{
   std::map<std::string,PGAttribute<float> >::const_iterator a = float_attributes.find("frequency");
   return a->second.get(p);
}

const Vector& ParticleGroup::translate( const size_t p ) const
{
   std::map<std::string,PGAttribute<Vector> >::const_iterator a = vector_attributes.find("translate");
   return a->second.get(p);
}

    
void ParticleGroup::set_attr( const std::string& nam, const size_t p, const int& value )
{
   int_attributes[nam].set(p, value);
}

void ParticleGroup::set_attr( const std::string& nam, const size_t p, const float& value )
{
   float_attributes[nam].set(p, value);
}

void ParticleGroup::set_attr( const std::string& nam, const size_t p, const Vector& value ) 
{
   vector_attributes[nam].set(p, value);
}

void ParticleGroup::set_attr( const std::string& nam, const size_t p, const Color& value ) 
{
   color_attributes[nam].set(p, value);
}


void ParticleGroup::set_id( const size_t p, const int& value )
{
   int_attributes["id"].set(p, value);
}

void ParticleGroup::set_pos( const size_t p, const Vector& value )
{
   vector_attributes["pos"].set(p, value);
}

void ParticleGroup::set_vel( const size_t p, const Vector& value )
{
   vector_attributes["vel"].set(p, value);
}

void ParticleGroup::set_ci( const size_t p, const Color& value )
{
   color_attributes["ci"].set(p, value);
}

void ParticleGroup::set_pscale( const size_t p, const float& value )
{
   float_attributes["pscale"].set(p, value);
}

void ParticleGroup::set_density( const size_t p, const float& value )
{
   float_attributes["density"].set(p, value);
}

void ParticleGroup::set_roughness( const size_t p, const float& value )
{
   float_attributes["roughness"].set(p, value);
}

void ParticleGroup::set_octaves( const size_t p, const float& value )
{
   float_attributes["octaves"].set(p, value);
}

void ParticleGroup::set_fjump( const size_t p, const float& value )
{
   float_attributes["fjump"].set(p, value);
}

void ParticleGroup::set_frequency( const size_t p, const float& value )
{
   float_attributes["frequency"].set(p, value);
}

void ParticleGroup::set_translate( const size_t p, const Vector& value )
{
   vector_attributes["translate"].set(p, value);
}


std::vector<std::string> ParticleGroup::show_int_attrs() const
{
   std::vector<std::string> keys;
   for( std::map<std::string, PGAttribute<int> >::const_iterator i = int_attributes.begin(); i != int_attributes.end(); i++ )
   {
      keys.push_back( i->first ); 
   }
   return keys;
}

std::vector<std::string> ParticleGroup::show_float_attrs() const
{
   std::vector<std::string> keys;
   for( std::map<std::string, PGAttribute<float> >::const_iterator i = float_attributes.begin(); i != float_attributes.end(); i++ )
   {
      keys.push_back( i->first ); 
   }
   return keys;
}


std::vector<std::string> ParticleGroup::show_vector_attrs() const
{
   std::vector<std::string> keys;
   for( std::map<std::string, PGAttribute<Vector> >::const_iterator i = vector_attributes.begin(); i != vector_attributes.end(); i++ )
   {
      keys.push_back( i->first ); 
   }
   return keys;
}

std::vector<std::string> ParticleGroup::show_color_attrs() const
{
   std::vector<std::string> keys;
   for( std::map<std::string, PGAttribute<Color> >::const_iterator i = color_attributes.begin(); i != color_attributes.end(); i++ )
   {
      keys.push_back( i->first ); 
   }
   return keys;
}

std::vector<std::string> ParticleGroup::show_all_attrs() const
{
   std::vector<std::string> keys;
   for( std::map<std::string, PGAttribute<int> >::const_iterator i = int_attributes.begin(); i != int_attributes.end(); i++ )
   {
      keys.push_back( i->first ); 
   }
   for( std::map<std::string, PGAttribute<float> >::const_iterator i = float_attributes.begin(); i != float_attributes.end(); i++ )
   {
      keys.push_back( i->first ); 
   }
   for( std::map<std::string, PGAttribute<Vector> >::const_iterator i = vector_attributes.begin(); i != vector_attributes.end(); i++ )
   {
      keys.push_back( i->first ); 
   }
   for( std::map<std::string, PGAttribute<Color> >::const_iterator i = color_attributes.begin(); i != color_attributes.end(); i++ )
   {
      keys.push_back( i->first ); 
   }
   return keys;
}

bool ParticleGroup::attr_exists( const std::string& nam ) const
{
   if( int_attributes.find(nam) != int_attributes.end() ){ return true; }
   if( float_attributes.find(nam) != float_attributes.end() ){ return true; }
   if( vector_attributes.find(nam) != vector_attributes.end() ){ return true; }
   if( color_attributes.find(nam) != color_attributes.end() ){ return true; }
   return false;
}


bool ParticleGroup::int_attr_exists( const std::string& nam ) const
{
   if( int_attributes.find(nam) != int_attributes.end() ){ return true; }
   return false;
}

bool ParticleGroup::float_attr_exists( const std::string& nam ) const
{
   if( float_attributes.find(nam) != float_attributes.end() ){ return true; }
   return false;
}

bool ParticleGroup::vector_attr_exists( const std::string& nam ) const
{
   if( vector_attributes.find(nam) != vector_attributes.end() ){ return true; }
   return false;
}

bool ParticleGroup::color_attr_exists( const std::string& nam ) const
{
   if( color_attributes.find(nam) != color_attributes.end() ){ return true; }
   return false;
}



void ParticleGroup::merge( const ParticleGroup& g )
{
   size_t next_particle = nb_particles();
   for( std::map< std::string, PGAttribute<int> >::const_iterator a = g.int_attributes.begin(); a != g.int_attributes.end(); a++ )
   {
      create_attr( a->second.attr_name(), a->second.default_value() );
   }
   for( std::map< std::string, PGAttribute<float> >::const_iterator a = g.float_attributes.begin(); a != g.float_attributes.end(); a++ )
   {
      create_attr( a->second.attr_name(), a->second.default_value() );
   }
   for( std::map< std::string, PGAttribute<Vector> >::const_iterator a = g.vector_attributes.begin(); a != g.vector_attributes.end(); a++ )
   {
      create_attr( a->second.attr_name(), a->second.default_value() );
   }
   for( std::map< std::string, PGAttribute<Color> >::const_iterator a = g.color_attributes.begin(); a != g.color_attributes.end(); a++ )
   {
      create_attr( a->second.attr_name(), a->second.default_value() );
   }
   add_particles( g.nb_particles() );
   for( std::map< std::string, PGAttribute<int> >::const_iterator a = g.int_attributes.begin(); a != g.int_attributes.end(); a++ )
   {
      for(size_t i=0;i<g.nb_particles();i++ )
      {
         set_attr( a->first, i+next_particle, g.get_int_attr( a->first, i ) );
      } 
   }
   for( std::map< std::string, PGAttribute<float> >::const_iterator a = g.float_attributes.begin(); a != g.float_attributes.end(); a++ )
   {
      for(size_t i=0;i<g.nb_particles();i++ )
      {
         set_attr( a->first, i+next_particle, g.get_float_attr( a->first, i ) );
      } 
   }
   for( std::map< std::string, PGAttribute<Vector> >::const_iterator a = g.vector_attributes.begin(); a != g.vector_attributes.end(); a++ )
   {
      for(size_t i=0;i<g.nb_particles();i++ )
      {
         set_attr( a->first, i+next_particle, g.get_vector_attr( a->first, i ) );
      } 
   }
   for( std::map< std::string, PGAttribute<Color> >::const_iterator a = g.color_attributes.begin(); a != g.color_attributes.end(); a++ )
   {
      for(size_t i=0;i<g.nb_particles();i++ )
      {
         set_attr( a->first, i+next_particle, g.get_color_attr( a->first, i ) );
      } 
   }
}

char* PointCloud::__str__()
     {
       static char typeLabel[2048];
       std::string lbl = "PointCloud";
       size_t lbllength = lbl.size();
       if( lbllength > 2047 ){ lbllength = 2047; }
       lbllength = lbl.copy( typeLabel, lbllength);
       typeLabel[lbllength] = '\0';
       return typeLabel;
    }


char* PointCloud::__doc__() 
     {
       static char docLabel[2048];
       std::string lbl = "PointCloud is a collection of points with attributes";
       size_t lbllength = lbl.size();
       if( lbllength > 2047 ){ lbllength = 2047; }
       lbllength = lbl.copy( docLabel, lbllength);
       docLabel[lbllength] = '\0';
       return docLabel;
    }



PointCloud::PointCloud() :  std::shared_ptr<ParticleGroup>() {}
PointCloud::PointCloud( ParticleGroup* f ) :  std::shared_ptr<ParticleGroup>( f ) {}
PointCloud::~PointCloud(){}

PointCloud PointCloud::operator+( const PointCloud& e2 )
{
   (*this)->merge( *e2 );
   return *this;
}



