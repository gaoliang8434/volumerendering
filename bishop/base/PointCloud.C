
#include "PointCloud.h"

#include <iostream>
#include <fstream>

using namespace lux;

PointCloud lux::makePointCloud( int nb )
{
   PointCloud pc( new ParticleGroup() );
   pc->add_particles(nb);
   return pc;
}



void lux::create_attr( PointCloud& p, const std::string& nam, const int& def ){ p->create_attr( nam, def ); }
void lux::create_attr( PointCloud& p, const std::string& nam, const float& def ){ p->create_attr( nam, def ); }
void lux::create_attr( PointCloud& p, const std::string& nam, const Vector& def ){ p->create_attr( nam, def ); }
void lux::create_attr( PointCloud& p, const std::string& nam, const Color& def ){ p->create_attr( nam, def ); }

size_t lux::add_particle( PointCloud& p ) { return p->add_particle(); }
size_t lux::add_particles( PointCloud& p, const size_t nb ) { return p->add_particles(nb); }
size_t lux::nb_particles( PointCloud& p ) { return p->nb_particles(); }

int lux::get_int_attr( PointCloud& p, const std::string& nam, const size_t i ) { return p->get_int_attr( nam, i ); }
float lux::get_float_attr( PointCloud& p, const std::string& nam, const size_t i ) { return p->get_float_attr( nam, i ); } 
Vector lux::get_vector_attr( PointCloud& p, const std::string& nam, const size_t i )  { return p->get_vector_attr( nam, i ); }
Color lux::get_color_attr( PointCloud& p, const std::string& nam, const size_t i )  { return p->get_color_attr( nam, i ); }
    
int lux::id(  PointCloud& p, const size_t i ){ return p->id(i); }
Vector lux::pos(  PointCloud& p, const size_t i ){ return p->pos(i); }
Vector lux::vel(  PointCloud& p, const size_t i ){ return p->vel(i); }
Color lux::ci(  PointCloud& p, const size_t i ){ return p->ci(i); }
float lux::pscale(  PointCloud& p, const size_t i ){ return p->pscale(i); }
float lux::density(  PointCloud& p, const size_t i ){ return p->pscale(i); }
float lux::roughness(  PointCloud& p, const size_t i ){ return p->pscale(i); }
float lux::octaves(  PointCloud& p, const size_t i ){ return p->pscale(i); }
float lux::fjump(  PointCloud& p, const size_t i ){ return p->pscale(i); }
float lux::frequency(  PointCloud& p, const size_t i ){ return p->pscale(i); }
Vector lux::translate(  PointCloud& p, const size_t i ){ return p->vel(i); }
    
void lux::set_attr( PointCloud& p, const std::string& nam, const size_t i, const int& value ){ p->set_attr( nam, i, value ); }
void lux::set_attr( PointCloud& p, const std::string& nam, const size_t i, const float& value ) { p->set_attr( nam, i, value ); }
void lux::set_attr( PointCloud& p, const std::string& nam, const size_t i, const Vector& value ) { p->set_attr( nam, i, value ); }
void lux::set_attr( PointCloud& p, const std::string& nam, const size_t i, const Color& value ) { p->set_attr( nam, i, value ); }

void lux::set_id(  PointCloud& p, const size_t i, const int& value ){ p->set_id(i,value); }
void lux::set_pos(  PointCloud& p, const size_t i, const Vector& value ){ p->set_pos(i,value); }
void lux::set_vel(  PointCloud& p, const size_t i, const Vector& value ){ p->set_vel(i,value); }
void lux::set_ci(  PointCloud& p, const size_t i, const Color& value ){ p->set_ci(i,value); }
void lux::set_pscale(  PointCloud& p, const size_t i, const float& value ){ p->set_pscale(i,value); }
void lux::set_density(  PointCloud& p, const size_t i, const float& value ){ p->set_density(i,value); }
void lux::set_roughness(  PointCloud& p, const size_t i, const float& value ){ p->set_roughness(i,value); }
void lux::set_octaves(  PointCloud& p, const size_t i, const float& value ){ p->set_octaves(i,value); }
void lux::set_fjump(  PointCloud& p, const size_t i, const float& value ){ p->set_fjump(i,value); }
void lux::set_frequency(  PointCloud& p, const size_t i, const float& value ){ p->set_frequency(i,value); }
void lux::set_translate(  PointCloud& p, const size_t i, const Vector& value ){ p->set_translate(i,value); }

std::vector<std::string> lux::show_int_attrs( PointCloud& p){ return p->show_int_attrs(); }
std::vector<std::string> lux::show_float_attrs( PointCloud& p){ return p->show_float_attrs(); }
std::vector<std::string> lux::show_vector_attrs( PointCloud& p){ return p->show_vector_attrs(); }
std::vector<std::string> lux::show_color_attrs( PointCloud& p){ return p->show_color_attrs(); }
std::vector<std::string> lux::show_all_attrs( PointCloud& p){ return p->show_all_attrs(); }

int lux::attr_exists(  PointCloud& p, const std::string& nam ){ if( p->attr_exists(nam) ){ return 1; } else { return 0; } }
int lux::int_attr_exists(  PointCloud& p, const std::string& nam ){ if( p->int_attr_exists(nam) ){ return 1; } else { return 0; } }
int lux::float_attr_exists(  PointCloud& p, const std::string& nam ){ if( p->float_attr_exists(nam) ){ return 1; } else { return 0; } }
int lux::vector_attr_exists(  PointCloud& p, const std::string& nam ){ if( p->vector_attr_exists(nam) ){ return 1; } else { return 0; } }
int lux::color_attr_exists(  PointCloud& p, const std::string& nam ){ if( p->color_attr_exists(nam) ){ return 1; } else { return 0; } }

void lux::merge(  PointCloud& p, const PointCloud& g ) { p->merge(*g); }



Vector lux::llc(const PointCloud& p )
{
   Vector value = p->pos(0);
   for( size_t i=1;i<p->nb_particles();i++)
   {
      const Vector& P = p->pos(i);
      if( P.X() < value.X() ){ value[0] = P.X(); }
      if( P.Y() < value.Y() ){ value[1] = P.Y(); }
      if( P.Z() < value.Z() ){ value[2] = P.Z(); }
   }
   return value;
}
Vector lux::urc(const PointCloud& p )
{
   Vector value = p->pos(0);
   for( size_t i=1;i<p->nb_particles();i++)
   {
      const Vector& P = p->pos(i);
      if( P.X() > value.X() ){ value[0] = P.X(); }
      if( P.Y() > value.Y() ){ value[1] = P.Y(); }
      if( P.Z() > value.Z() ){ value[2] = P.Z(); }
   }
   return value;
}





PointCloud lux::makeWispCloud( int nb )
{
   PointCloud pg =  PointCloud(new ParticleGroup());
   std::map<std::string,int> intattrs;
   std::map<std::string,float> floatattrs;
   std::map<std::string,Vector> vectorattrs;
   std::map<std::string,Color> colorattrs;
 
// int attributes
    intattrs["wispchildren"] = 0;

// float attributes
    floatattrs["opacity"] = 1.0;
    floatattrs["pscale"] = 1.0;
    floatattrs["octaves"] = 1.0;
    floatattrs["roughness"] = 0.5;
    floatattrs["frequency"] = 1.0;
    floatattrs["fjump"] = 2.0;
    floatattrs["offset"] = 0.0;
    floatattrs["shutter"] = 0.0;
    floatattrs["framerate"] = 24.0;
    floatattrs["lifetime"] = 1.0;
    floatattrs["age"] = 1.0;
    floatattrs["woctaves"] = 1.0;
    floatattrs["wroughness"] = 0.5;
    floatattrs["wfreq"] = 1.0;
    floatattrs["wfjump"] = 2.0;
    floatattrs["woffset"] = 0.0;
    floatattrs["correlation"] = 0.0;
    floatattrs["wclump"] = 1.0/3.0;
    floatattrs["wscale"] = 1.0;

// Vector attributes
    vectorattrs["translate"]  = Vector(0,0,0);
    vectorattrs["wtranslate"] = Vector(0,0,0);
    vectorattrs["normal"]     = Vector(0,0,1);
    vectorattrs["right"]      = Vector(1,0,0);
    vectorattrs["up"]         = Vector(0,1,0);
    vectorattrs["accel"]      = Vector(0,0,0);

// Color attributes

// Populate point cloud
   for( std::map<std::string,int>::const_iterator i=intattrs.begin(); i != intattrs.end(); i++ )
   {
      pg->create_attr( i->first, i->second );
   }
   for( std::map<std::string,float>::const_iterator i=floatattrs.begin(); i != floatattrs.end(); i++ )
   {
      pg->create_attr( i->first, i->second );
   }
   for( std::map<std::string,Vector>::const_iterator i=vectorattrs.begin(); i != vectorattrs.end(); i++ )
   {
      pg->create_attr( i->first, i->second );
   }
   for( std::map<std::string,Color>::const_iterator i=colorattrs.begin(); i != colorattrs.end(); i++ )
   {
      pg->create_attr( i->first, i->second );
   }
   pg->add_particles(nb);
   return pg;
}







PointCloud lux::makeNoiseCloud( int nb )
{
   PointCloud pg =  PointCloud(new ParticleGroup());
   std::map<std::string,int> intattrs;
   std::map<std::string,float> floatattrs;
   std::map<std::string,Vector> vectorattrs;
   std::map<std::string,Color> colorattrs;


 
// int attributes

//    intattrs["wispchildren"] = 0;


// float attributes

    floatattrs["opacity"] = 1.0;
    floatattrs["pscale"] = 1.0;
    floatattrs["octaves"] = 1.0;
    floatattrs["roughness"] = 0.5;
    floatattrs["freq"] = 1.0;
    floatattrs["fjump"] = 2.0;
    floatattrs["offset"] = 0.0;
    floatattrs["shutter"] = 0.0;
    floatattrs["framerate"] = 24.0;
    floatattrs["lifetime"] = 1.0;
    floatattrs["age"] = 1.0;
    floatattrs["fade"] = 1.0;
    floatattrs["amplitude"] = 1.0;


// Vector attributes

    vectorattrs["translate"]  = Vector(0,0,0);
    vectorattrs["accel"]      = Vector(0,0,0);


// Color attributes



// Populate point cloud
   for( std::map<std::string,int>::const_iterator i=intattrs.begin(); i != intattrs.end(); i++ )
   {
      pg->create_attr( i->first, i->second );
   }
   for( std::map<std::string,float>::const_iterator i=floatattrs.begin(); i != floatattrs.end(); i++ )
   {
      pg->create_attr( i->first, i->second );
   }
   for( std::map<std::string,Vector>::const_iterator i=vectorattrs.begin(); i != vectorattrs.end(); i++ )
   {
      pg->create_attr( i->first, i->second );
   }
   for( std::map<std::string,Color>::const_iterator i=colorattrs.begin(); i != colorattrs.end(); i++ )
   {
      pg->create_attr( i->first, i->second );
   }
   pg->add_particles(nb);
   return pg;
}





PointCloud lux::makePyroCloud( int nb )
{
   PointCloud pg =  PointCloud(new ParticleGroup());
   std::map<std::string,int> intattrs;
   std::map<std::string,float> floatattrs;
   std::map<std::string,Vector> vectorattrs;
   std::map<std::string,Color> colorattrs;


 
// int attributes

//    intattrs["wispchildren"] = 0;


// float attributes

    floatattrs["opacity"] = 1.0;
    floatattrs["pscale"] = 1.0;
    floatattrs["octaves"] = 1.0;
    floatattrs["roughness"] = 0.5;
    floatattrs["freq"] = 1.0;
    floatattrs["fjump"] = 2.0;
    floatattrs["offset"] = 0.0;
    floatattrs["shutter"] = 0.0;
    floatattrs["framerate"] = 24.0;
    floatattrs["lifetime"] = 1.0;
    floatattrs["age"] = 1.0;
    floatattrs["gamma"] = 1.0;


// Vector attributes

    vectorattrs["translate"]  = Vector(0,0,0);
    vectorattrs["accel"]      = Vector(0,0,0);


// Color attributes



// Populate point cloud
   for( std::map<std::string,int>::const_iterator i=intattrs.begin(); i != intattrs.end(); i++ )
   {
      pg->create_attr( i->first, i->second );
   }
   for( std::map<std::string,float>::const_iterator i=floatattrs.begin(); i != floatattrs.end(); i++ )
   {
      pg->create_attr( i->first, i->second );
   }
   for( std::map<std::string,Vector>::const_iterator i=vectorattrs.begin(); i != vectorattrs.end(); i++ )
   {
      pg->create_attr( i->first, i->second );
   }
   for( std::map<std::string,Color>::const_iterator i=colorattrs.begin(); i != colorattrs.end(); i++ )
   {
      pg->create_attr( i->first, i->second );
   }
   pg->add_particles(nb);
   return pg;
}

void lux::writeObj( const std::string& fname, const PointCloud& pc )
{
   std::ofstream file( fname.c_str() );
   if( !file ){ return; }
   file << "# OBJ file generated by bishop\n";
   file << "# vertices: " << pc->nb_particles() << "\n";
   for( size_t i=0;i<pc->nb_particles();i++ )
   {
      const Vector& v = pc->pos(i);
      file << "v " << v.X() << " " << v.Y() << " " << v.Z() << std::endl;
   }
  file.close();
}
