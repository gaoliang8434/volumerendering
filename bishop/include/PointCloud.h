#ifndef ____POINTCLOUD_H____
#define ____POINTCLOUD_H____

#include "ParticleGroup.h"
#include <string>
#include <cstring>


namespace lux
{

PointCloud makePointCloud( int nb ); 

void create_attr( PointCloud& p, const std::string& nam, const int& def );
void create_attr( PointCloud& p, const std::string& nam, const float& def );
void create_attr( PointCloud& p, const std::string& nam, const Vector& def );
void create_attr( PointCloud& p, const std::string& nam, const Color& def );

size_t add_particle( PointCloud& p );
size_t add_particles( PointCloud& p, const size_t nb );
size_t nb_particles( PointCloud& p );

int get_int_attr( PointCloud& p, const std::string& nam, const size_t i );
float get_float_attr( PointCloud& p, const std::string& nam, const size_t i );
Vector get_vector_attr( PointCloud& p, const std::string& nam, const size_t i );
Color get_color_attr( PointCloud& p, const std::string& nam, const size_t i );
    
int id(  PointCloud& p, const size_t i );
Vector pos(  PointCloud& p, const size_t i );
Vector vel(  PointCloud& p, const size_t i );
Color ci(  PointCloud& p, const size_t i );
float pscale(  PointCloud& p, const size_t i );
float density(  PointCloud& p, const size_t i );
float roughness(  PointCloud& p, const size_t i );
float octaves(  PointCloud& p, const size_t i );
float fjump(  PointCloud& p, const size_t i );
float frequency(  PointCloud& p, const size_t i );
Vector translate(  PointCloud& p, const size_t i );
    
void set_attr(  PointCloud& p, const std::string& nam, const size_t i, const int& value ); 
void set_attr(  PointCloud& p, const std::string& nam, const size_t i, const float& value ); 
void set_attr(  PointCloud& p, const std::string& nam, const size_t i, const Vector& value ); 
void set_attr(  PointCloud& p, const std::string& nam, const size_t i, const Color& value ); 

void set_id(  PointCloud& p, const size_t i, const int& value );
void set_pos(  PointCloud& p, const size_t i, const Vector& value );
void set_vel(  PointCloud& p, const size_t i, const Vector& value );
void set_ci(  PointCloud& p, const size_t i, const Color& value );
void set_pscale(  PointCloud& p, const size_t i, const float& value );
void set_density(  PointCloud& p, const size_t i, const float& value );
void set_roughness(  PointCloud& p, const size_t i, const float& value );
void set_octaves(  PointCloud& p, const size_t i, const float& value );
void set_fjump(  PointCloud& p, const size_t i, const float& value );
void set_frequency(  PointCloud& p, const size_t i, const float& value );
void set_translate(  PointCloud& p, const size_t i, const Vector& value );

std::vector<std::string> show_int_attrs( PointCloud& p);
std::vector<std::string> show_float_attrs( PointCloud& p);
std::vector<std::string> show_vector_attrs( PointCloud& p);
std::vector<std::string> show_color_attrs( PointCloud& p);
std::vector<std::string> show_all_attrs( PointCloud& p);

int attr_exists(  PointCloud& p, const std::string& nam );
int int_attr_exists(  PointCloud& p, const std::string& nam );
int float_attr_exists(  PointCloud& p, const std::string& nam );
int vector_attr_exists(  PointCloud& p, const std::string& nam );
int color_attr_exists(  PointCloud& p, const std::string& nam );

void merge(  PointCloud& p, const PointCloud& g );

Vector llc(const PointCloud& p );
Vector urc(const PointCloud& p );


PointCloud makeWispCloud( int nb );
PointCloud makeNoiseCloud( int nb );
PointCloud makePyroCloud( int nb );

void writeObj( const std::string& nam, const PointCloud& pc );



}
#endif
