
#ifndef ____PARTICLEGROUP_H____
#define ____PARTICLEGROUP_H____

//
//  ParticleGroup.h
//
//  Generic particles with attributes
//


#include "Vector.h"
#include "Color.h"
#include <vector>
#include <cstring>
#include <map>
#include <memory>
#include <string>


namespace lux
{

template<typename T>
class PGAttribute
{
  public:

    PGAttribute() : name("unknown") {}
    PGAttribute( const std::string& nam, const T& def ) : name(nam), defVal(def) {}
   ~PGAttribute(){}

    const size_t size() const { return data.size(); }
    const bool empty() const { return data.empty(); }
    void set(size_t i, const T& value ) { data[i] = value; }
    const T& get(size_t i ) const { return data[i]; }
    void expand_to( size_t n )
    {
       size_t current = data.size();
       for( size_t i=current;i<n;i++ )
       {
          data.push_back( defVal );
       }
    }
    void clear() { data.clear(); }
    const std::string& attr_name() const { return name; }
    const T& default_value() const { return defVal; }
    typename std::vector<T>::const_iterator cbegin() const { return data.begin(); }
    typename std::vector<T>::const_iterator cend() const { return data.end(); }
    typename std::vector<T>::iterator begin() { return data.begin(); }
    typename std::vector<T>::iterator end() { return data.end(); }


  private:
    std::vector<T> data;
    std::string name;
    T defVal;
};


class ParticleGroup
{
  public:

    ParticleGroup();
   ~ParticleGroup();


    void create_attr( const std::string& nam, const int& def );
    void create_attr( const std::string& nam, const float& def );
    void create_attr( const std::string& nam, const Vector& def );
    void create_attr( const std::string& nam, const Color& def );

    const size_t add_particle();
    const size_t add_particles( const size_t nb );
    size_t nb_particles() const;

    const int& get_int_attr( const std::string& nam, const size_t p ) const;
    const float& get_float_attr( const std::string& nam, const size_t p ) const;
    const Vector& get_vector_attr( const std::string& nam, const size_t p ) const;
    const Color& get_color_attr( const std::string& nam, const size_t p ) const;
    
    const int& id( const size_t p ) const;
    const Vector& pos( const size_t p ) const;
    const Vector& vel( const size_t p ) const;
    const Color& ci( const size_t p ) const;
    const float& pscale( const size_t p ) const;
    const float& density( const size_t p ) const;
    const float& roughness( const size_t p ) const;
    const float& octaves( const size_t p ) const;
    const float& fjump( const size_t p ) const;
    const float& frequency( const size_t p ) const;
    const Vector& translate( const size_t p ) const;
    
    void set_attr( const std::string& nam, const size_t p, const int& value ); 
    void set_attr( const std::string& nam, const size_t p, const float& value ); 
    void set_attr( const std::string& nam, const size_t p, const Vector& value ); 
    void set_attr( const std::string& nam, const size_t p, const Color& value ); 

    void set_id( const size_t p, const int& value );
    void set_pos( const size_t p, const Vector& value );
    void set_vel( const size_t p, const Vector& value );
    void set_ci( const size_t p, const Color& value );
    void set_pscale( const size_t p, const float& value );
    void set_density( const size_t p, const float& value );
    void set_roughness( const size_t p, const float& value );
    void set_octaves( const size_t p, const float& value );
    void set_fjump( const size_t p, const float& value );
    void set_frequency( const size_t p, const float& value );
    void set_translate( const size_t p, const Vector& value );

    std::vector<std::string> show_int_attrs() const;
    std::vector<std::string> show_float_attrs() const;
    std::vector<std::string> show_vector_attrs() const;
    std::vector<std::string> show_color_attrs() const;
    std::vector<std::string> show_all_attrs() const;

    bool attr_exists( const std::string& nam ) const;
    bool int_attr_exists( const std::string& nam ) const;
    bool float_attr_exists( const std::string& nam ) const;
    bool vector_attr_exists( const std::string& nam ) const;
    bool color_attr_exists( const std::string& nam ) const;

    void merge( const ParticleGroup& g );


  private:

// attribute collections
    std::map< std::string, PGAttribute<int>  > int_attributes;
    std::map< std::string, PGAttribute<float>  > float_attributes;
    std::map< std::string, PGAttribute<Vector> > vector_attributes;  
    std::map< std::string, PGAttribute<Color> > color_attributes;  



};



typedef std::shared_ptr<ParticleGroup> PointCloudBase;

class PointCloud : public PointCloudBase
{
  public:

    PointCloud();
    PointCloud( ParticleGroup* f );
   ~PointCloud();


     char* __str__();
     char* __doc__();

   // Merge point clouds
   PointCloud operator+( const PointCloud& e2 );

};


}
#endif
