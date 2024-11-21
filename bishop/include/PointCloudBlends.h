#ifndef ____POINTCLOUDBLENDS_H____
#define ____POINTCLOUDBLENDS_H____

#include "Vector.h"
#include "Color.h"

namespace lux
{

class PointCloudBlend
{
  public:
    PointCloudBlend(){}
    virtual ~PointCloudBlendMethod(){}

    virtual int blend( const int v1, const int v2 ) const = 0;
    virtual int blend( const int v1, const float d1, const int v2, const float d2 ) const = 0;

    virtual float blend( const float v1, const float v2 ) const = 0;
    virtual float blend( const float v1, const float d1, const float v2, const float d2 ) const = 0;

    virtual Vector blend( const Vector& v1, const Vector& v2 ) const = 0;
    virtual Vector blend( const Vector& v1, const float d1, const Vector& v2, const float d2 ) const = 0;

    virtual Color blend( const Color& v1, const Color& v2 ) const = 0;
    virtual Color blend( const Color& v1, const float d1, const Color& v2, const float d2 ) const = 0;
};


class PointCloudAddBlend : public PointCloudBlendMethod
{
  public:
    PointCloudAddBlend(){}
   ~PointCloudAddBlend(){}

    int blend( const int v1, const int v2 ) const { return v1+v2; }
    int blend( const int v1, const float d1, const int v2, const float d2 ) const { return (int)((v1*d1 + v2*d2)/(d1+d2)); }

    float blend( const float v1, const float v2 ) const { return v1+v2; }
    float blend( const float v1, const float d1, const float v2, const float d2 ) const { return (v1*d1 + v2*d2)/(d1+d2); }

    Vector blend( const Vector& v1, const Vector& v2 ) const  { return v1+v2; }
    Vector blend( const Vector& v1, const float d1, const Vector& v2, const float d2 ) const { return (v1*d1 + v2*d2)/(d1+d2); }

    Color blend( const Color& v1, const Color& v2 ) const { return v1+v2; }
    Color blend( const Color& v1, const float d1, const Color& v2, const float d2 ) const { return (v1*d1 + v2*d2)/(d1+d2); }
};

class PointCloudMaxBlend : public PointCloudBlendMethod
{
  public:
    PointCloudMaxBlend(){}
   ~PointCloudMaxBlend(){}

    int blend( const int v1, const int v2 ) const { return (v1>v2) ? v1 : v2; }
    int blend( const int v1, const float d1, const int v2, const float d2 ) const { return (d1>d2)  ? v1 : v2; }

    float blend( const float v1, const float v2 ) const { return (v1>v2) ? v1 : v2; }
    float blend( const float v1, const float d1, const float v2, const float d2 ) const { return (d1>d2) ? v1 : v2; }

    Vector blend( const Vector& v1, const Vector& v2 ) const  { return (v1>v2) ? v1 : v2; }
    Vector blend( const Vector& v1, const float d1, const Vector& v2, const float d2 ) const { return (d1>d2) ? v1 : v2; }

    Color blend( const Color& v1, const Color& v2 ) const { return (v1>v2) ? v1 : v2; }
    Color blend( const Color& v1, const float d1, const Color& v2, const float d2 ) const { return (d1>d2) ? v1 : v2; }
};



}
#endif
