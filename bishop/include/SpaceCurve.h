
#ifndef __SPACECURVE_H__
#define __SPACECURVE_H__

#include "Vector.h"
#include "Matrix.h"
#include <vector>
#include <map>
#include <string>
#include <memory>
#include <cstring>


namespace lux
{

class CurveFS
{
  public:
  
    CurveFS(const double Qmax, const double Qmin);
   ~CurveFS();
    
    virtual const Vector eval( const float q ) const ;
    virtual const Vector grad( const float q ) const ;

    virtual const Vector fsT( const float q ) const ;
    virtual const Vector fsN( const float q ) const ;
    virtual const Vector fsB( const float q ) const ;

    virtual double speed    ( const float q ) const ;
    virtual double curvature( const float q ) const ;
    virtual double torsion  ( const float q ) const ;

    double qMax() const;
    double qMin() const;
    
    virtual std::string typelabel() { return "unknown"; }
    virtual std::string documentation() { return "No documentation"; }

  protected:

    double qmax, qmin;

};

typedef std::shared_ptr< CurveFS > SpaceCurveBase;

class SpaceCurve : public SpaceCurveBase
{
  public:

    SpaceCurve();
    SpaceCurve( CurveFS* f );
   ~SpaceCurve();

     char* __str__() 
     {
       static char typeLabel[2048];
       std::string lbl = (*this)->typelabel();
       size_t lbllength = lbl.size();
       if( lbllength > 2047 ){ lbllength = 2047; }
       lbllength = lbl.copy( typeLabel, lbllength);
       typeLabel[lbllength] = '\0';
       return typeLabel;
    }


     char* __doc__() 
     {
       static char docLabel[2048];
       std::string lbl = (*this)->documentation();
       size_t lbllength = lbl.size();
       if( lbllength > 2047 ){ lbllength = 2047; }
       lbllength = lbl.copy( docLabel, lbllength);
       docLabel[lbllength] = '\0';
       return docLabel;
    }

};



}

#endif
