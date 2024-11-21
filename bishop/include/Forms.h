//*******************************************************************
//
//   Forms.h
//
// 3D Forms class in the namespace lux
//
//
//
//*******************************************************************

#ifndef __LUX_FORMS_H__
#define __LUX_FORMS_H__

#include <cmath>
#include <cstdio>

#include "Vector.h"

namespace lux
{

//! Form is a 3D form class
class Form
{
  public:

   Form()
     { 
        zeroform = 0;
        oneform[0] = oneform[1] = oneform[2] = 0; 
        twoform[0] = twoform[1] = twoform[2] = 0;
	threeform = 0;
     }

   Form(const Form& v)
   { 
      zeroform = v.zeroform;
      oneform = v.oneform; 
      twoform = v.twoform;
      threeform = v.threeform;
   }
   
   Form(const double a, const Vector& b, const Vector& c, const double d)
   {
        zeroform = a;
        oneform = b; 
        twoform = c;
	threeform = d;
   }

   ~Form(){}

   //!  Set all components
   void set( const float vx, const Vector& vy, const Vector& vz, const float vw )
   {
      zeroform = vx;
      oneform = vy; 
      twoform = vz;
      threeform = vw;
   }

   //! Add two forms together
   const Form operator+        (const Form& v) const 
   { 
      return Form( zeroform + v.zeroform, oneform + v.oneform, twoform + v.twoform, threeform + v.threeform ); 
   }
  
   //! Subtract one form from another
   const Form operator-        (const Form& v) const
   { 
      return Form(  zeroform - v.zeroform, oneform - v.oneform, twoform - v.twoform, threeform - v.threeform  ); 
   }

   //! Unary minus
   friend const Form operator- (const Form& v)
   { return Form(-v.zeroform,-v.oneform,-v.twoform, -v.threeform); }

   //! Multiplication of a constant with a form
   friend const Form operator* (const double w, const Form& v)
   { return v*w; }
	  
   //! Multiplication of a form with a constant
   const Form operator*        (const double v) const
   { return Form(zeroform*v, oneform*v, twoform*v, threeform*v); }

   const Form operator/        (const double v) const
   { return Form(zeroform/v, oneform/v, twoform/v, threeform/v); }

   //! wedge product
   const Form operator^        (const Form& v) const 
   { return   Form(zeroform*v.zeroform, 
		   (zeroform*v.oneform) + (oneform*v.zeroform),
		   (zeroform*v.twoform) + (twoform*v.zeroform) + (oneform^v.oneform),
		   (zeroform*v.threeform) + (threeform*v.zeroform) + (oneform*v.twoform) - (twoform*v.oneform)
		   ); 
   }


   const Form star() const { return Form( threeform, twoform, oneform, zeroform ); }

   Form& operator=       (const Form& v)
   { zeroform = v.zeroform; oneform = v.oneform; twoform = v.twoform; threeform = v.threeform; return *this; }
  
   Form& operator+=      (const Form& v)
   { zeroform += v.zeroform; oneform += v.oneform; twoform += v.twoform; threeform += v.threeform; return *this; }
  
   Form& operator-=      (const Form& v)
   { zeroform -= v.zeroform; oneform -= v.oneform; twoform -= v.twoform; threeform -= v.threeform; return *this; }
  
   Form& operator*=      (const double v)
   { zeroform *= v; oneform *= v; twoform *= v; threeform *= v; return *this; }
  
   Form& operator/=      (const double v)
   { zeroform /= v; oneform /= v; twoform /= v; threeform /= v; return *this; }
  

   const double& zero() const { return zeroform; } 
         double& zero()       { return zeroform; }

   const Vector& one() const { return oneform; } 
         Vector& one()       { return oneform; }

   const Vector& two() const { return twoform; } 
         Vector& two()       { return twoform; }

   const double& three() const { return threeform; } 
         double& three()       { return threeform; }


//  Comparisons

   const bool operator==         (const Form& v) const
       { return ( zeroform==v.zeroform && oneform==v.oneform && twoform==v.twoform && threeform==v.threeform ); }
  
   const bool operator!=         (const Form& v) const
       { return ( zeroform!=v.zeroform && oneform!=v.oneform && twoform!=v.twoform && threeform!=v.threeform ); }
  

   char *__str__() {
       static char tmp[1024];
       std::sprintf(tmp,"Form[%g,(%g,%g,%g),(%g,%g,%g),%g]", zeroform, oneform[0],oneform[1],oneform[2], twoform[0], twoform[1],twoform[2], threeform);
       return tmp;
   }

  private:
  double zeroform;
  Vector oneform;
  Vector twoform;
  double threeform;
};

const Form wedge( const Form& f1, const Form& f2 );
const Form star( const Form& f );

}



#endif
