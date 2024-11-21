
#ifndef ____HIGHORDERINTERPOLATOR_H____
#define ____HIGHORDERINTERPOLATOR_H____


#include <vector>

namespace lux {

class HighOrderInterpolator
{
  public:

    HighOrderInterpolator(){}
   ~HighOrderInterpolator(){}


    void weights( const double epsilon, int order, std::vector<double>& waits ) const;


  private:

    static std::vector<  std::vector< std::vector<double> >   > tables;


};

}

#endif
