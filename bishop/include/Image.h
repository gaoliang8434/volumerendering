
#ifndef __IMAGE_H__
#define __IMAGE_H__

#include <vector>
#include <cstddef>

using namespace std;

namespace lux
{

typedef std::vector<float> Pixel;


class Image
{
  public:


    //! Constructor for an empty image
    Image() :
       width  (0),
       height (0),
       depth  (0)
    {}

   ~Image(){}

    //! Set the size and depth of the image.  Default is a depth of 3 (rgb).
    void reset(  int w, int h, int d=3 )
    {
       width = w;
       height = h;
       depth = d;
       Pixel pixel;
       pixel.resize(depth);
       for( size_t i=0;i<(size_t)depth;i++ ){ pixel[i] = 0.0; }
       data.resize(width*height);
       for( size_t i=0;i<data.size();i++ ){ data[i] = pixel; }
    }
   
    //! Returns a const reference to a pixel value
    const float& value( int x, int y, int c ) const { return data[ index(x,y) ][(size_t)c]; }

    //! Returns a read-write reference to a pixel value
          float& value( int x, int y, int c ) { return data[ index(x,y) ][(size_t)c]; }

    //! Retrieve a pixel with all of its components by coordiantes 
    std::vector<float>& pixel(int x, int y ) { return data[ index(x,y) ]; }

    //! Retrieve a pixel with all of its components by index
    std::vector<float>& pixel(int ind ) { return data[ ind ]; }
   
    //! Get horizontal dimension size
    const int Width() const { return width; }

    //! Get vertical dimension dimension size
    const int Height() const { return height; }

    //! Get number of channels / image depth
    const int Depth() const { return depth; }

    //! Calculate index of a particular pixel from coordinates
    const size_t index( int x, int y ) const { return (size_t) ( x + width*y ) ; }

    //! Interpolates a channel
    const float interpolatedValue( float x, float y, int c ) const; 
    //! Interpolates a pixel 
    std::vector<float> interpolatedPixel( float x, float y ) const;

    //! Constant pointer to pixel buffer for quick iteration in Python
    const std::vector< std::vector<float> >* pixels() const { return &data; }

  protected:

    int width, height, depth;
    std::vector< std::vector<float> > data;
    

    void interpolationCoefficients( float x, float y, 
                                   float& wx, float& wwx,
				   float& wy, float& wwy,
				   int& ix, int& iix,
				   int& iy, int& iiy 
				 ) const;
};

void setPixel( Image& img, int x, int y, std::vector<float>& value );
void setPixel( Image& img, int x, int y, float r, float g, float b, float a );

}
#endif



