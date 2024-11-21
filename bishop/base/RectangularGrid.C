#include <cfloat>
#include <sstream>
#include "RectangularGrid.h"
using namespace lux;

#include <iostream>
using namespace std;

void RectangularGrid::init( int nnnx, int nnny, int nnnz, 
               double Lx, double Ly, double Lz,
	       const Vector& Origin        )
{
   nX = nnnx;
   nY = nnny;
   nZ = nnnz;

   lX = Lx;
   lY = Ly;
   lZ = Lz;

   origin = Origin;
   topleft = origin + Vector(lX, lY, lZ );

   dX = lX/nX;
   dY = lY/nY;
   dZ = lZ/nZ;
}


void RectangularGrid::init( const Vector&llc, const Vector& urc, const Vector& res )
{
   dX = res[0];
   dY = res[1];
   dZ = res[2];

   origin = llc;

   topleft = urc;

   lX = urc[0]-llc[0];
   lY = urc[1]-llc[1];
   lZ = urc[2]-llc[2];

   nX = (int)(lX/dX);
   nY = (int)(lY/dY);
   nZ = (int)(lZ/dZ);

   dX = lX/nX;
   dY = lY/nY;
   dZ = lZ/nZ;
}

void RectangularGrid::reset( double Lx, double Ly, double Lz, const Vector& Origin )
{
   lX = Lx;
   lY = Ly;
   lZ = Lz;

   origin = Origin;
   topleft = origin + Vector(lX, lY, lZ );

   dX = lX/nX;
   dY = lY/nY;
   dZ = lZ/nZ;
}


const Vector RectangularGrid::evalP( int i, int j, int k ) const
{
   return origin + Vector( i*dX, j*dY, k*dZ );
}


const bool RectangularGrid::isInside( const Vector& P ) const
{
   Vector X = transform(P);
   double test = X[0]/lX;
   if( test < 0 || test >= 1.0 ){ return false; }
   test = X[1]/lY;
   if( test < 0 || test >= 1.0 ){ return false; }
   test = X[2]/lZ;
   if( test < 0 || test >= 1.0 ){ return false; }
   return true;
}




int RectangularGrid::index( int i, int j, int k ) const { return ( i + nX*( j + nY*k) ); }
int RectangularGrid::index( const Vector& P ) const
{
   int ii,jj,kk;
   if(getGridIndex( P, ii,jj,kk ))
   {
      return index( ii, jj, kk );
   }
   return -1;
}


void RectangularGrid::triple( const int indx, int& i, int& j, int& k ) const
{
   k = indx / (nX*nY);
   j = (indx - k*nX*nY)/nX;
   i = indx - nX*(j+k*nY);
}




const float fade( const float t )
{
   return  ( t * t * t * ( t * ( t * 6.0 - 15.0 ) + 10.0 ) ); 
}



void linearInterpolation( float x, float dx, int nx, int& i, int& ii, float& w, float& ww, bool isperiodic )
{
    float r = x/dx;
    i  =  r;
    if( !isperiodic )
    {
    if( i >= 0 && i < nx )
    {
       ii = i + 1;
       w = r-i;
       ww = 1.0 - w;
       if(  ii >= nx )
       {
	  ii = nx-1;
       }
    }
    else
    {
       i = ii = 0;
       w = ww = 0;
    }
    }
    else
    {
       w = r-i;
       while( i < 0 ){ i += nx; }
       while( i >= nx ){ i -= nx; }

       ii = i+1;
       ww = 1.0 - w;
       if( ii >= nx ){ ii -= nx; }
    }
}

void RectangularGrid::getLinearInterpolation( const Vector& P,  
                                 int& ix, int& iix, float& wx, float& wwx, 
				 int& iy, int& iiy,  float& wy, float& wwy,
				 int& iz, int& iiz,  float& wz, float& wwz ) const
{
   Vector X = transform( P );
   linearInterpolation( X[0], dX, nX, ix, iix, wx, wwx, isPeriodicX );
   linearInterpolation( X[1], dY, nY, iy, iiy, wy, wwy, isPeriodicY );
   linearInterpolation( X[2], dZ, nZ, iz, iiz, wz, wwz, isPeriodicZ );

   //wx = fade(wx);
   //wy = fade(wy);
   //wz = fade(wz);
   //wwx = 1.0-wx;
   //wwy = 1.0-wy;
   //wwz = 1.0-wz;
}




void RectangularGrid::getHighOrderInterpolation( const Vector& P, 
                                    vector<int>& indices_x, vector<int>& indices_y, vector<int>& indices_z, 
                                    vector<double>& weights_x, vector<double>&
                                    weights_y, vector<double>& weights_z )
                                    const
{
   indices_x.clear();
   indices_y.clear();
   indices_z.clear();

   weights_x.clear();
   weights_y.clear();
   weights_z.clear();

   Vector X = transform( P );

   double rx = X[0]/dX;
   if( rx < 0.0 ){ return; }
   int ix = (int)rx;
   if( ix >= nX ){ return; }
   double eps_x = rx -(double)ix;
   int order_x = ho_order;
   if( ix + order_x >= nX ){ order_x = nX - ix - 1; }
   if( ix - order_x + 1 < 0  ){ order_x = ix + 1; }
   ho_interp.weights( eps_x, order_x, weights_x );
   for( int i = ix - order_x + 1; i<= ix+order_x; i++ ){ indices_x.push_back(i); }

   double ry = X[1]/dY;
   if( ry < 0.0 ){ return; }
   int iy = (int)ry;
   if( iy >= nY ){ return; }
   double eps_y = ry -(double)iy;
   int order_y = ho_order;
   if( iy + order_y >= nY ){ order_y = nY - iy - 1; }
   if( iy - order_y + 1 < 0  ){ order_y = iy + 1; }
   ho_interp.weights( eps_y, order_y, weights_y );
   for( int i = iy - order_y + 1; i<= iy+order_y; i++ ){ indices_y.push_back(i); }

   double rz = X[2]/dZ;
   if( rz < 0.0 ){ return; }
   int iz = (int)rz;
   if( iz >= nZ ){ return; }
   double eps_z = rz -(double)iz;
   int order_z = ho_order;
   if( iz + order_z >= nZ ){ order_z = nZ - iz - 1; }
   if( iz - order_z + 1 < 0  ){ order_z = iz + 1; }
   ho_interp.weights( eps_z, order_z, weights_z );
   for( int i = iz - order_z + 1; i<= iz+order_z; i++ ){ indices_z.push_back(i); }

}

















const bool RectangularGrid::isInside( int i, int j, int k ) const
{
   if( i<0 ){ return false; }
   if( j<0 ){ return false; }
   if( k<0 ){ return false; }
   if( i>=nX ){ return false; }
   if( j>=nY ){ return false; }
   if( k>=nZ ){ return false; }
   return true;
}






const bool RectangularGrid::getGridIndex( const Vector& P, int& ix, int& iy, int& iz ) const
{
   const Vector lower = transform( P );
   ix =  (int)(lower[0]/dX);
   iy =  (int)(lower[1]/dY);
   iz =  (int)(lower[2]/dZ);
   if( ix < 0 ){ return false; }
   if( iy < 0 ){ return false; }
   if( iz < 0 ){ return false; }
   if( ix >= nX ){ return false; }
   if( iy >= nY ){ return false; }
   if( iz >= nZ ){ return false; }
   return true;
}

const bool RectangularGrid::getBox( const Vector& ll, const Vector& uu, int& ixl, int& ixu, int& iyl, int& iyu, int& izl, int& izu ) const
{
   const Vector lower = transform( ll );
   const Vector upper = transform( uu );
   ixl =  (int)(lower[0]/dX - 0.0);
   ixu =  (int)(upper[0]/dX + 1.0);
   iyl =  (int)(lower[1]/dY - 0.0);
   iyu =  (int)(upper[1]/dY + 1.0);
   izl =  (int)(lower[2]/dZ - 0.0);
   izu =  (int)(upper[2]/dZ + 1.0);

   // tests
   if( !isPeriodicX )
   {
      if( ixu < 0 ){ return false; }
      if( ixl >= nX ){ return false; }
      if( ixl < 0 ){ ixl = 0; }
      if( ixu >= nX ){ ixu = nX-1; }
   }

   if( !isPeriodicY )
   {
      if( iyu < 0 ){ return false; }
      if( iyl >= nY ){ return false; }
      if( iyl < 0 ){ iyl = 0; }
      if( iyu >= nY ){ iyu = nY-1; }
   }

   if( !isPeriodicZ )
   {
      if( izu < 0 ){ return false; }
      if( izl >= nZ ){ return false; }
      if( izl < 0 ){ izl = 0; }
      if( izu >= nZ ){ izu = nZ-1; }
   }
   return true;
}




RectangularGrid lux::Union( const RectangularGrid& g1, const RectangularGrid& g2 )
{
   Vector llc1 = g1.llc();
   Vector ullc = g2.llc();
   ullc[0] = ( ullc[0] <= llc1[0] ) ? ullc[0] : llc1[0];
   ullc[1] = ( ullc[1] <= llc1[1] ) ? ullc[1] : llc1[1];
   ullc[2] = ( ullc[2] <= llc1[2] ) ? ullc[2] : llc1[2];

   Vector urc1 = g1.urc();
   Vector uurc = g2.llc();
   uurc[0] = ( uurc[0] <= urc1[0] ) ? uurc[0] : urc1[0];
   uurc[1] = ( uurc[1] <= urc1[1] ) ? uurc[1] : urc1[1];
   uurc[2] = ( uurc[2] <= urc1[2] ) ? uurc[2] : urc1[2];

   double dx = ( g1.dx() < g2.dx() ) ? g1.dx() : g2.dx();
   double dy = ( g1.dy() < g2.dy() ) ? g1.dy() : g2.dy();
   double dz = ( g1.dz() < g2.dz() ) ? g1.dz() : g2.dz();

   RectangularGrid out;
   out.init( ullc, uurc, Vector( dx,dy,dz ) );
   return out;
}




























void SparseMapRectangularGrid::init(  
               float dx, float dy, float dz,
	       const Vector& Origin        )
{
   origin = Origin;
   dX = dx;
   dY = dy;
   dZ = dz;
}




const Vector SparseMapRectangularGrid::evalP( int i, int j, int k ) const
{
   return origin + Vector( (i)*dX, (j)*dY, (k)*dZ );
}


void sparseLinearInterpolation( float x, float dx, int& i, int& ii, float& w, float& ww )
{
    float r = x/dx;
    i  =  r;
    ii = i + 1;
    w = r-i;
    ww = 1.0 - w;
}

void SparseMapRectangularGrid::getLinearInterpolation( const Vector& P,  
                                 int& ix, int& iix, float& wx, float& wwx, 
				 int& iy, int& iiy,  float& wy, float& wwy,
				 int& iz, int& iiz,  float& wz, float& wwz ) const
{
   Vector X = transform( P );
   sparseLinearInterpolation( X[0], dX, ix, iix, wx, wwx );
   sparseLinearInterpolation( X[1], dY, iy, iiy, wy, wwy );
   sparseLinearInterpolation( X[2], dZ, iz, iiz, wz, wwz );
}

const bool SparseMapRectangularGrid::getGridIndex( const Vector& P, int& ix, int& iy, int& iz ) const
{
   const Vector lower = transform( P );
   ix =  (int)(lower[0]/dX + 0.5);
   iy =  (int)(lower[1]/dY + 0.5);
   iz =  (int)(lower[2]/dZ + 0.5);
   return true;
}

const bool SparseMapRectangularGrid::getBox( const Vector& ll, const Vector& uu, int& ixl, int& ixu, int& iyl, int& iyu, int& izl, int& izu ) const
{
   const Vector lower = transform( ll );
   const Vector upper = transform( uu );
   ixl =  (int)(lower[0]/dX - 0.5);
   ixu =  (int)(upper[0]/dX + 0.5);
   iyl =  (int)(lower[1]/dY - 0.5);
   iyu =  (int)(upper[1]/dY + 0.5);
   izl =  (int)(lower[2]/dZ - 0.5);
   izu =  (int)(upper[2]/dZ + 0.5);
   return true;
}





GridBox::GridBox( RectangularGrid* f ) : GridBoxBase(f) {}

GridBox lux::makeGridBox( const Vector& llc, const Vector& urc, const Vector& dx )
{
   RectangularGrid* rg = new RectangularGrid();
   rg->init( llc, urc, dx );
   return GridBox(rg);
}


GridBox lux::makeGridBox( const RectangularGrid& rg )
{
   return lux::makeGridBox( rg.llc(), rg.urc(), Vector( rg.dx(), rg.dy(), rg.dz() ) );
}


GridBox::~GridBox(){}

const Vector GridBox::evalP( int i, int j, int k ) const { return (*this)->evalP(i,j,k); }
int GridBox::index( const Vector& P ) const { return (*this)->index(P); }

GridBox GridBox::operator+=( const GridBox& e2 )
{
   Vector l1 = (*this)->llc();
   Vector l2 = e2->llc();
   l1[0] = ( l1[0] <= l2[0] ) ? l1[0] : l2[0];
   l1[1] = ( l1[1] <= l2[1] ) ? l1[1] : l2[1];
   l1[2] = ( l1[2] <= l2[2] ) ? l1[2] : l2[2];

   Vector L1 = (*this)->urc();
   Vector L2 = e2->urc();
   L1[0] = ( L1[0] >= L2[0] ) ? L1[0] : L2[0];
   L1[1] = ( L1[1] >= L2[1] ) ? L1[1] : L2[1];
   L1[2] = ( L1[2] >= L2[2] ) ? L1[2] : L2[2];

   Vector D1( (*this)->dx(), (*this)->dy(), (*this)->dz() );
   Vector D2( e2->dx(), e2->dy(), e2->dz() );
   D1[0] = ( D1[0] <= D2[0] ) ? D1[0] : D2[0];
   D1[1] = ( D1[1] <= D2[1] ) ? D1[1] : D2[1];
   D1[2] = ( D1[2] <= D2[2] ) ? D1[2] : D2[2];
   
   RectangularGrid* rg = new RectangularGrid();
   rg->init( l1, L1, D1 );

   return GridBox( rg );
}


const double lux::dx( const GridBox& gb ){ return gb->dx(); }
const double lux::dy( const GridBox& gb ){ return gb->dy(); }
const double lux::dz( const GridBox& gb ){ return gb->dz(); }
const int lux::nx( const GridBox& gb ){ return gb->nx(); }
const int lux::ny( const GridBox& gb ){ return gb->ny(); }
const int lux::nz( const GridBox& gb ){ return gb->nz(); }
const double lux::Lx( const GridBox& gb ){ return gb->Lx(); }
const double lux::Ly( const GridBox& gb ){ return gb->Ly(); }
const double lux::Lz( const GridBox& gb ){ return gb->Lz(); }
const Vector lux::llc( const GridBox& gb ) { return gb->llc(); }
const Vector lux::urc( const GridBox& gb ) { return gb->urc(); }
void lux::setPeriodic( GridBox& gb ) { gb->setPeriodic(); }
void lux::setInterpOrder( GridBox& gb, int order ) { gb->setInterpolationOrder(order); }
int lux::getInterpOrder( const GridBox& gb ) { return gb->getInterpolationOrder(); }
AABB lux::getAABB( const GridBox& gb ) { return AABB( gb->llc(), gb->urc() ); }

template <typename T> 
std::string tostr(const T& t) { std::stringstream os; os<<t; return os.str(); }


char * GridBox::__str__()
{
   static char docLabel[2048];
   std::string lbl = "Points: " + tostr(nx(*this)) + " X " + tostr(ny(*this)) + " X " + tostr(nz(*this));
   lbl = lbl + "  cellsize: "  + tostr(dx(*this)) + " X " + tostr(dy(*this)) + " X " + tostr(dz(*this));
   lbl = lbl + "  length: "  + tostr(Lx(*this)) + " X " + tostr(Ly(*this)) + " X " + tostr(Lz(*this));
   lbl = lbl + "  llc: "  + std::string(llc(*this).__str__());
   lbl = lbl + "  urc: "  + std::string(urc(*this).__str__());
   size_t lbllength = lbl.size();
   if( lbllength > 2047 ){ lbllength = 2047; }
   lbllength = lbl.copy( docLabel, lbllength);
   docLabel[lbllength] = '\0';
   return docLabel;
}
