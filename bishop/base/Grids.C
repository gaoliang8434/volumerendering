

#include "Grids.h"
#include "Stamp.h"
#include "LinearAlgebra.h"
#include "Fields.h"

#include <omp.h>

using namespace std;
using namespace lux;





ScalarGrid lux::makeScalarGrid( const RectangularGrid& rg, float defValue ) { return makeGrid( rg, defValue); }

ScalarGrid lux::makeGrid( const RectangularGrid& rg, float defValue )
{
   SGrid<float>* f = new SGrid<float>();
   f->setDefVal( defValue );
   f->init( rg.nx(), rg.ny(), rg.nz(), rg.Lx(), rg.Ly(), rg.Lz(), rg.llc() );
   f->setInterpolationOrder( rg.getInterpolationOrder() );
   if (rg.periodicX() ) {f->setPeriodicX();}
   if (rg.periodicY() ) {f->setPeriodicY();}
   if (rg.periodicZ() ) {f->setPeriodicZ();}
   return ScalarGrid(f);
}

ScalarGrid lux::makeGrid( const GridBox& rg, float defValue )
{
   SGrid<float>* f = new SGrid<float>();
   f->setDefVal( defValue );
   f->init( rg->nx(), rg->ny(), rg->nz(), rg->Lx(), rg->Ly(), rg->Lz(), rg->llc() );
   f->setInterpolationOrder( rg->getInterpolationOrder() );
   if (rg->periodicX() ) {f->setPeriodicX();}
   if (rg->periodicY() ) {f->setPeriodicY();}
   if (rg->periodicZ() ) {f->setPeriodicZ();}
   return ScalarGrid(f);
}

ScalarGrid lux::makeGrid( const GridBox& rg, float defValue, int psize )
{
   SGrid<float>* f = new SGrid<float>(psize);
   f->setDefVal( defValue );
   f->init( rg->nx(), rg->ny(), rg->nz(), rg->Lx(), rg->Ly(), rg->Lz(), rg->llc() );
   f->setInterpolationOrder( rg->getInterpolationOrder() );
   if (rg->periodicX() ) {f->setPeriodicX();}
   if (rg->periodicY() ) {f->setPeriodicY();}
   if (rg->periodicZ() ) {f->setPeriodicZ();}
   return ScalarGrid(f);
}


VectorGrid lux::makeVectorGrid( const RectangularGrid& rg, Vector defValue ) { return makeGrid( rg, defValue); }

VectorGrid lux::makeGrid( const RectangularGrid& rg, Vector defValue )
{
   SGrid<Vector>* f = new SGrid<Vector>();
   f->setDefVal( defValue );
   f->init( rg.nx(), rg.ny(), rg.nz(), rg.Lx(), rg.Ly(), rg.Lz(), rg.llc() );
   f->setInterpolationOrder( rg.getInterpolationOrder() );
   if (rg.periodicX() ) {f->setPeriodicX();}
   if (rg.periodicY() ) {f->setPeriodicY();}
   if (rg.periodicZ() ) {f->setPeriodicZ();}
   return VectorGrid(f);
}




VectorGrid lux::makeGrid( const GridBox& rg, Vector defValue )
{
   SGrid<Vector>* f = new SGrid<Vector>();
   f->setDefVal( defValue );
   f->init( rg->nx(), rg->ny(), rg->nz(), rg->Lx(), rg->Ly(), rg->Lz(), rg->llc() );
   f->setInterpolationOrder( rg->getInterpolationOrder() );
   if (rg->periodicX() ) {f->setPeriodicX();}
   if (rg->periodicY() ) {f->setPeriodicY();}
   if (rg->periodicZ() ) {f->setPeriodicZ();}
   return VectorGrid(f);
}

VectorGrid lux::makeGrid( const GridBox& rg, Vector defValue, int psize )
{
   SGrid<Vector>* f = new SGrid<Vector>(psize);
   f->setDefVal( defValue );
   f->init( rg->nx(), rg->ny(), rg->nz(), rg->Lx(), rg->Ly(), rg->Lz(), rg->llc() );
   f->setInterpolationOrder( rg->getInterpolationOrder() );
   if (rg->periodicX() ) {f->setPeriodicX();}
   if (rg->periodicY() ) {f->setPeriodicY();}
   if (rg->periodicZ() ) {f->setPeriodicZ();}
   return VectorGrid(f);
}






ColorGrid lux::makeColorGrid( const RectangularGrid& rg, Color defValue ) { return makeGrid( rg, defValue); }

ColorGrid lux::makeGrid( const RectangularGrid& rg, Color defValue )
{
   SGrid<Color>* f = new SGrid<Color>();
   f->setDefVal( defValue );
   f->init( rg.nx(), rg.ny(), rg.nz(), rg.Lx(), rg.Ly(), rg.Lz(), rg.llc() );
   f->setInterpolationOrder( rg.getInterpolationOrder() );
   if (rg.periodicX() ) {f->setPeriodicX();}
   if (rg.periodicY() ) {f->setPeriodicY();}
   if (rg.periodicZ() ) {f->setPeriodicZ();}
   return ColorGrid(f);
}

ColorGrid lux::makeGrid( const GridBox& rg, Color defValue )
{
   SGrid<Color>* f = new SGrid<Color>();
   f->setDefVal( defValue );
   f->init( rg->nx(), rg->ny(), rg->nz(), rg->Lx(), rg->Ly(), rg->Lz(), rg->llc() );
   f->setInterpolationOrder( rg->getInterpolationOrder() );
   if (rg->periodicX() ) {f->setPeriodicX();}
   if (rg->periodicY() ) {f->setPeriodicY();}
   if (rg->periodicZ() ) {f->setPeriodicZ();}
   return ColorGrid(f);
}

ColorGrid lux::makeGrid( const GridBox& rg, Color defValue, int psize )
{
   SGrid<Color>* f = new SGrid<Color>(psize);
   f->setDefVal( defValue );
   f->init( rg->nx(), rg->ny(), rg->nz(), rg->Lx(), rg->Ly(), rg->Lz(), rg->llc() );
   f->setInterpolationOrder( rg->getInterpolationOrder() );
   if (rg->periodicX() ) {f->setPeriodicX();}
   if (rg->periodicY() ) {f->setPeriodicY();}
   if (rg->periodicZ() ) {f->setPeriodicZ();}
   return ColorGrid(f);
}



MatrixGrid lux::makeMatrixGrid( const RectangularGrid& rg, Matrix defValue ) { return makeGrid( rg, defValue); }

MatrixGrid lux::makeGrid( const RectangularGrid& rg, Matrix defValue )
{
   SGrid<Matrix>* f = new SGrid<Matrix>();
   f->setDefVal( defValue );
   f->init( rg.nx(), rg.ny(), rg.nz(), rg.Lx(), rg.Ly(), rg.Lz(), rg.llc() );
   f->setInterpolationOrder( rg.getInterpolationOrder() );
   if (rg.periodicX() ) {f->setPeriodicX();}
   if (rg.periodicY() ) {f->setPeriodicY();}
   if (rg.periodicZ() ) {f->setPeriodicZ();}
   return MatrixGrid(f);
}

MatrixGrid lux::makeGrid( const GridBox& rg, Matrix defValue )
{
   SGrid<Matrix>* f = new SGrid<Matrix>();
   f->setDefVal( defValue );
   f->init( rg->nx(), rg->ny(), rg->nz(), rg->Lx(), rg->Ly(), rg->Lz(), rg->llc() );
   f->setInterpolationOrder( rg->getInterpolationOrder() );
   if (rg->periodicX() ) {f->setPeriodicX();}
   if (rg->periodicY() ) {f->setPeriodicY();}
   if (rg->periodicZ() ) {f->setPeriodicZ();}
   return MatrixGrid(f);
}


MatrixGrid lux::makeGrid( const GridBox& rg, Matrix defValue, int psize )
{
   SGrid<Matrix>* f = new SGrid<Matrix>(psize);
   f->setDefVal( defValue );
   f->init( rg->nx(), rg->ny(), rg->nz(), rg->Lx(), rg->Ly(), rg->Lz(), rg->llc() );
   f->setInterpolationOrder( rg->getInterpolationOrder() );
   if (rg->periodicX() ) {f->setPeriodicX();}
   if (rg->periodicY() ) {f->setPeriodicY();}
   if (rg->periodicZ() ) {f->setPeriodicZ();}
   return MatrixGrid(f);
}







ScalarGrid::ScalarGrid( SGrid<float>* f ) : ScalarGridBase(f) {}


ScalarGrid::~ScalarGrid(){}


ScalarGrid ScalarGrid::operator+=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) +  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

ScalarGrid ScalarGrid::operator-=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) -  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

ScalarGrid ScalarGrid::operator-()
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  -( (*this)->get(i,j,k) ) );
	 }
      }
   }
   return *this;
}

ScalarGrid ScalarGrid::operator*=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) *  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

ScalarGrid ScalarGrid::operator/=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) /  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}









VectorGrid::VectorGrid( SGrid<Vector>* f ) : VectorGridBase(f) {}


VectorGrid::~VectorGrid(){}


VectorGrid VectorGrid::operator+=( const VectorField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) +  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

VectorGrid VectorGrid::operator-=( const VectorField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) -  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

VectorGrid VectorGrid::operator-()
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) *  (-1.0)  );
	 }
      }
   }
   return *this;
}

VectorGrid VectorGrid::operator*=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) *  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

VectorGrid VectorGrid::operator/=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) /  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}











MatrixGrid::MatrixGrid( SGrid<Matrix>* f ) : MatrixGridBase(f) {}


MatrixGrid::~MatrixGrid(){}



MatrixGrid MatrixGrid::operator+=( const MatrixField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) +  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

MatrixGrid MatrixGrid::operator-=( const MatrixField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) -  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

MatrixGrid MatrixGrid::operator-()
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k)*(-1.0)  );
	 }
      }
   }
   return *this;
}

MatrixGrid MatrixGrid::operator*=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) *  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

MatrixGrid MatrixGrid::operator/=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) /  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}









ColorGrid::ColorGrid( SGrid<Color>* f ) : ColorGridBase(f) {}


ColorGrid::~ColorGrid(){}


ColorGrid ColorGrid::operator+=( const ColorField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) +  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

ColorGrid ColorGrid::operator-=( const ColorField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) -  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

ColorGrid ColorGrid::operator-()
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k)*(-1.0)  );
	 }
      }
   }
   return *this;
}

ColorGrid ColorGrid::operator*=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) *  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

ColorGrid ColorGrid::operator/=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) /  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}





void lux::readScalarGrid(ScalarGrid& grid, string fname ){ readGrid( grid, fname ); }

void lux::readGrid( ScalarGrid& grid, string fname)
{
   ifstream out;
   out.open( fname.c_str() );
   int nx;
   int ny;
   int nz;

   int i, j, k;

   float x;
   float y;
   float z;

   float Lx;
   float Ly;
   float Lz;

   float value;

   out.read( (char*)&nx, sizeof(nx) );
   out.read( (char*)&ny, sizeof(ny) );
   out.read( (char*)&nz, sizeof(nz) );
   out.read( (char*)&Lx, sizeof(Lx) );
   out.read( (char*)&Ly, sizeof(Ly) );
   out.read( (char*)&Lz, sizeof(Lz) );

   out.read( (char*)&x, sizeof(x) );
   out.read( (char*)&y, sizeof(y) );
   out.read( (char*)&z, sizeof(z) );
   out.read( (char*)&value, sizeof(value) );

   grid->setDefVal(value);
   grid->init( nx, ny, nz, Lx, Ly, Lz, Vector(x,y,z) );

   ProgressMeter metr( nx*ny*nz, "Read ScalarGrid");
   while(out.read((char *)&i, sizeof(int)))
   {
       out.read((char *)&j, sizeof(int));
       out.read((char *)&k, sizeof(int));
       out.read((char *)&value, sizeof(float));
       grid->set(i, j, k, value);
       metr.update();
   }
}

void lux::writeScalarGrid(ScalarGrid& grid, string fname ){ writeGrid( grid, fname ); }

void lux::writeGrid(ScalarGrid& grid, string fname )
{
   ofstream out;
   out.open( fname.c_str() );
   int nx = grid->nx();
   int ny = grid->ny();
   int nz = grid->nz();

   Vector llc = grid->llc();

   float x = llc[0];
   float y = llc[1];
   float z = llc[2];

   float Lx = grid->Lx();
   float Ly = grid->Ly();
   float Lz = grid->Lz();

   out.write( (char*)&nx, sizeof(nx) );
   out.write( (char*)&ny, sizeof(ny) );
   out.write( (char*)&nz, sizeof(nz) );

   out.write( (char*)&Lx, sizeof(Lx) );
   out.write( (char*)&Ly, sizeof(Ly) );
   out.write( (char*)&Lz, sizeof(Lz) );
   out.write( (char*)&x, sizeof(x) );
   out.write( (char*)&y, sizeof(y) );
   out.write( (char*)&z, sizeof(z) );

   float temp = grid->getDefVal();
   out.write( (char*)&temp, sizeof(float) );

   ProgressMeter metr( nx*ny*nz, "Write ScalarGrid");
   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
            temp = grid->get(i,j,k);
            if(temp != grid->getDefVal()){
               out.write((char*)&i, sizeof(int));
               out.write((char*)&j, sizeof(int));
               out.write((char*)&k, sizeof(int));
               out.write((char*)&temp, sizeof(float));
            }
	    metr.update();
         }
      }
   }
   out.close();
}














void lux::readVectorGrid(VectorGrid& grid, string fname ){ readGrid( grid, fname ); }

void lux::readGrid( VectorGrid& grid, string fname)
{
   ifstream out;
   out.open( fname.c_str() );
   int nx;
   int ny;
   int nz;

   int i, j, k;

   float x;
   float y;
   float z;

   float Lx;
   float Ly;
   float Lz;

   Vector value;

   out.read( (char*)&nx, sizeof(nx) );
   out.read( (char*)&ny, sizeof(ny) );
   out.read( (char*)&nz, sizeof(nz) );
   out.read( (char*)&Lx, sizeof(Lx) );
   out.read( (char*)&Ly, sizeof(Ly) );
   out.read( (char*)&Lz, sizeof(Lz) );

   out.read( (char*)&x, sizeof(x) );
   out.read( (char*)&y, sizeof(y) );
   out.read( (char*)&z, sizeof(z) );
   out.read( (char*)&value, sizeof(value) );

   grid->setDefVal(value);
   grid->init( nx, ny, nz, Lx, Ly, Lz, Vector(x,y,z) );

   ProgressMeter metr( nx*ny*nz, "Read VectorGrid");
   while(out.read((char *)&i, sizeof(int)))
   {
       out.read((char *)&j, sizeof(int));
       out.read((char *)&k, sizeof(int));
       out.read((char *)&value, sizeof(float));
       grid->set(i, j, k, value);
       metr.update();
   }
}

void lux::writeVectorGrid(VectorGrid& grid, string fname ){ writeGrid( grid, fname ); }
void lux::writeGrid(VectorGrid& grid, string fname )
{
   ofstream out;
   out.open( fname.c_str() );
   int nx = grid->nx();
   int ny = grid->ny();
   int nz = grid->nz();


   Vector llc = grid->llc();

   float x = llc[0];
   float y = llc[1];
   float z = llc[2];

   float Lx = grid->Lx();
   float Ly = grid->Ly();
   float Lz = grid->Lz();

   out.write( (char*)&nx, sizeof(nx) );
   out.write( (char*)&ny, sizeof(ny) );
   out.write( (char*)&nz, sizeof(nz) );

   out.write( (char*)&Lx, sizeof(Lx) );
   out.write( (char*)&Ly, sizeof(Ly) );
   out.write( (char*)&Lz, sizeof(Lz) );
   out.write( (char*)&x, sizeof(x) );
   out.write( (char*)&y, sizeof(y) );
   out.write( (char*)&z, sizeof(z) );

   Vector temp = grid->getDefVal();
   out.write( (char*)&temp[0], sizeof(float) );
   out.write( (char*)&temp[1], sizeof(float) );
   out.write( (char*)&temp[2], sizeof(float) );

   ProgressMeter metr( nx*ny*nz, "Write VectorGrid");
   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
            temp = grid->get(i,j,k);
            if(temp != grid->getDefVal())
	    {
               out.write((char*)&i, sizeof(int));
               out.write((char*)&j, sizeof(int));
               out.write((char*)&k, sizeof(int));
               out.write((char*)&temp[0], sizeof(double));
               out.write((char*)&temp[1], sizeof(double));
               out.write((char*)&temp[2], sizeof(double));
            }
	    metr.update();
         }
      }
   }
   out.close();
}











void lux::readColorGrid(ColorGrid& grid, string fname ){ readGrid( grid, fname ); }

void lux::readGrid( ColorGrid& grid, string fname)
{
   ifstream out;
   out.open( fname.c_str() );
   int nx;
   int ny;
   int nz;

   int i, j, k;

   float x;
   float y;
   float z;

   float Lx;
   float Ly;
   float Lz;

   Color value;

   out.read( (char*)&nx, sizeof(nx) );
   out.read( (char*)&ny, sizeof(ny) );
   out.read( (char*)&nz, sizeof(nz) );
   out.read( (char*)&Lx, sizeof(Lx) );
   out.read( (char*)&Ly, sizeof(Ly) );
   out.read( (char*)&Lz, sizeof(Lz) );

   out.read( (char*)&x, sizeof(x) );
   out.read( (char*)&y, sizeof(y) );
   out.read( (char*)&z, sizeof(z) );
   out.read( (char*)&value, sizeof(value) );

   grid->setDefVal(value);
   grid->init( nx, ny, nz, Lx, Ly, Lz, Vector(x,y,z) );

   ProgressMeter metr( nx*ny*nz, "Read ColorGrid");
   while(out.read((char *)&i, sizeof(int)))
   {
       out.read((char *)&j, sizeof(int));
       out.read((char *)&k, sizeof(int));
       out.read((char *)&value, sizeof(value));
       grid->set(i, j, k, value);
       metr.update();
   }
}

void lux::writeColorGrid(ColorGrid& grid, string fname ){ writeGrid( grid, fname ); }
void lux::writeGrid(ColorGrid& grid, string fname )
{
   ofstream out;
   out.open( fname.c_str() );
   int nx = grid->nx();
   int ny = grid->ny();
   int nz = grid->nz();

   Vector llc = grid->llc();

   float x = llc[0];
   float y = llc[1];
   float z = llc[2];

   float Lx = grid->Lx();
   float Ly = grid->Ly();
   float Lz = grid->Lz();

   out.write( (char*)&nx, sizeof(nx) );
   out.write( (char*)&ny, sizeof(ny) );
   out.write( (char*)&nz, sizeof(nz) );

   out.write( (char*)&Lx, sizeof(Lx) );
   out.write( (char*)&Ly, sizeof(Ly) );
   out.write( (char*)&Lz, sizeof(Lz) );
   out.write( (char*)&x, sizeof(x) );
   out.write( (char*)&y, sizeof(y) );
   out.write( (char*)&z, sizeof(z) );

   Color temp = grid->getDefVal();
   out.write( (char*)&temp[0], sizeof(double) );
   out.write( (char*)&temp[1], sizeof(double) );
   out.write( (char*)&temp[2], sizeof(double) );
   out.write( (char*)&temp[3], sizeof(double) );

   ProgressMeter metr( nx*ny*nz, "Write ColorGrid");
   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
            temp = grid->get(i,j,k);
            if(temp != grid->getDefVal()){
               out.write((char*)&i, sizeof(int));
               out.write((char*)&j, sizeof(int));
               out.write((char*)&k, sizeof(int));
               out.write((char*)&temp[0], sizeof(double));
               out.write((char*)&temp[1], sizeof(double));
               out.write((char*)&temp[2], sizeof(double));
               out.write((char*)&temp[3], sizeof(double));
            }
	    metr.update();
         }
      }
   }
}











void lux::Blur( ScalarGrid& g )
{
   ScalarGrid temp = makeGrid( *(g.get()), g->getDefVal() );
   ProgressMeter meter( (g->nx())*(g->ny())*(g->nz()) * 2, "Blur" );
   for( int k=0;k<g->nz();k++ )
   {
      for( int j=0;j<g->ny();j++ )
      {
         for( int i=0;i<g->nx();i++ )
         {
	    temp->set( i,j,k, g->get(i,j,k) );
	    meter.update();
         }
      }
   }


   for( int k=0;k<g->nz();k++ )
   {
      int kmin = k-1;
      if( kmin < 0 ){ kmin = k; }
      int kmax = k+1;
      if( kmax >= g->nz() ){ kmax = k;}
      for( int j=0;j<g->ny();j++ )
      {
         int jmin = j-1;
         if( jmin < 0 ){ jmin = j; }
         int jmax = j+1;
         if( jmax >= g->ny() ){ jmax = j;}
         for( int i=0;i<g->nx();i++ )
         {
            int imin = i-1;
            if( imin < 0 ){ imin = i; }
            int imax = i+1;
            if( imax >= g->nx() ){ imax = i;}
	    float sum = 0;
	    int nb = 0;
	    for( int kk=kmin;kk<=kmax;kk++ )
	    {
	       for( int jj=jmin;jj<=jmax;jj++ )
	       {
	          for( int ii=imin;ii<=imax;ii++ )
	          {
		     sum += temp->get( ii, jj, kk );
		     ++nb;
	          }
	       }
	    }
	    sum = sum/nb;
	    g->set(i,j,k, sum);
	    meter.update();
         }
      }
   }
}




void lux::Blur( VectorGrid& g )
{
   VectorGrid temp = makeGrid( *(g.get()), g->getDefVal() );
   ProgressMeter meter( (g->nx())*(g->ny())*(g->nz()) * 2, "Blur" );
   for( int k=0;k<g->nz();k++ )
   {
      for( int j=0;j<g->ny();j++ )
      {
         for( int i=0;i<g->nx();i++ )
         {
	    temp->set( i,j,k, g->get(i,j,k) );
	    meter.update();
         }
      }
   }


   for( int k=0;k<g->nz();k++ )
   {
      int kmin = k-1;
      if( kmin < 0 ){ kmin = k; }
      int kmax = k+1;
      if( kmax >= g->nz() ){ kmax = k;}
      for( int j=0;j<g->ny();j++ )
      {
         int jmin = j-1;
         if( jmin < 0 ){ jmin = j; }
         int jmax = j+1;
         if( jmax >= g->ny() ){ jmax = j;}
         for( int i=0;i<g->nx();i++ )
         {
            int imin = i-1;
            if( imin < 0 ){ imin = i; }
            int imax = i+1;
            if( imax >= g->nx() ){ imax = i;}
	    Vector sum;
	    int nb = 0;
	    for( int kk=kmin;kk<=kmax;kk++ )
	    {
	       for( int jj=jmin;jj<=jmax;jj++ )
	       {
	          for( int ii=imin;ii<=imax;ii++ )
	          {
		     sum += temp->get( ii, jj, kk );
		     ++nb;
	          }
	       }
	    }
	    sum = sum/nb;
	    g->set(i,j,k, sum);
	    meter.update();
         }
      }
   }
}




void lux::Blur( ColorGrid& g )
{
   ColorGrid temp = makeGrid( *(g.get()), g->getDefVal() );
   ProgressMeter meter( (g->nx())*(g->ny())*(g->nz()) * 2, "Blur" );
   for( int k=0;k<g->nz();k++ )
   {
      for( int j=0;j<g->ny();j++ )
      {
         for( int i=0;i<g->nx();i++ )
         {
	    temp->set( i,j,k, g->get(i,j,k) );
	    meter.update();
         }
      }
   }


   for( int k=0;k<g->nz();k++ )
   {
      int kmin = k-1;
      if( kmin < 0 ){ kmin = k; }
      int kmax = k+1;
      if( kmax >= g->nz() ){ kmax = k;}
      for( int j=0;j<g->ny();j++ )
      {
         int jmin = j-1;
         if( jmin < 0 ){ jmin = j; }
         int jmax = j+1;
         if( jmax >= g->ny() ){ jmax = j;}
         for( int i=0;i<g->nx();i++ )
         {
            int imin = i-1;
            if( imin < 0 ){ imin = i; }
            int imax = i+1;
            if( imax >= g->nx() ){ imax = i;}
	    Color sum;
	    int nb = 0;
	    for( int kk=kmin;kk<=kmax;kk++ )
	    {
	       for( int jj=jmin;jj<=jmax;jj++ )
	       {
	          for( int ii=imin;ii<=imax;ii++ )
	          {
		     sum += temp->get( ii, jj, kk );
		     ++nb;
	          }
	       }
	    }
	    sum = sum/nb;
	    g->set(i,j,k, sum);
	    meter.update();
         }
      }
   }
}






void lux::stampNoise( ScalarGrid& g, Noise* noise, const Vector& center, const float radius )
{
   int nx = g->nx();
   int ny = g->ny();
   int nz = g->nz();

   ProgressMeter meter( nx*ny*nz, "stampNoise" );

   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
	    Vector P = g->evalP(i,j,k);
	    float distance = (P-center).magnitude();
	    if( distance <= radius )
	    {
	       float nv = noise->eval( P );
	       if( nv > g->get(i,j,k) ) { g->set(i,j,k,nv); }
	    }
	    meter.update();
         }
      }
   }
}



void lux::stampNoise( ScalarGrid& g, Noise* noise, const Vector& center, const float radius, const Vector& vel, const Vector& accel, const float timestep, const int seed )
{
   int nx = g->nx();
   int ny = g->ny();
   int nz = g->nz();

   float dx = g->dx();
   dx = ( dx < g->dy() ) ? dx : g->dy();
   dx = ( dx < g->dz() ) ? dx : g->dz();

   float approxDistance = timestep * vel.magnitude() + 0.5*timestep*timestep*accel.magnitude();
   int nbsamples = 1 + (int)(  approxDistance/dx  );



   UniformPRN prn;
   Noise_t parms;
   parms.seed = seed + 646383;
   prn.setParameters(parms);

   int ii,jj,kk;

   ProgressMeter meter( nx*ny*nz, "StampNoise" );
   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
	    Vector P = g->evalP(i,j,k);
	    float distance = (P-center).magnitude();
	    if( distance <= radius )
	    {
	       float nv = noise->eval( P ) / nbsamples;

	       // accumulate along the streak
	       if( nv > 0 )
	       {
               for( int q=0;q<nbsamples;q++ )
	       {
	          float local = prn.eval() * timestep;
		  Vector PP = P + local*vel + (0.5*local*local)*accel;
                  if( g->getGridIndex( PP, ii,jj,kk ) ) { g->set(ii,jj,kk, g->get(ii,jj,kk) + nv ); }
	       }
	       }
	    }
	    meter.update();
         }
      }
   }
}

void lux::stampNoise( ScalarGrid& g, Noise* noise, const ParticleGroupA& particles )
{
   Noise_t parms;
   for( size_t p=0;p<particles.size();p++ )
   {
      const Vector& center = particles[p].P();
      const float& radius = particles[p].pscale();
      parms.octaves = particles[p].octaves();
      parms.roughness = particles[p].roughness();
      parms.frequency = particles[p].freq();
      parms.fjump = particles[p].fjump();
      parms.offset = particles[p].offset();
      parms.translate = particles[p].translate();
      parms.translate = particles[p].translate();
      noise->setParameters( parms );
      if( (particles[p].v().magnitude() == 0 && particles[p].accel().magnitude() == 0) || particles[p].shutter() == 0 )
      {
         stampNoise( g, noise, center, radius );
      }
      else
      {
         float timestep = particles[p].shutter()/24.0;
         stampNoise( g, noise, center, radius, particles[p].v(), particles[p].accel(), timestep, particles[p].id() );
      }
   }
}



void lux::stampPointWisps( ScalarGrid& g, const ParticleGroupA& particles )
{

   PointWispWanderer pww;
   int metertotal = 0;

   float dx = g->dx();
   if( dx > g->dy() ){ dx = g->dy(); }
   if( dx > g->dz() ){ dx = g->dz(); }

   ProgressMeter meter( particles.size(), "StampPointWisps");
   for( size_t p=0;p<particles.size();p++ )
   {
      pww.reset(particles[p]);
      int nb = particles[p].nbWisps();
      Vector V = particles[p].v();
      Vector A = particles[p].accel();
      float opacity = particles[p].opacity();
      float shutter = particles[p].shutter();
      int id = particles[p].id();
      const float timestep = shutter/24.0;
      for( int i=0;i<nb;i++ )
      {
            const Vector P = pww.step();
            // now stamp a dot into the grid
	    stampBlurredWisps( g, P, timestep, V, A, opacity, id );

      }
      meter.update();
   }
}


void lux::stampBlurredWisps( ScalarGrid& g, const Vector& P, const float timestep, const Vector& velocity, const Vector& acceleration, const float opacity, int seed )
{
   float dx = g->dx();
   dx = ( dx < g->dy() ) ? dx : g->dy();
   dx = ( dx < g->dz() ) ? dx : g->dz();

   float approxDistance = timestep * velocity.magnitude() + 0.5*timestep*timestep*acceleration.magnitude();
   int nbsamples = 1 + (int)(approxDistance/dx);

   UniformPRN rn;
   if( nbsamples > 1 )
   {
      Noise_t parms;
      parms.seed = seed;
      rn.setParameters(parms);
   }

   int ix,iy,iz;

   float gvalue = opacity/nbsamples;
   for( int q=0;q<nbsamples;q++ )
   {
      float local = rn.eval() * timestep;
      if( q==0 ){ local = 0; }

      Vector position = P + local * velocity + (0.5*local*local)*acceleration;

      if( g->getGridIndex( position, ix,iy,iz ) )
      {
         Vector P000 = g->evalP(ix,iy,iz);
         position -= P000;

         float wx = position[0]/g->dx();
         float wy = position[1]/g->dy();
         float wz = position[2]/g->dz();

         float vx = 1.0-wx;
         float vy = 1.0-wy;
         float vz = 1.0-wz;

         float w000 = vx*vy*vz;
         float w100 = wx*vy*vz;
         float w010 = vx*wy*vz;
         float w001 = vx*vy*wz;
         float w110 = wx*wy*vz;
         float w011 = vx*wy*wz;
         float w101 = wx*vy*wz;
         float w111 = wx*wy*wz;

         bool tx = ( ix < g->nx()-1 );
         bool ty = ( iy < g->ny()-1 );
         bool tz = ( iz < g->nz()-1 );

         g->set(ix,iy,iz, g->get(ix,iy,iz) + gvalue*w000 );
         if( tx ){ g->set(ix+1,iy,iz, g->get(ix+1,iy,iz) + gvalue*w100 ); }
         if( ty ){ g->set(ix,iy+1,iz, g->get(ix,iy+1,iz) + gvalue*w010 ); }
         if( tz ){ g->set(ix,iy,iz+1, g->get(ix,iy,iz+1) + gvalue*w001 ); }
         if( tx && ty ){ g->set(ix+1,iy+1,iz, g->get(ix+1,iy+1,iz) + gvalue*w110 ); }
         if( tx && tz ){ g->set(ix+1,iy,iz+1, g->get(ix+1,iy,iz+1) + gvalue*w101 ); }
         if( ty && tz ){ g->set(ix,iy+1,iz+1, g->get(ix,iy+1,iz+1) + gvalue*w011 ); }
         if( tx && ty && tz ){ g->set(ix+1,iy+1,iz+1, g->get(ix+1,iy+1,iz+1) + gvalue*w111 ); }
      }
   }
}

void lux::stampBlurredWisps( ScalarGrid& g, ColorGrid& cg, const Vector& P, const float timestep, const Vector& velocity, const Vector& acceleration, const float opacity, const Color& cd, int seed )
{
   float dx = g->dx();
   dx = ( dx < g->dy() ) ? dx : g->dy();
   dx = ( dx < g->dz() ) ? dx : g->dz();

   float approxDistance = timestep * velocity.magnitude() + 0.5*timestep*timestep*acceleration.magnitude();
   int nbsamples = (1 + (int)(  approxDistance/dx  ))*10;

   UniformPRN rn;
   if( nbsamples > 1 )
   {
      Noise_t parms;
      parms.seed = seed;
      rn.setParameters(parms);
   }

   int ix,iy,iz;

   float gvalue = opacity/nbsamples;
   Color cvalue = cd/nbsamples;
   for( int q=0;q<nbsamples;q++ )
   {
      float local = rn.eval() * timestep;
      if( q==0 ){ local = 0; }

      Vector position = P + local * velocity + (0.5*local*local)*acceleration;

      if( g->getGridIndex( position, ix,iy,iz ) )
      {
         Vector P000 = g->evalP(ix,iy,iz);
         position -= P000;

         float wx = position[0]/g->dx();
         float wy = position[1]/g->dy();
         float wz = position[2]/g->dz();

         float vx = 1.0-wx;
         float vy = 1.0-wy;
         float vz = 1.0-wz;

         float w000 = vx*vy*vz;
         float w100 = wx*vy*vz;
         float w010 = vx*wy*vz;
         float w001 = vx*vy*wz;
         float w110 = wx*wy*vz;
         float w011 = vx*wy*wz;
         float w101 = wx*vy*wz;
         float w111 = wx*wy*wz;

         bool tx = ( ix < g->nx()-1 );
         bool ty = ( iy < g->ny()-1 );
         bool tz = ( iz < g->nz()-1 );

         g->set(ix,iy,iz, g->get(ix,iy,iz) + gvalue*w000 );
         cg->set(ix,iy,iz, cg->get(ix,iy,iz) + cvalue*w000 );
         if( tx ){ g->set(ix+1,iy,iz, g->get(ix+1,iy,iz) + gvalue*w100 ); cg->set(ix+1,iy,iz, cg->get(ix+1,iy,iz) + cvalue*w100 );    }
         if( ty ){ g->set(ix,iy+1,iz, g->get(ix,iy+1,iz) + gvalue*w010 ); cg->set(ix,iy+1,iz, cg->get(ix,iy+1,iz) + cvalue*w010 ); }
         if( tz ){ g->set(ix,iy,iz+1, g->get(ix,iy,iz+1) + gvalue*w001 ); }
         if( tx && ty ){ g->set(ix+1,iy+1,iz, g->get(ix+1,iy+1,iz) + gvalue*w110 ); }
         if( tx && tz ){ g->set(ix+1,iy,iz+1, g->get(ix+1,iy,iz+1) + gvalue*w101 ); }
         if( ty && tz ){ g->set(ix,iy+1,iz+1, g->get(ix,iy+1,iz+1) + gvalue*w011 ); }
         if( tx && ty && tz ){ g->set(ix+1,iy+1,iz+1, g->get(ix+1,iy+1,iz+1) + gvalue*w111 ); }
      }
   }
}



void lux::stamp( ScalarGrid& g, const ScalarField& field, const int nbsamples )
{
   int nx = g->nx();
   int ny = g->ny();
   int nz = g->nz();
   float dx = g->dx();
   float dy = g->dy();
   float dz = g->dz();

   UniformPRN rng;
   Noise_t p;
   p.seed = 474738;
   rng.setParameters(p);

   ProgressMeter meter( nz, "StampField" );

   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
	    Vector X = g->evalP(i,j,k);
	    float value = field->eval( X );
	    for( int s=1;s<nbsamples;s++ )
	    {
	       Vector P = X + Vector( dx*(rng.eval()), dy*(rng.eval()), dz*(rng.eval()) );
	       value += field->eval( P );
	    }
	    if( nbsamples > 1 ){ value /= nbsamples; }
	    g->set(i,j,k, value );
         }
      }
      meter.update();
   }
}



void lux::stamp( VectorGrid& g, const VectorField& field, const int nbsamples )
{
   int nx = g->nx();
   int ny = g->ny();
   int nz = g->nz();

   ProgressMeter meter( nz, "StampField" );

   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
	    g->set(i,j,k, field->eval( g->evalP(i,j,k) ) );
         }
      }
      meter.update();
   }
}



void lux::stamp( ColorGrid& g, const ColorField& field, const int nbsamples )
{
   int nx = g->nx();
   int ny = g->ny();
   int nz = g->nz();

   ProgressMeter meter( nz, "StampField" );

   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
	    g->set(i,j,k, field->eval( g->evalP(i,j,k) ) );
         }
      }
      meter.update();
   }
}


void lux::stamp( MatrixGrid& g, const MatrixField& field, const int nbsamples )
{
   int nx = g->nx();
   int ny = g->ny();
   int nz = g->nz();

   ProgressMeter meter( nz, "StampField" );

   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
         for( int i=0;i<nx;i++ )
         {
	    g->set(i,j,k, field->eval( g->evalP(i,j,k) ) );
         }
      }
      meter.update();
   }
}






void lux::stampParticles( ScalarGrid& g, const ParticleGroupA& particles, const float timestep )
{
   float dx = g->dx();
   dx = ( dx < g->dy() ) ? dx : g->dy();
   dx = ( dx < g->dz() ) ? dx : g->dz();

   UniformPRN rn;

   for( size_t p=0;p<particles.size();p++ )
   {
      const Vector& V = particles[p].v();
      const Vector& A = particles[p].accel();
      const Vector& P = particles[p].P();
      float approxDistance = timestep * V.magnitude() + 0.5*timestep*timestep*A.magnitude();
      int nbsamples = 1 + (int)(  approxDistance/dx  );

      if( nbsamples > 1 )
      {
         Noise_t parms;
         parms.seed = particles[p].id();
         rn.setParameters(parms);
      }

      float dt = timestep/nbsamples;

      int ix,iy,iz;

      float gvalue = particles[p].opacity()/nbsamples;
      for( int q=0;q<nbsamples;q++ )
      {
         float local = q * dt;
         Vector position = P + local*V + (0.5*local*local)*A;
         if( g->getGridIndex( position, ix,iy,iz ) )
         {
         Vector P000 = g->evalP(ix,iy,iz);
         position -= P000;

         float wx = position[0]/g->dx();
         float wy = position[1]/g->dy();
         float wz = position[2]/g->dz();

         float vx = 1.0-wx;
         float vy = 1.0-wy;
         float vz = 1.0-wz;

         float w000 = vx*vy*vz;
         float w100 = wx*vy*vz;
         float w010 = vx*wy*vz;
         float w001 = vx*vy*wz;
         float w110 = wx*wy*vz;
         float w011 = vx*wy*wz;
         float w101 = wx*vy*wz;
         float w111 = wx*wy*wz;

         bool tx = ( ix < g->nx()-1 );
         bool ty = ( iy < g->ny()-1 );
         bool tz = ( iz < g->nz()-1 );

         g->set(ix,iy,iz, g->get(ix,iy,iz) + gvalue*w000 );
         if( tx ){ g->set(ix+1,iy,iz, g->get(ix+1,iy,iz) + gvalue*w100 ); }
         if( ty ){ g->set(ix,iy+1,iz, g->get(ix,iy+1,iz) + gvalue*w010 ); }
         if( tz ){ g->set(ix,iy,iz+1, g->get(ix,iy,iz+1) + gvalue*w001 ); }
         if( tx && ty ){ g->set(ix+1,iy+1,iz, g->get(ix+1,iy+1,iz) + gvalue*w110 ); }
         if( tx && tz ){ g->set(ix+1,iy,iz+1, g->get(ix+1,iy,iz+1) + gvalue*w101 ); }
         if( ty && tz ){ g->set(ix,iy+1,iz+1, g->get(ix,iy+1,iz+1) + gvalue*w011 ); }
         if( tx && ty && tz ){ g->set(ix+1,iy+1,iz+1, g->get(ix+1,iy+1,iz+1) + gvalue*w111 ); }
         }
      }
   }
}



int lux::getNbAvailablePartitions( const ScalarGrid& g ){ return g->NbPartitions(); }
int lux::getNbUsedPartitions( const ScalarGrid& g ) { return g->Size(); }
int lux::getPartitionSize( const ScalarGrid& g ) { return g->blockSize(); }




void lux::GreenConvolve( const MatrixField& m, VectorGrid& result )
{
   int nx = result->nx();
   int ny = result->ny();
   int nz = result->nz();
   double scalefactor = result->dx() * result->dy() * result->dz()/( 4.0*M_PI );
   ProgressMeter metr( nz, "GreenConvolve" );
   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
#pragma omp parallel for
         for( int i=0;i<nx;i++ )
	 {
	    Vector X = result->evalP(i,j,k);
	    Vector accum;
            for( int kk=0;kk<nz;kk++ )
            {
               for( int jj=0;jj<ny;jj++ )
               {
                  for( int ii=0;ii<nx;ii++ )
	          {
		     Vector Y = result->evalP(ii,jj,kk);
		     Vector D = X-Y;
		     double Dmagsq = D*D;
		     if( Dmagsq > 0.0 )
		     {
		        Matrix M = m->eval(Y);
			Vector daccum = D*M;
			daccum /= std::pow( Dmagsq, 1.5 );
			accum += daccum;
		     }
	          }
	       }
            }
            result->set(i,j,k,accum*scalefactor);
	 }
      }
      metr.update();
   }
}


void lux::GreenSurfaceConvolve( const VectorField& XX, VectorGrid& result )
{
   int nx = result->nx();
   int ny = result->ny();
   int nz = result->nz();
   double scalefactorx = result->dy() * result->dz()/( 4.0*M_PI );
   double scalefactory = result->dx() * result->dz()/( 4.0*M_PI );
   double scalefactorz = result->dx() * result->dy()/( 4.0*M_PI );
   ProgressMeter metr( nz, "GreenSurfaceConvolve" );
   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
#pragma omp parallel for
         for( int i=0;i<nx;i++ )
	 {
	    Vector X = result->evalP(i,j,k);
	    Vector accum;
            for( int kk=0;kk<nz;kk++ )
            {
               for( int jj=0;jj<ny;jj++ )
               {
                  int ii=0;
	          {
		     Vector Y = result->evalP(ii,jj,kk);
		     Vector D = X-Y;
		     double Dmagsq = D*D;
		     if( Dmagsq > 0.0 )
		     {
		        Vector XXX = XX->eval(Y);
			Vector daccum = D[0]*XXX;
			daccum /= Dmagsq*Dmagsq;
			accum += daccum*scalefactorx;
		     }
	          }
	       }
            }
	    for( int kk=0;kk<nz;kk++ )
            {
               for( int jj=0;jj<ny;jj++ )
               {
                  int ii=nx-1;
	          {
		     Vector Y = result->evalP(ii,jj,kk);
		     Vector D = X-Y;
		     double Dmagsq = D*D;
		     if( Dmagsq > 0.0 )
		     {
		        Vector XXX = XX->eval(Y);
			Vector daccum = -D[0]*XXX;
			daccum /= Dmagsq*Dmagsq;
			accum += daccum*scalefactorx;
		     }
	          }
	       }
            }

	    int kk=0;
            {
               for( int jj=0;jj<ny;jj++ )
               {
                  for( int ii=0;ii<nx;ii++ )
	          {
		     Vector Y = result->evalP(ii,jj,kk);
		     Vector D = X-Y;
		     double Dmagsq = D*D;
		     if( Dmagsq > 0.0 )
		     {
		        Vector XXX = XX->eval(Y);
			Vector daccum = D[2]*XXX;
			daccum /= Dmagsq*Dmagsq;
			accum += daccum*scalefactorz;
		     }
	          }
	       }
            }
	    kk=nz-1;
            {
               for( int jj=0;jj<ny;jj++ )
               {
                  for( int ii=0;ii<nx;ii++ )
	          {
		     Vector Y = result->evalP(ii,jj,kk);
		     Vector D = X-Y;
		     double Dmagsq = D*D;
		     if( Dmagsq > 0.0 )
		     {
		        Vector XXX = XX->eval(Y);
			Vector daccum = -D[2]*XXX;
			daccum /= Dmagsq*Dmagsq;
			accum += daccum*scalefactorz;
		     }
	          }
	       }
            }

	    for( int kk=0;kk<nz;kk++ )
            {
               int jj=0;
               {
                  for( int ii=0;ii<nx;ii++ )
	          {
		     Vector Y = result->evalP(ii,jj,kk);
		     Vector D = X-Y;
		     double Dmagsq = D*D;
		     if( Dmagsq > 0.0 )
		     {
		        Vector XXX = XX->eval(Y);
			Vector daccum = D[1]*XXX;
			daccum /= Dmagsq*Dmagsq;
			accum += daccum*scalefactory;
		     }
	          }
	       }
            }
	    for( int kk=0;kk<nz;kk++ )
            {
               int jj=ny-1;
               {
                  for( int ii=0;ii<nx;ii++ )
	          {
		     Vector Y = result->evalP(ii,jj,kk);
		     Vector D = X-Y;
		     double Dmagsq = D*D;
		     if( Dmagsq > 0.0 )
		     {
		        Vector XXX = XX->eval(Y);
			Vector daccum = -D[1]*XXX;
			daccum /= Dmagsq*Dmagsq;
			accum += daccum*scalefactory;
		     }
	          }
	       }
            }
            result->set(i,j,k,accum);
	 }
      }
      metr.update();
   }
}












ScalarFrustumGrid lux::makeScalarFrustumGrid( const FrustumGrid& rg, float defValue ) { return makeFrustumGrid( rg, defValue); }

ScalarFrustumGrid lux::makeFrustumGrid( const FrustumGrid& rg, float defValue )
{
   ScalarFrustumGrid f( new FSGrid<float>() );
   f->setDefVal( defValue );
   f->init( rg.nx(), rg.ny(), rg.nz(), rg.camera() );
   return f;
}



ScalarFrustumGrid lux::makeFrustumGrid( const FrustumBox& rg, float defValue )
{
   ScalarFrustumGrid f( new FSGrid<float>() );
   f->setDefVal( defValue );
   f->init( rg->nx(), rg->ny(), rg->nz(), rg->camera() );
   return f;
}


ScalarFrustumGrid lux::makeFrustumGrid( const FrustumBox& rg, float defValue, int psize )
{
   ScalarFrustumGrid f( new FSGrid<float>(psize) );
   f->setDefVal( defValue );
   f->init( rg->nx(), rg->ny(), rg->nz(), rg->camera() );
   return f;
}


VectorFrustumGrid lux::makeVectorFrustumGrid( const FrustumGrid& rg, Vector defValue ) { return makeFrustumGrid( rg, defValue); }

VectorFrustumGrid lux::makeFrustumGrid( const FrustumGrid& rg, Vector defValue )
{
   FSGrid<Vector>* f = new FSGrid<Vector>();
   f->setDefVal( defValue );
   f->init( rg.nx(), rg.ny(), rg.nz(), rg.camera() );
   return VectorFrustumGrid(f);
}


VectorFrustumGrid lux::makeFrustumGrid( const FrustumBox& rg, Vector defValue )
{
   FSGrid<Vector>* f = new FSGrid<Vector>();
   f->setDefVal( defValue );
   f->init( rg->nx(), rg->ny(), rg->nz(), rg->camera() );
   return VectorFrustumGrid(f);
}



VectorFrustumGrid lux::makeFrustumGrid( const FrustumBox& rg, Vector defValue, int psize )
{
   FSGrid<Vector>* f = new FSGrid<Vector>(psize);
   f->setDefVal( defValue );
   f->init( rg->nx(), rg->ny(), rg->nz(), rg->camera() );
   return VectorFrustumGrid(f);
}


// Aliases so client code does not explicitly need to name the type in the
// function name.
ScalarFrustumGrid lux::makeGrid( const FrustumBox& rg, float defValue, int partitionSize )
{
  return makeFrustumGrid(rg, defValue, partitionSize);
}

VectorFrustumGrid lux::makeGrid( const FrustumBox& rg, Vector defValue, int partitionSize )
{
  return makeFrustumGrid(rg, defValue, partitionSize);
}

ColorFrustumGrid lux::makeGrid( const FrustumBox& rg, Color defValue, int partitionSize )
{
  return makeFrustumGrid(rg, defValue, partitionSize);
}

MatrixFrustumGrid lux::makeGrid( const FrustumBox& rg, Matrix defValue, int partitionSize )
{
  return makeFrustumGrid(rg, defValue, partitionSize);
}



ColorFrustumGrid lux::makeColorFrustumGrid( const FrustumGrid& rg, Color defValue ) { return makeFrustumGrid( rg, defValue); }

ColorFrustumGrid lux::makeFrustumGrid( const FrustumGrid& rg, Color defValue )
{
   FSGrid<Color>* f = new FSGrid<Color>();
   f->setDefVal( defValue );
   f->init( rg.nx(), rg.ny(), rg.nz(), rg.camera() );
   return ColorFrustumGrid(f);
}

ColorFrustumGrid lux::makeFrustumGrid( const FrustumBox& rg, Color defValue )
{
   FSGrid<Color>* f = new FSGrid<Color>();
   f->setDefVal( defValue );
   f->init( rg->nx(), rg->ny(), rg->nz(), rg->camera() );
   return ColorFrustumGrid(f);
}


ColorFrustumGrid lux::makeFrustumGrid( const FrustumBox& rg, Color defValue, int psize )
{
   FSGrid<Color>* f = new FSGrid<Color>(psize);
   f->setDefVal( defValue );
   f->init( rg->nx(), rg->ny(), rg->nz(), rg->camera() );
   return ColorFrustumGrid(f);
}



MatrixFrustumGrid lux::makeMatrixFrustumGrid( const FrustumGrid& rg, Matrix defValue ) { return makeFrustumGrid( rg, defValue); }

MatrixFrustumGrid lux::makeFrustumGrid( const FrustumGrid& rg, Matrix defValue )
{
   FSGrid<Matrix>* f = new FSGrid<Matrix>();
   f->setDefVal( defValue );
   f->init( rg.nx(), rg.ny(), rg.nz(), rg.camera() );
   return MatrixFrustumGrid(f);
}



MatrixFrustumGrid lux::makeFrustumGrid( const FrustumBox& rg, Matrix defValue )
{
   FSGrid<Matrix>* f = new FSGrid<Matrix>();
   f->setDefVal( defValue );
   f->init( rg->nx(), rg->ny(), rg->nz(), rg->camera() );
   return MatrixFrustumGrid(f);
}



MatrixFrustumGrid lux::makeFrustumGrid( const FrustumBox& rg, Matrix defValue, int psize )
{
   FSGrid<Matrix>* f = new FSGrid<Matrix>(psize);
   f->setDefVal( defValue );
   f->init( rg->nx(), rg->ny(), rg->nz(), rg->camera() );
   return MatrixFrustumGrid(f);
}




ScalarFrustumGrid::ScalarFrustumGrid( FSGrid<float>* f ) : ScalarFrustumGridBase(f) {}


ScalarFrustumGrid::~ScalarFrustumGrid(){}


ScalarFrustumGrid ScalarFrustumGrid::operator+=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) +  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

ScalarFrustumGrid ScalarFrustumGrid::operator-=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) -  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

ScalarFrustumGrid ScalarFrustumGrid::operator-()
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  -( (*this)->get(i,j,k) ) );
	 }
      }
   }
   return *this;
}

ScalarFrustumGrid ScalarFrustumGrid::operator*=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) *  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

ScalarFrustumGrid ScalarFrustumGrid::operator/=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) /  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}












VectorFrustumGrid::VectorFrustumGrid( FSGrid<Vector>* f ) : VectorFrustumGridBase(f) {}


VectorFrustumGrid::~VectorFrustumGrid(){}


VectorFrustumGrid VectorFrustumGrid::operator+=( const VectorField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) +  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

VectorFrustumGrid VectorFrustumGrid::operator-=( const VectorField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) -  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

VectorFrustumGrid VectorFrustumGrid::operator-()
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) *  (-1.0)  );
	 }
      }
   }
   return *this;
}

VectorFrustumGrid VectorFrustumGrid::operator*=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) *  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

VectorFrustumGrid VectorFrustumGrid::operator/=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) /  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}











MatrixFrustumGrid::MatrixFrustumGrid( FSGrid<Matrix>* f ) : MatrixFrustumGridBase(f) {}


MatrixFrustumGrid::~MatrixFrustumGrid(){}



MatrixFrustumGrid MatrixFrustumGrid::operator+=( const MatrixField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) +  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

MatrixFrustumGrid MatrixFrustumGrid::operator-=( const MatrixField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) -  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

MatrixFrustumGrid MatrixFrustumGrid::operator-()
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k)*(-1.0)  );
	 }
      }
   }
   return *this;
}

MatrixFrustumGrid MatrixFrustumGrid::operator*=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) *  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

MatrixFrustumGrid MatrixFrustumGrid::operator/=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) /  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}









ColorFrustumGrid::ColorFrustumGrid( FSGrid<Color>* f ) : ColorFrustumGridBase(f) {}


ColorFrustumGrid::~ColorFrustumGrid(){}


ColorFrustumGrid ColorFrustumGrid::operator+=( const ColorField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) +  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

ColorFrustumGrid ColorFrustumGrid::operator-=( const ColorField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) -  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

ColorFrustumGrid ColorFrustumGrid::operator-()
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k)*(-1.0)  );
	 }
      }
   }
   return *this;
}

ColorFrustumGrid ColorFrustumGrid::operator*=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) *  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}

ColorFrustumGrid ColorFrustumGrid::operator/=( const ScalarField& e2 )
{
   for(int k=0;k<(*this)->nz();k++)
   {
      for(int j=0;j<(*this)->ny();j++)
      {
         for(int i=0;i<(*this)->nx();i++)
	 {
	    (*this)->set(i,j,k,  (*this)->get(i,j,k) /  e2->eval( (*this)->evalP(i,j,k) )  );
	 }
      }
   }
   return *this;
}



void lux::gridStatistics( const ScalarGrid& g, double& mean, double& stddev, double& max, double& min )
{
   mean = stddev = 0.0;
   max = g->get(0,0,0);
   min = g->get(0,0,0);
   for( int k=0;k<g->nz();k++ )
   {
      for( int j=0;j<g->ny();j++ )
      {
         for( int i=0;i<g->nx();i++ )
         {
            float data = g->get(i,j,k);
            mean += data;
            stddev += data*data;
            if (max < data) { max = data; }
            if (min > data) { min = data; }
         }
      }
   }
   mean /= (g->nx())*(g->ny())*(g->nz());
   stddev /= (g->nx())*(g->ny())*(g->nz());
   stddev -= mean*mean;
}

void lux::gridStatistics( const VectorGrid& g, Vector& mean, double& stddev, Vector& max, Vector& min )
{
   mean = Vector(0,0,0);
   stddev = 0.0;
   max = g->get(0,0,0);
   min = g->get(0,0,0);
   for( int k=0;k<g->nz();k++ )
   {
      for( int j=0;j<g->ny();j++ )
      {
         for( int i=0;i<g->nx();i++ )
         {
            Vector data = g->get(i,j,k);
            mean += data;
            stddev += data*data;
            if (max < data) { max = data; }
            if (min > data) { min = data; }
         }
      }
   }
   mean /= (g->nx())*(g->ny())*(g->nz());
   stddev /= (g->nx())*(g->ny())*(g->nz());
   stddev -= mean*mean;
}




void lux::stamp( ScalarGrid& g, const Vector& pos, const float pscale, const float value )
{
   float dx = g->dx();
   float dy = g->dy();
   float dz = g->dz();
   float dr = dx;
   if( dr > dy ){ dr = dy; }
   if( dr > dz ){ dr = dz; }
   int ixl, ixu, iyl, iyu, izl, izu;
   Vector llc = pos - Vector(pscale, pscale, pscale);
   Vector urc = pos + Vector(pscale, pscale, pscale);
   if( g->getBox( llc, urc, ixl, ixu, iyl, iyu, izl, izu ) )
   {
      for( int k=izl;k<=izu;k++ )
      {
         for( int j=iyl;j<=iyu;j++ )
         {
            for( int i=ixl;i<=ixu;i++ )
            {
	       Vector X = g->evalP(i,j,k) - pos;
               float current = g->get(i,j,k);
               if( X.magnitude() <= pscale ) 
               {
                  float fraction = (pscale-X.magnitude())/dr;
                  if( fraction > 1.0 ){ fraction = 1.0; }
                  g->set(i,j,k,value*fraction + current);
               }
            }
         }
      }
   }
}

void lux::stamp( VectorGrid& g, const Vector& pos, const float pscale, const Vector& value )
{
   float dx = g->dx();
   float dy = g->dy();
   float dz = g->dz();
   float dr = dx;
   if( dr > dy ){ dr = dy; }
   if( dr > dz ){ dr = dz; }
   int ixl, ixu, iyl, iyu, izl, izu;
   Vector llc = pos - Vector(pscale, pscale, pscale);
   Vector urc = pos + Vector(pscale, pscale, pscale);
   if( g->getBox( llc, urc, ixl, ixu, iyl, iyu, izl, izu ) )
   {
      for( int k=izl;k<=izu;k++ )
      {
         for( int j=iyl;j<=iyu;j++ )
         {
            for( int i=ixl;i<=ixu;i++ )
            {
	       Vector X = g->evalP(i,j,k) - pos;
               Vector current = g->get(i,j,k);
               if( X.magnitude() <= pscale ) 
               {
                  float fraction = (pscale-X.magnitude())/dr;
                  if( fraction > 1.0 ){ fraction = 1.0; }
                  g->set(i,j,k,value*fraction + current);
               }
            }
         }
      }
   }
}
 
void lux::stamp( ColorGrid& g, const Vector& pos, const float pscale, const Color& value )
{
   float dx = g->dx();
   float dy = g->dy();
   float dz = g->dz();
   float dr = dx;
   if( dr > dy ){ dr = dy; }
   if( dr > dz ){ dr = dz; }
   int ixl, ixu, iyl, iyu, izl, izu;
   Vector llc = pos - Vector(pscale, pscale, pscale);
   Vector urc = pos + Vector(pscale, pscale, pscale);
   if( g->getBox( llc, urc, ixl, ixu, iyl, iyu, izl, izu ) )
   {
      for( int k=izl;k<=izu;k++ )
      {
         for( int j=iyl;j<=iyu;j++ )
         {
            for( int i=ixl;i<=ixu;i++ )
            {
	       Vector X = g->evalP(i,j,k) - pos;
               Color current = g->get(i,j,k);
               if( X.magnitude() <= pscale ) 
               {
                  float fraction = (pscale-X.magnitude())/dr;
                  if( fraction > 1.0 ){ fraction = 1.0; }
                  g->set(i,j,k,value*fraction + current);
               }
            }
         }
      }
   }
}






VectorGrid lux::optimumVelocityFromGrad( const MatrixField& m, const GridBox& gb, const int nb_iterations )
{
   VectorGrid u = makeGrid( gb, Vector(0,0,0) );
   int nx = u->nx();
   int ny = u->ny();
   int nz = u->nz();
   float dx = u->dx();
   float dy = u->dy();
   float dz = u->dz();


   VectorGrid mdiv = makeGrid( gb, Vector(0,0,0) );
   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
#pragma omp parallel for
         for( int i=0;i<nx;i++ )
         {
            Matrix mmx =( m->eval( u->evalP(i+1,j,k) ) - m->eval( u->evalP(i-1,j,k) ) )/dx;
            Matrix mmy =( m->eval( u->evalP(i,j+1,k) ) - m->eval( u->evalP(i,j-1,k) ) )/dy;
            Matrix mmz =( m->eval( u->evalP(i,j,k+1) ) - m->eval( u->evalP(i,j,k-1) ) )/dz;
            float vx = mmx.Get(0,0) + mmy.Get(1,0) + mmz.Get(2,0);
            float vy = mmx.Get(0,1) + mmy.Get(1,1) + mmz.Get(2,1);
            float vz = mmx.Get(0,2) + mmy.Get(1,2) + mmz.Get(2,2);
            mdiv->set(i,j,k,Vector(vx,vy,vz)); 
         }
      }
   }

   float dd = 2.0/(dx*dx) + 2.0/(dy*dy) + 2.0/(dz*dz);
   for( int loop=0;loop<nb_iterations;loop++ )
   {
      for( int k=1;k<nz-1;k++ )
      {
         for( int j=1;j<ny-1;j++ )
         {
#pragma omp parallel for
            for( int i=1;i<nx-1;i++ )
            {
               Vector V = ( u->get(i-1,j,k) + u->get(i+1,j,k) )/(dx*dx)
                        + ( u->get(i,j-1,k) + u->get(i,j+1,k) )/(dy*dy)
                        + ( u->get(i,j,k-1) + u->get(i,j,k+1) )/(dz*dz);
               V -= mdiv->get(i,j,k); 
               V /= dd;
               u->set(i,j,k,V);
            }
         }
      }
   }
    
   return u;
}


VectorField lux::gradientStretchCMFromGrad( const MatrixField& m, const GridBox& gb, const int nb_iterations, const float T, const int nbgs  )
{
   VectorGrid velocity = optimumVelocityFromGrad( m, gb, nb_iterations );
   VectorGrid Xg = makeGrid( gb, Vector(0,0,0) );
   const double dt = T/nbgs;
   int nx = Xg->nx();
   int ny = Xg->ny();
   int nz = Xg->nz();
   for( int k=0;k<nz;k++ )
   {
      for( int j=0;j<ny;j++ )
      {
#pragma omp parallel for
         for( int i=0;i<nx;i++ )
         {
            Vector P = Xg->evalP(i,j,k);
            Vector u = velocity->get(i,j,k);
            Matrix gu = m->eval(P);
            Matrix GX = exp(-dt*gu );
            Vector X = P;
            Matrix integral;
            Matrix expintegral = dt*sinch(dt*gu);
            for( int it=0;it<nbgs;it++ )
            {
               integral +=  expintegral;
               X = P - u*integral;
               gu = m->eval(X);
               GX = exp( -dt*gu );
               expintegral = expintegral * GX;
            }
            Xg->set(i,j,k,X-P);
         }
      }
   }
   return identity() + gridded(Xg);
}






/*
void lux::readScalarOGrid( const string& fname, ScalarOGrid& g ) { g->read( fname ); }
void lux::readVectorOGrid( const string& fname, VectorOGrid& g ) { g->read( fname ); }
void lux::readColorOGrid( const string& fname, ColorOGrid& g ) { g->read( fname ); }

void lux::writeScalarOGrid( const string& fname, const ScalarOGrid& g ) { g->write( fname ); }
void lux::writeVectorOGrid( const string& fname, const VectorOGrid& g ) { g->write( fname ); }
void lux::writeColorOGrid( const string& fname, const ColorOGrid& g ) { g->write( fname ); }

ScalarOGrid lux::makeOGrid( const GridBox& rg, float defValue )
{
   OVDBGrid<float>* f = new OVDBGrid<float>();
   f->setDefVal( defValue );
   f->init( rg->Lx(), rg->Ly(), rg->Lz(), rg->llc() );
   f->setInterpolationOrder( rg->getInterpolationOrder() );
   return ScalarOGrid(f);
}

VectorOGrid lux::makeOGrid( const GridBox& rg, Vector defValue )
{
   OVDBGrid<Vector>* f = new OVDBGrid<Vector>();
   f->setDefVal( defValue );
   f->init( rg->Lx(), rg->Ly(), rg->Lz(), rg->llc() );
   f->setInterpolationOrder( rg->getInterpolationOrder() );
   return VectorOGrid(f);
}

ColorOGrid lux::makeOGrid( const GridBox& rg, Color defValue );
{
   OVDBGrid<Color>* f = new OVDBGrid<Color>();
   f->setDefVal( defValue );
   f->init( rg->Lx(), rg->Ly(), rg->Lz(), rg->llc() );
   f->setInterpolationOrder( rg->getInterpolationOrder() );
   return ColorOGrid(f);
}

MatrixOGrid lux::makeOGrid( const GridBox& rg, Matrix defValue );
{
   OVDBGrid<Matrix>* f = new OVDBGrid<Matrix>();
   f->setDefVal( defValue );
   f->init( rg->Lx(), rg->Ly(), rg->Lz(), rg->llc() );
   f->setInterpolationOrder( rg->getInterpolationOrder() );
   return MatrixOGrid(f);
}



void lux::stamp( ScalarOGrid& grid, const ScalarField& field, const int nbsamples ) {}
void lux::stamp( VectorOGrid& grid, const VectorField& field, const int nbsamples ) {}
void lux::stamp( ColorOGrid& grid, const ColorField& field, const int nbsamples ) {}
void lux::stamp( MatrixOGrid& grid, const MatrixField& field, const int nbsamples ) {}

*/




IntervalSet lux::getIntervalSet( const ScalarGrid& g, int maxlvl, int minobj )
{
   IntervalSet is = makeIntervalSet( g->llc(), g->urc(), 0, maxlvl, minobj );
   int nbblocks = g->NbPartitions();
   int aabbcount = 0;
   for( int i=0;i<nbblocks;i++ )
   {
      if( g->goodBlock(i) )
      {
         int i0,j0,k0,i1,j1,k1;
         g->blockBounds( i, i0, j0, k0, i1, j1, k1 );
         Vector p = g->evalP(i0,j0,k0);
         double xmin = p[0];
         double ymin = p[1];
         double zmin = p[2];
         double xmax = p[0];
         double ymax = p[1];
         double zmax = p[2];
         p = g->evalP(i0,j0,k1);
         xmin = ( xmin <= p[0] ) ? xmin : p[0];
         ymin = ( ymin <= p[1] ) ? ymin : p[1];
         zmin = ( zmin <= p[2] ) ? zmin : p[2];
         xmax = ( xmax >= p[0] ) ? xmax : p[0];
         ymax = ( ymax >= p[1] ) ? ymax : p[1];
         zmax = ( zmax >= p[2] ) ? zmax : p[2];
         p = g->evalP(i0,j1,k0);
         xmin = ( xmin <= p[0] ) ? xmin : p[0];
         ymin = ( ymin <= p[1] ) ? ymin : p[1];
         zmin = ( zmin <= p[2] ) ? zmin : p[2];
         xmax = ( xmax >= p[0] ) ? xmax : p[0];
         ymax = ( ymax >= p[1] ) ? ymax : p[1];
         zmax = ( zmax >= p[2] ) ? zmax : p[2];
         p = g->evalP(i1,j0,k0);
         xmin = ( xmin <= p[0] ) ? xmin : p[0];
         ymin = ( ymin <= p[1] ) ? ymin : p[1];
         zmin = ( zmin <= p[2] ) ? zmin : p[2];
         xmax = ( xmax >= p[0] ) ? xmax : p[0];
         ymax = ( ymax >= p[1] ) ? ymax : p[1];
         zmax = ( zmax >= p[2] ) ? zmax : p[2];
         p = g->evalP(i1,j0,k1);
         xmin = ( xmin <= p[0] ) ? xmin : p[0];
         ymin = ( ymin <= p[1] ) ? ymin : p[1];
         zmin = ( zmin <= p[2] ) ? zmin : p[2];
         xmax = ( xmax >= p[0] ) ? xmax : p[0];
         ymax = ( ymax >= p[1] ) ? ymax : p[1];
         zmax = ( zmax >= p[2] ) ? zmax : p[2];
         p = g->evalP(i1,j1,k0);
         xmin = ( xmin <= p[0] ) ? xmin : p[0];
         ymin = ( ymin <= p[1] ) ? ymin : p[1];
         zmin = ( zmin <= p[2] ) ? zmin : p[2];
         xmax = ( xmax >= p[0] ) ? xmax : p[0];
         ymax = ( ymax >= p[1] ) ? ymax : p[1];
         zmax = ( zmax >= p[2] ) ? zmax : p[2];
         p = g->evalP(i0,j1,k1);
         xmin = ( xmin <= p[0] ) ? xmin : p[0];
         ymin = ( ymin <= p[1] ) ? ymin : p[1];
         zmin = ( zmin <= p[2] ) ? zmin : p[2];
         xmax = ( xmax >= p[0] ) ? xmax : p[0];
         ymax = ( ymax >= p[1] ) ? ymax : p[1];
         zmax = ( zmax >= p[2] ) ? zmax : p[2];
         p = g->evalP(i1,j1,k1);
         xmin = ( xmin <= p[0] ) ? xmin : p[0];
         ymin = ( ymin <= p[1] ) ? ymin : p[1];
         zmin = ( zmin <= p[2] ) ? zmin : p[2];
         xmax = ( xmax >= p[0] ) ? xmax : p[0];
         ymax = ( ymax >= p[1] ) ? ymax : p[1];
         zmax = ( zmax >= p[2] ) ? zmax : p[2];
         Vector llc = Vector(xmin,ymin,zmin);
         Vector urc = Vector(xmax,ymax,zmax);
         AABB aabb(llc,urc);
         is->addObject( aabb );  
         aabbcount++;       
      }
   }
   is->Divide();
   cout << "Number of objects in intervalset: " << aabbcount << endl;
   return is;   
}

IntervalSet lux::getIntervalSet( const ScalarFrustumGrid& g, int maxlvl, int minobj )
{
   IntervalSet is = makeIntervalSet( g->llc(), g->urc(), 0, maxlvl, minobj );
   int nbblocks = g->NbPartitions();
   for( int i=0;i<nbblocks;i++ )
   {
      if( g->goodBlock(i) )
      {
         int i0,j0,k0,i1,j1,k1;
         g->blockBounds( i, i0, j0, k0, i1, j1, k1 );
         Vector p = g->evalP(i0,j0,k0);
         double xmin = p[0];
         double ymin = p[1];
         double zmin = p[2];
         double xmax = p[0];
         double ymax = p[1];
         double zmax = p[2];
         p = g->evalP(i0,j0,k1);
         xmin = ( xmin <= p[0] ) ? xmin : p[0];
         ymin = ( ymin <= p[1] ) ? ymin : p[1];
         zmin = ( zmin <= p[2] ) ? zmin : p[2];
         xmax = ( xmax >= p[0] ) ? xmax : p[0];
         ymax = ( ymax >= p[1] ) ? ymax : p[1];
         zmax = ( zmax >= p[2] ) ? zmax : p[2];
         p = g->evalP(i0,j1,k0);
         xmin = ( xmin <= p[0] ) ? xmin : p[0];
         ymin = ( ymin <= p[1] ) ? ymin : p[1];
         zmin = ( zmin <= p[2] ) ? zmin : p[2];
         xmax = ( xmax >= p[0] ) ? xmax : p[0];
         ymax = ( ymax >= p[1] ) ? ymax : p[1];
         zmax = ( zmax >= p[2] ) ? zmax : p[2];
         p = g->evalP(i1,j0,k0);
         xmin = ( xmin <= p[0] ) ? xmin : p[0];
         ymin = ( ymin <= p[1] ) ? ymin : p[1];
         zmin = ( zmin <= p[2] ) ? zmin : p[2];
         xmax = ( xmax >= p[0] ) ? xmax : p[0];
         ymax = ( ymax >= p[1] ) ? ymax : p[1];
         zmax = ( zmax >= p[2] ) ? zmax : p[2];
         p = g->evalP(i1,j0,k1);
         xmin = ( xmin <= p[0] ) ? xmin : p[0];
         ymin = ( ymin <= p[1] ) ? ymin : p[1];
         zmin = ( zmin <= p[2] ) ? zmin : p[2];
         xmax = ( xmax >= p[0] ) ? xmax : p[0];
         ymax = ( ymax >= p[1] ) ? ymax : p[1];
         zmax = ( zmax >= p[2] ) ? zmax : p[2];
         p = g->evalP(i1,j1,k0);
         xmin = ( xmin <= p[0] ) ? xmin : p[0];
         ymin = ( ymin <= p[1] ) ? ymin : p[1];
         zmin = ( zmin <= p[2] ) ? zmin : p[2];
         xmax = ( xmax >= p[0] ) ? xmax : p[0];
         ymax = ( ymax >= p[1] ) ? ymax : p[1];
         zmax = ( zmax >= p[2] ) ? zmax : p[2];
         p = g->evalP(i0,j1,k1);
         xmin = ( xmin <= p[0] ) ? xmin : p[0];
         ymin = ( ymin <= p[1] ) ? ymin : p[1];
         zmin = ( zmin <= p[2] ) ? zmin : p[2];
         xmax = ( xmax >= p[0] ) ? xmax : p[0];
         ymax = ( ymax >= p[1] ) ? ymax : p[1];
         zmax = ( zmax >= p[2] ) ? zmax : p[2];
         p = g->evalP(i1,j1,k1);
         xmin = ( xmin <= p[0] ) ? xmin : p[0];
         ymin = ( ymin <= p[1] ) ? ymin : p[1];
         zmin = ( zmin <= p[2] ) ? zmin : p[2];
         xmax = ( xmax >= p[0] ) ? xmax : p[0];
         ymax = ( ymax >= p[1] ) ? ymax : p[1];
         zmax = ( zmax >= p[2] ) ? zmax : p[2];
         Vector llc = Vector(xmin,ymin,zmin);
         Vector urc = Vector(xmax,ymax,zmax);
         AABB aabb(llc,urc);
         is->addObject( aabb );         
      }
   }
   is->Divide();
   cout << "Number of frustum objects in intervalset: " << is->nbObjects() << endl;
   return is;   
}



double lux::gridMean( const ScalarGrid& g )
{
   double mean = 0.0;
   ProgressMeter pm( g->nz(), "Grid Mean" );
   for( int k=0;k<g->nz();k++ )
   {
      for( int j=0;j<g->ny();j++ )
      {
         for( int i=0;i<g->nx();i++ )
         {
            float data = g->get(i,j,k);
            mean += data;
         }
      }
      pm.update();
   }
   double norm = g->nx();
   norm *= g->ny();
   norm *= g->nz();
   mean /= norm;
   return mean;
}

double lux::gridStdDev( const ScalarGrid& g )
{
   double mean = 0.0;
   double stddev = 0.0;
   ProgressMeter pm( g->nz(), "Grid StdDev" );
   for( int k=0;k<g->nz();k++ )
   {
      for( int j=0;j<g->ny();j++ )
      {
         for( int i=0;i<g->nx();i++ )
         {
            float data = g->get(i,j,k);
            mean += data;
            stddev += data*data;
         }
      }
      pm.update();
   }
   double norm = g->nx();
   norm *= g->ny();
   norm *= g->nz();
   mean /= norm;
   stddev /= norm;
   stddev -= mean*mean;
   return stddev;
}

double lux::gridMax( const ScalarGrid& g )
{
   double max = g->get(0,0,0);
   ProgressMeter pm( g->nz(), "Grid Max" );
   for( int k=0;k<g->nz();k++ )
   {
      for( int j=0;j<g->ny();j++ )
      {
         for( int i=0;i<g->nx();i++ )
         {
            float data = g->get(i,j,k);
            if (max < data) { max = data; }
         }
      }
      pm.update();
   }
   return max;
}

double lux::gridMin( const ScalarGrid& g )
{
   double min = g->get(0,0,0);
   ProgressMeter pm( g->nz(), "Grid Min" );
   for( int k=0;k<g->nz();k++ )
   {
      for( int j=0;j<g->ny();j++ )
      {
         for( int i=0;i<g->nx();i++ )
         {
            float data = g->get(i,j,k);
            if (min > data) { min = data; }
         }
      }
      pm.update();
   }
   return min;
}
