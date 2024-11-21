
#include "LevelSetGeometry.h"
#include "ProgressMeter.h"



using namespace lux;


const bool lux::ClosestPerpDistance( const Triangle& t, const Vector& P, double eps, double& dist)
{
   Vector bary = t.barycentric( P );
   dist = ClosestDistance( t, P );
   if( bary.Z() < 0.0 )
   {
      dist = -dist;
   }
   return true;
}




void lux::GenerateLevelSet( const TriangleGeometry& geom, const int nbIterations, const double slop, const double narrow_band_expansion, ScalarGrid& ls )
{
//  Pass 0: Initialize grid with default distant value
   ls->clear();
   float maxdistance = -sqrt(  pow(ls->Lx(),2) + pow(ls->Ly(),2) + pow(ls->Lz(),2) )*10.0;
   ls->setDefVal( maxdistance );
   int gridsize = (ls->nx())*(ls->ny())*(ls->nz());
//  Pass 1: rasterize triangles into the grid
   std::map<int,int> adjacents;
   std::map<int,int> hits;
   int* closest_face = new int[gridsize];
   for( int i=0;i<gridsize;i++ ){ closest_face[i] = -1; }
   ProgressMeter * pm = new ProgressMeter( geom.nbFaces(), "Level Set Pass 1" );
   for( int f=0;f<(int)geom.nbFaces();f++ )
   {
      int e1,e2,e3;
      geom.getFace( f, e1, e2, e3 );
      const Vector& v1 = geom.getVertex(e1);
      const Vector& v2 = geom.getVertex(e2);
      const Vector& v3 = geom.getVertex(e3);
      // bounding box for triangle
      double xmax = v1.X();
      xmax = ( xmax > v2.X() ) ? xmax : v2.X();
      xmax = ( xmax > v3.X() ) ? xmax : v3.X();
      double ymax = v1.Y();
      ymax = ( ymax > v2.Y() ) ? ymax : v2.Y();
      ymax = ( ymax > v3.Y() ) ? ymax : v3.Y();
      double zmax = v1.Z();
      zmax = ( zmax > v2.Z() ) ? zmax : v2.Z();
      zmax = ( zmax > v3.Z() ) ? zmax : v3.Z();
      double xmin = v1.X();
      xmin = ( xmin < v2.X() ) ? xmin : v2.X();
      xmin = ( xmin < v3.X() ) ? xmin : v3.X();
      double ymin = v1.Y();
      ymin = ( ymin < v2.Y() ) ? ymin : v2.Y();
      ymin = ( ymin < v3.Y() ) ? ymin : v3.Y();
      double zmin = v1.Z();
      zmin = ( zmin < v2.Z() ) ? zmin : v2.Z();
      zmin = ( zmin < v3.Z() ) ? zmin : v3.Z();

      // Add cells on each side to bounding box

      xmin -= ls->dx()*narrow_band_expansion;
      ymin -= ls->dy()*narrow_band_expansion;
      zmin -= ls->dz()*narrow_band_expansion;
      xmax += ls->dx()*narrow_band_expansion;
      ymax += ls->dy()*narrow_band_expansion;
      zmax += ls->dz()*narrow_band_expansion;


      // convert bounding box to range of grid indices
      int imin,jmin,kmin,imax,jmax,kmax;
      if( ls->getBox( Vector(xmin,ymin,zmin), Vector(xmax,ymax,zmax), imin, imax, jmin, jmax, kmin, kmax ) )
      {
         Triangle t( v1,v2,v3 );
         for( int k=kmin;k<=kmax;k++ )
         {
            for( int j=jmin;j<=jmax;j++ )
            {
               for( int i=imin;i<=imax;i++ )
               {
                  int indx = ls->index(i,j,k);
                  Vector P = ls->evalP(i,j,k);
                  double dist;
                  if( ClosestPerpDistance( t, P, slop, dist ) )
                  {
                     if( closest_face[indx] < 0 )
                     {
                        ls->set(i,j,k,-dist);
                        closest_face[indx] = f; 
                        hits[indx] = 1;
                     }
                     else
                     {
                        double current_value = ls->get(i,j,k);
                        if( fabs(dist) < fabs(current_value) )
                        {
                           ls->set(i,j,k,-dist);
                           closest_face[indx] = f;
                           hits[indx] = 1;
                        }
                     }
                  }
               }
            }
         }
      }
      pm->update();
   }
   delete pm;

   adjacents.clear(); 
   for( std::map<int,int>::iterator p=hits.begin();p!=hits.end();p++)
   {
      int i,j,k;
      ls->triple(p->first,i,j,k);
      // mark points around this one as adjacent.

      if(i>0)
      {
         int iindx = ls->index(i-1,j,k);
         if( closest_face[iindx] < 0 ) { adjacents[iindx] = 1; }
      }
      if(i<ls->nx()-1)
      {
         int iindx = ls->index(i+1,j,k);
         if( closest_face[iindx] < 0 ) { adjacents[iindx] = 1; }
      }
      if(j>0)
      {
         int iindx = ls->index(i,j-1,k);
         if( closest_face[iindx] < 0 ) { adjacents[iindx] = 1; }
      }
      if(j<ls->ny()-1)
      {
         int iindx = ls->index(i,j+1,k);
         if( closest_face[iindx] < 0 ) { adjacents[iindx] = 1; }
      }
      if(k>0)
      {
         int iindx = ls->index(i,j,k-1);
         if( closest_face[iindx] < 0 ) { adjacents[iindx] = 1; }
      }
      if(k<ls->nz()-1)
      {
         int iindx = ls->index(i,j,k+1);
         if( closest_face[iindx] < 0 ) { adjacents[iindx] = 1; }
      }
   }

//  Loop:
//  Pass 2: Move out to adjacent grid points and update them

   pm = new ProgressMeter( nbIterations, "Level Set Pass 2" );
   for( int iter=0;iter<nbIterations;iter++ )
   {
     hits.clear();
     for( std::map<int,int>::iterator p=adjacents.begin();p!=adjacents.end();p++)
     {
      int indx = p->first;
      if( closest_face[indx] < 0 )
      {
         int i,j,k;
         ls->triple(indx,i,j,k);
         Vector P = ls->evalP(i,j,k);
         // check the grid points around this one to find the list of closest faces for them.
         int imin = i-1;
         if( imin < 0 ){ imin=0; }
         int imax = i+1;
         if( imax >= ls->nx() ){ imax=ls->nx()-1; }
         int jmin = j-1;
         if( jmin < 0 ){ jmin=0; }
         int jmax = j+1;
         if( jmax >= ls->ny() ){ jmax=ls->ny()-1; }
         int kmin = k-1;
         if( kmin < 0 ){ kmin=0; }
         int kmax = k+1;
         if( kmax >= ls->nz() ){ kmax=ls->nz()-1; }
         for( int kk=kmin;kk<=kmax;kk++ )
         {
            for( int jj=jmin;jj<=jmax;jj++ )
            {
               for( int ii=imin;ii<=imax;ii++ )
               {
                  int iindx = ls->index(ii,jj,kk);
                  if( iindx != indx )
                  {
                     if( closest_face[iindx] >= 0 )
                     {
                        // check distance from this face
                        int e1,e2,e3;
                        geom.getFace( closest_face[iindx], e1, e2, e3 );
                        const Vector& v1 = geom.getVertex(e1);
                        const Vector& v2 = geom.getVertex(e2);
                        const Vector& v3 = geom.getVertex(e3);
                        const Triangle t(v1,v2,v3);
                        double dist;
                        if( ClosestPerpDistance( t, P, slop, dist ) )
                        {
                           float current_value = ls->get(i,j,k);
                           if( fabs(dist) < fabs(current_value) )
                           {
                              ls->set(i,j,k,-dist);
                              closest_face[indx] = closest_face[iindx];
                              hits[indx] = 1;
                           }
                        }
                     }
                  }
               }
            }
         }
      }
     }

     adjacents.clear(); 
     for( std::map<int,int>::iterator p=hits.begin();p!=hits.end();p++)
     {
         int i,j,k;
         ls->triple(p->first,i,j,k);
         // mark points around this one as adjacent.

         if(i>0)
         {
            int iindx = ls->index(i-1,j,k);
            if( closest_face[iindx] < 0 ) { adjacents[iindx] = 1; }
         }
         if(i<ls->nx()-1)
         {
            int iindx = ls->index(i+1,j,k);
            if( closest_face[iindx] < 0 ) { adjacents[iindx] = 1; }
         }
         if(j>0)
         {
            int iindx = ls->index(i,j-1,k);
            if( closest_face[iindx] < 0 ) { adjacents[iindx] = 1; }
         }
         if(j<ls->ny()-1)
         {
            int iindx = ls->index(i,j+1,k);
            if( closest_face[iindx] < 0 ) { adjacents[iindx] = 1; }
         }
         if(k>0)
         {
            int iindx = ls->index(i,j,k-1);
            if( closest_face[iindx] < 0 ) { adjacents[iindx] = 1; }
         }
         if(k<ls->nz()-1)
         {
            int iindx = ls->index(i,j,k+1);
            if( closest_face[iindx] < 0 ) { adjacents[iindx] = 1; }
         }
     }
     pm->update();
   }
   delete pm;

/*
//  Loop:
//  Pass 3: Fill in points that have not been touched.

   cout << "Level Set Pass 3\n";
   while( adjacents.size() > 0 )
   {
     cout << " Size: " << adjacents.size() << "\r" << flush;
     hits.clear();
     for( std::map<int,int>::iterator p=adjacents.begin();p!=adjacents.end();p++)
     {
      int indx = p->first;
      if( closest_face[indx] < 0 )
      {
         int i,j,k;
         ls->triple(indx,i,j,k);
         Vector P = ls->evalP(i,j,k);
         // check the grid points around this one to find the list of closest faces for them.
         int imin = i-1;
         if( imin < 0 ){ imin=0; }
         int imax = i+1;
         if( imax >= ls->nx() ){ imax=ls->nx()-1; }
         int jmin = j-1;
         if( jmin < 0 ){ jmin=0; }
         int jmax = j+1;
         if( jmax >= ls->ny() ){ jmax=ls->ny()-1; }
         int kmin = k-1;
         if( kmin < 0 ){ kmin=0; }
         int kmax = k+1;
         if( kmax >= ls->nz() ){ kmax=ls->nz()-1; }
         float value = maxdistance;
         int face = -1;
         for( int kk=kmin;kk<=kmax;kk++ )
         {
            for( int jj=jmin;jj<=jmax;jj++ )
            {
               for( int ii=imin;ii<=imax;ii++ )
               {
                  if( closest_face[indx] < 0 )
                  {
                     int iindx = ls->index(ii,jj,kk);
                     if( iindx != indx )
                     {
                        if( closest_face[iindx] >= 0 )
                        {
                           float current_value = ls->get(ii,jj,kk);
                           if( fabs(current_value) < fabs(value) )
                           { 
                              value = current_value; 
                              face = closest_face[iindx];
                           }
                        }
                     }
                  }
               }
            }
         }
         if( value > 0 ){ls->set(i,j,k,value);}
         closest_face[indx] = face;
         // only need to update inside values.  Outside
         // values are defaulted.
         if( value >= 0 ){ hits[indx] = 1; }
      }
     }

     adjacents.clear(); 
     for( std::map<int,int>::iterator p=hits.begin();p!=hits.end();p++)
     {
         int i,j,k;
         ls->triple(p->first,i,j,k);
         // mark points around this one as adjacent.
         int imin = i-1;
         if( imin < 0 ){ imin=0; }
         int imax = i+1;
         if( imax >= ls->nx() ){ imax=ls->nx()-1; }
         int jmin = j-1;
         if( jmin < 0 ){ jmin=0; }
         int jmax = j+1;
         if( jmax >= ls->ny() ){ jmax=ls->ny()-1; }
         int kmin = k-1;
         if( kmin < 0 ){ kmin=0; }
         int kmax = k+1;
         if( kmax >= ls->nz() ){ kmax=ls->nz()-1; }
         for( int kk=kmin;kk<=kmax;kk++ )
         {
            for( int jj=jmin;jj<=jmax;jj++ )
            {
               for( int ii=imin;ii<=imax;ii++ )
               {
                  int iindx = ls->index(ii,jj,kk);
                  if( closest_face[iindx] < 0 )
                  {
                     adjacents[iindx] = 1;
                  }
               }
            }
         }
     }
   }
   cout << endl;
*/


   delete[] closest_face;
}

