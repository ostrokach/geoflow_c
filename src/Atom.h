#ifndef __Atom_h_
#define __Atom_h_

#include <iostream>
#include <fstream>
#include <ostream>
#include <vector>
#include <map>
#include <iterator>
#include <sstream>
#include <algorithm>
#include <climits>
#include <stdio.h>

#include "Mat.h"
#include "ComData.h"

using namespace std;

class Atom
{
   private:

      double p_x, p_y, p_z;
      double p_radius;
      double p_pqr;
      double p_ljepsilon;

      void setRadius( const int ffmodel, double r )
      { 
         if (ffmodel == 1 )
            p_radius = ( r < 1e-6 ) ? 1.21 : r; 
         else
            p_radius = r;
      }

   public:

      //
      //  empty constructor
      //
      Atom( );
      Atom( const int ffmodel, double x, double y, double z, double r, double pqr);
      Atom( const int ffmodel, double x, double y, double z, double r, double pqr, double e );
      //~Path() { };

      //
      //  copy constructor
      // 
      Atom( const Atom& A ) ;
      Atom( const Atom* A ) ;

      //
      //  operator=
      //
      //Atom& operator=( const Atom& A ) ;
      //ostream& operator<<( ostream& os, const Atom& A );

      double x() const { return p_x; }
      double y() const { return p_y; }
      double z() const { return p_z; }
      double r() const { return p_radius; }
      double pqr() const { return p_pqr; }
      double epsilon() const { return p_ljepsilon; }

      //
      //  for debugging
      //
      void print() const;
};

class AtomList
{
   private:

      vector< Atom > p_atomList;

   public:

      AtomList();

      AtomList( string xyzr_file, const double radexp, const int ffmodel );

      AtomList( const AtomList& AL ) { p_atomList = AL.p_atomList ; }

      AtomList( const AtomList* AL ) { p_atomList = AL->p_atomList ; }

      unsigned int size() const { return p_atomList.size(); }

      const Atom& get( unsigned int i ) const { return p_atomList[i]; }

      void changeChargeDistribution 
         ( Mat<>& charget, Mat<>& corlocqt, Mat< size_t>& loc_qt,
           const ComData& comData ) const;

      void print() const;
};

#endif

      

