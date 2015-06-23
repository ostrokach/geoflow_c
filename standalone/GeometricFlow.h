#ifndef __GeometricFlow_h_
#define __GeometricFlow_h_

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
#include <valarray>

#include "Atom.h"
#include "Mat.h"
#include "ComData.h"

#include <Eigen/Sparse>


using namespace std;


class GeometricFlow
{
   private:

      //
      //  Input Parameters
      //
      //int nmol = 17,	
      double p_expervalue;
      double p_pres_i;
      double p_gama_i;
      int    p_npiter;
      int    p_ngiter;
      double p_tauval;
      double p_prob;
      int    p_ffmodel;
      double p_sigmas;
      double p_epsilonw;
      int    p_vdwdispersion;
      double p_extvalue;
      // 0,		// iprec
      // 10,		// istep
      int    p_iadi;
      double p_alpha;
      // 1,		// IPBIN
      double p_tol;
      double p_tottf;
      double p_dcel;
      int    p_maxstep;
      double p_epsilons;
      double p_epsilonp;
      int    p_radexp;
      double p_crevalue;
      // 0,		// idacsl
      double p_density;
      double p_foo;  // MAGIC_FOO

      ComData p_comdata;

      int p_lj_iosetar;
      int p_lj_iosetaa;
      int p_lj_iwca;


      /*
      double p_comdata_deltax;
      double p_comdata_deltay;
      double p_comdata_deltaz;
      double p_comdata_dcel;
      double p_comdata_xleft;
      double p_comdata_yleft;
      double p_comdata_zleft;
      double p_comdata_xright;
      double p_comdata_yright;
      double p_comdata_zright;
      double p_comdata_nx;
      double p_comdata_ny;
      double p_comdata_nz;
      double p_comdata_pi;
      */

      double left(const valarray<double>& pr, double h, double ev)
      {
         return floor( (pr - ev).min()/h ) * h - ev;
      }

      double right(const valarray<double>& pr, double h, double ev)
      {
         return ceil( (pr + ev).max()/h ) * h + ev;
      }

      double dot (double x, double y, double z) { return x*x + y*y + z*z; }


      void domainInitialization( const AtomList& atomlist );

      void yhsurface( const AtomList& atomList, 
         double tott, double dt, Mat<>& phitotx, Mat<>& surfu, int iloop,
         double& area, double& volume, double& attint, double alpha, int iadi,
         int igfin, double roro, double conms );

      void potIntegral(double rcfactor, size_t natm,
            valarray<double>& atom_x, valarray<double>& atom_y,
            valarray<double>& atom_z, valarray<double>& seta12,
            valarray<double>& seta6, valarray<double>& epsilon,
            valarray<double>& sigma, Mat<>& g, Mat<>& potr, Mat<>& pota);
      
      double volumeIntegration(const Mat<>& f);

      void upwinding(double dt, int nt, 
                              Mat<>& g, Mat<>& su, Mat<>& phitotx);

      void initial(double xl, double yl, double zl, int n_atom,
            const std::valarray<double>& atom_x, const std::valarray<double>& atom_y,
            const std::valarray<double>& atom_z, const std::valarray<double>& atom_r,
            Mat<>& g, Mat<>& phi);

      void normalizeSurfuAndEps (Mat<>& surfu, Mat<>& eps);

      void computeSoleng(double& soleng, 
                   Mat<>& phi, Mat<>& charget, Mat<size_t>& loc_qt);

      void seteqb( Mat<>& bg, const AtomList& al, const Mat<>& charget,
            const Mat<>& corlocqt);

      double qb(size_t i,size_t j,size_t k, const AtomList& al,
            const Mat<>& charget, const Mat<>& corlocqt );
      
      double qbboundary( double x, double y, double z, const AtomList& al );
      
      double qbinterior(double x, double y, double z, 
            const Mat<>& charget, const Mat<>& corlocqt);

      void pbsolver(Mat<>& eps, Mat<>& phi, Mat<>& bgf, double tol, int iter);

   public:

      //
      //  empty constructor
      //
      GeometricFlow( );
      //~Path() { };

      //
      //  copy constructor
      // 
      GeometricFlow( const GeometricFlow& gf ) ;
      GeometricFlow( const GeometricFlow* gf ) ;

      //
      //  operator=
      //
      //GeometricFlow& operator=( const GeometricFlow& gf ) ;

      //
      //  Set and Get methods
      //
      void setExperValue( double exper_val ) { p_expervalue = exper_val; }

      void setPressure( double pres_i ) { p_pres_i = pres_i; }
      
      void setGama( double gama_i ) { p_gama_i = gama_i; }
      
      void setFFModel( int ffmodel ) { p_ffmodel = ffmodel; }

      void setVDWDispersion( int vdwdispersion ) 
            { p_vdwdispersion = vdwdispersion; }

      void setExtValue( double extvalue ) { p_extvalue = extvalue; }
      
      void setEpsilonS( double epsilons ) { p_epsilons = epsilons; }
      
      void setEpsilonP( double epsilonp ) { p_epsilonp = epsilonp; }

      // uncomment if needed:
      //void setNPiter( int npiter ) { p_npiter = npiter; }
      //void setNGiter( int ngiter ) { p_ngiter = ngiter; }
      //void setTauVal( double tauval ) { p_tauval = p_tauval; }
      //void setProb( double prob ) { p_prob = prob; }

      double getExperValue() { return p_expervalue; }

      double getPressure() { return p_pres_i; }
      
      double getGama() { return p_gama_i; }
      
      int getFFModel() { return p_ffmodel; }

      double getVDWDispersion() { return p_vdwdispersion; }

      double getExtValue() { return p_extvalue; }
      
      double getEpsilonS() { return p_epsilons; }
      
      double getEpsilonP() { return p_epsilonp; }

      double getRadExp() { return p_radexp; }

      //
      //  for debugging
      //
      void printAllParams();

      void run( const AtomList& atomList );
      


};

#endif

      

