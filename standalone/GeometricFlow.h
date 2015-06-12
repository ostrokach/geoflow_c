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

      //
      //  Out Parameters
      //

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

      //
      //  for debugging
      //
      void printAllParams();
      


};

#endif

      

