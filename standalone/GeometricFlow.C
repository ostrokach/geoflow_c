#include "GeometricFlow.h"

using namespace std;

GeometricFlow::GeometricFlow( )
{
   p_pres_i = 0.008;
   p_gama_i = 0.0001;
   p_npiter = 1;
   p_ngiter = 1;

   // probe radius for creating the solvent accessible surface
   p_tauval = 1.40;

   p_prob = 0.0;

   p_ffmodel = 1;  // FFMODEL: 1 for ZAP-9/AM1-BCCv1; 2 for OPLS/AA

   // SIGMAS: Angstrom (radius of water molecule based on LJ parameter sigma)
   p_sigmas = 1.5828;
   // EPSILONW:  epsilon parameter of O (kcal/mol) of water molecule.
   p_epsilonw = 0.1554;
   // VDWDISPERSION:  1(on) or 0 (off)- previously called REPULSIVE
   p_vdwdispersion = 0;
   // EXTVALUE:  (distance atom surface and box boundary)
   p_extvalue = 1.90;
     // iprec
     	// istep
   // IPREC: flag to indicate the usage of preconditioner iprec =1 (yes); 0 (no)
   p_iadi = 0;
   // ALPHA: weight of previous solution to change the next solution in geometry flow
   p_alpha = 0.50;
     // IPBIN  //IPBIN: start guess for PB 1; inherit '0' - not used?
   p_tol = 1e-4;
   p_tottf = 3.5;
   p_dcel = 0.25;
   p_maxstep = 20;
   p_epsilons = 80.00;
   p_epsilonp = 1.5;
   p_radexp = 1;
   p_crevalue = 0.01;
     // idacsl //idacsl: 0 for solvation force calculation; 1 or accuracy test
   p_density = 0.03346;
}

//~GeometricFlow() { };

//
//  copy constructor
// 
GeometricFlow::GeometricFlow( const GeometricFlow& gf )
{
   p_expervalue = gf.p_expervalue;
   p_pres_i = gf.p_pres_i;
   p_gama_i = gf.p_gama_i;
   p_npiter = gf.p_npiter;
   p_ngiter = gf.p_ngiter;
   p_tauval = gf.p_tauval;
   p_prob = gf.p_prob;
   p_ffmodel = gf.p_ffmodel;
   p_sigmas = gf.p_sigmas;
   p_epsilonw = gf.p_epsilonw;
   p_vdwdispersion = gf.p_vdwdispersion;
   p_extvalue = gf.p_extvalue;
   p_iadi = gf.p_iadi;
   p_alpha = gf.p_alpha;
   p_tol = gf.p_tol;
   p_tottf = gf.p_tottf;
   p_dcel = gf.p_dcel;
   p_maxstep = gf.p_maxstep;
   p_epsilons = gf.p_epsilons;
   p_epsilonp = gf.p_epsilonp;
   p_radexp = gf.p_radexp;
   p_crevalue = gf.p_crevalue;
   p_density = gf.p_density;
}

GeometricFlow::GeometricFlow( const GeometricFlow* gf )
{
   p_expervalue = gf->p_expervalue;
   p_pres_i = gf->p_pres_i;
   p_gama_i = gf->p_gama_i;
   p_npiter = gf->p_npiter;
   p_ngiter = gf->p_ngiter;
   p_tauval = gf->p_tauval;
   p_prob = gf->p_prob;
   p_ffmodel = gf->p_ffmodel;
   p_sigmas = gf->p_sigmas;
   p_epsilonw = gf->p_epsilonw;
   p_vdwdispersion = gf->p_vdwdispersion;
   p_extvalue = gf->p_extvalue;
   p_iadi = gf->p_iadi;
   p_alpha = gf->p_alpha;
   p_tol = gf->p_tol;
   p_tottf = gf->p_tottf;
   p_dcel = gf->p_dcel;
   p_maxstep = gf->p_maxstep;
   p_epsilons = gf->p_epsilons;
   p_epsilonp = gf->p_epsilonp;
   p_radexp = gf->p_radexp;
   p_crevalue = gf->p_crevalue;
   p_density = gf->p_density;
}

//
//  operator=
//
//GeometricFlow& GeometricFlow::operator=( const GeometricFlow& gf )
//{
//}

void printAllParams()
{
}
