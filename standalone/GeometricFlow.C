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
   p_maxstep = 20;
   p_epsilons = 80.00;
   p_epsilonp = 1.5;
   p_radexp = 1;
   p_crevalue = 0.01;
     // idacsl //idacsl: 0 for solvation force calculation; 1 or accuracy test
   p_density = 0.03346;

   p_dcel = 0.25;  // distance per cell - grid spacing
   p_comdata.init( p_dcel );
   //p_comdata_deltax = p_dcel;
   //p_comdata_deltay = p_dcel;
   //p_comdata_deltaz = p_dcel;
   //p_comdata_pi = acos(-1.0);  //Pi

   // http://ccom.ucsd.edu/~mholst/pubs/dist/Hols94d.pdf (page 12)
   p_foo = 332.06364;  // TODO: WTH did this come from???

}

//~GeometricFlow() { };

/*
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
   p_maxstep = gf.p_maxstep;
   p_epsilons = gf.p_epsilons;
   p_epsilonp = gf.p_epsilonp;
   p_radexp = gf.p_radexp;
   p_crevalue = gf.p_crevalue;
   p_density = gf.p_density;

   p_dcel = gf.p_dcel;
   p_comdata_deltax = gf.p_comdata_deltax;
   p_comdata_deltax = gf.p_comdata_deltay;
   p_comdata_deltax = gf.p_comdata_deltaz;
   p_comdata_pi = gf.p_comdata_pi;
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
   p_maxstep = gf->p_maxstep;
   p_epsilons = gf->p_epsilons;
   p_epsilonp = gf->p_epsilonp;
   p_radexp = gf->p_radexp;
   p_crevalue = gf->p_crevalue;
   p_density = gf->p_density;

   p_dcel = gf->p_dcel;
   p_comdata_deltax = gf->p_comdata_deltax;
   p_comdata_deltax = gf->p_comdata_deltay;
   p_comdata_deltax = gf->p_comdata_deltaz;
   p_comdata_pi = gf->p_comdata_pi;
}
*/

//
//  operator=
//
//GeometricFlow& GeometricFlow::operator=( const GeometricFlow& gf )
//{
//}

void printAllParams()
{
}

void GeometricFlow::setup( const AtomList& atomList )
{
   //
   //  initialize the domain
   //
   domainInitialization( atomList );
	std::cout << "dimensions:\t" << p_comdata.nx() << " " <<
      p_comdata.ny() << " " << p_comdata.nz() << std::endl;
	std::cout << "grid spacing:\t" << p_comdata.deltax() << " " <<
      p_comdata.deltay() << " " << p_comdata.deltaz() << std::endl;

  
   //
   //  assign charge distributions
   //
   unsigned int natm = atomList.size();
   Mat<> charget(natm, 8);
	Mat<> corlocqt(natm, 8, 3);
	Mat<size_t> loc_qt(natm,8,3);
   atomList.changeChargeDistribution( charget, corlocqt, loc_qt, p_comdata );

   //
   //  compute nx, ny, nz
   //
   /*
	vector< double > xc( p_comdata.nx() );
	for (int i=0; i < p_comdata.nx(); i++) {
		xc[i] = p_comdata.xleft() + (i-1) * p_dcel;
	}

	vector< double > yc( p_comdata.ny() );
	for (int i=0; i < p_comdata.ny(); i++) {
		yc[i] = p_comdata.yleft() + (i-1) * p_dcel;
	}

	vector< double > zc( p_comdata.nz() );
	for (int i=0; i < p_comdata.nz(); i++) {
		zc[i] = p_comdata.zleft() + (i-1) * p_dcel;
	}
   */


   //
   //  setup phi
   //
	Mat<> phi( p_comdata.nx(), p_comdata.ny(), p_comdata.nz(),
              p_comdata.deltax(), p_comdata.deltay(), p_comdata.deltaz() ), 
         phix( phi ), 
         phivoc( phi ), 
         surfu( phi ), 
         eps( phi ),
			bg( phi );

   vector<double> solv(p_maxstep+1);  // keep track of the solvation energies
	double diffEnergy = 100;

	int iloop = 0; 
	double tott = 0.0;
	double elec = 0.0, area = 0.0, volume = 0.0, attint = 0.0;
   double tpb = 0.0;  // time step calculation for total pb
   int iterf = 0, itert = 0; // iteration num for first iteration and total
   double potcoe = 1.0 / p_gama_i;
   double lj_roro = p_density / p_gama_i;
   double lj_conms = p_pres_i / p_gama_i;

   //
   // iteration coupling surface generation and poisson solver
   //
   while ( (iloop < p_maxstep) && (diffEnergy > p_crevalue) ) 
   {
      iloop++;
      double deltat = 0; //this is wrong for adi...
      if (!p_iadi) {
         deltat =
            pow(p_comdata.deltax()*p_comdata.deltay()*p_comdata.deltaz(), 2.0/3.0)/4.5;
      }
      std::cout << "deltat = " << deltat << std::endl;

      double totnow = p_tottf - iloop + 1;
      if (totnow > 1.0) {
         tott = totnow;
      } else {
         tott = 1.0;
      }

		area = volume = attint = 0.0;
//		yhsurface(xyzr, ljepsilon, natm, tott, deltat, phix, surfu, iloop, area,
//				    volume, attint, alpha, iadi, igfin);
//		normalizeSurfuAndEps(surfu, eps, epsilons, epsilonp);

		if (iloop == 1) {
//			seteqb(bg, xyzr, pqr, charget, corlocqt, epsilons);
		}

		int iter = 1000;
		double fpb, titer = 0.0;
//		pbsolver(eps, phi, bg, tol, iter);
		if (iloop == 1) {
			fpb = titer;
			iterf = iter;
		}
		tpb = tpb + titer;
		itert += iter;

      //
      //  call the PB Solver!
      //
		eps = p_epsilonp;
		if (iloop == 1) {
//			pbsolver(eps, phivoc, bg, tol, iter);
		}

		for (size_t ix = 2; ix <= p_comdata.nx() - 1; ix++) {
			for (size_t iy = 2; iy <= p_comdata.ny() - 1; iy++) {
				for (size_t iz = 2; iz <= p_comdata.nz() - 1; iz++) {
					double phixx = phi(ix+1,iy,iz) -
                  phi(ix-1,iy,iz)/(2.0*p_comdata.deltax());
					double phixy = phi(ix,iy+1,iz) -
                  phi(ix,iy-1,iz)/(2.0*p_comdata.deltay());
					double phixz = phi(ix,iy,iz+1) -
                  phi(ix,iy,iz-1)/(2.0*p_comdata.deltaz());

					phix(ix,iy,iz) = 0.5 * (p_epsilons - p_epsilonp) * (phixx *
							phixx +	phixy * phixy + phixz * phixz) * potcoe;
				}
			}
		}

		// solvation
		cout << "iloop = " << iloop << std::endl;
		double soleng1, soleng2;
		soleng1 = soleng2 = 0.0;
//		computeSoleng(soleng1, phi, charget, loc_qt);
//		computeSoleng(soleng2, phivoc, charget, loc_qt);
		std::cout << "soleng1 = " << soleng1 << std::endl;
		std::cout << "soleng2 = " << soleng2 << std::endl;
		//std::cout << "soleng2 is too small!!" << std::endl;  // why is this here?
		elec = (soleng1 - soleng2) * p_foo;
		solv[iloop - 1] = elec + p_gama_i * (area + volume * lj_conms + attint *
				lj_roro);
		if (iloop > 1) {
//			diffEnergy = fabs((solv[iloop - 1] - solv[iloop - 2]));
		}
      
      // print the solvation energies by loop index; want to only print
      // the last two energies???  this prints them all.
      for (int i = 0; i < iloop; i++ )
      {
         std::cout << "solv[" << i << "] = " << solv[i] << std::endl;
      }
		std::cout << "diffEnergy = " << diffEnergy << std::endl;
	}



	double sumpot = area + volume * lj_conms + attint * lj_roro;
	double nonpolarSolvation = sumpot*p_gama_i;
	double totalSolvation = nonpolarSolvation + elec;
   cout << "totalSolv:\t" << totalSolvation << "\t";
   cout << "nonpolar: " << nonpolarSolvation << "\t";
   cout << "electro: " << elecSolvation << "\n" << std::endl;

}

void GeometricFlow::domainInitialization( const AtomList& atomList )
{
   unsigned int natm = atomList.size();
    std::valarray<double> atom_x(natm), 
                          atom_y(natm),
                          atom_z(natm), 
                          atom_r(natm);

	for(size_t i = 0; i < natm; ++i) 
   {
		atom_x[i] = atomList.get(i).x(); //xyzr[i][0];
		atom_y[i] = atomList.get(i).y(); //xyzr[i][1];
		atom_z[i] = atomList.get(i).z(); //xyzr[i][2];
		atom_r[i] = atomList.get(i).r(); //xyzr[i][3];
	}

	double xleft = left(atom_x - atom_r, p_comdata.deltax(), p_extvalue);
	double yleft = left(atom_y - atom_r, p_comdata.deltay(), p_extvalue);
	double zleft = left(atom_z - atom_r, p_comdata.deltaz(), p_extvalue);

	double xright = right(atom_x + atom_r, p_comdata.deltax(), p_extvalue);
	double yright = right(atom_y + atom_r, p_comdata.deltay(), p_extvalue);
	double zright = right(atom_z + atom_r, p_comdata.deltaz(), p_extvalue);

	int nx = (xright - xleft)/p_comdata.deltax() + 1;
	int ny = (yright - yleft)/p_comdata.deltay() + 1;
	int nz = (zright - zleft)/p_comdata.deltaz() + 1;

	xright = xleft + p_comdata.deltax()*(nx - 1);
	yright = yleft + p_comdata.deltay()*(ny - 1);
	zright = zleft + p_comdata.deltaz()*(nz - 1);

	//set the stupid globals...
   p_comdata.setBounds( xleft, xright, 
                        yleft, yright,
                        zleft, zright,
                        nx, ny, nz );

}


//		yhsurface(xyzr, ljepsilon, natm, tott, deltat, phix, surfu, iloop, area,
//				    volume, attint, alpha, iadi, igfin);
void GeometricFlow::yhsurface( const AtomList& atomList,
     double* ljepsilon, 
		double tott, double dt, Mat<>& phitotx, Mat<>& surfu, int iloop,
		double& area, double& volume, double& attint, double alpha, int iadi,
		int igfin)
{
   const int natm = atomList.size();
	size_t nx = p_comdata.nx(), ny = p_comdata.ny(), nz = p_comdata.nz();
	double xl = p_comdata.xleft(), yl = p_comdata.yleft(), zl = p_comdata.zleft();
	std::valarray<double> atom_x(natm), atom_y(natm), atom_z(natm), atom_r(natm);
	for (size_t i = 0; i < natm; ++i) {
		atom_x[i] = xyzr[i][0];
		atom_y[i] = xyzr[i][1];
		atom_z[i] = xyzr[i][2];
		atom_r[i] = xyzr[i][3];
	}

	Mat<> su(surfu), g(surfu);
	initial(xl, yl, zl, natm, atom_x,atom_y, atom_z, atom_r, g, su);
	if (iloop > 1 && igfin == 1)
		su = surfu;

	double rcfactor = (p_ffmodel == 1) ? 1.0 : pow(2.0, 1.0/6.0);
	std::valarray<double> sigma(atom_r);
	std::valarray<double> seta12(natm), seta6(natm), epsilon(natm);

	if (p_ffmodel == 1) {
		for (size_t i = 0; i < natm; ++i) {
			sigma[i] = atom_r[i] + p_sigmas;
			if (p_vdwdispersion != 0) {
				double se = sigma[i]/(atom_r[i] + p_prob);
				epsilon[i] = pow( pow(se, 12.0) - 2.0*pow(se, 6.0) , -1.0);
			}
			seta12[i] = lj.iosetar * p_vdwdispersion * epsilon[i];
			seta6[i] = 2.0*lj.iosetaa * p_vdwdispersion * epsilon[i];
		}
	} else {
		for (size_t i = 0; i < natm; ++i) {
			sigma[i] = sqrt(4.0* atom_r[i] * p_sigmas);
			if (p_vdwdispersion != 0) {
				epsilon[i] = sqrt(ljepsilon[i] * p_epsilonw);
				seta12[i] = 4.0*epsilon[i];
				seta6[i] = 4.0*epsilon[i];
			}
		}
	}

	Mat<> potr(nx,ny,nz), pota(nx,ny,nz);
	potIntegral(rcfactor, natm, atom_x, atom_y, atom_z, seta12, seta6,
			epsilon, sigma, g, potr, pota);

	if (lj.iwca == 1)
		potr = 0;

	for (size_t i = 0; i < phitotx.size(); ++i) {
		phitotx[i] = -lj.conms - phitotx[i] + lj.roro*(potr[i] + pota[i]);
	}

	if (iadi == 0 || iloop > 1) {
		int nt = ceil(tott/dt) + 1;
		upwinding(dt, nt, g, su, phitotx);
	} else {
		std::cerr << "ADI not implemented..." << std::endl;
		exit(1);
	}

	if (iloop > 1) {
		for (size_t i = 0; i < surfu.size(); ++i) {
		   surfu[i] = surfu[i]*alpha + su[i]*(1.0 - alpha);
		}
		su = surfu;
	} else {
		surfu = su;
	}

	volume = volumeIntegration(su);
	std::cout << "volume = " << volume << std::endl;

	Mat<> fintegr(nx,ny,nz);
	double dx = comdata.deltax, dy = comdata.deltay, dz = comdata.deltaz;
	for (size_t x = 2; x < nx; ++x) {
		for (size_t y = 2; y < ny; ++y) {
			for (size_t z = 2; z < nz; ++z) {
				double sux = su(x+1,y,z) - su(x-1,y,z);
				double suy = su(x,y+1,z) - su(x,y-1,z);
				double suz = su(x,y,z+1) - su(x,y,z-1);
				fintegr(x,y,z) = sqrt(dot(sux/(2.0*dx), suy/(2.0*dy),
						suz/(2.0*dz)));
			}
		}
	}

	area = volumeIntegration(fintegr);
	std::cout << "area = " << area << std::endl;

	potIntegral(rcfactor, natm, atom_x, atom_y, atom_z, seta12, seta6,
			epsilon, sigma, g, potr, pota);

	if (lj.iwca == 1) {
		for (size_t i = 0; i < fintegr.size(); ++i) {
			fintegr[i] = pota[i]*(1e3 - su[i]);
		}
	} else {
		for (size_t i = 0; i < fintegr.size(); ++i) {
			fintegr[i] = (pota[i] + potr[i])*(1e3 - su[i]);
		}
	}

	attint = volumeIntegration(fintegr);
	std::cout << "attint = " << attint << std::endl;
}
