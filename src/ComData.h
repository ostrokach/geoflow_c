#ifndef __ComData_h_
#define __ComData_h_

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

class ComData
{
   private:

      // initialize these in the constructor 
      double p_deltax ;
      double p_deltay ;
      double p_deltaz ;
      double p_pi;  //Pi

      // set these later
      double p_xleft;
      double p_xright;
      double p_yleft;
      double p_yright;
      double p_zleft;
      double p_zright;
      double p_nx;
      double p_ny;
      double p_nz;


   public:

      //
      //  empty constructor
      //
      ComData()
      {
         p_deltax = 1;
         p_deltay = 1;
         p_deltaz = 1;

         p_pi = acos(-1.0);  //Pi
      }

      void init( double dcel )
      {
         p_deltax = dcel;
         p_deltay = dcel;
         p_deltaz = dcel;
      }

      double pi() const { return p_pi; }

      double deltax() const { return p_deltax; }
      double deltay() const { return p_deltay; }
      double deltaz() const { return p_deltaz; }

      void setBounds( double xleft, double xright,
                      double yleft, double yright,
                      double zleft, double zright,
                      double nx, double ny, double nz )
      {
         p_xleft = xleft; p_xright = xright;
         p_yleft = yleft; p_yright = yright;
         p_zleft = zleft; p_zright = zright;
         p_nx = nx; p_ny = ny; p_nz = nz;
      }


      double nx() const { return p_nx; }
      double ny() const { return p_ny; }
      double nz() const { return p_nz; }

      double xleft() const { return p_xleft; }
      double yleft() const { return p_yleft; }
      double zleft() const { return p_zleft; }
      double xright() const { return p_xright; }
      double yright() const { return p_yright; }
      double zright() const { return p_zright; }

      double xvalue(size_t i) const { return (i - 1)*p_deltax + p_xleft; }
      double yvalue(size_t i) const { return (i - 1)*p_deltay + p_yleft; }
      double zvalue(size_t i) const { return (i - 1)*p_deltaz + p_zleft; }

      size_t inverx(double x) const { return size_t( (x - p_xleft)/p_deltax ) + 1; }
      size_t invery(double y) const { return size_t( (y - p_yleft)/p_deltay ) + 1; }
      size_t inverz(double z) const { return size_t( (z - p_zleft)/p_deltaz ) + 1; }

      void print() const
      {
         cout << "nx, ny, nz: " << p_nx 
            << ", " << p_ny << ", " << p_nz << endl;
         cout << "xleft, xright: " << p_xleft << ", " << p_xright <<
            endl;
         cout << "yleft, yright: " << p_yleft << ", " << p_yright <<
            endl;
         cout << "zleft, zright: " << p_zleft << ", " << p_zright <<
            endl;
         cout << "deltax, deltay, deltaz: " << p_deltax 
            << ", " << p_deltay << ", " << p_deltaz << endl;
         //cout << "dcel, pi: " << comdata.dcel << ", " << comdata.pi << std::endl;
      }
  
};


#endif

      

