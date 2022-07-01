#include <cmath>

#include "make_z_table.h"

using namespace std;

//--------------- comoving distance from redshift ---------------
double Dcom(const double Om_0,const double redshift){

    // initialize free parameters
    double zmax = redshift;
    double zmin = 0.;
    int Ndz = 10000;
    double dz = (zmax - zmin)/double(Ndz);

    double c = 2.99792458E5;  //speed fo light in km/s

    double OL = 1.-Om_0;

    double integral = 0;

    double z_lo, z_hi, E_lo, E_hi;

    //integrate
    for(int k = 0; k < Ndz; k++){

            z_lo = zmin + (k)*dz;
            z_hi = zmin + (k+1)*dz;
            E_lo = sqrt(Om_0*pow(1+z_lo,3) + OL);
            E_hi = sqrt(Om_0*pow(1+z_hi,3) + OL);

            integral = integral + 0.5*(1./E_lo + 1./E_hi);
    }

    return integral*c*dz/100.;  // comov distance in Mpc/h
}

//--------------------------------------------------------------------------
void make_table(vector <double> &z_vec, vector <double> &Dc_vec){

    double Om = 0.25;  // MICE matter density

    // make vectors of redshifts and comoving distances
    // this needs to be done only once before looping over te halos

    double zmin = 0.;
    double zmax = 2.;
    int Ndz = 1000;
    double dz = (zmax - zmin)/double(Ndz);

    double z, Dout;

    z_vec.resize(Ndz);
    Dc_vec.resize(Ndz);

    for(int i=0; i<Ndz; i++){
        z = zmin + (i+0.5)*dz;
        Dout = Dcom(Om,z);
        //cout << Dout <<"\t"<< z <<endl;
        z_vec.push_back(z);
        Dc_vec.push_back(Dout);
    }
}
