#include <vector>
#include <cmath>

#include "position_to_redshift.h"

using namespace std;

//--------------- comoving distance from redshift ---------------
double Dcom(const double Om_0,const double redshift){

    //initialize free parameters
    double zmax = redshift;
    double zmin = 0.;
    int Ndz = 10000;
    double dz = (zmax - zmin)/double(Ndz);

    double c = 2.99792458E5; //speed fo light in km/s

    double OL = 1.-Om_0;

    double integral = 0;

    //integrate
    for(int k = 0; k < Ndz; k++){

            double z_lo = zmin + (k)*dz;
            double z_hi = zmin + (k+1)*dz;
            double E_lo = sqrt(Om_0*pow(1+z_lo, 3) + OL);
            double E_hi = sqrt(Om_0*pow(1+z_hi, 3) + OL);

            integral = integral + 0.5*(1./E_lo + 1./E_hi);
    }

    return integral*c*dz/100.;//comov distance in Mpc/h
}
//---------------------------------------------------------------


//------- function for searching nearest bin in array for input value ----------
int get_low_neighbor_bin(const vector<double>& val_vec, const double val_interpol){

    short step = 0;

    //search values below and above input val
    int bin1 = 0;
    int bin3 =  val_vec.size()-1.;
    int bin2 = int((bin1 + bin3)/2.);

    while (bin1 < bin3-1) {
        if(val_vec[bin2] < val_interpol){

            bin1 = bin2;
        }
        else{
            bin3 = bin2;
        }

        bin2 = int((bin1 + bin3)/2.);

        step++;
    }

    return bin3-1;
}
//---------------------------------------------------------------



//------ function for getting for y at x from interpolation of input x and y vectors -------
double get_vec_interpol(const vector<double>& x_vec, const vector<double>& y_vec, const double x_interpol){

    int bin1 = get_low_neighbor_bin(x_vec, x_interpol);
    int bin2 = bin1 + 1;

    //interpolate between bin below and above input x
    double xi_interpol = (y_vec[bin2] - y_vec[bin1])/(x_vec[bin2] - x_vec[bin1]) * (x_interpol - x_vec[bin1]) + y_vec[bin1];// = a*x + b;

    return xi_interpol;
}
//-----------------------------------------------------------------------------------------

//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
// int main(int argc, char **argv){
//
//
//     double Om = 0.25;//MICE matter density
//
//
//     // make vectors of redshifts and comoving distances
//     // this needs to be done only once before looping over te halos
//
//     double zmin = 0.;
//     double zmax = 2.;
//     int Ndz = 1000;
//     double dz = (zmax - zmin)/double(Ndz);
//
//     vector<double> z_vec, Dc_vec;
//
//     for(int i=0; i<Ndz; i++){
//         double z = zmin + (i+0.5)*dz;
//         z_vec.push_back(z);
//         Dc_vec.push_back(Dcom(Om,z));
//     }
//
//
//     // some comoving galaxy positions from MICE in [Mpc/h] and corresponding true redshifts for testing
//     //double pos_x = 717.8, pos_y= 1444.6, pos_z = 1803.7; //z_cgal = 1.01901
//     double pos_x = 870.9, pos_y= 267.3, pos_z = 800.7; //z_cgal = 0.44484
//     //double pos_x = 785.3, pos_y= 156.3, pos_z = 100.3; //z_cgal = 0.28544
//     //-> redshifts from cosmohub are well recovered :-)
//
//     // interpolate redshift at comoving distance
//     // this can be done once for each halo, assuming that all halo particles have the same redhsift
//     double Dc = sqrt(pow(pos_x,2) + pow(pos_y,2) + pow(pos_z,2));//comoving distance
//     double z = get_vec_interpol(Dc_vec, z_vec, Dc); // redshift
//     cout << Dc <<"\t"<< z <<endl;
//
//     return 0;
// }
