//g++ position_to_redshift.cpp -O3 -o position_to_redshift
#include "pos_to_z.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cmath>

using namespace std;

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
    int bin2 = bin1+1;
    
    //interpolate between bin below and above input x
    double xi_interpol = (y_vec[bin2] - y_vec[bin1])/(x_vec[bin2] - x_vec[bin1]) * (x_interpol - x_vec[bin1]) + y_vec[bin1];// = a*x + b;
    
    return xi_interpol;
}
//-----------------------------------------------------------------------------------------
void get_z(const double pos_x, const double pos_y, const double pos_z,
            const vector <double> z_vec, const vector <double> Dc_vec, double *zhalo){
                
    // some comoving galaxy positions from MICE in [Mpc/h] and corresponding true redshifts for testing
    double Dc = sqrt(pow(pos_x/1000.,2) + pow(pos_y/1000.,2) + pow(pos_z/1000.,2));//comoving distance
    double z = get_vec_interpol(Dc_vec, z_vec, Dc); // redshift
    
    *zhalo += z;
    
}
