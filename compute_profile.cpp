#include "compute_profile.h"

#include <cmath>
#include <iostream>

using namespace std;

void ro_r(const vector <float> x, const vector <float> y, const vector <float> z,
        const int nrings, const float max_distance, vector <double> &ro,
        const float a, const float b, const float c){

    double mp = 2.927e10; //1 particle mass [M_sun/h]
    double pi = 3.141592653589793;
    float rin = 0.;
    float s = c/a;
    float q = b/a;

    float V; //Volumen de la cascara
    float rsq_in; 
    float rsq_out; 

    float step;
    step = float(max_distance) / float(nrings);

    for (int i = 0; i < nrings; i++){       

        float a_in = rin/pow(q*s,1./3.);
        float b_in = a_in*q;
        float c_in = a_in*s;

        float a_out = (rin+step)/pow(q*s,1./3.);
        float b_out = a_out*q;
        float c_out = a_out*s;

        float npart = 0

        for (int j = 0; i < npart; j++){

            rsq_in  = (x[j]*x[j])/(a_in*a_in) + (y[j]*y[j])/(b_in*b_in) + (z[j]*z[j])/(c_in*c_in);
            rsq_out = (x[j]*x[j])/(a_out*a_out) + (y[j]*y[j])/(b_out*b_out) + (z[j]*z[j])/(c_out*c_out);

            if(rsq_in <= 1. && rsq_out < 1.){                    
                npart += 1;
            }
        }
         
        V = (4./3.) * pi * (pow((rin + step)/1.e3, 3) - pow(rin/1.e3, 3)); //In units of Mpc3/h3
        ro[i] = (mp*npart) / V; //In units of (M_sun h2)/Mpc3
        
        rin += step

    }
}


void Sigma_r(const vector <float> x, const vector <float> y,
        const int nrings, const float max_distance, vector <double> &Sigma,
        const float a, const float b){

    double mp = 2.927e10; //1 particle mass [M_sun/h]
    double pi = 3.141592653589793;
    float rin = 0.;
    float q = b/a;

    float A; //Volumen de la cascara
    float rsq_in; 
    float rsq_out; 

    float step;
    step = float(max_distance) / float(nrings);

    for (int i = 0; i < nrings; i++){       

        float a_in = rin/sqrt(q);
        float b_in = a_in*q;

        float a_out = (rin+step)/sqrt(q);
        float b_out = a_out*q;

        float npart = 0

        for (int j = 0; i < npart; j++){

            rsq_in  = (x[j]*x[j])/(a_in*a_in) + (y[j]*y[j])/(b_in*b_in);
            rsq_out = (x[j]*x[j])/(a_out*a_out) + (y[j]*y[j])/(b_out*b_out);

            if(rsq_in <= 1. && rsq_out < 1.){                    
                npart += 1;
            }
        }
         
        A = pi * (pow((rin + step)/1.e3, 2) - pow(rin/1.e3, 2)); //In units of Mpc3/h3
        ro[i] = (mp*npart) / A; //In units of (M_sun h2)/Mpc3

        rin += step

    }
}
