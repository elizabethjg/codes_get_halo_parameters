#include <cmath>
#include <iostream>

#include "compute_profile.h"

using namespace std;

void ro_r(const vector <float> x, const vector <float> y, const vector <float> z, \
        const double a_t, const int nrings, const float max_distance, \
        vector <double> &R, vector <double> &ro, \
        const float a, const float b, const float c){

    double mp = 2.927e10; //1 particle mass [M_sun/h]
    // double pi = 3.141592653589793;
    float rin = 10.;
    float s = c/a;
    float q = b/a;

    float V; //Volumen de la cascara
    float rsq_in;
    float rsq_out;

    float a_in, b_in, c_in;
    float a_out, b_out, c_out;
    float npart;

    int Npart = x.size();
    int i = 0;

    float step;
    //step = (0.7*a_t*float(max_distance) - rin) / float(nrings);
    step = (1000. - rin) / float(nrings);
    float RMAX = 0.7*a_t*float(max_distance);

    if(RMAX < 220.){
        RMAX = 220.;
    }

    //for (int i = 0; i < nrings; i++){
    while((rin + step) < RMAX && i < nrings){

        a_in = rin/pow(q*s,1./3.);
        b_in = a_in*q;
        c_in = a_in*s;

        a_out = (rin+step)/pow(q*s,1./3.);
        b_out = a_out*q;
        c_out = a_out*s;

        npart = 0;

        for (int j = 0; j < Npart; j++){

            rsq_in  = pow(a_t, 2) * (pow(x[j], 2)/pow(a_in, 2) + pow(y[j], 2)/pow(b_in, 2) + pow(z[j], 2)/pow(c_in, 2));
            rsq_out = pow(a_t, 2) * (pow(x[j], 2)/pow(a_out, 2) + pow(y[j], 2)/pow(b_out, 2) + pow(z[j], 2)/pow(c_out, 2));

            if(rsq_in >= 1. && rsq_out < 1.){
                npart += 1;
            }
        }


        V = (4./3.) * M_PI * (pow((rin + step)/1.e3, 3) - pow(rin/1.e3, 3)); //In units of Mpc3/h3
        ro[i] = (mp*npart) / V; //In units of (M_sun h2)/Mpc3
        R[i]  = (rin + 0.5*step); //In units kpc/h

        rin += step;
        i += 1;

    }
}


void Sigma_r(const vector <float> x, const vector <float> y, const double a_t,
        const int nrings, const float max_distance, vector <double> &R,
        vector <double> &Sigma, const float a, const float b){

    double mp = 2.927e10; //1 particle mass [M_sun/h]
    // double pi = 3.141592653589793;
    float rin = 10.;
    float q = b/a;

    int Npart = x.size();

    float A; //Volumen de la cascara
    float rsq_in;
    float rsq_out;

    float a_in, b_in;
    float a_out, b_out;

    float npart;

    float step;
    //step = (0.7*a_t*float(max_distance) - rin) / float(nrings);
    step = (1000. - rin) / float(nrings);
    float RMAX = 0.7*a_t*float(max_distance);

    if(RMAX < 220.){
        RMAX = 220.;
    }

    int i = 0;

    //for (int i = 0; i < nrings; i++){
    while((rin + step) < RMAX && i < nrings){

        a_in = rin/sqrt(q);
        b_in = a_in*q;

        a_out = (rin+step)/sqrt(q);
        b_out = a_out*q;

        npart = 0;

        for (int j = 0; j < Npart; j++){

            rsq_in  = pow(a_t, 2) * (pow(x[j], 2)/pow(a_in, 2) + pow(y[j], 2)/pow(b_in, 2));
            rsq_out = pow(a_t, 2) * (pow(x[j], 2)/pow(a_out, 2) + pow(y[j], 2)/pow(b_out, 2));

            if(rsq_in >= 1. && rsq_out < 1.){
                npart += 1;
            }
        }

        A = M_PI * (pow((rin + step)/1.e3, 2) - pow(rin/1.e3, 2)); //In units of Mpc3/h3
        Sigma[i] = (mp*npart) / A; //In units of (M_sun h2)/Mpc3
        R[i]  = (rin + 0.5*step); //In units kpc/h

        rin += step;
        i += 1;

    }
}
