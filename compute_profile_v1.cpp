#include "compute_profile.h"

#include <cmath>
#include <iostream>

using namespace std;

void ro_r(const vector <float> x, const vector <float> y, const vector <float> z,
        const int nrings, const float max_distance, vector <double> &ro,
        const float a, const float b, const float c){

    double mp = 2.927e10; //1 particle mass [M_sun/h]
    double pi = 3.141592653589793;
    float rin;

    float V; //Volumen de la cascara
    float rsq; 
    int npart = x.size(), idx = -1;

    float ring_width;
    ring_width = float(max_distance) / float(nrings);

    for (int i = 0; i < npart; i++){


        if(a == 1. && b == 1. && c == 1.){
            rsq = sqrt((x[i]*x[i]) + (y[i]*y[i]) + (z[i]*z[i]));
        }
        else{
            rsq = pow(((x[i]*x[i])*((c*b)/a) + (y[i]*y[i])*((a*c)/b) + (z[i]*z[i])*((a*b)/c)),1./3.);
        }
            
        idx = (int)(rsq/ring_width);

        if (idx <= (nrings-1)){
            ro[idx] += 1;
        }


	}

    ring_width = ring_width/1000.; //Change units from kpc to Mpc

    for (int i = 0; i < nrings; i++){

        rin = ring_width * i;
        V = (4./3.) * pi * (pow((rin + ring_width), 3) - pow(rin, 3)); //In units of Mpc3/h3
        ro[i] = (mp*ro[i]) / V; //In units of (M_sun h2)/Mpc3

    }
}


void Sigma_r(const vector <float> x, const vector <float> y,
        const int nrings, const float max_distance, vector <double> &Sigma,
        const float a, const float b){

    double mp = 2.927e10; //1 particle mass [M_sun/h]
    double pi = 3.141592653589793;
    float rin;

    float A; //Volumen de la cascara
    float Rsq; 
    int npart = x.size(), idx = -1;

    float ring_width;
    ring_width = float(max_distance) / float(nrings);

    for (int i = 0; i < npart; i++){

        Rsq = sqrt((x[i]*x[i])*(b/a) + (y[i]*y[i])*(a/b));

        idx = (int)(Rsq/ring_width);

        if (idx <= (nrings-1)){
            Sigma[idx] += 1;
        }

	}

    ring_width = ring_width/1000.; //Change units from kpc to Mpc

    for (int i = 0; i < nrings; i++){

        rin = ring_width * i;
        A = pi * (pow((rin + ring_width), 2) - pow(rin, 2)); //In units of Mpc2
        Sigma[i] = (mp*Sigma[i]) / A; //In units of M_sun/(h*Mpc2)

    }
}
