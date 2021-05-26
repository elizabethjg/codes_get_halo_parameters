#include "compute_profile.h"

#include <cmath>
#include <iostream>

using namespace std;

void ro_r(const vector <float> x, const vector <float> y, const vector <float> z,
        const int nrings, const float max_distance, vector <double> &ro){

    double mp = 2.927e10; //1 particle mass [M_sun/h]
    double pi = 3.141592653589793;
    float rin;

    float V; //Volumen de la cascara
    float rsq; 
    int npart = x.size(), idx = -1;

    float ring_width;
    ring_width = float(max_distance) / float(nrings);

    for (int i = 0; i < npart; i++){

        //rsq = pow(x_part[i], 2) + pow(y_part[i], 2) + pow(z_part[i], 2)
        rsq = sqrt(x[i]*x[i] + y[i]*y[i]+ z[i]*z[i]);

        idx = (int)(rsq/ring_width);

        /* Este if hace falta para el caso en que la distancia sea
         * exactamente igual a la máxima distancia posible porque se
         * suma partícula en cada anillo si r1 <= r < r2.*/
        if (idx >= nrings){
            idx = nrings - 1;
        }

        ro[idx] += 1;
	}

    /* Chequeo del número de partículas por anillo y del total.*/
    int total = 0;
    for (int i = 0; i < nrings; i++){
        total += ro[i];
    }

    if (total != npart){
		cout << (npart -total) << "particle(s) are missing" << endl;
    }

    ring_width = ring_width/1000.; //Change units from kpc to Mpc

    for (int i = 0; i < nrings; i++){

        rin = ring_width * i;
        V = (4./3.) * pi * (pow((rin + ring_width), 3) - pow(rin, 3)); //In units of Mpc3
        ro[i] = (mp*ro[i]) / V; //In units of M_sun/(h*Mpc3)

    }
}


void Sigma_r(const vector <float> x, const vector <float> y, 
        const int nrings, const float max_distance, vector <double> &Sigma){

    double mp = 2.927e10; //1 particle mass [M_sun/h]
    double pi = 3.141592653589793;
    float rin;

    float A; //Volumen de la cascara
    float Rsq; 
    int npart = x.size(), idx = -1;

    float ring_width;
    ring_width = float(max_distance) / float(nrings);

    for (int i = 0; i < npart; i++){

        //rsq = pow(x_part[i], 2) + pow(y_part[i], 2) + pow(z_part[i], 2)
        Rsq = sqrt(x[i]*x[i] + y[i]*y[i]);

        idx = (int)(Rsq/ring_width);

        /* Este if hace falta para el caso en que la distancia sea
         * exactamente igual a la máxima distancia posible porque se
         * suma partícula en cada anillo si r1 <= r < r2.*/
        if (idx >= nrings){
            idx = nrings - 1;
        }

        Sigma[idx] += 1;
	}

    /* Chequeo del número de partículas por anillo y del total.*/
    int total = 0;
    for (int i = 0; i < nrings; i++){
        total += Sigma[i];
    }

    if (total != npart){
		cout << (npart -total) << "particle(s) are missing" << endl;
    }

    ring_width = ring_width/1000.; //Change units from kpc to Mpc

    for (int i = 0; i < nrings; i++){

        rin = ring_width * i;
        A = pi * (pow((rin + ring_width), 2) - pow(rin, 2)); //In units of Mpc2
        Sigma[i] = (mp*Sigma[i]) / A; //In units of M_sun/(h*Mpc2)

    }
}
