#include "compute_profile.h"
#include <cmath>


void ro_r(const vector <float> x, const vector <float> y, const vector <float> z,
        const int nrings, const float max_distance, vector <float> &ro){

    double mp = 2.927e10; //1 particle mass [M_sun/h]

    float V; //Volumen de la cascara
    float rsq; //Volumen de la cascara
    int npart = x_part.size();

    const float ring_width;
    ring_width = float(max_distance) / float(nrings);

    for (i = 0; i < npart; i++){

        //rsq = pow(x_part[i], 2) + pow(y_part[i], 2) + pow(z_part[i], 2)
        rsq = sqrt(x[i]*x[i] + y[i]*y[i]+ z[i]*z[i])

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
    for (i = 0; i < nrings; i++){
        total += ro[i];
    }

    if (total != npart){
		printf("%d particle(s) are missing\n", npart - total);
    }

    ring_width = ring_width*1000. //Change units from kpc to pc

    for (i = 0; i < nrings; i++){

        rin = ring_width * i
        V = (4./3.) * pi * (pow((rin + ring_witdth), 3) - pow(rin, 3)) //In units of pc3
        ro[i] = (mp*ro[i]) / V //In units of M_sun/(h*pc3)

    }
}
