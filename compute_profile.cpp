#define _USE_MATH_DEFINES

#include <cmath>


void ro_r(const vector <float> x_part, const vector <float> y_part, const vector <float> z_part,
        const int NRINGS, const float MAX_DISTANCE, vector <float> ro){

    float V;
    int NPART = x_part.size();

    const float ring_width = MAX_DISTANCE / NRINGS;

    for (i = 0; i < NPART; i++){

        rsq = pow(x_part[i], 2) + pow(y_part[i], 2) + pow(z_part[i], 2)

        idx = (int)(rsq/ring_width);

        /* Este if hace falta para el caso en que la distancia sea
         * exactamente igual a la máxima distancia posible porque se
         * suma partícula en cada anillo si r1 <= r < r2.*/
        if (idx >= NRINGS){
            idx = NRINGS - 1;
        }

        ro[idx] += 1;
	}

    /* Chequeo del número de partículas por anillo y del total.*/
    int total = 0;
    for (i = 0; i < NRINGS; i++){
        total += ro[i];
    }

    if (total != NPART){
		printf("%d particle(s) are missing\n", NPART - total);
    }

    for (i = 0; i < NRINGS; i++){

        rin = ring_width * i
        V = (4./3.) * pi * (pow((rin + ring_witdth), 3) - pow(rin, 3))
        ro[i] = ri[i] / V

    }
}
