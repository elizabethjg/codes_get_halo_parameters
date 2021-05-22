#include "halo_energy.h"
#include <cmath>

void halo_energy(const vector <float> x, const vector <float> y,
                const vector <float> z, const vector <float> vx,
                const vector <float> vy, const vector <float> vz,
                double *Epot_halo, double *Ekin_halo){

    double mp = 2.927e10; //particle mass [M_sun/h]
    double G = 6.67384e-11; // [m3 kg-1 s-2]
    double Msun = 1.9891e30; // [kg/M_sun]
    double h = 0.7;
    double Mpc = 3.08567758e22; // [m/Mpc]
    double fEkin = 0.5*h*Msun*1.e-40; //for energy unit conversion
    double fEpot = ((Msun*Msun*h*G*1.e-6) / Mpc) * 1.e-40;
    double Ekin_part, Epot_part = 0.;

    float xi, yi, zi;
    float dxi, dyi, dzi;

    int np = x.size();

    #pragma omp parallel for \
        private(xi,yi,zi,dxi,dyi,dzi)
    for (int j = 0; j < np; j++) {

        xi=x[j]; yi=y[j]; zi=z[j];

        // in units of (kg*km^2/s^2)*10^40
        Ekin_part = fEkin*mp*(vx[j]*vx[j] + vy[j]*vy[j] + vz[j]*vz[j]);

        for (int k = 0; k < np; k++) {

            if(k != j){

                // conversion from kpc/h to mpc/h
                dxi = (xi - x[k]) / 1000.;
                dyi = (yi - y[k]) / 1000.;
                dzi = (zi - z[k]) / 1000.;

                // in units of (kg*m^2/s^2)*10^40
                Epot_part = Epot_part + fEpot*mp*mp/sqrt(dxi*dxi + dyi*dyi + dzi*dzi);

            }

        }

        *Epot_halo=*Epot_halo + Epot_part;
        *Ekin_halo=*Ekin_halo + Ekin_part;


    }
    *Epot_halo = Epot_part;
    *Ekin_halo = Ekin_part;
}
            //-----------------------------------------------------------------------------------------
