#include "halo_energy.h"
#include <cmath>

void halo_energy(const vector <float> x, const vector <float> y,
                const vector <float> z, const vector <float> vx,
                const vector <float> vy, const vector <float> vz,
                const double a_t, const double mp,
                double *Epot_halo, double *Ekin_halo){

    double G = 6.67384; //1e-11 [m3 kg-1 s-2]
    double Msun = 1.989; //1.e30 [kg/M_sun]
    double h = 0.7;
    double Mpc = 3.08567758; //1e22 [m/Mpc]
    double fEkin = ((0.5*Msun)/h)* 1.e-4; //for energy unit conversion
    double fEpot = ((Msun*Msun*G) / (h*Mpc)) * 1.e-3;
    double Ekin_part, Epot_part, Ekin_acc = 0., Epot_acc=0.;

    float xi, yi, zi;
    float dxi, dyi, dzi;

    int np = x.size();

    #pragma omp parallel for \
        private(xi,yi,zi,dxi,dyi,dzi,Ekin_part,Epot_part) \
        reduction(+:Ekin_acc,Epot_acc)
    for (int j = 0; j < (np-1); j++) {

        xi=x[j]; yi=y[j]; zi=z[j];

        // in units of (kg*m^2/s^2)*10^30
        Ekin_part = fEkin*mp*(vx[j]*vx[j] + vy[j]*vy[j] + vz[j]*vz[j]);
        Epot_part = 0.;

        for (int k = (j+1); k < np; k++) {

            if(k != j){

                // conversion from kpc/h to mpc/h
                dxi = (xi - x[k]) / 1000.;
                dyi = (yi - y[k]) / 1000.;
                dzi = (zi - z[k]) / 1000.;


                // in units of (kg*m^2/s^2)*10^30
                Epot_part += fEpot*mp*mp/(a_t * sqrt((dxi*dxi + dyi*dyi + dzi*dzi)+pow(50./1000.,2)));

            }

        }

        Epot_acc += Epot_part;
        Ekin_acc += Ekin_part;

    }
    *Epot_halo += -1.*Epot_acc; 
    *Ekin_halo += Ekin_acc;
}
            //-----------------------------------------------------------------------------------------
