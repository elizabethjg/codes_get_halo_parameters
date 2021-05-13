#include <stdlib.h>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

 //calculate potential an kinetic energy of every particle
 
      
void halo_energy(const vector <float> x, const vector <float> y, const vector <float> z, const vector <float> vx, const vector <float> vy, const vector <float> vz, float Epot_halo, float Ekin_halo){

        double mp = 2.927e10; //particle mass [M_sun/h]
        double G = 6.67384e-11; // [m3 kg-1 s-2]
        double Msun = 1.9891e30; // [kg/M_sun]
        double h = 0.7;
        double Mpc = 3.08567758e22; // [m/Mpc]
        double fEkin = 0.5*h*Msun*1.e-40; //for energy unit conversion
        double fEpot = ((Msun*Msun*h*G)/Mpc)*1.e-40;



       for (j = 0; j < np; j++) {
          
           xi=x[j]; yi=y[j]; zi=z[j];
          
           Ekin_part = fEkin*mp*(vx[j]*vx[j] + vy[j]*vy[j] + vz[j]*vz[j]); // in units of (kg*m^2/s^2)*10^40
          
           Epot_part = 0;
          
           for (k = 0; k < np; k++) {
          
               if(k != j){

                   dxi = xi -x[k]; dxi = dxi/1000.; // conversion from kpc/h to mpc/h
                   dyi = yi -y[k]; dyi = dyi/1000.;
                   dzi = zi -z[k]; dzi = dzi/1000.;

                   Epot_part = Epot_part + fEpot*mp*mp/sqrt(dxi*dxi + dyi*dyi + dzi*dzi);// in units of (kg*m^2/s^2)*10^40
                                                    
               }
          
           }
          
           Epot_halo=Epot_halo + Epot_part;
           Ekin_halo=Ekin_halo + Ekin_part;
          
       }
   }
            //-----------------------------------------------------------------------------------------
