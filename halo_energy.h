#ifndef HALO_ENERGY_H
#define HALO_ENERGY_H

#include <vector>

using namespace std;

 //calculate potential an kinetic energy of every particle

 void halo_energy(const vector <float> x, const vector <float> y,
                 const vector <float> z, const vector <float> vx,
                 const vector <float> vy, const vector <float> vz,
                 double *Epot_halo, double *Ekin_halo);

#endif
