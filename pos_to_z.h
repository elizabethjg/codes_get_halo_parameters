#ifndef POS_TO_Z_H
#define POS_TO_Z_H

#include <vector>

using namespace std;

//--------------- Computes profiles ---------------
void get_z(const double pos_x, const double pos_y, const double pos_z,
           const vector <double> z_vec, const vector <double> Dc_vec, double *z);

#endif
