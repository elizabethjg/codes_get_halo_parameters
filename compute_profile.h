#ifndef PROFILE_H
#define PROFILE_H

#include <vector>

using namespace std;

//--------------- Computes profiles ---------------
void ro_r(const vector <float> x, const vector <float> y, const vector <float> z,
        const int nrings, const float max_distance, vector <double> &R,
        vector <double> &ro, const float a, const float b, const float c);

void Sigma_r(const vector <float> x, const vector <float> y,
        const int nrings, const float max_distance, vector <double> &R,
        vector <double> &Sigma, const float a, const float b);

#endif
