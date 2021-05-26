#ifndef PROFILE_H
#define PROFILE_H

#include <vector>

using namespace std;

//--------------- recentering coordinates---------------
void ro_r(const vector <float> x, const vector <float> y, const vector <float> z,
        const int nrings, const float max_distance, vector <double> &ro);

#endif
