#ifndef COORDINATES_H
#define COORDINATES_H

#include <vector>

using namespace std;

void save_coordinates(const int ihalo,
        const vector <float> x, const vector <float> y, const vector <float> z,
        const vector <float> x_rot, const vector <float> y_rot, const vector <float> z_rot,
        const vector <float> x2d, const vector <float> y2d,
        const vector <float> x2d_rot, const vector <float> y2d_rot);
#endif
