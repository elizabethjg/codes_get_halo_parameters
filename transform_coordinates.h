#ifndef TRANSFORM_COORDINATES_H
#define TRANSFORM_COORDINATES_H

#include <vector>

using namespace std;

 // Transform particle coordinates according to the elliptical axis


void transform_coordinates(const vector <float> x, const vector <float> y, const vector <float> z,
        vector <float> &x_rot, vector <float> &y_rot, vector <float> &z_rot,
        const vector <float> x2d, const vector <float> y2d,
        vector <float> &x2d_rot, vector <float> &y2d_rot,
        const float *a3D, const float *b3D, const float *c3D, 
        const float *a2D, const float *b2D);

#endif
