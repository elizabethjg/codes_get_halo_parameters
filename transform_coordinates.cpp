#include <cmath>

#include "transform_coordinates.h"

using namespace std;

void transform_coordinates(const vector <float> x, const vector <float> y, const vector <float> z,
        vector <float> &x_rot, vector <float> &y_rot, vector <float> &z_rot,
        const vector <float> x2d, const vector <float> y2d,
        vector <float> &x2d_rot, vector <float> &y2d_rot,
        const float *a3D, const float *b3D, const float *c3D,
        const float *a2D, const float *b2D){

    int np = x.size();

    float xi = 0., yi = 0., zi = 0.;
    float x2di = 0., y2di = 0.;

    for (int k = 0; k < np; k++) {

        x2di = (a2D[0]*x2d[k])+(a2D[1]*y2d[k]);
        y2di = (b2D[0]*x2d[k])+(b2D[1]*y2d[k]);

        xi = (a3D[0]*x[k])+(a3D[1]*y[k])+(a3D[2]*z[k]);
        yi = (b3D[0]*x[k])+(b3D[1]*y[k])+(b3D[2]*z[k]);
        zi = (c3D[0]*x[k])+(c3D[1]*y[k])+(c3D[2]*z[k]);

        x2d_rot.push_back(x2di);
        y2d_rot.push_back(y2di);

        x_rot.push_back(xi);
        y_rot.push_back(yi);
        z_rot.push_back(zi);


    }

}
