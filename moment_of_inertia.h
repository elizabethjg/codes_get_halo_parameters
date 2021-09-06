#ifndef MOMENT_OF_INERTIA_H
#define MOMENT_OF_INERTIA_H

#include <vector>
#include <string>

using namespace std;

void ini_MI_2D(const vector <float> x_part, const vector <float> y_part,\
                const double a_t, double MI[2*2], const string type);

void ini_MI_3D(const vector <float> x_part, const vector <float> y_part,\
                const vector <float> z_part, const double a_t, \
                double MI[3*3], const string type);

#endif
