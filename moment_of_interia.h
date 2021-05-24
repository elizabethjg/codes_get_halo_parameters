#ifndef MOMENT_OF_INTERIA_H
#define MOMENT_OF_INTERIA_H

#include <vector>

void ini_MI_2D(const vector <float> x_part, const vector <float> y_part,\
                double MI[2*2], const string type);

void ini_MI_3D(const vector <float> x_part, const vector <float> y_part,\
                const vector <float> z_part, double MI[3*3], const string type);

#endif
