#ifndef CALCULATE_SHAPES_H
#define CALCULATE_SHAPES_H

#include <vector>

using namespace std;

void calculate_2d_shapes(const vector <float> x_part_proj, const vector <float> y_part_proj, \
                        const double a_t, float *a2D, float *b2D, float *a2Dr, float *b2Dr, \
                        float *a2Dr_abs, float *b2Dr_abs, float *a2D_abs, float *b2D_abs);

void calculate_3d_shapes(const vector <float> x_part, const vector <float> y_part, \
                const vector <float> z_part, const double a_t, \
                float *a3D, float *b3D, float *c3D, \
                float *a3Dr, float *b3Dr, float *c3Dr, \
                float *a3D_abs, float *b3D_abs, float *c3D_abs, \
                float *a3Dr_abs, float *b3Dr_abs, float *c3Dr_abs);

#endif
