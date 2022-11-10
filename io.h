#ifndef IO_H
#define IO_H

#include <vector>
#include <fstream>

using namespace std;

void print_header(ofstream &outdata);

void print_profile_header(ofstream &outdata_pro);

void save_output(ofstream &outdata, int ihalo, int Npart, float mass, \
                float xc_fof, float yc_fof, float zc_fof, double z_h, float r_max, \
                float xc, float yc, float zc, \
                float vxc, float vyc, float vzc, \
                double *J, double EKin, double EPot, \
                float a2D_abs, float b2D_abs, float a2Dr_abs, float b2Dr_abs, \
                float *a2D, float *b2D, float *a2Dr, float *b2Dr, \
                float a3Dr_abs, float b3Dr_abs, float c3Dr_abs, \
                float a3D_abs, float b3D_abs, float c3D_abs, \
                float *a3D, float *b3D, float *c3D, float *a3Dr, float *b3Dr, float *c3Dr, \
                float a2D_abs_it, float b2D_abs_it, float a2Dr_abs_it, float b2Dr_abs_it, \
                float *a2D_it, float *b2D_it, float *a2Dr_it, float *b2Dr_it, \
                float a3Dr_abs_it, float b3Dr_abs_it, float c3Dr_abs_it, \
                float a3D_abs_it, float b3D_abs_it, float c3D_abs_it, \
                float *a3D_it, float *b3D_it, float *c3D_it, float *a3Dr_it, float *b3Dr_it, float *c3Dr_it
                );

void save_profile(ofstream &outdata_pro, vector <float> ro, int ihalo, float r_max);

#endif
