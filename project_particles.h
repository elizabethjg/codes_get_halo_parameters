#ifndef PROJECT_H
#define PROJECT_H

#include <vector>

using namespace std;

void project(const vector <float> x_part, const vector <float> y_part,
             const vector <float> z_part, const float xc, 
             const float yc, const float zc,
             vector <float> &x_part_proj, vector <float> &y_part_proj);

#endif
