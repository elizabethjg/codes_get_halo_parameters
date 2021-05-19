#include <stdlib.h>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

//--------------- recentering coordinates---------------
void recenter(const float xc_fof, const float yc_fof, const float zc_fof,
            vector <float> &x, vector <float> &y, vector <float> &z,
            float *xc_rc, float *yc_rc, float *zc_rc, float *r_max){
