#include "project_particles.h"
#include <cmath>

void project(const vector <float> x_part, const vector <float> y_part,
             const vector <float> z_part, const float xc,
             const float yc, const float zc,
             vector <float> &x_part_proj, vector <float> &y_part_proj){

    int Npart = x_part.size();
    //----------- project particles on tangential plain (perpendicular to observers line of sight) -----------
    //get ra & dec of center
    double ra_center = 0;

    if(yc > 0){
        ra_center = atan(xc/yc);
    }

    double dec_center = asin(zc/sqrt(xc*xc + yc*yc + zc*zc));

    //define normalized 3Dvectors which span tangential plain perpendicular to los vector pointing to halo
    //vector with constant latitude (no dependence on dec),
    double e1x = cos(ra_center);
    double e1y = -sin(ra_center);
    //double e1z =   0;//not needed

    //vector with constant longitude (perpendiculat to los and e1)
    double e2x = - sin(dec_center) * sin(ra_center);
    double e2y = - sin(dec_center) * cos(ra_center);
    double e2z = cos(dec_center);

    float xi, yi;

    //project halo particles on plain via scalar product

    for(int i = 0; i < Npart; i++){

        xi = e1x*x_part[i] + e1y*y_part[i];  // + e1z*z_part[i]
        yi = e2x*x_part[i] + e2y*y_part[i] + e2z*z_part[i];

        x_part_proj.push_back(xi);
        y_part_proj.push_back(yi);
    }
  }
    //-------------------------------------------------------------------------------------------------------
