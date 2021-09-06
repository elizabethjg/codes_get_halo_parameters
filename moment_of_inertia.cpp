#include "moment_of_inertia.h"

#include <cmath>

using namespace std;

//=================== 2D moment of interia ===================
void ini_MI_2D(const vector <float> x_part, const vector <float> y_part,\
                const double a_t, double MI[2*2], const string type){

    int Npart = x_part.size();

    for(int i = 0; i < 2*2; i++) MI[i] = 0;

    for (int i = 0; i < Npart; i++) {

        double rsq=1.;

        if(type=="reduced"){
            rsq = pow(x_part[i],2) + pow(y_part[i],2);
        }

        MI[0] += (pow(a_t,2) * x_part[i] * x_part[i] / rsq);
        MI[1] += (pow(a_t,2) * x_part[i] * y_part[i] / rsq);
                  
        MI[2] += (pow(a_t,2) * y_part[i] * x_part[i] / rsq);
        MI[3] += (pow(a_t,2) * y_part[i] * y_part[i] / rsq);
    }
}
//============================================================


//=================== 3D moment of interia ===================
void ini_MI_3D(const vector <float> x_part, const vector <float> y_part,\
                const vector <float> z_part, const double a_t, \
                double MI[3*3], const string type){

    int Npart = x_part.size();

    for(int i = 0; i < 3*3; i++) MI[i] = 0;


    //define moment of inertia MI and search max distance to center
    for (int i = 0; i < Npart; i++) {

        double rsq=1.;

        if(type=="reduced"){
            rsq = pow(x_part[i],2) + pow(y_part[i],2) + pow(z_part[i],2);
        }

        MI[0] += pow(at,2) * x_part[i] * x_part[i] / rsq;
        MI[1] += pow(at,2) * x_part[i] * y_part[i] / rsq;
        MI[2] += pow(at,2) * x_part[i] * z_part[i] / rsq;
                 
        MI[3] += pow(at,2) * y_part[i] * x_part[i] / rsq;
        MI[4] += pow(at,2) * y_part[i] * y_part[i] / rsq;
        MI[5] += pow(at,2) * y_part[i] * z_part[i] / rsq;
                 
        MI[6] += pow(at,2) * z_part[i] * x_part[i] / rsq;
        MI[7] += pow(at,2) * z_part[i] * y_part[i] / rsq;
        MI[8] += pow(at,2) * z_part[i] * z_part[i] / rsq;
    }
}
//============================================================
