#include "save_coordinates.h"

#include <cmath>
#include <iostream>

using namespace std;

void save_coordinates(const vector <float> x, const vector <float> y, const vector <float> z,
        const vector <float> x_rot, const vector <float> y_rot, const vector <float> z_rot,
        const vector <float> x2d, const vector <float> y2d,
        const vector <float> x2d_rot, const vector <float> y2d_rot){

    int np = x.size();


    //open output file to save profiles
    ofstream outdata_coords;
    string out_file_coords = "coords";
    outdata_coords.open(out_file_coords.c_str());

    //set format for output
    outdata_coords.setf(ios::fixed);
    outdata_coords.precision(3);


    outdata_coords <<
    
    ro[0] <<delim<<  ro[1] <<delim<<  ro[2] <<delim<<
    ro[3] <<delim<<  ro[4] <<delim<<  ro[5] <<delim<<
    ro[6] <<delim<<  ro[7] <<delim<<  ro[8] <<delim<<
    ro[9] <<delim<<  ro[10] <<delim<< ro[11] <<delim<<
    ro[12] <<delim<< ro[13] <<delim<< ro[14] <<delim<<
    Sigma[0] <<delim<<  Sigma[1] <<delim<<  Sigma[2] <<delim<<
    Sigma[3] <<delim<<  Sigma[4] <<delim<<  Sigma[5] <<delim<<
    Sigma[6] <<delim<<  Sigma[7] <<delim<<  Sigma[8] <<delim<<
    Sigma[9] <<delim<<  Sigma[10] <<delim<< Sigma[11] <<delim<<
    Sigma[12] <<delim<< Sigma[13] <<delim<< Sigma[14] <<

    endl;


}
