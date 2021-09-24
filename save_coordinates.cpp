#include "save_coordinates.h"

#include <cmath>
#include <iostream>
#include <string> 
#include <fstream>

using namespace std;

void save_coordinates(const int ihalo,
        const vector <float> x, const vector <float> y, const vector <float> z,
        const vector <float> x_rot, const vector <float> y_rot, const vector <float> z_rot,
        const vector <float> x2d, const vector <float> y2d,
        const vector <float> x2d_rot, const vector <float> y2d_rot){

    int np = x.size();
    std::string id = std::to_string(ihalo);

    string delim = "    ";

    //open output file to save profiles
    ofstream outdata_coords;
    string out_file_coords = "../catalogs/ind_halos_lM/coords"+id;
    outdata_coords.open(out_file_coords.c_str());

    //set format for output
    outdata_coords.setf(ios::fixed);
    outdata_coords.precision(3);


    for (int k = 0; k < np; k++) {
        outdata_coords <<
        x[k] <<delim<<  y[k] <<delim<<  z[k] <<delim<<
        x_rot[k] <<delim<<  y_rot[k] <<delim<<  z_rot[k] <<delim<<
        x2d[k] <<delim<<  y2d[k] <<delim<<  
        x2d_rot[k] <<delim<<  y2d_rot[k] <<  
        endl;
    }

    outdata_coords.close();

}
