//Kai Hoffmann, ICE-IEEC, 2016
//code measures shapes of dm particle distributions in mice-fof groups by diagonalizing the reduced moment of interia
//shapes and orientations are measured in 3d and in 2d on a plain perpendicular the observers line of sight


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>  // memcpy
#include <cmath>
// #include <gsl/gsl_eigen.h>
#include <chrono>

#include "halo_energy.h"
#include "recentering.h"
#include "compute_profile.h"
#include "calculate_shapes.h"
#include "project_particles.h"
// #include "save_coordinates.h"
#include "transform_coordinates.h"
#include "io.h"
#include "make_z_table.h"
#include "pos_to_z.h"
#include "read_sidm_simu.h"

#define PATH_PREFIX "/mnt/simulations/SIDM_simus/Lentes/V2/CDM/"

#define NRINGS 25

using namespace std;


//---------------------------------------------------------------------------------
//----------------------------------- main code -----------------------------------
int main(int argc, char **argv){

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
/*
    //eat input and output filenames
    string filename_input = argv[1];
    string filename_output = argv[2];

    //delimiter for output
    string delim = ",";

    //open output file for catalog
    ofstream outdata;
    string out_file_main = filename_output+"_main.csv";
    outdata.open(out_file_main.c_str());

    //open output file to save profiles
    ofstream outdata_pro;
    string out_file_pro = filename_output+"_pro.csv";
    outdata_pro.open(out_file_pro.c_str());

    //set format for output
    outdata.setf(ios::fixed);
    outdata.precision(3);
    outdata_pro.setf(ios::fixed);
    outdata_pro.precision(3);


    //-------------------- open input file --------------------
    ifstream indata;
    indata.open(filename_input.c_str(), ios::in|ios::binary);

    if (indata.is_open()){
        cout<<"# read "<<filename_input<<endl;
        cout<<"# write "<<filename_output<<endl;
    } else {
        cout<<"####### ERROR: can't open "<<filename_input<<" #######"<<endl;

        return 0;
    }
    //---------------------------------------------------------

    //----- read catalog properties from input file-----
    float limitmass=0, mass=0;
    int Nparttot=0, nhalos=0;
    int length = 4; // length of entries in bytes
    char buffer[2*length];

    indata.read(buffer, length);
    // total number of halos in the file
    indata.read(reinterpret_cast<char*>(&nhalos), length);

    // only halos larger that limitmass were stored
    indata.read(reinterpret_cast<char*>(&limitmass), length);

    //total number of particles in all halos
    indata.read(reinterpret_cast<char*>(&Nparttot), length);

    indata.read(buffer, 2*length);

    cout << ("###########################") << endl;
    cout << ("TOTAL NUMBER OF HALOS = %d", nhalos) << endl;
    cout << ("Limit mass = %.1f", limitmass) << endl;
    cout << ("Total number of particles = %d", Nparttot) << endl;
    cout << ("---------------------------") << endl;
    //cout << ("WARNING") << endl;
    //cout << ("Computation will be performed for halos") << endl;
    //cout << ("with more than 10 particles") << endl;
    //cout << ("---------------------------") << endl;

    // print headers
    print_header(outdata);
    print_profile_header(outdata_pro);

    // Make table to compute z_h
    vector<double> z_vec, Dc_vec;
    make_table(z_vec, Dc_vec);
*/

    string filename_output;
    string path_prefix;
    int nhalos = 0;
    double z_halo = 0.;

    if (argc <= 1){
        cout << "check the arguments" << endl;
        exit(1);
    }
    else if (argc == 5){
        path_prefix = argv[1];
        nhalos = atoi(argv[2]);
        z_halo = stod(argv[3]);
        filename_output = argv[4];

        cout << "### Used parameters ###" << endl;
        cout << "PATH_PREFIX :" << path_prefix << endl;
        cout << "NHALOS: " << nhalos << endl;
        cout << "ZHALOS: " << z_halo << endl;
        cout << "FILE OUTPUT: " << filename_output << endl;
    }

    //delimiter for output
    string delim = ",";

    //open output file for catalog
    ofstream outdata;
    string out_file_main = filename_output + "_main.csv";
    outdata.open(out_file_main.c_str());

    //open output file to save profiles
    ofstream outdata_pro;
    string out_file_pro = filename_output + "_pro.csv";
    outdata_pro.open(out_file_pro.c_str());

    //set format for output
    outdata.setf(ios::fixed);
    outdata.precision(3);
    outdata_pro.setf(ios::fixed);
    outdata_pro.precision(3);

    // print headers
    print_header(outdata);
    print_profile_header(outdata_pro);

    //double mp = 2.927; //1.e10 particle mass [M_sun/h]
    double mp = 0.013398587;  // 1.e10 particle mass [M_sun/h]

    string path;

    float avance = 0.02, mass = 0.;

    float xc_fof=0., yc_fof=0., zc_fof=0., vxc=0., vyc=0., vzc=0.;
    float lm = 0.;

    // Define variables that will save semi-axis vectors
    // and modulus
    float a2D[2], b2D[2], a2Dr[2], b2Dr[2];
    float a2Dr_abs, b2Dr_abs, a2D_abs, b2D_abs;

    float a3D[3], b3D[3], c3D[3];
    float a3Dr[3], b3Dr[3], c3Dr[3];
    float a3Dr_abs, b3Dr_abs, c3Dr_abs;
    float a3D_abs, b3D_abs, c3D_abs;

    double a_t = 0.;
    double EKin = 0.;
    double EPot = 0.;

    float r_max = 0;
    float xc = 0;
    float yc = 0;
    float zc = 0;

    double J[3] = {0., 0., 0.};

    vector <double> R(25), ro(25), ro_E(25);
    vector <double> Sigma(25), Sigma_E(25);

    int Npart = -1;

    vector<double> z_vec, Dc_vec;
    make_table(z_vec, Dc_vec);

    //--------------- begin loop over halos --------------------
    for (int ihalo = 0; ihalo < nhalos; ihalo++) {

        path = path_prefix + "halo_" + to_string(ihalo) + ".hdf5";
        cout << path << endl;

        if((float(ihalo)/float(nhalos)) > avance){

            avance = avance + 0.02;

        }

        //--------------------- read data ---------------------
        //read halo properties
        Npart = read_int_attr(path, "Npart");

        xc_fof = read_float_attr(path, "x0")* 1.e3;
        yc_fof = read_float_attr(path, "y0")* 1.e3;
        zc_fof = read_float_attr(path, "z0")* 1.e3;

        vxc = read_float_attr(path, "Vx_halo");
        vyc = read_float_attr(path, "Vy_halo");
        vzc = read_float_attr(path, "Vz_halo");

        mass = read_float_attr(path, "Mvir");

        // indata.read(buffer, 2*length);

        lm = log10(mass);

        float *x_part_arr = (float *) read_dataset(path, "X");
        vector<float> x_part(Npart);
        memcpy(&x_part[0], &x_part_arr[0], Npart * sizeof(float));
        free(x_part_arr);

        float *y_part_arr = (float *) read_dataset(path, "Y");
        vector<float> y_part(Npart);
        memcpy(&y_part[0], &y_part_arr[0], Npart * sizeof(float));
        free(y_part_arr);

        float *z_part_arr = (float *) read_dataset(path, "Z");
        vector<float> z_part(Npart);
        memcpy(&z_part[0], &z_part_arr[0], Npart * sizeof(float));
        free(z_part_arr);

        float *vx_part_arr = (float *) read_dataset(path, "vel x");
        vector<float> vx_part(Npart);
        memcpy(&vx_part[0], &vx_part_arr[0], Npart * sizeof(float));
        free(vx_part_arr);

        float *vy_part_arr = (float *) read_dataset(path, "vel y");
        vector<float> vy_part(Npart);
        memcpy(&vy_part[0], &vy_part_arr[0], Npart * sizeof(float));
        free(vy_part_arr);

        float *vz_part_arr = (float *) read_dataset(path, "vel z");
        vector<float> vz_part(Npart);
        memcpy(&vz_part[0], &vz_part_arr[0], Npart * sizeof(float));
        free(vz_part_arr);

        for (int i = 0; i < Npart; i++) {

            x_part[i] = (x_part[i]* 1.e3 - xc_fof);
            y_part[i] = (y_part[i]* 1.e3 - yc_fof);
            z_part[i] = (z_part[i]* 1.e3 - zc_fof);

        }

        //read particle velocities
        for (int i = 0; i < Npart; i++) {

            vx_part[i] = vx_part[i] - vxc;
            vy_part[i] = vy_part[i] - vyc;
            vz_part[i] = vz_part[i] - vzc;

        }

        // indata.read(buffer, 2*length);
        //-----------------------------------------------------


        if(Npart < 250000){

            // COMPUTE HALO REDSHIFT
            z_halo = 0;
            get_z(xc_fof, yc_fof, zc_fof, z_vec, Dc_vec, &z_halo);
            a_t = 1./(1.+ z_halo);

            // COMPUTE KINETIC AND POTENTIAL ENERGIES
            EKin = 0.;
            EPot = 0.;
            halo_energy(x_part, y_part, z_part, vx_part, vy_part, vz_part, a_t, mp, &EPot, &EKin);

            // RECENTER THE HALO
            r_max = 0.;
            xc = 0.;
            yc = 0.;
            zc = 0.;

            recenter(xc_fof, yc_fof, zc_fof, x_part, y_part, z_part, &xc, &yc, &zc, &r_max);

            //PROJECT POSITION OF PARTICLES
            vector <float> x_part_proj, y_part_proj;
            project(x_part, y_part, z_part, xc, yc, zc, x_part_proj, y_part_proj);

            // COMPUTE SEMI-AXIS USING INTERTIAL TENSOR

            calculate_2d_shapes(x_part_proj, y_part_proj, a_t, \
                                    a2D, b2D, a2Dr, b2Dr, \
                                    &a2Dr_abs, &b2Dr_abs, &a2D_abs, &b2D_abs);

            calculate_3d_shapes(x_part, y_part, z_part,\
                            a_t, a3D, b3D, c3D, \
                            a3Dr, b3Dr, c3Dr, \
                            &a3D_abs, &b3D_abs, &c3D_abs, \
                            &a3Dr_abs, &b3Dr_abs, &c3Dr_abs);

            //---------------------- angular momentum -------------------------------

            // fill J array with zeros
            memset(J, 0, 3 * sizeof(double));

            //sum cross-products J = r x p, where p = mv and m = 1 for each particle, i.e. Jtot=sum(part)
            for(int i = 0; i < Npart; i++){
                J[0] = J[0] + (a_t * y_part[i] * vz_part[i] - a_t * z_part[i] * vy_part[i]);
                J[1] = J[1] - (a_t * x_part[i] * vz_part[i] - a_t * z_part[i] * vx_part[i]);
                J[2] = J[2] + (a_t * x_part[i] * vy_part[i] - a_t * y_part[i] * vx_part[i]);
            }

            //specific angular momentum
            J[0] /= double(Npart);
            J[1] /= double(Npart);
            J[2] /= double(Npart);

            //-----------------------------------------------------------------------
            save_output(outdata, ihalo, Npart, mass, \
                        xc_fof, yc_fof, zc_fof, z_halo, r_max, \
                        xc, yc, zc, \
                        vxc, vyc, vzc, \
                        J, EKin, EPot, \
                        a2D_abs, b2D_abs, a2Dr_abs, b2Dr_abs, \
                        a2D, b2D, a2Dr,  b2Dr, \
                        a3Dr_abs, b3Dr_abs, c3Dr_abs, \
                        a3D_abs, b3D_abs, c3D_abs, \
                        a3D, b3D, c3D, a3Dr, b3Dr, c3Dr);

            // TRANSFORM COORDINATES
            vector <float> x_rot, y_rot, z_rot;
            vector <float> x2d_rot, y2d_rot;

            transform_coordinates(x_part,y_part,z_part, x_rot, y_rot, z_rot,
                                  x_part_proj, y_part_proj, x2d_rot, y2d_rot,
                                  a3D, b3D, c3D, a2D, b2D);

            // save_coordinates(ihalo, x_part, y_part, z_part, x_rot, y_rot, z_rot,
            //                 x_part_proj, y_part_proj, x2d_rot, y2d_rot);


            // COMPUTE DENSITY PROFILE

            outdata_pro << ihalo << delim << r_max << delim;

            // 3D profile
            ro_r(x_part, y_part, z_part, a_t, NRINGS, r_max, R, ro, 1., 1., 1.);

            for (int k = 0; k < NRINGS; k++) {
                    outdata_pro << R[k] << delim;
            }

            for (int k = 0; k < NRINGS; k++) {
                    outdata_pro << ro[k] << delim;
            }

            // 3D elliptical profile
            ro_r(x_rot, y_rot, z_rot, a_t, NRINGS, r_max, R, ro_E, a3D_abs, b3D_abs, c3D_abs);

            for (int k = 0; k < NRINGS; k++) {
                    outdata_pro << ro_E[k] << delim;
            }

            // 2D profile
            Sigma_r(x_part_proj, y_part_proj, a_t, NRINGS, r_max, R, Sigma, 1., 1.);

            for (int k = 0; k < NRINGS; k++) {
                    outdata_pro << Sigma[k] << delim;
            }

            // 2D elliptical profile
            Sigma_r(x2d_rot, y2d_rot, a_t, NRINGS, r_max, R, Sigma_E, a2D_abs, b2D_abs);

            for (int k = 0; k < NRINGS-1; k++) {
                    outdata_pro << Sigma_E[k] << delim;
            }

            outdata_pro << Sigma_E[NRINGS-1] << endl;


        }

    }//--------------- end loop over halos --------------------

    // indata.close();
    outdata.close();
    outdata_pro.close();

    string cmd = "bzip2 " + out_file_main;
    system(cmd.c_str());

    string cmd_pro = "bzip2 " + out_file_pro;
    system(cmd_pro.c_str());

    cout << (">\n");

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Total TIME = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

    return 0;

}
