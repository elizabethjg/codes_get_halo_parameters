//Kai Hoffmann, ICE-IEEC, 2016
//code measures shapes of dm particle distributions in mice-fof groups by diagonalizing the reduced moment of interia
//shapes and orientations are measured in 3d and in 2d on a plain perpendicular the observers line of sight


#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <gsl/gsl_eigen.h>
#include <chrono>

#include "halo_energy.h"
#include "recentering.h"
#include "compute_profile.h"
#include "calculate_shapes.h"
#include "project_particles.h"
//#include "save_coordinates.h"
#include "transform_coordinates.h"

using namespace std;


void print_header(ofstream &outdata);

void print_profile_header(ofstream &outdata);


void save_profile(ofstream &outdata_pro, vector <float> ro, int ihalo, float r_max);

void save_output(ofstream &outdata, int ihalo, int Npart, float mass, \
                float xc_fof, float yc_fof, float zc_fof, float r_max, \
                float xc, float yc, float zc, \
                float vxc, float vyc, float vzc, \
                double *J, double EKin, double EPot, \
                float a2D_abs, float b2D_abs, float a2Dr_abs, float b2Dr_abs, \
                float *a2D, float *b2D, float *a2Dr, float *b2Dr, \
                float a3Dr_abs, float b3Dr_abs, float c3Dr_abs, \
                float a3D_abs, float b3D_abs, float c3D_abs, \
                float *a3D, float *b3D, float *c3D, float *a3Dr, float *b3Dr, float *c3Dr);

//---------------------------------------------------------------------------------
//----------------------------------- main code -----------------------------------
int main(int argc, char **argv){

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    //eat input and output filenames
    string filename_input = argv[1];
    string filename_output = argv[2];

    //delimiter for output
    string delim = ",";

    //open output file for catalog
    ofstream outdata;
    outdata.open(filename_output.c_str());

    //open output file to save profiles
    ofstream outdata_pro;
    string out_file_pro = filename_output+"_profile";
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

    printf("########################### \n");
    printf("TOTAL NUMBER OF HALOS = %d \n", nhalos);
    printf("Limit mass = %.1f \n", limitmass);
    printf("Total number of particles = %d \n", Nparttot);
    printf("--------------------------- \n");
    printf("WARNING \n");
    printf("Computation will be performed for halos \n");
    printf("with more than 10 particles \n");
    printf("--------------------------- \n");

    // print headers
    print_header(outdata);
    print_profile_header(outdata_pro);

    float avance = 0.02;

    //--------------- begin loop over halos --------------------
    for (int ihalo = 0; ihalo < nhalos; ihalo++) {

        // Define variables that will save semi-axis vectors
        // and modulus
        float a2D[2], b2D[2], a2Dr[2], b2Dr[2];
        float a2Dr_abs, b2Dr_abs, a2D_abs, b2D_abs;

        float a3D[3], b3D[3], c3D[3];
        float a3Dr[3], b3Dr[3], c3Dr[3];
        float a3Dr_abs, b3Dr_abs, c3Dr_abs;
        float a3D_abs, b3D_abs, c3D_abs;

        if((float(ihalo)/float(nhalos))  > avance){

            printf("=");
            avance = avance + 0.02;

        }

        //--------------------- read data ---------------------
        //read halo properties
        int Npart = 0;
        unsigned int haloID=0;
        float xc_fof=0, yc_fof=0, zc_fof=0, vxc=0, vyc=0, vzc=0;

        //number of fof particles
        indata.read(reinterpret_cast<char*>(&Npart), length);
        //fof mass = particle mass * sNpart
        indata.read(reinterpret_cast<char*>(&mass), length);
        //x,y,z coordinates of fof center of mass in kpc/h
        indata.read(reinterpret_cast<char*>(&xc_fof), length);
        indata.read(reinterpret_cast<char*>(&yc_fof), length);
        indata.read(reinterpret_cast<char*>(&zc_fof), length);
        //x,y,z velocity components of fof center of mass in km/s
        indata.read(reinterpret_cast<char*>(&vxc), length);
        indata.read(reinterpret_cast<char*>(&vyc), length);
        indata.read(reinterpret_cast<char*>(&vzc), length);

        indata.read(buffer, 2*length);

        float lm = log10(mass);

        //read particle coordinates
        vector <float> x_part, y_part, z_part;

        for (int i = 0; i < Npart; i++) {
            float xi=0, yi=0, zi=0;

            indata.read(reinterpret_cast<char*>(&xi), length);
            indata.read(reinterpret_cast<char*>(&yi), length);
            indata.read(reinterpret_cast<char*>(&zi), length);

            xi = xi - xc_fof;
            yi = yi - yc_fof;
            zi = zi - zc_fof;

            x_part.push_back(xi);
            y_part.push_back(yi);
            z_part.push_back(zi);

        }

        indata.read(buffer, 2*length);

        //read particle velocities
        vector <float>  vx_part, vy_part, vz_part;
        for (int i = 0; i < Npart; i++) {
            float vxi=0, vyi=0, vzi=0;

            indata.read(reinterpret_cast<char*>(&vxi), length);
            indata.read(reinterpret_cast<char*>(&vyi), length);
            indata.read(reinterpret_cast<char*>(&vzi), length);

            vxi = vxi - vxc;
            vyi = vyi - vyc;
            vzi = vzi - vzc;

            vx_part.push_back(vxi);
            vy_part.push_back(vyi);
            vz_part.push_back(vzi);
        }

        indata.read(buffer, 2*length);
        //-----------------------------------------------------


        if(Npart > 3000.){


            // COMPUTE KENETIC AND POTENTIAL ENERGIES
            double EKin = 0;
            double EPot = 0;
            halo_energy(x_part, y_part, z_part, vx_part, vy_part, vz_part, &EPot, &EKin);

            // RECENTER THE HALO
            float r_max = 0;
            float xc = 0;
            float yc = 0;
            float zc = 0;

            recenter(xc_fof, yc_fof, zc_fof, x_part, y_part, z_part, &xc, &yc, &zc, &r_max);

            //PROJECT POSITION OF PARTICLES
            vector <float> x_part_proj, y_part_proj;
            project(x_part, y_part, z_part, xc, yc, zc, x_part_proj, y_part_proj);


            // COMPUTE SEMI-AXIS USING INTERTIAL TENSOR

            calculate_2d_shapes(x_part_proj, y_part_proj, \
                                    a2D, b2D, a2Dr, b2Dr, \
                                    &a2Dr_abs, &b2Dr_abs, &a2D_abs, &b2D_abs);

            calculate_3d_shapes(x_part, y_part, \
                            z_part, a3D, b3D, c3D, \
                            a3Dr, b3Dr, c3Dr, \
                            &a3D_abs, &b3D_abs, &c3D_abs, \
                            &a3Dr_abs, &b3Dr_abs, &c3Dr_abs);

            //---------------------- angular momentum -------------------------------
            //array for angular momentum vector
            double J[3];

            //initialize J
            for(int i = 0; i < 3; i++){ J[i] = 0; }

            //sum cross-products J = r x p, where p = mv and m = 1 for each particle, i.e. Jtot=sum(part)
            for(int i = 0; i < Npart; i++){
                J[0] = J[0] + (y_part[i] * vz_part[i] - z_part[i] * vy_part[i]);
                J[1] = J[1] - (x_part[i] * vz_part[i] - z_part[i] * vx_part[i]);
                J[2] = J[2] + (x_part[i] * vy_part[i] - y_part[i] * vx_part[i]);
            }

            //specific angular momentum
            for(int i = 0; i < 3; i++){ J[i] /= double(Npart); }
            //-----------------------------------------------------------------------

            // //------------------ print output -----------------------
            //
            // //axis ratios
            // //double q2D = b2D_abs / a2D_abs;
            // //double q2Dr = b2Dr_abs / a2Dr_abs;
            // //
            // //double q3D = b3D_abs / a3D_abs;
            // //double q3Dr = b3Dr_abs / a3Dr_abs;
            // //
            // //double s3D = c3D_abs / a3D_abs;
            // //double s3Dr = c3Dr_abs / a3Dr_abs;
            //
            // //cout <<
            // outdata <<
            //
            // ihalo <<delim<< Npart <<delim<< log10(mass) <<delim<<     //0,1,2
            //
            // //position
            // xc_fof <<delim<< yc_fof <<delim<< zc_fof <<delim<<        //3,4,5
            // xc     <<delim<< yc     <<delim<< zc <<delim<<            //6,7,8
            //
            // //max radius
            // r_max <<delim<<                                           //9
            //
            // //velocity
            // vxc <<delim<< vyc <<delim<< vzc <<delim<<                 //10,11,12
            //
            // //angular momentum
            // J[0] <<delim<< J[1] <<delim<< J[2] <<delim<<              //13,14,15
            //
            // //Energies
            // EKin <<delim<< EPot <<delim<<                             //16,17
            //
            // //2D
            // a2D_abs <<delim<< b2D_abs <<delim<<                       //18,19
            // a2D[0] <<delim<< a2D[1] <<delim<<                         //20,21
            // b2D[0] <<delim<< b2D[1] <<delim<<                         //22,23
            //
            // //2D (reduced)
            // a2Dr_abs <<delim<< b2Dr_abs <<delim<<                     //24,25
            // a2Dr[0] <<delim<< a2Dr[1] <<delim<<                       //26,27
            // b2Dr[0] <<delim<< b2Dr[1] <<delim<<                       //28,29
            //
            // //3D
            // a3D_abs <<delim<< b3D_abs <<delim<< c3D_abs <<delim<<     //30,31,32
            // a3D[0] <<delim<< a3D[1] <<delim<< a3D[2] <<delim<<        //33,34,35
            // b3D[0] <<delim<< b3D[1] <<delim<< b3D[2] <<delim<<        //36,37,38
            // c3D[0] <<delim<< c3D[1] <<delim<< c3D[2] <<delim<<        //39,40,41
            //
            // //3D (reduced)
            // a3Dr_abs <<delim<< b3Dr_abs <<delim<< c3Dr_abs <<delim<<  //42,43,44
            // a3Dr[0] <<delim<< a3Dr[1] <<delim<< a3Dr[2] <<delim<<     //45,46,47
            // b3Dr[0] <<delim<< b3Dr[1] <<delim<< b3Dr[2] <<delim<<     //48,49,50
            // c3Dr[0] <<delim<< c3Dr[1] <<delim<< c3Dr[2] <<            //51,52,53
            //
            // endl;
            // //-------------------------------------------------------
            save_output(outdata, ihalo, Npart, mass, \
                        xc_fof, yc_fof, zc_fof, r_max, \
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
                                  
            //save_coordinates(ihalo, x_part, y_part, z_part, x_rot, y_rot, z_rot,
            //                 x_part_proj, y_part_proj, x2d_rot, y2d_rot);


            // COMPUTE DENSITY PROFILE

            outdata_pro <<
            ihalo << delim << r_max << delim;
            
            int NRINGS = 15;
            
            // 3D profile
            vector <double> ro = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
            ro_r(x_part, y_part, z_part, NRINGS, r_max, ro, 1., 1., 1.);            
            for (int k = 0; k < NRINGS; k++) {
                    outdata_pro <<
                    ro[k] << delim;
            }

            // 3D elliptical profile
            vector <double> ro_E = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
            ro_r(x_rot, y_rot, z_rot, NRINGS, r_max, ro_E, a3D_abs, b3D_abs, c3D_abs);            
            for (int k = 0; k < NRINGS; k++) {
                    outdata_pro <<
                    ro_E[k] << delim; 
            }

            // 2D profile
            vector <double> Sigma = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
            Sigma_r(x_part_proj, y_part_proj, NRINGS, r_max, Sigma,1.,1.);
            for (int k = 0; k < NRINGS; k++) {
                    outdata_pro <<
                    Sigma[k] << delim;
            }
            
            // 2D elliptical profile
            vector <double> Sigma_E = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
            Sigma_r(x2d_rot, y2d_rot, NRINGS, r_max, Sigma_E, a2D_abs, b2D_abs);
            for (int k = 0; k < NRINGS-1; k++) {
                    outdata_pro <<
                    Sigma_E[k] << delim;
            }

            outdata_pro << Sigma_E[14] <<
            endl;


        }

    }//--------------- end loop over halos --------------------

    indata.close();
    outdata.close();
    outdata_pro.close();

    string cmd = "bzip2 " + filename_output;
    system(cmd.c_str());

    string cmd_pro = "bzip2 " + out_file_pro;
    system(cmd_pro.c_str());

    printf(">\n");

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::cout << "Total TIME = " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]" << std::endl;

    return 0;

}


void print_header(ofstream &outdata){

    //delimiter for output
    string delim = ",";

    //--------------------------------------------------
    //------------------ print header -----------------------

    // For output params
    outdata <<
    "Halo number" <<delim<< "Npart" <<delim<< "log10(mass)" <<delim<<

    //position
    "xc_fof" <<delim<< "yc_fof" <<delim<< "zc_fof" <<delim<<
    "xc_rc" <<delim<< "yc_rc" <<delim<< "zc_rc" <<delim<<

    //max radius
    "r_max" <<delim<<

    //velocity
    "vxc" <<delim<< "vyc" <<delim<< "vzc" <<delim<<

    //angular momentum
    "J0"  <<delim<< "J1" <<delim<< "J2" <<delim<<

    //Energies
    "EKin" <<delim<< "EPot" <<delim<<

    //2D
    "a2D_mod" <<delim<< "b2D_mod" <<delim<<
    "a2D_x  " <<delim<< "a2D_y  " <<delim<<
    "b2D_x  " <<delim<< "b2D_y  " <<delim<<

    //2D (reduced)
    "a2Dr_mod" <<delim<< "b2Dr_mod" <<delim<<
    "a2Dr_x  " <<delim<< "a2Dr_y  " <<delim<<
    "b2Dr_x  " <<delim<< "b2Dr_y  " <<delim<<

    //3D
    "a3D_mod" <<delim<< "b3D_mod" <<delim<< "c3D_mod" <<delim<<
    "a3D_x  " <<delim<< "a3D_y  " <<delim<< "a3D_z  " <<delim<<
    "b3D_x  " <<delim<< "b3D_y  " <<delim<< "b3D_z  " <<delim<<
    "c3D_x  "  <<delim<< "c3D_y  " <<delim<< "c3D_z  " <<delim<<

    //3D (reduced)
    "a3Dr_mod" <<delim<< "b3Dr_mod" <<delim<< "c3Dr_mod" <<delim<<
    "a3Dr_x  " <<delim<< "a3Dr_y  " <<delim<< "a3Dr_z  " <<delim<<
    "b3Dr_x  " <<delim<< "b3Dr_y  " <<delim<< "b3Dr_z  " <<delim<<
    "c3Dr_x  "  <<delim<< "c3Dr_y  " <<delim<< "c3Dr_z  " <<
    endl;
    //-------------------------------------------------------
}


void print_profile_header(ofstream &outdata_pro){

    //delimiter for output
    string delim = ",";

    // For output profile

    outdata_pro <<
    "Halo number" <<delim<< "r_max" <<delim<<
    "rho_0" <<delim<< "rho_1" <<delim<< "rho_2" <<delim<<
    "rho_3" <<delim<< "rho_4" <<delim<< "rho_5" <<delim<<
    "rho_6" <<delim<< "rho_7" <<delim<< "rho_8" <<delim<<
    "rho_9" <<delim<< "rho_10" <<delim<< "rho_11" <<delim<<
    "rho_12" <<delim<< "rho_13" <<delim<< "rho_14" <<delim<<
    "rhoE_0" <<delim<< "rhoE_1" <<delim<< "rhoE_2" <<delim<<
    "rhoE_3" <<delim<< "rhoE_4" <<delim<< "rhoE_5" <<delim<<
    "rhoE_6" <<delim<< "rhoE_7" <<delim<< "rhoE_8" <<delim<<
    "rhoE_9" <<delim<< "rhoE_10" <<delim<< "rhoE_11" <<delim<<
    "rhoE_12" <<delim<< "rhoE_13" <<delim<< "rhoE_14" <<delim<<
    "Sigma_0" <<delim<< "Sigma_1" <<delim<< "Sigma_2" <<delim<<
    "Sigma_3" <<delim<< "Sigma_4" <<delim<< "Sigma_5" <<delim<<
    "Sigma_6" <<delim<< "Sigma_7" <<delim<< "Sigma_8" <<delim<<
    "Sigma_9" <<delim<< "Sigma_10" <<delim<< "Sigma_11" <<delim<<
    "Sigma_12" <<delim<< "Sigma_13" <<delim<< "Sigma_14" <<delim<<
    "SigmaE_0" <<delim<< "SigmaE_1" <<delim<< "SigmaE_2" <<delim<<
    "SigmaE_3" <<delim<< "SigmaE_4" <<delim<< "SigmaE_5" <<delim<<
    "SigmaE_6" <<delim<< "SigmaE_7" <<delim<< "SigmaE_8" <<delim<<
    "SigmaE_9" <<delim<< "SigmaE_10" <<delim<< "SigmaE_11" <<delim<<
    "SigmaE_12" <<delim<< "SigmaE_13" <<delim<< "SigmaE_14" <<
    endl;

}

void save_output(ofstream &outdata, int ihalo, int Npart, float mass, \
                float xc_fof, float yc_fof, float zc_fof, float r_max, \
                float xc, float yc, float zc, \
                float vxc, float vyc, float vzc, \
                double *J, double EKin, double EPot, \
                float a2D_abs, float b2D_abs, float a2Dr_abs, float b2Dr_abs, \
                float *a2D, float *b2D, float *a2Dr, float *b2Dr, \
                float a3Dr_abs, float b3Dr_abs, float c3Dr_abs, \
                float a3D_abs, float b3D_abs, float c3D_abs, \
                float *a3D, float *b3D, float *c3D, float *a3Dr, float *b3Dr, float *c3Dr
                ){

    //delimiter for output
    string delim = ",";

    outdata <<

    ihalo <<delim<< Npart <<delim<< log10(mass) <<delim<<     //0,1,2

    //position
    xc_fof <<delim<< yc_fof <<delim<< zc_fof <<delim<<        //3,4,5
    xc     <<delim<< yc     <<delim<< zc <<delim<<            //6,7,8

    //max radius
    r_max <<delim<<                                           //9

    //velocity
    vxc <<delim<< vyc <<delim<< vzc <<delim<<                 //10,11,12

    //angular momentum
    J[0] <<delim<< J[1] <<delim<< J[2] <<delim<<              //13,14,15

    //Energies
    EKin <<delim<< EPot <<delim<<                             //16,17

    //2D
    a2D_abs <<delim<< b2D_abs <<delim<<                       //18,19
    a2D[0] <<delim<< a2D[1] <<delim<<                         //20,21
    b2D[0] <<delim<< b2D[1] <<delim<<                         //22,23

    //2D (reduced)
    a2Dr_abs <<delim<< b2Dr_abs <<delim<<                     //24,25
    a2Dr[0] <<delim<< a2Dr[1] <<delim<<                       //26,27
    b2Dr[0] <<delim<< b2Dr[1] <<delim<<                       //28,29

    //3D
    a3D_abs <<delim<< b3D_abs <<delim<< c3D_abs <<delim<<     //30,31,32
    a3D[0] <<delim<< a3D[1] <<delim<< a3D[2] <<delim<<        //33,34,35
    b3D[0] <<delim<< b3D[1] <<delim<< b3D[2] <<delim<<        //36,37,38
    c3D[0] <<delim<< c3D[1] <<delim<< c3D[2] <<delim<<        //39,40,41

    //3D (reduced)
    a3Dr_abs <<delim<< b3Dr_abs <<delim<< c3Dr_abs <<delim<<  //42,43,44
    a3Dr[0] <<delim<< a3Dr[1] <<delim<< a3Dr[2] <<delim<<     //45,46,47
    b3Dr[0] <<delim<< b3Dr[1] <<delim<< b3Dr[2] <<delim<<     //48,49,50
    c3Dr[0] <<delim<< c3Dr[1] <<delim<< c3Dr[2] <<            //51,52,53

    endl;
    //-------------------------------------------------------

}

void save_profile(ofstream &outdata_pro, vector <float> ro, int ihalo, float r_max){

    //delimiter for output
    string delim = ",";

    outdata_pro <<
    ihalo << delim<< r_max <<delim<<
    ro[0] <<delim<<  ro[1] <<delim<<  ro[2] <<delim<<
    ro[3] <<delim<<  ro[4] <<delim<<  ro[5] <<delim<<
    ro[6] <<delim<<  ro[7] <<delim<<  ro[8] <<delim<<
    ro[9] <<delim<<  ro[10] <<delim<< ro[11] <<delim<<
    ro[12] <<delim<< ro[13] <<delim<< ro[14] <<
    endl;

}
