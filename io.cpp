#include "io.h"
#include <cmath>

using namespace std;

void print_header(ofstream &outdata){

    //delimiter for output
    string delim = ",";

    //--------------------------------------------------
    //------------------ print header -----------------------

    // For output params
    outdata <<
    "Halo_number" <<delim<< "Npart" <<delim<< "lMfof" <<delim<<

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
