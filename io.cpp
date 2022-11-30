#include <cmath>

#include "io.h"

using namespace std;

void print_header(ofstream &outdata, int Np){

    //delimiter for output
    string delim = ",";
  	string ext[3] = { "xy", "zx", "yz" };

    //--------------------------------------------------
    //------------------ print header -----------------------

    // For output params
    outdata <<
    "column_halo_id" <<delim<< "Npart" <<delim<< "lgM" <<delim<<

    //position
    "xc" <<delim<< "yc" <<delim<< "zc" <<delim<<
    "xc_rc" <<delim<< "yc_rc" <<delim<< "zc_rc" <<delim<<

    //max radius
    "redshift" <<delim<<
    "r_max" <<delim<<

    //velocity
    "vxc" <<delim<< "vyc" <<delim<< "vzc" <<delim<<

    //angular momentum
    "Jx"  <<delim<< "Jy" <<delim<< "Jz" <<delim<<

    //Energies
    "EKin" <<delim<< "EPot";

     for(int iproj=0; iproj<Np; iproj++)
     {
       string tmp_delim = "_" + ext[iproj] + delim;
       
       outdata <<delim<<
       //2D
       "a2D"     <<tmp_delim<< "b2D"     <<tmp_delim<<
       "a2Dx"   <<tmp_delim<< "a2Dy"   <<tmp_delim<<
       "b2Dx"   <<tmp_delim<< "b2Dy"   <<tmp_delim<<

       //2D (reduced)
       "a2Dr"     <<tmp_delim<< "b2Dr"     <<tmp_delim<<
       "a2Drx"   <<tmp_delim<< "a2Dry"   <<tmp_delim<<
       "b2Drx"   <<tmp_delim<< "b2Dry"   <<tmp_delim<<

       //3D
       "a3D"     <<tmp_delim<< "b3D"     <<tmp_delim<< "c3D"     <<tmp_delim<<
       "a3Dx"   <<tmp_delim<< "a3Dy"   <<tmp_delim<< "a3Dz"   <<tmp_delim<<
       "b3Dx"   <<tmp_delim<< "b3Dy"   <<tmp_delim<< "b3Dz"   <<tmp_delim<<
       "c3Dx"   <<tmp_delim<< "c3Dy"   <<tmp_delim<< "c3Dz"   <<tmp_delim<<

       //3D (reduced)
       "a3Dr"     <<tmp_delim<< "b3Dr"     <<tmp_delim<< "c3Dr"     <<tmp_delim<<
       "a3Drx"   <<tmp_delim<< "a3Dry"   <<tmp_delim<< "a3Drz"   <<tmp_delim<<
       "b3Drx"   <<tmp_delim<< "b3Dry"   <<tmp_delim<< "b3Drz"   <<tmp_delim<<
       "c3Drx"   <<tmp_delim<< "c3Dry"   <<tmp_delim<< "c3Drz"   <<tmp_delim<<

       //2D
       "a2D_it"     <<tmp_delim<< "b2D_it"     <<tmp_delim<<
       "a2Dx_it"   <<tmp_delim<< "a2Dy_it"   <<tmp_delim<<
       "b2Dx_it"   <<tmp_delim<< "b2Dy_it"   <<tmp_delim<<

       //2D (reduced)
       "a2Dr_it"     <<tmp_delim<< "b2Dr_it"     <<tmp_delim<<
       "a2Drx_it"   <<tmp_delim<< "a2Dry_it"   <<tmp_delim<<
       "b2Drx_it"   <<tmp_delim<< "b2Dry_it"   <<tmp_delim<<

       //3D
       "a3D_it"    <<tmp_delim<< "b3D_it"    <<tmp_delim<< "c3D_it"    <<tmp_delim<<
       "a3Dx_it"   <<tmp_delim<< "a3Dy_it"   <<tmp_delim<< "a3Dz_it"   <<tmp_delim<<
       "b3Dx_it"   <<tmp_delim<< "b3Dy_it"   <<tmp_delim<< "b3Dz_it"   <<tmp_delim<<
       "c3Dx_it"   <<tmp_delim<< "c3Dy_it"   <<tmp_delim<< "c3Dz_it"   <<tmp_delim<<

       //3D (reduced)
       "a3Dr_it"    <<tmp_delim<< "b3Dr_it"    <<tmp_delim<< "c3Dr_it"    <<tmp_delim<<
       "a3Drx_it"   <<tmp_delim<< "a3Dry_it"   <<tmp_delim<< "a3Drz_it"   <<tmp_delim<<
       "b3Drx_it"   <<tmp_delim<< "b3Dry_it"   <<tmp_delim<< "b3Drz_it"   <<tmp_delim<<
       "c3Drx_it"   <<tmp_delim<< "c3Dry_it"   <<tmp_delim<< "c3Drz_it"   << "_" <<  ext[iproj];
    }    
    outdata << endl;
    //-------------------------------------------------------
}


void print_profile_header(ofstream &outdata_pro, int Np){

    //delimiter for output
    string delim = ",";
  	string ext[3] = { "xy", "zx", "yz" };

    // For output profile
    outdata_pro << "column_halo_id" << delim<< "r_max";

    for(int iproj=0; iproj<Np; iproj++)
    {
     	string tmp_delim = "_" + ext[iproj] + delim;
      
      outdata_pro <<delim<<
      "r_0" <<tmp_delim<< "r_1" <<tmp_delim<< "r_2" <<tmp_delim<<
      "r_3" <<tmp_delim<< "r_4" <<tmp_delim<< "r_5" <<tmp_delim<<
      "r_6" <<tmp_delim<< "r_7" <<tmp_delim<< "r_8" <<tmp_delim<<
      "r_9" <<tmp_delim<< "r_10" <<tmp_delim<< "r_11" <<tmp_delim<<
      "r_12" <<tmp_delim<< "r_13" <<tmp_delim<< "r_14" <<tmp_delim<<
      "r_15" <<tmp_delim<< "r_16" <<tmp_delim<< "r_17" <<tmp_delim<<
      "r_18" <<tmp_delim<< "r_19" <<tmp_delim<< "r_20" <<tmp_delim<<
      "r_21" <<tmp_delim<< "r_22" <<tmp_delim<< "r_23" <<tmp_delim<<
      "r_24" <<tmp_delim<<
      "rho_0" <<tmp_delim<< "rho_1" <<tmp_delim<< "rho_2" <<tmp_delim<<
      "rho_3" <<tmp_delim<< "rho_4" <<tmp_delim<< "rho_5" <<tmp_delim<<
      "rho_6" <<tmp_delim<< "rho_7" <<tmp_delim<< "rho_8" <<tmp_delim<<
      "rho_9" <<tmp_delim<< "rho_10" <<tmp_delim<< "rho_11" <<tmp_delim<<
      "rho_12" <<tmp_delim<< "rho_13" <<tmp_delim<< "rho_14" <<tmp_delim<<
      "rho_15" <<tmp_delim<< "rho_16" <<tmp_delim<< "rho_17" <<tmp_delim<<
      "rho_18" <<tmp_delim<< "rho_19" <<tmp_delim<< "rho_20" <<tmp_delim<<
      "rho_21" <<tmp_delim<< "rho_22" <<tmp_delim<< "rho_23" <<tmp_delim<<
      "rho_24" <<tmp_delim<<
      "rhoE_0" <<tmp_delim<< "rhoE_1" <<tmp_delim<< "rhoE_2" <<tmp_delim<<
      "rhoE_3" <<tmp_delim<< "rhoE_4" <<tmp_delim<< "rhoE_5" <<tmp_delim<<
      "rhoE_6" <<tmp_delim<< "rhoE_7" <<tmp_delim<< "rhoE_8" <<tmp_delim<<
      "rhoE_9" <<tmp_delim<< "rhoE_10" <<tmp_delim<< "rhoE_11" <<tmp_delim<<
      "rhoE_12" <<tmp_delim<< "rhoE_13" <<tmp_delim<< "rhoE_14" <<tmp_delim<<
      "rhoE_15" <<tmp_delim<< "rhoE_16" <<tmp_delim<< "rhoE_17" <<tmp_delim<<
      "rhoE_18" <<tmp_delim<< "rhoE_19" <<tmp_delim<< "rhoE_20" <<tmp_delim<<
      "rhoE_21" <<tmp_delim<< "rhoE_22" <<tmp_delim<< "rhoE_23" <<tmp_delim<<
      "rhoE_24" <<tmp_delim<<
      "Sigma_0" <<tmp_delim<< "Sigma_1" <<tmp_delim<< "Sigma_2" <<tmp_delim<<
      "Sigma_3" <<tmp_delim<< "Sigma_4" <<tmp_delim<< "Sigma_5" <<tmp_delim<<
      "Sigma_6" <<tmp_delim<< "Sigma_7" <<tmp_delim<< "Sigma_8" <<tmp_delim<<
      "Sigma_9" <<tmp_delim<< "Sigma_10" <<tmp_delim<< "Sigma_11" <<tmp_delim<<
      "Sigma_12" <<tmp_delim<< "Sigma_13" <<tmp_delim<< "Sigma_14" <<tmp_delim<<
      "Sigma_15" <<tmp_delim<< "Sigma_16" <<tmp_delim<< "Sigma_17" <<tmp_delim<<
      "Sigma_18" <<tmp_delim<< "Sigma_19" <<tmp_delim<< "Sigma_20" <<tmp_delim<<
      "Sigma_21" <<tmp_delim<< "Sigma_22" <<tmp_delim<< "Sigma_23" <<tmp_delim<<
      "Sigma_24" <<tmp_delim<<
      "SigmaE_0" <<tmp_delim<< "SigmaE_1" <<tmp_delim<< "SigmaE_2" <<tmp_delim<<
      "SigmaE_3" <<tmp_delim<< "SigmaE_4" <<tmp_delim<< "SigmaE_5" <<tmp_delim<<
      "SigmaE_6" <<tmp_delim<< "SigmaE_7" <<tmp_delim<< "SigmaE_8" <<tmp_delim<<
      "SigmaE_9" <<tmp_delim<< "SigmaE_10" <<tmp_delim<< "SigmaE_11" <<tmp_delim<<
      "SigmaE_12" <<tmp_delim<< "SigmaE_13" <<tmp_delim<< "SigmaE_14" <<tmp_delim<<
      "SigmaE_15" <<tmp_delim<< "SigmaE_16" <<tmp_delim<< "SigmaE_17" <<tmp_delim<<
      "SigmaE_18" <<tmp_delim<< "SigmaE_19" <<tmp_delim<< "SigmaE_20" <<tmp_delim<<
      "SigmaE_21" <<tmp_delim<< "SigmaE_22" <<tmp_delim<< "SigmaE_23" <<tmp_delim<<
      "SigmaE_24" << "_" <<  ext[iproj];
    }
    outdata_pro << endl;
}

void save_output(ofstream &outdata, int ihalo, int Npart, float mass, \
                float xc_fof, float yc_fof, float zc_fof, double z_h, float r_max, \
                float xc, float yc, float zc, \
                float vxc, float vyc, float vzc, \
                double *J, double EKin, double EPot
                ){

    //delimiter for output
    string delim = ",";

    outdata <<

    ihalo <<delim<< Npart <<delim<< log10(mass) <<delim<<     //0,1,2

    //position
    xc_fof <<delim<< yc_fof <<delim<< zc_fof <<delim<<        //3,4,5
    xc     <<delim<< yc     <<delim<< zc <<delim<<            //6,7,8

    //max radius
    z_h <<delim<<                                             //9
    r_max <<delim<<                                           //10

    //velocity
    vxc <<delim<< vyc <<delim<< vzc <<delim<<                 //11,12,13

    //angular momentum
    J[0] <<delim<< J[1] <<delim<< J[2] <<delim<<              //14,15,16

    //Energies
    EKin <<delim<< EPot <<delim;                             //17,18

    //-------------------------------------------------------

}

void save_output_iteration(ofstream &outdata, \
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
                ){

    //delimiter for output
    string delim = ",";

    outdata <<

    //2D
    a2D_abs <<delim<< b2D_abs <<delim<<                       //19,20
    a2D[0] <<delim<< a2D[1] <<delim<<                         //21,22
    b2D[0] <<delim<< b2D[1] <<delim<<                         //23,24

    //2D (reduced)
    a2Dr_abs <<delim<< b2Dr_abs <<delim<<                     //25,26
    a2Dr[0] <<delim<< a2Dr[1] <<delim<<                       //27,28
    b2Dr[0] <<delim<< b2Dr[1] <<delim<<                       //29,30

    //3D
    a3D_abs <<delim<< b3D_abs <<delim<< c3D_abs <<delim<<     //31,32,33
    a3D[0] <<delim<< a3D[1] <<delim<< a3D[2] <<delim<<        //34,35,36
    b3D[0] <<delim<< b3D[1] <<delim<< b3D[2] <<delim<<        //37,38,39
    c3D[0] <<delim<< c3D[1] <<delim<< c3D[2] <<delim<<        //40,41,42

    //3D (reduced)
    a3Dr_abs <<delim<< b3Dr_abs <<delim<< c3Dr_abs <<delim<<  //43,44,45
    a3Dr[0] <<delim<< a3Dr[1] <<delim<< a3Dr[2] <<delim<<     //46,47,48
    b3Dr[0] <<delim<< b3Dr[1] <<delim<< b3Dr[2] <<delim<<     //49,50,51
    c3Dr[0] <<delim<< c3Dr[1] <<delim<< c3Dr[2] <<delim<<     //52,53,54

    //2D iterative
    a2D_abs_it <<delim<< b2D_abs_it <<delim<<                       //19,20
    a2D_it[0] <<delim<< a2D_it[1] <<delim<<                         //21,22
    b2D_it[0] <<delim<< b2D_it[1] <<delim<<                         //23,24

    //2D iterative (reduced)
    a2Dr_abs_it <<delim<< b2Dr_abs_it <<delim<<                     //25,26
    a2Dr_it[0] <<delim<< a2Dr_it[1] <<delim<<                       //27,28
    b2Dr_it[0] <<delim<< b2Dr_it[1] <<delim<<                       //29,30

    //3D iterative
    a3D_abs_it <<delim<< b3D_abs_it <<delim<< c3D_abs_it <<delim<<     //31,32,33
    a3D_it[0] <<delim<< a3D_it[1] <<delim<< a3D_it[2] <<delim<<        //34,35,36
    b3D_it[0] <<delim<< b3D_it[1] <<delim<< b3D_it[2] <<delim<<        //37,38,39
    c3D_it[0] <<delim<< c3D_it[1] <<delim<< c3D_it[2] <<delim<<        //40,41,42

    //3D iterative (reduced)
    a3Dr_abs_it <<delim<< b3Dr_abs_it <<delim<< c3Dr_abs_it <<delim<<  //43,44,45
    a3Dr_it[0] <<delim<< a3Dr_it[1] <<delim<< a3Dr_it[2] <<delim<<     //46,47,48
    b3Dr_it[0] <<delim<< b3Dr_it[1] <<delim<< b3Dr_it[2] <<delim<<     //49,50,51
    c3Dr_it[0] <<delim<< c3Dr_it[1] <<delim<< c3Dr_it[2];              //52,53,54

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
