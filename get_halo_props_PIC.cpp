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


using namespace std;


//=================== 2D moment of interia ===================
void ini_MI_2D(const vector <float> x_part, const vector <float> y_part, double MI[2*2], const string type){
    
    int Npart = x_part.size();
    
    for(int i = 0; i < 2*2; i++) MI[i] = 0;

    for (int i = 0; i < Npart; i++) {
        
        double rsq=1.;
        
        if(type=="reduced"){
            rsq = pow(x_part[i],2) + pow(y_part[i],2);
        }
        
        MI[0] = MI[0] + (x_part[i] * x_part[i] / rsq);
        MI[1] = MI[1] + (x_part[i] * y_part[i] / rsq);

        MI[2] = MI[2] + (y_part[i] * x_part[i] / rsq);
        MI[3] = MI[3] + (y_part[i] * y_part[i] / rsq);            
    }
}
//============================================================


//=================== 3D moment of interia ===================
void ini_MI_3D(const vector <float> x_part, const vector <float> y_part, const vector <float> z_part, double MI[3*3], const string type){
    
    int Npart = x_part.size();

    for(int i = 0; i < 3*3; i++) MI[i] = 0;
    

    //define moment of inertia MI and search max distance to center
    for (int i = 0; i < Npart; i++) {

        double rsq=1.;
        
        if(type=="reduced"){
            rsq = pow(x_part[i],2) + pow(y_part[i],2) + pow(z_part[i],2);
        }
        
        MI[0] += x_part[i] * x_part[i] / rsq;
        MI[1] += x_part[i] * y_part[i] / rsq;
        MI[2] += x_part[i] * z_part[i] / rsq;

        MI[3] += y_part[i] * x_part[i] / rsq;
        MI[4] += y_part[i] * y_part[i] / rsq;
        MI[5] += y_part[i] * z_part[i] / rsq;

        MI[6] += z_part[i] * x_part[i] / rsq;
        MI[7] += z_part[i] * y_part[i] / rsq;
        MI[8] += z_part[i] * z_part[i] / rsq;
    }
}
//============================================================


//---------------------------------------------------------------------------------
//----------------------------------- main code -----------------------------------
int main(int argc, char **argv){


//eat input and output filenames
string filename_input = argv[1];
string filename_output = argv[2];

//delimiter for output
string delim = ",";


//----- gsl variables for eigenvalue computation ------
gsl_matrix_view M2D;
gsl_vector *eval2D;
gsl_matrix *evec2D;
gsl_eigen_symmv_workspace *w2D;
eval2D = gsl_vector_alloc(2);
evec2D = gsl_matrix_alloc (2, 2);
w2D = gsl_eigen_symmv_alloc (2);

gsl_matrix_view M3D;
gsl_vector *eval3D;
gsl_matrix *evec3D;
gsl_eigen_symmv_workspace *w3D;
eval3D = gsl_vector_alloc(3);
evec3D = gsl_matrix_alloc (3, 3);
w3D = gsl_eigen_symmv_alloc (3);
//-----------------------------------------------------



//open output file for catalog
ofstream outdata;
outdata.open(filename_output.c_str());


//set format for output
outdata.setf(ios::fixed);
outdata.precision(3);



//-------------------- open input file --------------------
ifstream indata;
indata.open(filename_input.c_str(), ios::in|ios::binary);
 
if (indata.is_open()){
    cout<<"# read "<<filename_input<<endl;
    cout<<"# write "<<filename_output<<endl;
}else{
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
indata.read(reinterpret_cast<char*>(&nhalos), length);// total number of halos in the file
indata.read(reinterpret_cast<char*>(&limitmass), length);// only halos larger that limitmass were stored
indata.read(reinterpret_cast<char*>(&Nparttot), length);//total number of particles in all halos

indata.read(buffer, 2*length);

printf("########################### \n");
printf("TOTAL NUMBER OF HALOS = %d \n", nhalos);
printf("Limit mass = %.1f \n", limitmass);
printf("Total number of particles = %d \n", Nparttot);
printf("--------------------------- \n");

//--------------------------------------------------



//------------------ print header -----------------------


outdata <<

//"#" <<

"Npart" <<delim<< "lgM" <<delim<<

//position
"xc" <<delim<< "yc" <<delim<< "zc" <<delim<<

//velocity
"vxc" <<delim<< "vyc" <<delim<< "vzc" <<delim<<

//angular momentum
"Jx" <<delim<< "Jy" <<delim<< "Jz" <<delim<<

//axis ratios, major, semi-major, minor axes

//2D
"q2D" <<delim<<
"a2Dx" <<delim<< "a2Dy" <<delim<<
"b2Dx" <<delim<< "b2Dy" <<delim<<

//2D (reduced)
"q2Dr" <<delim<<
"a2Drx" <<delim<< "a2Dry" <<delim<<
"b2Drx" <<delim<< "b2Dry" <<delim<<

//3D
"q3D" <<delim<< "s3D" <<delim<<
"a3Dx" <<delim<< "a3Dy" <<delim<< "a3Dz" <<delim<<
"c3Dx" <<delim<< "c3Dy" <<delim<< "c3Dz" <<delim<<

//3D (reduced)
"q3Dr" <<delim<< "s3Dr" <<delim<<
"a3Drx" <<delim<< "a3Dry" <<delim<< "a3Drz" <<delim<<
"c3Drx" <<delim<< "c3Dry" <<delim<< "c3Drz" <<

endl;
//-------------------------------------------------------
        

//--------------- begin loop over halos --------------------
//for (int ihalo = 0; ihalo < nhalos; ihalo++) {
for (int ihalo = 0; ihalo < 10; ihalo++) {
    
    
    //--------------------- read data ---------------------
    //read halo properties
    int Npart = 0;
    unsigned int haloID=0;
    float xc=0, yc=0, zc=0, vxc=0, vyc=0, vzc=0;

    //indata.read(reinterpret_cast<char*>(&haloID),   length); //fof ID
    indata.read(reinterpret_cast<char*>(&Npart),   length); //number of fof particles
    indata.read(reinterpret_cast<char*>(&mass), length); //fof mass = particle mass * sNpart
    indata.read(reinterpret_cast<char*>(&xc),   length); //x,y,z coordinates of fof center of mass in kpc/h
    indata.read(reinterpret_cast<char*>(&yc),   length);
    indata.read(reinterpret_cast<char*>(&zc),   length);
    indata.read(reinterpret_cast<char*>(&vxc),  length); //x,y,z velocity components of fof center of mass in km/s
    indata.read(reinterpret_cast<char*>(&vyc),  length);
    indata.read(reinterpret_cast<char*>(&vzc),  length);
    
    indata.read(buffer, 2*length);

    float lm = log10(mass);
    
    printf("Halo %d: %d haloID\n", ihalo, haloID);
    printf("%d particles\n", Npart);
    printf("log(M_fof) = %.1f\n", lm);
    printf("xc = %.1f\n", xc_fof);
    printf("yc = %.1f\n", yc_fof);
    printf("zc = %.1f\n", zc_fof);
    printf("vxc = %.1f\n", vxc);
    printf("vyc = %.1f\n", vyc);
    printf("vzc = %.1f\n", vzc);



    //read particle coordinates
    vector <float> x_part, y_part, z_part;
    for (int i = 0; i < Npart; i++) {
        float xi=0, yi=0, zi=0;

        indata.read(reinterpret_cast<char*>(&xi), length);
        indata.read(reinterpret_cast<char*>(&yi), length);
        indata.read(reinterpret_cast<char*>(&zi), length);
    
        xi = xi - xc;
        yi = yi - yc;
        zi = zi - zc;
        
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

    
    if(Npart>0){
            
        
        //----------- project particles on tangential plain (perpendicular to observers line of sight) -----------
        //get ra & dec of center
        double ra_center = 0;
        if(yc > 0){ra_center = atan(xc/yc);}
        double dec_center = asin(zc/sqrt(xc*xc + yc*yc + zc*zc));
        
        //define normalized 3Dvectors which span tangential plain perpendicular to los vector pointing to halo
        //vector with constant latitude (no dependence on dec),
        double e1x =   cos(ra_center);
        double e1y =   -sin(ra_center);
        //double e1z =   0;//not needed

        //vector with constant longitude (perpendiculat to los and e1)
        double e2x = - sin(dec_center) * sin(ra_center);
        double e2y = - sin(dec_center) * cos(ra_center);
        double e2z =   cos(dec_center);

        //project halo particles on plain via scalar product
        vector <float> x_part_proj, y_part_proj;
        for(int i = 0; i < Npart; i++){
            
            float xi = e1x*x_part[i] + e1y*y_part[i]; // + e1z*z_part[i]
            float yi = e2x*x_part[i] + e2y*y_part[i] + e2z*z_part[i];
            
            x_part_proj.push_back(xi);
            y_part_proj.push_back(yi);
        }
        //-------------------------------------------------------------------------------------------------------

        
        
        //--------------------------- 2D shapes ----------------------------------
        double MI_2D[2*2];

        
        //----- standard MI -----
        
        //initialize moment of inertia
        ini_MI_2D(x_part_proj, y_part_proj, MI_2D, "standard");
        
        //fill gsl matrix
        M2D = gsl_matrix_view_array (MI_2D, 2, 2);

        //get eigenvalues and eigenvectors  of MI
        gsl_eigen_symmv (&M2D.matrix, eval2D, evec2D, w2D);
        
        //sort eigenvalues and eigenvectors of MI in descending order
        gsl_eigen_symmv_sort (eval2D, evec2D, GSL_EIGEN_SORT_ABS_DESC);

        //get normalized eigenvectors from gsl
        float a2D[2];
        a2D[0] = gsl_matrix_get(evec2D,(0),(0));
        a2D[1] = gsl_matrix_get(evec2D,(1),(0));

        float b2D[2];
        b2D[0] = gsl_matrix_get(evec2D,(0),(1));
        b2D[1] = gsl_matrix_get(evec2D,(1),(1));
        
        //get eigenvalues from gsl
        float a2D_abs = sqrt(fabs(gsl_vector_get(eval2D,(0))));
        float b2D_abs = sqrt(fabs(gsl_vector_get(eval2D,(1))));
        
        
        //----- reduced MI -----
        
        //initialize moment of inertia 2DMi
        ini_MI_2D(x_part_proj, y_part_proj, MI_2D, "reduced");
        
        //fill gsl matrix
        M2D = gsl_matrix_view_array (MI_2D, 2, 2);

        //get eigenvalues and eigenvectors  of MI
        gsl_eigen_symmv (&M2D.matrix, eval2D, evec2D, w2D);
        
        //sort eigenvalues and eigenvectors of MI in descending order
        gsl_eigen_symmv_sort (eval2D, evec2D, GSL_EIGEN_SORT_ABS_DESC);

        //get normalized eigenvectors from gsl
        float a2Dr[2];
        a2Dr[0] = gsl_matrix_get(evec2D,(0),(0));
        a2Dr[1] = gsl_matrix_get(evec2D,(1),(0));
        
        float b2Dr[2];
        b2Dr[0] = gsl_matrix_get(evec2D,(0),(1));
        b2Dr[1] = gsl_matrix_get(evec2D,(1),(1));
        
        //get eigenvalues from gsl
        float a2Dr_abs = sqrt(fabs(gsl_vector_get(eval2D,(0))));
        float b2Dr_abs = sqrt(fabs(gsl_vector_get(eval2D,(1))));
        //------------------------------------------------------------------------
        
        
        
        
        //--------------------------- 3D shapes ----------------------------------
        double MI_3D[3*3];

        //----- standard MI -----
        
        //initialize moment of inertia
        ini_MI_3D(x_part, y_part, z_part, MI_3D, "standard");
        
        //get eigenvalues and eigenvectors of MI using gsl
        M3D = gsl_matrix_view_array (MI_3D, 3, 3);
        gsl_eigen_symmv (&M3D.matrix, eval3D, evec3D, w3D);

        //sort eigenvalues and eigenvectors of MI in descending order
        gsl_eigen_symmv_sort (eval3D, evec3D, GSL_EIGEN_SORT_ABS_DESC);
        
        //get normalized eigenvectors from gsl
        float a3D[3];
        a3D[0] = gsl_matrix_get(evec3D,(0),(0));
        a3D[1] = gsl_matrix_get(evec3D,(1),(0));
        a3D[2] = gsl_matrix_get(evec3D,(2),(0));

        float b3D[3];
        b3D[0] = gsl_matrix_get(evec3D,(0),(1));
        b3D[1] = gsl_matrix_get(evec3D,(1),(1));
        b3D[2] = gsl_matrix_get(evec3D,(2),(1));

        float c3D[3];
        c3D[0] = gsl_matrix_get(evec3D,(0),(2));
        c3D[1] = gsl_matrix_get(evec3D,(1),(2));
        c3D[2] = gsl_matrix_get(evec3D,(2),(2));

        //get eigenvalues from gsl
        float a3D_abs = sqrt(fabs(gsl_vector_get(eval3D,(0))));
        float b3D_abs = sqrt(fabs(gsl_vector_get(eval3D,(1))));
        float c3D_abs = sqrt(fabs(gsl_vector_get(eval3D,(2))));
        
        
        //----- reduced MI -----
        
        //initialize moment of inertia
        ini_MI_3D(x_part, y_part, z_part, MI_3D, "reduced");
        
        //get eigenvalues and eigenvectors of MI using gsl
        M3D = gsl_matrix_view_array (MI_3D, 3, 3);
        gsl_eigen_symmv (&M3D.matrix, eval3D, evec3D, w3D);

        //sort eigenvalues and eigenvectors of MI in descending order
        gsl_eigen_symmv_sort (eval3D, evec3D, GSL_EIGEN_SORT_ABS_DESC);
        
        //get normalized eigenvectors from gsl
        float a3Dr[3];
        a3Dr[0] = gsl_matrix_get(evec3D,(0),(0));
        a3Dr[1] = gsl_matrix_get(evec3D,(1),(0));
        a3Dr[2] = gsl_matrix_get(evec3D,(2),(0));

        float b3Dr[3];
        b3Dr[0] = gsl_matrix_get(evec3D,(0),(1));
        b3Dr[1] = gsl_matrix_get(evec3D,(1),(1));
        b3Dr[2] = gsl_matrix_get(evec3D,(2),(1));

        float c3Dr[3];
        c3Dr[0] = gsl_matrix_get(evec3D,(0),(2));
        c3Dr[1] = gsl_matrix_get(evec3D,(1),(2));
        c3Dr[2] = gsl_matrix_get(evec3D,(2),(2));

        //get eigenvalues from gsl
        float a3Dr_abs = sqrt(fabs(gsl_vector_get(eval3D,(0))));
        float b3Dr_abs = sqrt(fabs(gsl_vector_get(eval3D,(1))));
        float c3Dr_abs = sqrt(fabs(gsl_vector_get(eval3D,(2))));
        //------------------------------------------------------------------------

        
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
        
        
        
        //------------------ print output -----------------------

        //axis ratios
        double q2D = b2D_abs / a2D_abs;
        double q2Dr = b2Dr_abs / a2Dr_abs;

        double q3D = b3D_abs / a3D_abs;
        double q3Dr = b3Dr_abs / a3Dr_abs;

        double s3D = c3D_abs / a3D_abs;
        double s3Dr = c3Dr_abs / a3Dr_abs;
        
        
        //cout <<
        outdata <<
        
        Npart <<delim<< log10(mass) <<delim<<                   //1,2
        
        //position
        xc <<delim<< yc <<delim<< zc <<delim<<                  //3,4,5
        
        //velocity
        vxc <<delim<< vyc <<delim<< vzc <<delim<<               //6,7,8
        
        //angular momentum
        J[0] <<delim<< J[1] <<delim<< J[2] <<delim<<            //9,10,11
        
        //2D
        q2D <<delim<<                                           //12
        a2D[0] <<delim<< a2D[1] <<delim<<                       //13,14
        b2D[0] <<delim<< b2D[1] <<delim<<                       //15,16

        //2D (reduced)
        q2Dr <<delim<<                                          //17
        a2Dr[0] <<delim<< a2Dr[1] <<delim<<                     //18,19
        b2Dr[0] <<delim<< b2Dr[1] <<delim<<                     //20,21

        //3D
        q3D <<delim<< s3D <<delim<<                             //22, 23
        a3D[0] <<delim<< a3D[1] <<delim<< a3D[2] <<delim<<      //24,25,26
        c3D[0] <<delim<< c3D[1] <<delim<< c3D[2] <<delim<<      //27,28,29
        
        //3D (reduced)
        q3Dr <<delim<< s3Dr <<delim<<                           //30,31
        a3Dr[0] <<delim<< a3Dr[1] <<delim<< a3Dr[2] <<delim<<   //32,33,34
        c3Dr[0] <<delim<< c3Dr[1] <<delim<< c3Dr[2] <<          //35,36,37
        
        endl;
        //-------------------------------------------------------
        
    }    
    
    
}//--------------- end loop over halos --------------------

indata.close();
outdata.close();

string cmd = "bzip2 " + filename_output;
system(cmd.c_str());

return 0;
}
