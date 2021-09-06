#include "calculate_shapes.h"
#include "moment_of_inertia.h"

#include <cmath>
#include <gsl/gsl_eigen.h>

using namespace std;


void calculate_2d_shapes(const vector <float> x_part_proj, const vector <float> y_part_proj, \
                        const double a_t, float *a2D, float *b2D, float *a2Dr, float *b2Dr, \
                        float *a2Dr_abs, float *b2Dr_abs, float *a2D_abs, float *b2D_abs){

    //----- gsl variables for eigenvalue computation ------
    gsl_matrix_view M2D;
    gsl_vector *eval2D;
    gsl_matrix *evec2D;
    gsl_eigen_symmv_workspace *w2D;
    eval2D = gsl_vector_alloc(2);
    evec2D = gsl_matrix_alloc (2, 2);
    w2D = gsl_eigen_symmv_alloc (2);

    //--------------------------- 2D shapes ----------------------------------
    double MI_2D[2*2];


    //----- standard MI -----

    //initialize moment of inertia
    ini_MI_2D(x_part_proj, y_part_proj, a_t, MI_2D, "standard");

    //fill gsl matrix
    M2D = gsl_matrix_view_array (MI_2D, 2, 2);

    //get eigenvalues and eigenvectors  of MI
    gsl_eigen_symmv (&M2D.matrix, eval2D, evec2D, w2D);

    //sort eigenvalues and eigenvectors of MI in descending order
    gsl_eigen_symmv_sort (eval2D, evec2D, GSL_EIGEN_SORT_ABS_DESC);

    //get normalized eigenvectors from gsl
    a2D[0] = gsl_matrix_get(evec2D,(0),(0));
    a2D[1] = gsl_matrix_get(evec2D,(1),(0));

    b2D[0] = gsl_matrix_get(evec2D,(0),(1));
    b2D[1] = gsl_matrix_get(evec2D,(1),(1));

    //get eigenvalues from gsl
    *a2D_abs = sqrt(fabs(gsl_vector_get(eval2D,(0))));
    *b2D_abs = sqrt(fabs(gsl_vector_get(eval2D,(1))));

    //----- reduced MI -----

    //initialize moment of inertia 2DMi
    ini_MI_2D(x_part_proj, y_part_proj, a_t, MI_2D, "reduced");

    //fill gsl matrix
    M2D = gsl_matrix_view_array (MI_2D, 2, 2);

    //get eigenvalues and eigenvectors  of MI
    gsl_eigen_symmv (&M2D.matrix, eval2D, evec2D, w2D);

    //sort eigenvalues and eigenvectors of MI in descending order
    gsl_eigen_symmv_sort (eval2D, evec2D, GSL_EIGEN_SORT_ABS_DESC);

    //get normalized eigenvectors from gsl
    a2Dr[0] = gsl_matrix_get(evec2D,(0),(0));
    a2Dr[1] = gsl_matrix_get(evec2D,(1),(0));

    b2Dr[0] = gsl_matrix_get(evec2D,(0),(1));
    b2Dr[1] = gsl_matrix_get(evec2D,(1),(1));

    //get eigenvalues from gsl
    *a2Dr_abs = sqrt(fabs(gsl_vector_get(eval2D,(0))));
    *b2Dr_abs = sqrt(fabs(gsl_vector_get(eval2D,(1))));
    //------------------------------------------------------------------------
}

void calculate_3d_shapes(const vector <float> x_part, const vector <float> y_part, \
                const vector <float> z_part, const double a_t, \
                float *a3D, float *b3D, float *c3D, \
                float *a3Dr, float *b3Dr, float *c3Dr, \
                float *a3D_abs, float *b3D_abs, float *c3D_abs, \
                float *a3Dr_abs, float *b3Dr_abs, float *c3Dr_abs){

    //----- gsl variables for eigenvalue computation ------
    gsl_matrix_view M3D;
    gsl_vector *eval3D;
    gsl_matrix *evec3D;
    gsl_eigen_symmv_workspace *w3D;
    eval3D = gsl_vector_alloc(3);
    evec3D = gsl_matrix_alloc (3, 3);
    w3D = gsl_eigen_symmv_alloc (3);

    //--------------------------- 3D shapes ----------------------------------
    double MI_3D[3*3];

    //----- standard MI -----

    //initialize moment of inertia
    ini_MI_3D(x_part, y_part, z_part, a_t, MI_3D, "standard");

    //get eigenvalues and eigenvectors of MI using gsl
    M3D = gsl_matrix_view_array (MI_3D, 3, 3);
    gsl_eigen_symmv (&M3D.matrix, eval3D, evec3D, w3D);

    //sort eigenvalues and eigenvectors of MI in descending order
    gsl_eigen_symmv_sort (eval3D, evec3D, GSL_EIGEN_SORT_ABS_DESC);

    //get normalized eigenvectors from gsl
    a3D[0] = gsl_matrix_get(evec3D,(0),(0));
    a3D[1] = gsl_matrix_get(evec3D,(1),(0));
    a3D[2] = gsl_matrix_get(evec3D,(2),(0));

    b3D[0] = gsl_matrix_get(evec3D,(0),(1));
    b3D[1] = gsl_matrix_get(evec3D,(1),(1));
    b3D[2] = gsl_matrix_get(evec3D,(2),(1));

    c3D[0] = gsl_matrix_get(evec3D,(0),(2));
    c3D[1] = gsl_matrix_get(evec3D,(1),(2));
    c3D[2] = gsl_matrix_get(evec3D,(2),(2));

    //get eigenvalues from gsl
    *a3D_abs = sqrt(fabs(gsl_vector_get(eval3D,(0))));
    *b3D_abs = sqrt(fabs(gsl_vector_get(eval3D,(1))));
    *c3D_abs = sqrt(fabs(gsl_vector_get(eval3D,(2))));

    //----- reduced MI -----

    //initialize moment of inertia
    ini_MI_3D(x_part, y_part, z_part, a_t, MI_3D, "reduced");

    //get eigenvalues and eigenvectors of MI using gsl
    M3D = gsl_matrix_view_array (MI_3D, 3, 3);
    gsl_eigen_symmv (&M3D.matrix, eval3D, evec3D, w3D);

    //sort eigenvalues and eigenvectors of MI in descending order
    gsl_eigen_symmv_sort (eval3D, evec3D, GSL_EIGEN_SORT_ABS_DESC);

    //get normalized eigenvectors from gsl
    a3Dr[0] = gsl_matrix_get(evec3D,(0),(0));
    a3Dr[1] = gsl_matrix_get(evec3D,(1),(0));
    a3Dr[2] = gsl_matrix_get(evec3D,(2),(0));

    b3Dr[0] = gsl_matrix_get(evec3D,(0),(1));
    b3Dr[1] = gsl_matrix_get(evec3D,(1),(1));
    b3Dr[2] = gsl_matrix_get(evec3D,(2),(1));

    c3Dr[0] = gsl_matrix_get(evec3D,(0),(2));
    c3Dr[1] = gsl_matrix_get(evec3D,(1),(2));
    c3Dr[2] = gsl_matrix_get(evec3D,(2),(2));

    //get eigenvalues from gsl
    *a3Dr_abs = sqrt(fabs(gsl_vector_get(eval3D,(0))));
    *b3Dr_abs = sqrt(fabs(gsl_vector_get(eval3D,(1))));
    *c3Dr_abs = sqrt(fabs(gsl_vector_get(eval3D,(2))));
    //------------------------------------------------------------------------

}
