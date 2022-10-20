#include <cmath>
#include <gsl/gsl_eigen.h>

#include "calculate_shapes.h"
#include "moment_of_inertia.h"

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
    evec2D = gsl_matrix_alloc(2, 2);
    w2D = gsl_eigen_symmv_alloc(2);

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
    a2D[0] = gsl_matrix_get(evec2D, (0), (0));
    a2D[1] = gsl_matrix_get(evec2D, (1), (0));

    b2D[0] = gsl_matrix_get(evec2D, (0), (1));
    b2D[1] = gsl_matrix_get(evec2D, (1), (1));

    //get eigenvalues from gsl
    *a2D_abs = sqrt(fabs(gsl_vector_get(eval2D, (0))));
    *b2D_abs = sqrt(fabs(gsl_vector_get(eval2D, (1))));

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
    a2Dr[0] = gsl_matrix_get(evec2D, (0), (0));
    a2Dr[1] = gsl_matrix_get(evec2D, (1), (0));

    b2Dr[0] = gsl_matrix_get(evec2D, (0), (1));
    b2Dr[1] = gsl_matrix_get(evec2D, (1), (1));

    //get eigenvalues from gsl
    *a2Dr_abs = sqrt(fabs(gsl_vector_get(eval2D, (0))));
    *b2Dr_abs = sqrt(fabs(gsl_vector_get(eval2D, (1))));
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
    evec3D = gsl_matrix_alloc(3, 3);
    w3D = gsl_eigen_symmv_alloc(3);

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
    a3D[0] = gsl_matrix_get(evec3D, (0), (0));
    a3D[1] = gsl_matrix_get(evec3D, (1), (0));
    a3D[2] = gsl_matrix_get(evec3D, (2), (0));

    b3D[0] = gsl_matrix_get(evec3D, (0), (1));
    b3D[1] = gsl_matrix_get(evec3D, (1), (1));
    b3D[2] = gsl_matrix_get(evec3D, (2), (1));

    c3D[0] = gsl_matrix_get(evec3D, (0), (2));
    c3D[1] = gsl_matrix_get(evec3D, (1), (2));
    c3D[2] = gsl_matrix_get(evec3D, (2), (2));

    //get eigenvalues from gsl
    *a3D_abs = sqrt(fabs(gsl_vector_get(eval3D, (0))));
    *b3D_abs = sqrt(fabs(gsl_vector_get(eval3D, (1))));
    *c3D_abs = sqrt(fabs(gsl_vector_get(eval3D, (2))));

    //----- reduced MI -----

    //initialize moment of inertia
    ini_MI_3D(x_part, y_part, z_part, a_t, MI_3D, "reduced");

    //get eigenvalues and eigenvectors of MI using gsl
    M3D = gsl_matrix_view_array (MI_3D, 3, 3);
    gsl_eigen_symmv (&M3D.matrix, eval3D, evec3D, w3D);

    //sort eigenvalues and eigenvectors of MI in descending order
    gsl_eigen_symmv_sort (eval3D, evec3D, GSL_EIGEN_SORT_ABS_DESC);

    //get normalized eigenvectors from gsl
    a3Dr[0] = gsl_matrix_get(evec3D, (0), (0));
    a3Dr[1] = gsl_matrix_get(evec3D, (1), (0));
    a3Dr[2] = gsl_matrix_get(evec3D, (2), (0));

    b3Dr[0] = gsl_matrix_get(evec3D, (0), (1));
    b3Dr[1] = gsl_matrix_get(evec3D, (1), (1));
    b3Dr[2] = gsl_matrix_get(evec3D, (2), (1));

    c3Dr[0] = gsl_matrix_get(evec3D, (0), (2));
    c3Dr[1] = gsl_matrix_get(evec3D, (1), (2));
    c3Dr[2] = gsl_matrix_get(evec3D, (2), (2));

    //get eigenvalues from gsl
    *a3Dr_abs = sqrt(fabs(gsl_vector_get(eval3D, (0))));
    *b3Dr_abs = sqrt(fabs(gsl_vector_get(eval3D, (1))));
    *c3Dr_abs = sqrt(fabs(gsl_vector_get(eval3D, (2))));
    //------------------------------------------------------------------------

}

/*
// bound = 1 - dismiss
// SHAPE_ITERATIONS = 10
// FORCE_RES - preguntar Agus (buscar)
// r -> radio del halo
// Despreciar K < U

void calc_shape(struct halo *h, int64_t total_p, int64_t bound) {
  int64_t i,j,k,l,iter=SHAPE_ITERATIONS, analyze_p=0, a,b,c;
  float b_to_a, c_to_a, min_r = FORCE_RES*FORCE_RES;
  double mass_t[3][3], orth[3][3], eig[3]={0},  r=0, dr, dr2, weight=0;
  h->b_to_a = h->c_to_a = 0;
  memset(h->A, 0, sizeof(float)*3);

  if (!(h->r>0)) return;
  
  min_r *= 1e6 / (h->r*h->r); // fuerza escalada por el radio (epsilon/radio)**2
  
  
  // CUENTA EL NUMERO DE PARTICULAS QUE VA A USAR (despreciar)
  for (j=0; j<total_p; j++) {
    if (bound && (po[j].pe < po[j].ke)) continue;
    analyze_p++; 
  }
  
  if (analyze_p < 3 || !(h->r>0)) return;
  if (analyze_p < iter) iter = analyze_p;
  /////////////


  //INICIALIZA
  for (i=0; i<3; i++) {
    memset(orth[i], 0, sizeof(double)*3);
    orth[i][i] = 1;
    eig[i] = (h->r*h->r)*1e-6;
  }

  // EMPIEZA A ITERAR  
  for (i=0; i<iter; i++) {
    
    // inicializa mass_t (momento de inercia) 
    for (k=0; k<3; k++) memset(mass_t[k], 0, sizeof(double)*3);
    weight=0;
    
    // recorre las particulas
    for (j=0; j<total_p; j++) {
      if (bound && (po[j].pe < po[j].ke)) continue;
      r=0;
      for (k=0; k<3; k++) 
      {        
        for (dr=0, l=0; l<3; l++) 
        {
            dr += orth[k][l]*(po[j].pos[l]-h->pos[l]);
        }
        r += dr*dr/eig[k]; // (x-x0)**2/rhalo**2
      }
      
      if (r < min_r) r = min_r;
      
      if (!(r>0 && r<=1)) continue;
      
      double tw = (WEIGHTED_SHAPES) ? 1.0/r : 1.0; // chequea standard o pesado
      weight +=tw;
      
      // CONSTRUYE EL TENSOR      
      for (k=0; k<3; k++) {
        dr = po[j].pos[k]-h->pos[k];
        mass_t[k][k] += dr*dr*tw;
    
        for (l=0; l<k; l++) {
            dr2 = po[j].pos[l]-h->pos[l];
            mass_t[k][l] += dr2*dr*tw;
            mass_t[l][k] = mass_t[k][l];
        }
      }
    }

    if (!weight) return; // chequeo (dismiss)
    
    for (k=0; k<3; k++) for (l=0; l<3; l++) mass_t[k][l] /= (double)weight;
    
    // CALCULA AUTOVALORES Y AUTOVECTORES
    
    jacobi_decompose(mass_t, eig, orth);
    
    // ordena los autovalores
    a = 0; b = 1; c = 2;
    if (eig[1]>eig[0]) { b=0; a=1; }
    if (eig[2]>eig[b]) { c=b; b=2; }
    if (eig[b]>eig[a]) { int64_t t=a; a=b; b=t; }
    
    if (!eig[a] || !eig[b] || !eig[c]) return; // chequeo
    
    b_to_a = sqrt(eig[b]/eig[a]);
    c_to_a = sqrt(eig[c]/eig[a]);
    
    // TERMINA LA ITERACION SI a/b = 0.01 a/b_prev y los mismo para c/a
    
    if ((fabs(b_to_a-h->b_to_a) < 0.01*h->b_to_a) &&
	(fabs(c_to_a-h->c_to_a) < 0.01*h->c_to_a)) return;
    
    // Pisa los viejos
    
    h->b_to_a = (b_to_a > 0) ? b_to_a : 0;
    h->c_to_a = (c_to_a > 0) ? c_to_a : 0;
    
    // Redefine el radio del halo = sqrt(a)
    
    r = sqrt(eig[a]);
    
    // Guarda un autovector
    for (k=0; k<3; k++) {
      h->A[k] = 1e3*r*orth[a][k];
      eig[k] *= (h->r*h->r*1e-6)/(r*r);
    }
  }
}
*/
