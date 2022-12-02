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

    for(int k = 0; k < 10; k++)
      fprintf(stdout,"inside %f %f\n",x_part_proj[k],y_part_proj[k]);

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

    *a2D_abs = 0;
    *b2D_abs = 0;

    //get eigenvalues from gsl
    *a2D_abs = sqrt(fabs(gsl_vector_get(eval2D, (0))));
    *b2D_abs = sqrt(fabs(gsl_vector_get(eval2D, (1))));
    
    fprintf(stdout,"%f %f\n", *a2D_abs, *b2D_abs);

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
    
    fprintf(stdout,"%f %f\n", *a2Dr_abs, *b2Dr_abs);

    
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

#ifndef ROCKSTAR_SHAPE

void calculate_2d_shapes_iterative(const vector <float> x_part_proj, const vector <float> y_part_proj, \
                        const double a_t, float *a2D, float *b2D, float *a2Dr, float *b2Dr, \
                        float *a2Dr_abs, float *b2Dr_abs, float *a2D_abs, float *b2D_abs, const float Halo_R)
{
  int i,j,k,a,b;
  int WEIGHTED_SHAPES; 
  int Npart = x_part_proj.size();
  int iter = Npart < SHAPE_ITERATION ? Npart : SHAPE_ITERATION;
  vector<float> tmp_x_part_proj, tmp_y_part_proj, tmp_rell_proj;
  float b_to_a, prev_b_to_a;
  float rscale;
  double eig[2], orth[2][2];

  for(WEIGHTED_SHAPES=0; WEIGHTED_SHAPES<=1; WEIGHTED_SHAPES++)
  {
    /// Inicializa
    rscale      = Halo_R;
    b_to_a      = 0;
    prev_b_to_a = 0;
    for (i=0; i<2; i++) 
    {
      for(j=0; j<2; j++) 
        orth[i][j] = 0;
      orth[i][i] = 1;
      eig[i] = 1;
    }

    // EMPIEZA A ITERAR  
    for (i=0; i<iter; i++) 
    {    
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

      for(j=0;j<Npart;j++)
      {
        float xell = (orth[0][0] * a_t * x_part_proj[j]) + (orth[0][1] * a_t * y_part_proj[j]);
        float yell = (orth[1][0] * a_t * x_part_proj[j]) + (orth[1][1] * a_t * y_part_proj[j]);
        float rell = sqrt(pow(xell,2)/eig[0] + pow(yell,2)/eig[1]);

        if(rell/rscale <= 1)
        {
          tmp_x_part_proj.push_back(x_part_proj[j]);
          tmp_y_part_proj.push_back(y_part_proj[j]);
          if(WEIGHTED_SHAPES) //----- reduced MI -----
            tmp_rell_proj.push_back(rell);
          else                //----- standard MI -----
            tmp_rell_proj.push_back(1.0);
        }
      }

      if(tmp_rell_proj.size() == 0)
        break;
 
      if(WEIGHTED_SHAPES)  //----- reduced MI -----
      {
        //initialize moment of inertia 2DMi
        ini_MI_2D_iterative(tmp_x_part_proj, tmp_y_part_proj, tmp_rell_proj, a_t, MI_2D);

        //fill gsl matrix
        M2D = gsl_matrix_view_array (MI_2D, 2, 2);

        //get eigenvalues and eigenvectors  of MI
        gsl_eigen_symmv (&M2D.matrix, eval2D, evec2D, w2D);

        //sort eigenvalues and eigenvectors of MI in descending order
        gsl_eigen_symmv_sort (eval2D, evec2D, GSL_EIGEN_SORT_ABS_DESC);

        if(fabs(gsl_vector_get(eval2D, (0))) == 0 ||
           fabs(gsl_vector_get(eval2D, (1))) == 0)
          break;

        //get normalized eigenvectors from gsl
        a2Dr[0] = gsl_matrix_get(evec2D, (0), (0));
        a2Dr[1] = gsl_matrix_get(evec2D, (1), (0));

        b2Dr[0] = gsl_matrix_get(evec2D, (0), (1));
        b2Dr[1] = gsl_matrix_get(evec2D, (1), (1));

        //get eigenvalues from gsl
        *a2Dr_abs = sqrt(fabs(gsl_vector_get(eval2D, (0))));
        *b2Dr_abs = sqrt(fabs(gsl_vector_get(eval2D, (1))));
        //------------------------------------------------------------------------

        b_to_a = (*b2Dr_abs)/(*a2Dr_abs);

      }else{ //----- standard MI -----

        //initialize moment of inertia
        ini_MI_2D_iterative(tmp_x_part_proj, tmp_y_part_proj, tmp_rell_proj, a_t, MI_2D);

        //fill gsl matrix
        M2D = gsl_matrix_view_array (MI_2D, 2, 2);

        //get eigenvalues and eigenvectors  of MI
        gsl_eigen_symmv (&M2D.matrix, eval2D, evec2D, w2D);

        //sort eigenvalues and eigenvectors of MI in descending order
        gsl_eigen_symmv_sort (eval2D, evec2D, GSL_EIGEN_SORT_ABS_DESC);

        if(fabs(gsl_vector_get(eval2D, (0))) == 0 ||
           fabs(gsl_vector_get(eval2D, (1))) == 0)
          break;

        //get normalized eigenvectors from gsl
        a2D[0] = gsl_matrix_get(evec2D, (0), (0));
        a2D[1] = gsl_matrix_get(evec2D, (1), (0));

        b2D[0] = gsl_matrix_get(evec2D, (0), (1));
        b2D[1] = gsl_matrix_get(evec2D, (1), (1));

        //get eigenvalues from gsl
        *a2D_abs = sqrt(fabs(gsl_vector_get(eval2D, (0))));
        *b2D_abs = sqrt(fabs(gsl_vector_get(eval2D, (1))));

        b_to_a = (*b2D_abs)/(*a2D_abs);

      }

      tmp_x_part_proj.clear();
      tmp_y_part_proj.clear();
      tmp_rell_proj.clear();

      if((fabs(b_to_a-prev_b_to_a) < 0.01*prev_b_to_a)) break;

      prev_b_to_a = (b_to_a > 0) ? b_to_a : 0;

      eig[0] = 1; // Estan al cuadrado
      eig[1] = fabs(gsl_vector_get(eval2D, (1))) / fabs(gsl_vector_get(eval2D, (0)));

      // Guarda un autovector
      for (j=0; j<2; j++) 
      {
        if(WEIGHTED_SHAPES)
        {
          orth[0][j] = a2Dr[j];
          orth[1][j] = b2Dr[j];
        }else{
          orth[0][j] =  a2D[j];
          orth[1][j] =  b2D[j];
        }
      }

    } // Iteracion
  } // 
}

void calculate_3d_shapes_iterative(\
                const vector <float> x_part, const vector <float> y_part, const vector <float> z_part, \
                const double a_t, float *a3D, float *b3D, float *c3D, \
                float *a3Dr, float *b3Dr, float *c3Dr, \
                float *a3D_abs, float *b3D_abs, float *c3D_abs, \
                float *a3Dr_abs, float *b3Dr_abs, float *c3Dr_abs, const float Halo_R)
{
  int i, j;
  int WEIGHTED_SHAPES; 
  int Npart = x_part.size();
  int iter = Npart < SHAPE_ITERATION ? Npart : SHAPE_ITERATION;
  vector<float> tmp_x_part, tmp_y_part, tmp_z_part, tmp_rell;
  float b_to_a, prev_b_to_a;
  float c_to_a, prev_c_to_a;
  float rscale;
  double eig[3], orth[3][3];

  for(WEIGHTED_SHAPES=0; WEIGHTED_SHAPES<=1; WEIGHTED_SHAPES++)
  {
    /// Inicializa
    rscale      = Halo_R;
    b_to_a      = c_to_a = 0;
    prev_b_to_a = prev_c_to_a = 0;
    for (i=0; i<3; i++) 
    {
      for(j=0; j<3; j++) 
        orth[i][j] = 0;
      orth[i][i] = 1;
      eig[i] = 1;
    }

    // EMPIEZA A ITERAR  
    for (i=0; i<iter; i++) 
    {    
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

      //------------------------ Select Particles ------------------------------

      for(j=0;j<Npart;j++)
      {
        float xell = (orth[0][0] * a_t * x_part[j]) + (orth[0][1] * a_t * y_part[j]) + (orth[0][2] * a_t * z_part[j]);
        float yell = (orth[1][0] * a_t * x_part[j]) + (orth[1][1] * a_t * y_part[j]) + (orth[1][2] * a_t * z_part[j]);
        float zell = (orth[2][0] * a_t * x_part[j]) + (orth[2][1] * a_t * y_part[j]) + (orth[2][2] * a_t * z_part[j]);
        
        float rell = sqrt(pow(xell,2)/eig[0] + pow(yell,2)/eig[1] + pow(zell,2)/eig[2]);

        if(rell/rscale <= 1)
        {
          tmp_x_part.push_back(x_part[j]);
          tmp_y_part.push_back(y_part[j]);
          tmp_z_part.push_back(z_part[j]);
          if(WEIGHTED_SHAPES) //----- reduced MI -----
            tmp_rell.push_back(rell);
          else                //----- standard MI -----
            tmp_rell.push_back(1.0);
        }
      }

      if(tmp_rell.size() == 0)
        break;
 
      if(WEIGHTED_SHAPES)  //----- reduced MI -----
      {
        //initialize moment of inertia
        ini_MI_3D_iterative(tmp_x_part, tmp_y_part, tmp_z_part, tmp_rell, a_t, MI_3D);

        //get eigenvalues and eigenvectors of MI using gsl
        M3D = gsl_matrix_view_array (MI_3D, 3, 3);
        gsl_eigen_symmv (&M3D.matrix, eval3D, evec3D, w3D);

        //sort eigenvalues and eigenvectors of MI in descending order
        gsl_eigen_symmv_sort (eval3D, evec3D, GSL_EIGEN_SORT_ABS_DESC);

        if(fabs(gsl_vector_get(eval3D, (0))) == 0 ||
           fabs(gsl_vector_get(eval3D, (1))) == 0 ||
           fabs(gsl_vector_get(eval3D, (2))) == 0 )
          break;

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

        b_to_a = (*b3Dr_abs)/(*a3Dr_abs);
        c_to_a = (*c3Dr_abs)/(*a3Dr_abs);

      }else{ //----- standard MI -----

        //initialize moment of inertia
        ini_MI_3D_iterative(tmp_x_part, tmp_y_part, tmp_z_part, tmp_rell, a_t, MI_3D);

        //get eigenvalues and eigenvectors of MI using gsl
        M3D = gsl_matrix_view_array (MI_3D, 3, 3);
        gsl_eigen_symmv (&M3D.matrix, eval3D, evec3D, w3D);

        //sort eigenvalues and eigenvectors of MI in descending order
        gsl_eigen_symmv_sort (eval3D, evec3D, GSL_EIGEN_SORT_ABS_DESC);

        if(fabs(gsl_vector_get(eval3D, (0))) == 0 ||
           fabs(gsl_vector_get(eval3D, (1))) == 0 ||
           fabs(gsl_vector_get(eval3D, (2))) == 0 )
          break;

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

        b_to_a = (*b3D_abs)/(*a3D_abs);
        c_to_a = (*c3D_abs)/(*a3D_abs);
      }

      tmp_x_part.clear();
      tmp_y_part.clear();
      tmp_z_part.clear();
      tmp_rell.clear();

      if((fabs(b_to_a-prev_b_to_a) < 0.01*prev_b_to_a) && (fabs(c_to_a-prev_c_to_a) < 0.01*prev_c_to_a)) break;

      prev_b_to_a = (b_to_a > 0) ? b_to_a : 0;
      prev_c_to_a = (c_to_a > 0) ? c_to_a : 0;

      eig[0] = 1; // Estan al cuadrado
      eig[1] = fabs(gsl_vector_get(eval3D, (1))) / fabs(gsl_vector_get(eval3D, (0)));
      eig[2] = fabs(gsl_vector_get(eval3D, (2))) / fabs(gsl_vector_get(eval3D, (0)));

      // Guarda un autovector
      for (j=0; j<3; j++) 
      {
        if(WEIGHTED_SHAPES)
        {
          orth[0][j] = a3Dr[j];
          orth[1][j] = b3Dr[j];
          orth[2][j] = c3Dr[j];
        }else{
          orth[0][j] =  a3D[j];
          orth[1][j] =  b3D[j];
          orth[2][j] =  c3D[j];
        }
      }

    } // Iteracion
  } // 

}

#else

static void jacobi_decompose_2D(double cov_matrix[][NUM_PARAMS_2D], double *eigenvalues, double orth_matrix[][NUM_PARAMS_2D]) 
{
  const int NUM_PARAMS = NUM_PARAMS_2D;
  int i,j,k,l;
  int max_col[NUM_PARAMS];
  int changed[NUM_PARAMS]; //To keep track of eigenvalues changing
  int n = NUM_PARAMS;
  int max_row;
  int state = 0;
  double max, c, s, t, u, a;
  int count = 0;

  //Set up the maximum value cache
  for (i=0; i<n; i++) {
    max=0;
    max_col[i] = i+1;
    for (j=i+1; j<n; j++) {
      if (fabs(cov_matrix[i][j])>max) {
        max_col[i] = j;
      	max = fabs(cov_matrix[i][j]);
      }
    }
  }

  //Set up the orthogonal matrix as the identity:
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      orth_matrix[i][j] = (i==j) ? 1 : 0;

  //Copy the diagonal values to the eigenvalue array:
  for (i=0; i<n; i++) {
    eigenvalues[i] = cov_matrix[i][i];
    changed[i] = 1;
    for (j=0; j<n; j++)
      if (j!=i && cov_matrix[i][j]) break;
    if (j==n) {
      state--;
      changed[i] = 0;
    }
  }

  //Sweep time: iterate until the eigenvalues stop changing.
  state += n;
  while (state) {
    count++;
    //Find the largest nonzero element in the matrix:
    max = fabs(cov_matrix[0][max_col[0]]);
    max_row = 0;
    for (i=1; i<n-1; i++) {
      if (fabs(cov_matrix[i][max_col[i]])>max) {
	max = fabs(cov_matrix[i][max_col[i]]);
	max_row = i;
      }
    }
    k = max_row; l = max_col[k];
    max = cov_matrix[k][l];
    if (max==0) break;

    //Calculate the Jacobi rotation matrix
    //Tan 2phi = 2S_kl / (S_kk - S_ll)
    a = (eigenvalues[l] - eigenvalues[k]) * 0.5;
    t = fabsl(a) + sqrtl(max*max + a*a);
    s = sqrtl(max*max + t*t);
    c = t/s;
    s = max/s;
    t = max*max / t;
    if (a<0) {
      s = -s; t = -t;
    }

    //Update eigenvalues
#define UPDATE(x,y)							\
    a = eigenvalues[x];							\
    eigenvalues[x] += y;						\
    if (changed[x] && (fabs(y) < fabs(eigenvalues[x]))) { /*Eigenvalue didn't change*/ \
      changed[x]=0;							\
      state--;								\
    } else if (!changed[x] && (fabs(y) > fabs(eigenvalues[x]))) { /*Egval did change*/ \
      changed[x] = 1;							\
      state++;								\
    }
      
    UPDATE(k, -t);
    UPDATE(l, t);

    //Update covariance matrix:
    cov_matrix[k][l] = 0;

#define ROTATE(m,w,x,y,z,r)  /*Perform a Jacobi rotation*/	\
      t = m[w][x]; u = m[y][z];					\
      m[w][x] = t*c - s*u;					\
      m[y][z] = s*t + c*u;				        \
      if (r) {							\
	if (fabs(m[w][x])>fabs(m[w][max_col[w]]))		\
	  max_col[w] = x;					\
	if (y < NUM_PARAMS && fabs(m[y][z])>fabs(m[y][max_col[y]])) \
	  max_col[y] = z;					\
      }
    

    for (i=0; i<k; i++) { ROTATE(cov_matrix,i,k,i,l,1); }
    for (i=k+1; i<l; i++) { ROTATE(cov_matrix,k,i,i,l,1); }
    for (i=l+1; i<n; i++) { ROTATE(cov_matrix,k,i,l,i,1); }
    for (i=0; i<n; i++) { ROTATE(orth_matrix,k,i,l,i,0); }
    
#undef ROTATE
#undef UPDATE
  }

  // Restore the matrix to its original form:
  for (k=0; k<n-1; k++) {
    for (l=k+1; l<n; l++) {
      cov_matrix[k][l] = cov_matrix[l][k];
    }
  }
}

static void jacobi_decompose_3D(double cov_matrix[][NUM_PARAMS_3D], double *eigenvalues, double orth_matrix[][NUM_PARAMS_3D]) 
{
  const int NUM_PARAMS = NUM_PARAMS_3D;
  int i,j,k,l;
  int max_col[NUM_PARAMS];
  int changed[NUM_PARAMS]; //To keep track of eigenvalues changing
  int n = NUM_PARAMS;
  int max_row;
  int state = 0;
  double max, c, s, t, u, a;
  int count = 0;

  //Set up the maximum value cache
  for (i=0; i<n; i++) {
    max=0;
    max_col[i] = i+1;
    for (j=i+1; j<n; j++) {
      if (fabs(cov_matrix[i][j])>max) {
        max_col[i] = j;
      	max = fabs(cov_matrix[i][j]);
      }
    }
  }

  //Set up the orthogonal matrix as the identity:
  for (i=0; i<n; i++)
    for (j=0; j<n; j++)
      orth_matrix[i][j] = (i==j) ? 1 : 0;

  //Copy the diagonal values to the eigenvalue array:
  for (i=0; i<n; i++) {
    eigenvalues[i] = cov_matrix[i][i];
    changed[i] = 1;
    for (j=0; j<n; j++)
      if (j!=i && cov_matrix[i][j]) break;
    if (j==n) {
      state--;
      changed[i] = 0;
    }
  }

  //Sweep time: iterate until the eigenvalues stop changing.
  state += n;
  while (state) {
    count++;
    //Find the largest nonzero element in the matrix:
    max = fabs(cov_matrix[0][max_col[0]]);
    max_row = 0;
    for (i=1; i<n-1; i++) {
      if (fabs(cov_matrix[i][max_col[i]])>max) {
	max = fabs(cov_matrix[i][max_col[i]]);
	max_row = i;
      }
    }
    k = max_row; l = max_col[k];
    max = cov_matrix[k][l];
    if (max==0) break;

    //Calculate the Jacobi rotation matrix
    //Tan 2phi = 2S_kl / (S_kk - S_ll)
    a = (eigenvalues[l] - eigenvalues[k]) * 0.5;
    t = fabsl(a) + sqrtl(max*max + a*a);
    s = sqrtl(max*max + t*t);
    c = t/s;
    s = max/s;
    t = max*max / t;
    if (a<0) {
      s = -s; t = -t;
    }

    //Update eigenvalues
#define UPDATE(x,y)							\
    a = eigenvalues[x];							\
    eigenvalues[x] += y;						\
    if (changed[x] && (fabs(y) < fabs(eigenvalues[x]))) { /*Eigenvalue didn't change*/ \
      changed[x]=0;							\
      state--;								\
    } else if (!changed[x] && (fabs(y) > fabs(eigenvalues[x]))) { /*Egval did change*/ \
      changed[x] = 1;							\
      state++;								\
    }
      
    UPDATE(k, -t);
    UPDATE(l, t);

    //Update covariance matrix:
    cov_matrix[k][l] = 0;

#define ROTATE(m,w,x,y,z,r)  /*Perform a Jacobi rotation*/	\
      t = m[w][x]; u = m[y][z];					\
      m[w][x] = t*c - s*u;					\
      m[y][z] = s*t + c*u;				        \
      if (r) {							\
	if (fabs(m[w][x])>fabs(m[w][max_col[w]]))		\
	  max_col[w] = x;					\
	if (y < NUM_PARAMS && fabs(m[y][z])>fabs(m[y][max_col[y]])) \
	  max_col[y] = z;					\
      }
    

    for (i=0; i<k; i++) { ROTATE(cov_matrix,i,k,i,l,1); }
    for (i=k+1; i<l; i++) { ROTATE(cov_matrix,k,i,i,l,1); }
    for (i=l+1; i<n; i++) { ROTATE(cov_matrix,k,i,l,i,1); }
    for (i=0; i<n; i++) { ROTATE(orth_matrix,k,i,l,i,0); }
    
#undef ROTATE
#undef UPDATE
  }

  // Restore the matrix to its original form:
  for (k=0; k<n-1; k++) {
    for (l=k+1; l<n; l++) {
      cov_matrix[k][l] = cov_matrix[l][k];
    }
  }
}

void calculate_2d_shapes_iterative(const vector <float> x_part_proj, const vector <float> y_part_proj, \
                        const double a_t, float *a2D, float *b2D, float *a2Dr, float *b2Dr, \
                        float *a2Dr_abs, float *b2Dr_abs, float *a2D_abs, float *b2D_abs, const float Halo_R)
{
  int Npart = x_part_proj.size();
  int iter = Npart < SHAPE_ITERATION ? Npart : SHAPE_ITERATION;
  int WEIGHTED_SHAPES; 
  int i,j,k,a,b;
  float b_to_a, prev_b_to_a;
  double r, weight, dr;
  double eig[2]={0};
  double mass_t[2][2], orth[2][2];

  for(WEIGHTED_SHAPES=0; WEIGHTED_SHAPES<=1; WEIGHTED_SHAPES++)
  {
    //INICIALIZA
    b_to_a      = 0;
    prev_b_to_a = 0;
    for (i=0; i<2; i++) 
    {
      for(j=0; j<2; j++) 
        orth[i][j] = 0;
      orth[i][i] = 1;
      eig[i] = (Halo_R*Halo_R);
    }

    // EMPIEZA A ITERAR  
    for (i=0; i<iter; i++) 
    {    
      // inicializa mass_t (momento de inercia) 
      weight=0;
      for (k=0; k<2; k++) 
        for(j=0; j<2; j++) 
          mass_t[k][j] = 0;
  
      // recorre las particulas
      for (j=0; j<Npart; j++) 
      {
        r = 0;
        for (k=0; k<2; k++) 
        {
          dr = (orth[k][0] * a_t * x_part_proj[j]) + (orth[k][1] * a_t * y_part_proj[j]);
          r += dr*dr/eig[k]; // (x-x0)**2/rhalo**2
        }      
        
        if (!(r>0 && r<=1)) continue;
        
        double tw = (WEIGHTED_SHAPES) ? 1.0/r : 1.0; // chequea standard o pesado

        mass_t[0][0] += pow(a_t,2) * x_part_proj[j] * x_part_proj[j] * tw;
        mass_t[0][1] += pow(a_t,2) * x_part_proj[j] * y_part_proj[j] * tw;

        mass_t[1][0] += pow(a_t,2) * y_part_proj[j] * x_part_proj[j] * tw;
        mass_t[1][1] += pow(a_t,2) * y_part_proj[j] * y_part_proj[j] * tw;
        
        weight +=tw;
      } // Cierro el numero de particulas

      for(k=0; k<2; k++) 
        for(j=0; j<2; j++) 
          mass_t[k][j] /= (double)weight;
      
      // CALCULA AUTOVALORES Y AUTOVECTORES
      jacobi_decompose_2D(mass_t, eig, orth);
      
      // ordena los autovalores
      a = 0; b = 1;
      if (eig[b]>eig[a]) { b=0; a=1; }
      
      // Chequeo si algun autovalor es nulo
      if(!eig[a] || !eig[b]) break;
      
      b_to_a = sqrt(eig[b]/eig[a]);
      
      // TERMINA LA ITERACION SI a/b = 0.01 a/b_prev y los mismo para c/a
      if ((fabs(b_to_a-prev_b_to_a) < 0.01*prev_b_to_a)) break;
      
      prev_b_to_a = (b_to_a > 0) ? b_to_a : 0;
      
      r = sqrt(eig[a]);
      
      // Guarda un autovector
      for (k=0; k<2; k++) 
      {
        if(WEIGHTED_SHAPES)
        {
          a2Dr[k] = orth[a][k];
          b2Dr[k] = orth[b][k];
        }else{
          a2D[k]  = orth[a][k];
          b2D[k]  = orth[b][k];
        }
        eig[k] *= (Halo_R*Halo_R)/(r*r);
      }

      if(WEIGHTED_SHAPES)
      {
        *a2D_abs  = sqrt(eig[a]);
        *b2D_abs  = sqrt(eig[b]);
      }else{
        *a2Dr_abs = sqrt(eig[a]);
        *b2Dr_abs = sqrt(eig[b]);
      }
    } // Cierra la iteracion
  }
}

void calculate_3d_shapes_iterative(\
                const vector <float> x_part, const vector <float> y_part, const vector <float> z_part, \
                const double a_t, float *a3D, float *b3D, float *c3D, \
                float *a3Dr, float *b3Dr, float *c3Dr, \
                float *a3D_abs, float *b3D_abs, float *c3D_abs, \
                float *a3Dr_abs, float *b3Dr_abs, float *c3Dr_abs, const float Halo_R)
{
  int Npart = x_part.size();
  int iter = Npart < SHAPE_ITERATION ? Npart : SHAPE_ITERATION;
  int WEIGHTED_SHAPES; 
  int i,j,k,a,b,c;
  float b_to_a, c_to_a, prev_b_to_a, prev_c_to_a;
  double r, weight, dr;
  double eig[3]={0};
  double mass_t[3][3], orth[3][3];

  for(WEIGHTED_SHAPES=0; WEIGHTED_SHAPES<=1; WEIGHTED_SHAPES++)
  {
    //INICIALIZA
    b_to_a      = c_to_a = 0;
    prev_b_to_a = prev_c_to_a = 0;
    for (i=0; i<3; i++) 
    {
      for(j=0; j<3; j++) 
        orth[i][j] = 0;
      orth[i][i] = 1;
      eig[i] = (Halo_R*Halo_R);
    }

    // EMPIEZA A ITERAR  
    for (i=0; i<iter; i++) 
    {    
      // inicializa mass_t (momento de inercia) 
      weight=0;
      for (k=0; k<3; k++) 
        for(j=0; j<3; j++) 
          mass_t[k][j] = 0;
  
      // recorre las particulas
      for (j=0; j<Npart; j++) 
      {
        r = 0;
        for (k=0; k<3; k++) 
        {
          dr = (orth[k][0] * a_t * x_part[j]) + (orth[k][1] * a_t * y_part[j]) + (orth[k][2] * a_t * z_part[j]);
          r += dr*dr/eig[k]; // (x-x0)**2/rhalo**2
        }      
        
        if (!(r>0 && r<=1)) continue;
        
        double tw = (WEIGHTED_SHAPES) ? 1.0/r : 1.0; // chequea standard o pesado

        mass_t[0][0] += pow(a_t,2) * x_part[j] * x_part[j] * tw;
        mass_t[0][1] += pow(a_t,2) * x_part[j] * y_part[j] * tw;
        mass_t[0][2] += pow(a_t,2) * x_part[j] * z_part[j] * tw;

        mass_t[1][0] += pow(a_t,2) * y_part[j] * x_part[j] * tw;
        mass_t[1][1] += pow(a_t,2) * y_part[j] * y_part[j] * tw;
        mass_t[1][2] += pow(a_t,2) * y_part[j] * z_part[j] * tw;

        mass_t[2][0] += pow(a_t,2) * z_part[j] * x_part[j] * tw;
        mass_t[2][1] += pow(a_t,2) * z_part[j] * y_part[j] * tw;
        mass_t[2][2] += pow(a_t,2) * z_part[j] * z_part[j] * tw;
        
        weight +=tw;
      } // Cierro el numero de particulas

      for(k=0; k<3; k++) 
        for(j=0; j<3; j++) 
          mass_t[k][j] /= (double)weight;
      
      // CALCULA AUTOVALORES Y AUTOVECTORES
      jacobi_decompose_3D(mass_t, eig, orth);
     
      // ordena los autovalores
      a = 0; b = 1; c = 2;
      if (eig[1]>eig[0]) { b=0; a=1; }
      if (eig[2]>eig[b]) { c=b; b=2; }
      if (eig[b]>eig[a]) { int64_t t=a; a=b; b=t; }
      
      // Chequeo si algun autovalor es nulo
      if(!eig[a] || !eig[b] || !eig[c]) break;
      
      b_to_a = sqrt(eig[b]/eig[a]);
      c_to_a = sqrt(eig[c]/eig[a]);
      
      // TERMINA LA ITERACION SI a/b = 0.01 a/b_prev y los mismo para c/a
      if ((fabs(b_to_a-prev_b_to_a) < 0.01*prev_b_to_a) && (fabs(c_to_a-prev_c_to_a) < 0.01*prev_c_to_a)) break;
      
      prev_b_to_a = (b_to_a > 0) ? b_to_a : 0;
      prev_c_to_a = (c_to_a > 0) ? c_to_a : 0;
      
      r = sqrt(eig[a]);
      
      // Guarda un autovector
      for (k=0; k<3; k++) 
      {
        if(WEIGHTED_SHAPES)
        {
          a3Dr[k] = orth[a][k];
          b3Dr[k] = orth[b][k];
          c3Dr[k] = orth[c][k];
        }else{
          a3D[k]  = orth[a][k];
          b3D[k]  = orth[b][k];
          c3D[k]  = orth[c][k];
        }
        eig[k] *= (Halo_R*Halo_R)/(r*r);
      }

      if(WEIGHTED_SHAPES)
      {
        *a3D_abs  = sqrt(eig[a]);
        *b3D_abs  = sqrt(eig[b]);
        *c3D_abs  = sqrt(eig[c]);
      }else{
        *a3Dr_abs = sqrt(eig[a]);
        *b3Dr_abs = sqrt(eig[b]);
        *c3Dr_abs = sqrt(eig[c]);
      }
    } // Cierra la iteracion
  }
}

#endif
