//------------------------------------------------------------------------------
#include <stdio.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
// gcc gsl.c -I/usr/include/gsl -lgsl -lcblas
//------------------------------------------------------------------------------
int main()
{
  //void cst();  cst();
  void eigen();  eigen();
}
//------------------------------------------------------------------------------
void cst(){
  double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
  double au = GSL_CONST_MKSA_ASTRONOMICAL_UNIT;
  double minutes = GSL_CONST_MKSA_MINUTE;
  double r_earth = 1.00 * au;
  double r_mars = 1.52 * au;
  double t_min, t_max;
  t_min = (r_mars - r_earth) / c;
  t_max = (r_mars + r_earth) / c;
  printf ("light travel time from Earth to Mars:\n");
  printf ("minimum = %.1f minutes\n", t_min / minutes);
  printf ("maximum = %.1f minutes\n", t_max / minutes);
}
//------------------------------------------------------------------------------
void eigen ()
{
  double data[] = { 1.0  , 1/2.0, 1/3.0, 1/4.0,
                    1/2.0, 1/3.0, 1/4.0, 1/5.0,
                    1/3.0, 1/4.0, 1/5.0, 1/6.0,
                    1/4.0, 1/5.0, 1/6.0, 1/7.0 };

  gsl_matrix_view m  = gsl_matrix_view_array (data, 4, 4);
  gsl_vector *eval = gsl_vector_alloc (4);
  gsl_matrix *evec = gsl_matrix_alloc (4, 4);

  gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc (4);

  gsl_eigen_symmv (&m.matrix, eval, evec, w);
  gsl_eigen_symmv_free (w);
  gsl_eigen_symmv_sort (eval, evec,
                        GSL_EIGEN_SORT_ABS_ASC);
  {
    int i;

    for (i = 0; i < 4; i++)
      {
        double eval_i
           = gsl_vector_get (eval, i);
        gsl_vector_view evec_i
           = gsl_matrix_column (evec, i);

        printf ("eigenvalue = %g\n", eval_i);
        printf ("eigenvector = \n");
        gsl_vector_fprintf (stdout,
                            &evec_i.vector, "%g");
      }
  }
  gsl_vector_free (eval);
  gsl_matrix_free (evec);
}
//------------------------------------------------------------------------------
