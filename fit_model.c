#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

#include "allvars.h"
#include "define.h"
#include "proto.h"


#ifdef GENERAL_NFW
int dens_f(const gsl_vector *x, void *data_now, gsl_vector *f)
{
  int n = ((data *) data_now)->n;
  double *y = ((data *) data_now)->y;
  double *sigma = ((data *) data_now)->sigma;
  
  double log_rhos = gsl_vector_get (x, 0);
  double log_rs = gsl_vector_get (x, 1);
  double alpha = gsl_vector_get (x, 2);
  
  int i;
  
  for(i=0; i<n; i++)
    {
	double log_r = profile[i][0];  /* log10(r) */
	double log_rhoi = log_rhos - alpha * (log_r - log_rs) - (3 - alpha) * (log10(pow(10, log_r) + pow(10, log_rs)) - log_rs);
	
	gsl_vector_set(f, i, (log_rhoi - y[i])/sigma[i]);
		
	}	

  return GSL_SUCCESS;
}        /* end dens_f */

int dens_df (const gsl_vector *x, void *data_now, gsl_matrix *J)
{
  int n = ((data *) data_now)->n;
  double *sigma = ((data *) data_now)->sigma;
  
  double log_rhos = gsl_vector_get (x, 0);
  double log_rs = gsl_vector_get(x, 1);
  double alpha = gsl_vector_get (x, 2);
  
  int i;
  
  for(i = 0; i<n ; i++)
    {
	double log_r = profile[i][0];
	double s = sigma[i];
	
	gsl_matrix_set(J, i, 0, 1.0 / s);
	gsl_matrix_set(J, i, 1,  (3-(3-alpha)*pow(10, log_rs)/(pow(10, log_r)+pow(10, log_rs))) / s);
	gsl_matrix_set(J, i, 2,  (log10(pow(10, log_r) + pow(10, log_rs)) - log_r) / s); 					
	}
	
   return GSL_SUCCESS;	
}       /* end dens_df */	


int dens_fdf(const gsl_vector *x, void *data_now, gsl_vector *f, gsl_matrix *J)
{
  dens_f(x, data_now, f);
  dens_df(x, data_now, J);
  
  return GSL_SUCCESS;
	
}      /* end dens_fdf */

#endif

#ifdef BURKERT

int dens_f(const gsl_vector *x, void *data_now, gsl_vector *f)
{
  int n = ((data *) data_now)->n;
  double *y = ((data *) data_now)->y;
  double *sigma = ((data *) data_now)->sigma;
  
  double log_rhos = gsl_vector_get (x, 0);
  double log_rs = gsl_vector_get (x, 1);
  
  int i;
  
  for(i=0; i<n; i++)
    {
	double log_r = profile[i][0];  /* log10(r) */
	double log_rhoi = log_rhos + 3 * log_rs - log10(pow(10, log_r)+pow(10, log_rs)) - log10(pow(10, 2*log_r)+pow(10, 2*log_rs));
	
	gsl_vector_set(f, i, (log_rhoi - y[i])/sigma[i]);
		
	}	

  return GSL_SUCCESS;
}        /* end dens_f */

int dens_df (const gsl_vector *x, void *data_now, gsl_matrix *J)
{
  int n = ((data *) data_now)->n;
  double *sigma = ((data *) data_now)->sigma;
  
  double log_rhos = gsl_vector_get (x, 0);
  double log_rs = gsl_vector_get(x, 1);
  
  int i;
  
  for(i = 0; i<n ; i++)
    {
	double log_r = profile[i][0];
	double s = sigma[i];
	
	gsl_matrix_set(J, i, 0, 1.0 / s);
	gsl_matrix_set(J, i, 1,  (3 - 1/(pow(10, log_r - log_rs)+1) - 2/(pow(10, 2*(log_r - log_rs)) + 1)) / s);					
	}
	
   return GSL_SUCCESS;	
}       /* end dens_df */	


int dens_fdf(const gsl_vector *x, void *data_now, gsl_vector *f, gsl_matrix *J)
{
  dens_f(x, data_now, f);
  dens_df(x, data_now, J);
  
  return GSL_SUCCESS;
	
}      /* end dens_fdf */


#endif


#ifdef BURKERT_NFW

int dens_f(const gsl_vector *x, void *data_now, gsl_vector *f)
{
  int n = ((data *) data_now)->n;
  double *y = ((data *) data_now)->y;
  double *sigma = ((data *) data_now)->sigma;
  
  double log_rhos = gsl_vector_get (x, 0);
  double log_rs = gsl_vector_get (x, 1);
  
  int i;
  
  for(i=0; i<n; i++)
    {
	double log_r = profile[i][0];  /* log10(r) */
	double log_rhoi = log_rhos - 3 * log10(pow(10, log_r)+pow(10, log_rs)) + 3 * log_rs;
	
	gsl_vector_set(f, i, (log_rhoi - y[i])/sigma[i]);
		
	}	

  return GSL_SUCCESS;
}        /* end dens_f */

int dens_df (const gsl_vector *x, void *data_now, gsl_matrix *J)
{
  int n = ((data *) data_now)->n;
  double *sigma = ((data *) data_now)->sigma;
  
  double log_rhos = gsl_vector_get (x, 0);
  double log_rs = gsl_vector_get(x, 1);
  
  int i;
  
  for(i = 0; i<n ; i++)
    {
	double log_r = profile[i][0];
	double s = sigma[i];
	
	gsl_matrix_set(J, i, 0, 1.0 / s);
	gsl_matrix_set(J, i, 1,  (-3 * (pow(10, log_rs))/(pow(10, log_rs)+pow(10, log_r)) + 3) / s);					
	}
	
   return GSL_SUCCESS;	
}       /* end dens_df */	


int dens_fdf(const gsl_vector *x, void *data_now, gsl_vector *f, gsl_matrix *J)
{
  dens_f(x, data_now, f);
  dens_df(x, data_now, J);
  
  return GSL_SUCCESS;
	
}      /* end dens_fdf */


#endif
