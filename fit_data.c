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
void fit_data(int m)  /* m==0, fit the total density, m==1, fit the daughter density */
{ 

  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  unsigned int i, iter = 0;
  const int n = BIN_NUM;
  const int p = 3;
  gsl_multifit_function_fdf f;
  double x_init[3] = {GUESS_LOG_RHOS, GUESS_LOG_RS, GUESS_ALPHA};
	
  gsl_matrix *covar = gsl_matrix_alloc(p, p);
  double y[BIN_NUM], sigma[BIN_NUM];
  data d = {n, y, sigma};
  gsl_vector_view x = gsl_vector_view_array(x_init, p);
   
  char file_name[500];
  FILE *fp;
  
  f.f = &dens_f;
  f.df = &dens_df;
  f.fdf= &dens_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

/* this part prepares the data */  
  for(i = 0; i < n; i++)
    {
	if(m==0)	
	  {
	  y[i] = profile[i][1];   /* this is the total density */  
	  sigma[i] = 1;           /* unweighted */
	  }	  
	else
	  {
	  y[i] = profile[i][4];   /* this is the daughter particle density */  
	  sigma[i] = 1;           /* unweighted */
	  }	  	
	}

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);

  if(m==0)
    sprintf(file_name,"%s%s%s%03d%s", OUTPUT_FOLDER, OUTPUT_FILE_BASE, "_fit_t_", file_name_num,".txt");
  else
    sprintf(file_name,"%s%s%s%03d%s", OUTPUT_FOLDER, OUTPUT_FILE_BASE, "_fit_d_", file_name_num,".txt");
    
  fp = fopen(file_name, "w");
  
  if(m==0)
    {
    fprintf(fp, "%s\n", "This is the fitting log for total density.");
    printf("%s\n", "This is the fitting log for total density.");
    }
  else
    {
    fprintf(fp, "%s\n", "This is the fitting log for daughter density.");
    printf("%s\n", "This is the fitting log for daughter density.");
    }
    
  fprintf(fp, "iter: %3u x0 = %15.8f, x1 = %15.8f, x2 = %15.8f |f(x)| = %g\n",
          iter, 
          gsl_vector_get(s->x, 0),
          gsl_vector_get(s->x, 1),
          gsl_vector_get(s->x, 2),
          gsl_blas_dnrm2(s->f));
          
  printf("iter: %3u x0 = %15.8f, x1 = %15.8f, x2 = %15.8f |f(x)| = %g\n",
          iter, 
          gsl_vector_get(s->x, 0),
          gsl_vector_get(s->x, 1),
          gsl_vector_get(s->x, 2),
          gsl_blas_dnrm2(s->f));          
          
  do
    {
	 iter++;	
	 status = gsl_multifit_fdfsolver_iterate(s);
	 
     fprintf(fp, "iter: %3u x0 = %15.8f, x1 = %15.8f, x2 = %15.8f |f(x)| = %g\n",
             iter, 
             gsl_vector_get(s->x, 0),
             gsl_vector_get(s->x, 1),
             gsl_vector_get(s->x, 2),
             gsl_blas_dnrm2(s->f));	
             
     printf("iter: %3u x0 = %15.8f, x1 = %15.8f, x2 = %15.8f |f(x)| = %g\n",
             iter, 
             gsl_vector_get(s->x, 0),
             gsl_vector_get(s->x, 1),
             gsl_vector_get(s->x, 2),
             gsl_blas_dnrm2(s->f));	             
             
     if(status)        	
	   break;
	   
	 status = gsl_multifit_test_delta(s->dx, s->x, 0.0, 1e-4);  
		
	}while(status == GSL_CONTINUE && iter < 500);
	
	gsl_multifit_covar(s->J, 0.0, covar);
	
    double chi = gsl_blas_dnrm2(s->f);
    double dof = n - p;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof));
    double log_rhos = gsl_vector_get(s->x, 0);
    double log_rs = gsl_vector_get(s->x, 1);
    double alpha = gsl_vector_get(s->x, 2);
    
    fprintf(fp, "chisq/dof = %g\n", pow(chi, 2) / dof);
    fprintf(fp, "log_rhos = %.6f +/- %.6f\n", log_rhos, c*sqrt(gsl_matrix_get(covar, 0, 0)));
	fprintf(fp, "log_rs   = %.6f +/- %.6f\n", log_rs, c*sqrt(gsl_matrix_get(covar, 1, 1)));
	fprintf(fp, "alpha    = %.6f +/- %.6f\n", alpha, c*sqrt(gsl_matrix_get(covar, 2, 2)));	
	fprintf(fp, "status = %s\n", gsl_strerror(status));
    fprintf(fp, "%s\n", "The fitting formula (in log scale) is:");
	fprintf(fp, "%s%lf%s%lf%s%lf%s%lf%s%lf%s%lf%s\n", "(", log_rhos, ")-(",alpha,")*(x-(",log_rs,"))-(3-(",alpha,"))*(log10(10^x+10^(",log_rs,"))-(",log_rs,"))" );
	
	printf("chisq/dof = %g\n", pow(chi, 2) / dof);
    printf("log_rhos = %.6f +/- %.6f\n", gsl_vector_get(s->x, 0), c*sqrt(gsl_matrix_get(covar, 0, 0)));
	printf("log_rs   = %.6f +/- %.6f\n", gsl_vector_get(s->x, 1), c*sqrt(gsl_matrix_get(covar, 1, 1)));
	printf("alpha    = %.6f +/- %.6f\n", gsl_vector_get(s->x, 2), c*sqrt(gsl_matrix_get(covar, 2, 2)));	
	printf("status = %s\n", gsl_strerror(status));
    printf("%s\n", "The fitting formula (in log scale) is:");
	printf("%s%lf%s%lf%s%lf%s%lf%s%lf%s%lf%s\n", "(", log_rhos, ")-(",alpha,")*(x-(",log_rs,"))-(3-(",alpha,"))*(log10(10^x+10^(",log_rs,"))-(",log_rs,"))" );	
	
	gsl_multifit_fdfsolver_free(s);
	gsl_matrix_free(covar);
	fclose(fp);
	
}    /* end fit_data */


#endif

#ifdef BURKERT

void fit_data(int m)  /* m==0, fit the total density, m==1, fit the daughter density */
{ 

  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  unsigned int i, iter = 0;
  const int n = BIN_NUM;
  const int p = 2;   /* two parameters here */
  gsl_multifit_function_fdf f;
  double x_init[2] = {GUESS_LOG_RHOS, GUESS_LOG_RS};
	
  gsl_matrix *covar = gsl_matrix_alloc(p, p);
  double y[BIN_NUM], sigma[BIN_NUM];
  data d = {n, y, sigma};
  gsl_vector_view x = gsl_vector_view_array(x_init, p);
   
  char file_name[500];
  FILE *fp;
  
  f.f = &dens_f;
  f.df = &dens_df;
  f.fdf= &dens_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

/* this part prepares the data */  
  for(i = 0; i < n; i++)
    {
	if(m==0)	
	  {
	  y[i] = profile[i][1];   /* this is the total density */  
	  sigma[i] = 1;           /* unweighted */
	  }	  
	else
	  {
	  y[i] = profile[i][4];   /* this is the daughter particle density */  
	  sigma[i] = 1;           /* unweighted */
	  }	  	
	}

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);

  if(m==0)
    sprintf(file_name,"%s%s%s%03d%s", OUTPUT_FOLDER, OUTPUT_FILE_BASE, "fit_t_", file_name_num,".txt");
  else
    sprintf(file_name,"%s%s%s%03d%s", OUTPUT_FOLDER, OUTPUT_FILE_BASE, "fit_d_", file_name_num,".txt");
    
  fp = fopen(file_name, "w");
  
  if(m==0)
    {
    fprintf(fp, "%s\n", "This is the fitting log for total density.");
    printf("%s\n", "This is the fitting log for total density.");
    }
  else
    {
    fprintf(fp, "%s\n", "This is the fitting log for daughter density.");
    printf("%s\n", "This is the fitting log for daughter density.");
    }
    
  fprintf(fp, "iter: %3u x0 = %15.8f, x1 = %15.8f, |f(x)| = %g\n",
          iter, 
          gsl_vector_get(s->x, 0),
          gsl_vector_get(s->x, 1),
          gsl_blas_dnrm2(s->f));
          
  printf("iter: %3u x0 = %15.8f, x1 = %15.8f, |f(x)| = %g\n",
          iter, 
          gsl_vector_get(s->x, 0),
          gsl_vector_get(s->x, 1),
          gsl_blas_dnrm2(s->f));          
          
  do
    {
	 iter++;	
	 status = gsl_multifit_fdfsolver_iterate(s);
	 
     fprintf(fp, "iter: %3u x0 = %15.8f, x1 = %15.8f, |f(x)| = %g\n",
             iter, 
             gsl_vector_get(s->x, 0),
             gsl_vector_get(s->x, 1),
             gsl_blas_dnrm2(s->f));	
             
     printf("iter: %3u x0 = %15.8f, x1 = %15.8f, |f(x)| = %g\n",
             iter, 
             gsl_vector_get(s->x, 0),
             gsl_vector_get(s->x, 1),
             gsl_blas_dnrm2(s->f));	             
             
     if(status)        	
	   break;
	   
	 status = gsl_multifit_test_delta(s->dx, s->x, 0.0, 1e-4);  
		
	}while(status == GSL_CONTINUE && iter < 500);
	
	gsl_multifit_covar(s->J, 0.0, covar);
	
    double chi = gsl_blas_dnrm2(s->f);
    double dof = n - p;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof));
    double log_rhos = gsl_vector_get(s->x, 0);
    double log_rs = gsl_vector_get(s->x, 1);
    
    fprintf(fp, "chisq/dof = %g\n", pow(chi, 2) / dof);
    fprintf(fp, "log_rhos = %.6f +/- %.6f\n", log_rhos, c*sqrt(gsl_matrix_get(covar, 0, 0)));
	fprintf(fp, "log_rs   = %.6f +/- %.6f\n", log_rs, c*sqrt(gsl_matrix_get(covar, 1, 1)));
	fprintf(fp, "status = %s\n", gsl_strerror(status));
	fprintf(fp, "%s\n", "The fitting formula (in log scale) is:");
	fprintf(fp, "%s%lf%s%lf%s%lf%s%lf%s\n", "(", log_rhos, ")+3*(", log_rs,")-log10(10^x+10^(", log_rs, "))-log10(10^(2*x)+10^(2*", log_rs, "))");
	
	
	printf("chisq/dof = %g\n", pow(chi, 2) / dof);
    printf("log_rhos = %.6f +/- %.6f\n", gsl_vector_get(s->x, 0), c*sqrt(gsl_matrix_get(covar, 0, 0)));
	printf("log_rs   = %.6f +/- %.6f\n", gsl_vector_get(s->x, 1), c*sqrt(gsl_matrix_get(covar, 1, 1)));
	printf("status = %s\n", gsl_strerror(status));
	printf("%s\n", "The fitting formula (in log scale) is:");
	printf("%s%lf%s%lf%s%lf%s%lf%s\n", "(", log_rhos, ")+3*(", log_rs,")-log10(10^x+10^(", log_rs, "))-log10(10^(2*x)+10^(2*", log_rs, "))");
	
	gsl_multifit_fdfsolver_free(s);
	gsl_matrix_free(covar);
	fclose(fp);
	
}    /* end fit_data */



#endif


#ifdef BURKERT_NFW

void fit_data(int m)  /* m==0, fit the total density, m==1, fit the daughter density */
{ 

  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  int status;
  unsigned int i, iter = 0;
  const int n = BIN_NUM;
  const int p = 2;   /* two parameters here */
  gsl_multifit_function_fdf f;
  double x_init[2] = {GUESS_LOG_RHOS, GUESS_LOG_RS};
	
  gsl_matrix *covar = gsl_matrix_alloc(p, p);
  double y[BIN_NUM], sigma[BIN_NUM];
  data d = {n, y, sigma};
  gsl_vector_view x = gsl_vector_view_array(x_init, p);
   
  char file_name[500];
  FILE *fp;
  
  f.f = &dens_f;
  f.df = &dens_df;
  f.fdf= &dens_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

/* this part prepares the data */  
  for(i = 0; i < n; i++)
    {
	if(m==0)	
	  {
	  y[i] = profile[i][1];   /* this is the total density */  
	  sigma[i] = 1;           /* unweighted */
	  }	  
	else
	  {
	  y[i] = profile[i][4];   /* this is the daughter particle density */  
	  sigma[i] = 1;           /* unweighted */
	  }	  	
	}

  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc (T, n, p);
  gsl_multifit_fdfsolver_set (s, &f, &x.vector);

  if(m==0)
    sprintf(file_name,"%s%s%s%03d%s", OUTPUT_FOLDER, OUTPUT_FILE_BASE, "fit_t_", file_name_num,".txt");
  else
    sprintf(file_name,"%s%s%s%03d%s", OUTPUT_FOLDER, OUTPUT_FILE_BASE, "fit_d_", file_name_num,".txt");
    
  fp = fopen(file_name, "w");
  
  if(m==0)
    {
    fprintf(fp, "%s\n", "This is the fitting log for total density.");
    printf("%s\n", "This is the fitting log for total density.");
    }
  else
    {
    fprintf(fp, "%s\n", "This is the fitting log for daughter density.");
    printf("%s\n", "This is the fitting log for daughter density.");
    }
    
  fprintf(fp, "iter: %3u x0 = %15.8f, x1 = %15.8f, |f(x)| = %g\n",
          iter, 
          gsl_vector_get(s->x, 0),
          gsl_vector_get(s->x, 1),
          gsl_blas_dnrm2(s->f));
          
  printf("iter: %3u x0 = %15.8f, x1 = %15.8f, |f(x)| = %g\n",
          iter, 
          gsl_vector_get(s->x, 0),
          gsl_vector_get(s->x, 1),
          gsl_blas_dnrm2(s->f));          
          
  do
    {
	 iter++;	
	 status = gsl_multifit_fdfsolver_iterate(s);
	 
     fprintf(fp, "iter: %3u x0 = %15.8f, x1 = %15.8f, |f(x)| = %g\n",
             iter, 
             gsl_vector_get(s->x, 0),
             gsl_vector_get(s->x, 1),
             gsl_blas_dnrm2(s->f));	
             
     printf("iter: %3u x0 = %15.8f, x1 = %15.8f, |f(x)| = %g\n",
             iter, 
             gsl_vector_get(s->x, 0),
             gsl_vector_get(s->x, 1),
             gsl_blas_dnrm2(s->f));	             
             
     if(status)        	
	   break;
	   
	 status = gsl_multifit_test_delta(s->dx, s->x, 0.0, 1e-4);  
		
	}while(status == GSL_CONTINUE && iter < 500);
	
	gsl_multifit_covar(s->J, 0.0, covar);
	
    double chi = gsl_blas_dnrm2(s->f);
    double dof = n - p;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof));
    double log_rhos = gsl_vector_get(s->x, 0);
    double log_rs = gsl_vector_get(s->x, 1);
    
    fprintf(fp, "chisq/dof = %g\n", pow(chi, 2) / dof);
    fprintf(fp, "log_rhos = %.6f +/- %.6f\n", log_rhos, c*sqrt(gsl_matrix_get(covar, 0, 0)));
	fprintf(fp, "log_rs   = %.6f +/- %.6f\n", log_rs, c*sqrt(gsl_matrix_get(covar, 1, 1)));
	fprintf(fp, "status = %s\n", gsl_strerror(status));
	fprintf(fp, "%s\n", "The fitting formula (in log scale) is:");
	fprintf(fp, "%s%lf%s%lf%s%lf%s\n", "(", log_rhos, ")-3*log10(10^x+10^(",log_rs,"))+3*(", log_rs, ")");
	
	
	printf("chisq/dof = %g\n", pow(chi, 2) / dof);
    printf("log_rhos = %.6f +/- %.6f\n", gsl_vector_get(s->x, 0), c*sqrt(gsl_matrix_get(covar, 0, 0)));
	printf("log_rs   = %.6f +/- %.6f\n", gsl_vector_get(s->x, 1), c*sqrt(gsl_matrix_get(covar, 1, 1)));
	printf("status = %s\n", gsl_strerror(status));
	printf("%s\n", "The fitting formula (in log scale) is:");
	printf("%s%lf%s%lf%s%lf%s\n", "(", log_rhos, ")-3*log10(10^x+10^(",log_rs,"))+3*(", log_rs, ")");
	
	gsl_multifit_fdfsolver_free(s);
	gsl_matrix_free(covar);
	fclose(fp);
	
}    /* end fit_data */

#endif
















