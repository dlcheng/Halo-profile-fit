#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>


#include "allvars.h"
#include "define.h"

void init_boundary();
void re_set_data();

void write_out_put_file();

void process_all_files();
int process_single_file();

void calculate_density();
void update_data(float *pos, float mass, unsigned int id, float *vel);
void read_the_head();
void locate_pos_block();
void locate_mass_block();
void locate_id_block();
#ifdef VEL_PROFILE
void locate_vel_block();
#endif

void fit_data(int m);

int dens_f(const gsl_vector *x, void *data_now, gsl_vector *f);
int dens_df (const gsl_vector *x, void *data_now, gsl_matrix *J);
int dens_fdf(const gsl_vector *x, void *data_now, gsl_vector *f, gsl_matrix *J);
