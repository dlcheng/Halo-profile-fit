#include <stdio.h>
#include "allvars.h"
#include "define.h"
#include "proto.h"


double profile[BIN_NUM][6]; /* profile[][0] = log10(radius)
                                profile[][1] = rho
                                profile[][2] = total part num in the bin
                                profile[][3] = new born particle, if DECAY_DARK_MATTER 
                             */
double bin_rt_bd[BIN_NUM];  /* right boundary of each bin */

#ifdef VEL_PROFILE
double vel_profile[VEL_BIN_END - VEL_BIN_START + 1][VEL_BIN_NUM];
#endif

int each_file_num;
int file_name_num;
gadget_head gh;

FILE *fp_1;        
FILE *fp_2;                  /* point to the mass block, if needed */
FILE *fp_3;                  /* point to the ID block */
#ifdef VEL_PROFILE
FILE *fp_4;                   /* point to the VEL block */
#endif

#ifdef DECAY_DARK_MATTER
double total_mass;
double decay_mass;
double decay_fraction;
#endif

double min_cut_r;            /* the lower limit cutoff radius */
