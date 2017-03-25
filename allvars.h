#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include "define.h"
#include "proto.h"

typedef struct
{  
  unsigned int Npart[6];
  double Massarr[6];
  double Time;
  double Redshift;
  int FlagSfr;
  int FlagFeedback;
  int Nall[6];
  int  FlagCooling;
  int NumFiles;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int FlagAge;
  int FlagMetals;
  int NallHW[6];
  int Flag_entr_ics;
#ifndef DECAY_DARK_MATTER
  char unused[60];         
#else
  unsigned int Old_num;
  double Old_soften_length;
  char fill[48];
#endif
} gadget_head; 

typedef struct {
  int n;
  double *y;
  double *sigma;	
} data;



extern double profile[BIN_NUM][6]; /* profile[][0] = log10(radius)
                                        profile[][1] = rho/log_rho
                                        profile[][2] = total part num in the bin
                                        profile[][3] = new born particle, if defined DECAY_DARK_MATTER 
                                        profile[][4] = density of the new born particle, if defined DECAY_DARK_MATTER
                                        profile[][5] = fraction of the new born particle, if defined DECAY_DARK_MATTERE
                                     */
extern double bin_rt_bd[BIN_NUM];  /* right boundary of each bin */

#ifdef VEL_PROFILE
extern double vel_profile[VEL_BIN_END - VEL_BIN_START + 1][VEL_BIN_NUM];
#endif
extern int each_file_num;
extern int file_name_num;
extern gadget_head gh;

extern FILE *fp_1;        
extern FILE *fp_2;                  /* point to the mass block, if needed */
extern FILE *fp_3;                  /* point to the ID block */
#ifdef VEL_PROFILE
extern FILE *fp_4;
#endif

#ifdef DECAY_DARK_MATTER
extern double total_mass;
extern double decay_mass;
extern double decay_fraction;
#endif

extern double min_cut_r;            /* the lower limit cutoff radius */

#endif
