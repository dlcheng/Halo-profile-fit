#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "allvars.h"
#include "define.h"
#include "proto.h"

void init_boundary()
{
  int i;
  double inter_temp, bin_temp, bin_min;
  
  min_cut_r = ABANDON_R;
  
#ifndef IGNORE_SL_CUT  

#ifdef DECAY_DARK_MATTER  
  bin_min = gh.Old_soften_length * FT_BIN_FACTOR;
#else
  bin_min = SOFTEN_LENGTH * FT_BIN_FACTOR;
#endif

  if(BIN_NUM * bin_min >= DET_RADIUS)
    {
      printf("Too much bins, can't proceed without cutting the SOFTEN_LENGTH!\n");
      exit(0);
    }	

  bin_rt_bd[0] = bin_min;

  for(i = 1; i < BIN_NUM; i++)
    {
#ifdef LOG_SCALE
  inter_temp = log10(DET_RADIUS / bin_rt_bd[i - 1]) / ((double) BIN_NUM - i);
  bin_temp = pow(10.0, log10(bin_rt_bd[i -1]) + inter_temp);	

  if((bin_temp - bin_rt_bd[i - 1]) >= bin_min)
     bin_rt_bd[i] = bin_temp;
  else
     bin_rt_bd[i] = bin_rt_bd[i - 1] + bin_min;	    
#else
  inter_temp = (DET_RADIUS - bin_rt_bd[i - 1]) / ((double) BIN_NUM - i);
  bin_temp = bin_rt_bd[i - 1] + inter_temp;	

  if((bin_temp - bin_rt_bd[i - 1]) >= bin_min)
     bin_rt_bd[i] = bin_temp;
  else
     bin_rt_bd[i] = bin_rt_bd[i - 1] + bin_min;	   
#endif
    }
    
#else
  double log_bin_temp = (log10(DET_RADIUS) - log10(min_cut_r)) / (double) BIN_NUM;  /* the low cut is 10^-3 * DET_RADIUS */
  for(i = 0 ; i < BIN_NUM; i++)
    bin_rt_bd[i] = pow(10, log_bin_temp * (i + 1)+log10(min_cut_r));
#endif        
}

void re_set_data()
{
 int i, j;

#ifdef DECAY_DARK_MATTER 
 total_mass = 0;
 decay_mass = 0;
#endif

 for(i = 0; i < BIN_NUM; i++)
   for(j = 0; j < 6; j++)  /* update the profile[] from 4 to 6 */
     profile[i][j] = 0;

#ifdef VEL_PROFILE
 for(i = 0; i < VEL_BIN_END - VEL_BIN_START + 1; i++)
   for(j = 0; j < VEL_BIN_NUM; j++)
     vel_profile[i][j] = 0;
#endif
}
