/*
  This code analyzes the density profile and velocity profile of a shperical halo. 
  Special care has been made to avoid cutting the softening length.

Density profile format:
#ifdef LOG_SCALE
  column 1    column 2    column 3          colum 4       column 5    column 6    column 7
  log10(r)   log10(rho) log10(rho_decay) decay_fraction   Num_inside   Old_part    New_part(ifdef DECAY_DARK_MATTER)

#else
  column 1    column 2    column 3          colum 4       column 5    column 6    column 7
     r          rho       rho_decay      decay_fraction   Num_inside   Old_part    New_part(ifdef DECAY_DARK_MATTER)
#endif    


Velocity profile format:
    column 1                 column 2
 velcocity bin                 mass 
*/

#include <math.h>
#include <stdio.h>

#include "allvars.h"
#include "define.h"
#include "proto.h"


int main()
{   
  process_all_files();

  return 0;
}                /* end main */

