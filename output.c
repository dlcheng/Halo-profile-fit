#include <math.h>
#include <stdio.h>

#include "allvars.h"
#include "define.h"
#include "proto.h"

void write_out_put_file()
{
  FILE *fp;
  char file_name[100];
  int i;

#ifdef VEL_PROFILE
  FILE *fp_v;
  int j;
  double bd_l, bd_r;
  double v_interval;
#endif

#ifdef DECAY_DARK_MATTER
  printf("decayed mass fraction = %.6e\n", decay_fraction);
#endif
    
  sprintf(file_name, "%s%s%03d%s", OUTPUT_FOLDER, OUTPUT_FILE_BASE, file_name_num,".txt");
  fp = fopen(file_name, "w");

  for(i = 0; i < BIN_NUM; i++)
    {
       
      if(profile[i][2] >= RESO_NUM)
        {
#ifdef DECAY_DARK_MATTER 
          fprintf(fp, "%.6e	%.6e	%.6e	%.6e	%.0f	%.0f	%.0f\n", profile[i][0], profile[i][1], 
                                                                             profile[i][4], profile[i][5],
                                                                             profile[i][2], profile[i][2] - profile[i][3], 
                                                                             profile[i][3]);
#else
          fprintf(fp, "%.6e	%.6e	%.0f\n", profile[i][0], profile[i][1], profile[i][2]);     
#endif      
        }
    }

  fclose(fp);

#ifdef VEL_PROFILE
  v_interval = (VEL_MAX - VEL_MIN) / (double) VEL_BIN_NUM;  

  for(i = 0; i < VEL_BIN_END - VEL_BIN_START + 1; i++)
    {
      if(i + VEL_BIN_START - 2 < 0)
        bd_l = 0;
      else
        bd_l = bin_rt_bd[i + VEL_BIN_START - 2];

      bd_r = bin_rt_bd[i + VEL_BIN_START - 1];

      sprintf(file_name,"%s%s%03d%s%g%s%g%s%s", 
                         OUTPUT_FOLDER, "halo_vel_profile_", file_name_num, "_r(" ,bd_l ,"," ,bd_r ,")" ,".txt" );

      fp_v = fopen(file_name, "w");

      for(j = 0; j < VEL_BIN_NUM; j++)
        fprintf(fp_v, "%.4e	%.4e\n", VEL_MIN + ((double) j + 0.5) * v_interval, vel_profile[i][j]);
      
      fclose(fp_v);
    }
#endif
}                                                           /* end write_out_put_file */
