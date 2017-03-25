#include <math.h>
#include <stdio.h>

#include "allvars.h"
#include "define.h"
#include "proto.h"

void process_all_files()
{
  int i;

  for(i = 0; i < TOTAL_FILE_NUM; i++)
    { 
      re_set_data();
      file_name_num = FILE_START_NUM + i * FILE_DIS;

      if(process_single_file())
        {
        write_out_put_file();
        fit_data(0);
#ifdef	DECAY_DARK_MATTER
        fit_data(1);
#endif        
	    }
      else
        printf("\nFound no file(s) in the base of %s%03d\nskip it.\n\n", INPUT_FILE_BASE, file_name_num);
    }
}                            /* end process_all_files */

int process_single_file()
{
  char file_name_test[300];
  char file_name[300];
  unsigned int i, id, j;
  float pos[3], mass, vel[3];

  sprintf(file_name_test, "%s%s%03d", INPUT_FOLDER, INPUT_FILE_BASE, file_name_num);   
    
  if((fp_1 = fopen(file_name_test, "rb")) != NULL)
    {
      printf("\n\n\nFound single-file snapshot %s%03d.\n", INPUT_FILE_BASE, file_name_num);
      each_file_num = 1;
      read_the_head();
      fclose(fp_1);
    }
  else
    { 
      sprintf(file_name_test, "%s%s%03d%s%d", INPUT_FOLDER, INPUT_FILE_BASE, file_name_num,".", 0);   
	
      if((fp_1 = fopen(file_name_test, "rb")) != NULL)
        read_the_head();
      else 
        return 0;
	
      each_file_num = gh.NumFiles;  
      printf("\n\n\nFound multi-files snapshot %s%03d.(0-%d)\n", INPUT_FILE_BASE, file_name_num, each_file_num - 1);        
      fclose(fp_1);
    }

#ifdef DECAY_DARK_MATTER
  printf("gh.Old_num = %d, gh.Old_soften_length = %.4f\n", gh.Old_num, gh.Old_soften_length);
#endif

  init_boundary();

  for(i = 0; i < each_file_num; i++)
    {
      if(each_file_num > 1)
        { 
          printf("processing %s%03d%s%d . . .\n", INPUT_FILE_BASE, file_name_num, ".", i);
          sprintf(file_name, "%s%s%03d%s%d", INPUT_FOLDER, INPUT_FILE_BASE, file_name_num,".", i);
        }
      else
        {
          printf("processing %s%03d . . .\n", INPUT_FILE_BASE, file_name_num);
          sprintf(file_name, "%s%s%03d", INPUT_FOLDER, INPUT_FILE_BASE, file_name_num);          
        }
        
      fp_1 = fopen(file_name, "rb");     /* locate the pos block */
      read_the_head();
      locate_pos_block();
      
      if(gh.Massarr[1] == 0)             /* need to read the mass block */
        {
          fp_2 = fopen(file_name, "rb");
          locate_mass_block();
        }

      fp_3 = fopen(file_name, "rb");    /* locate the id block */
      locate_id_block();        

#ifdef VEL_PROFILE
      fp_4 = fopen(file_name, "rb");    /* locacte the VEL block */
      locate_vel_block();
#endif 

      for(j = 0; j < gh.Npart[1]; j++)
        {
          fread(pos, 3 * sizeof(float), 1, fp_1);

          if(gh.Massarr[1] == 0)
            fread(&mass, sizeof(float), 1, fp_2);
          else
            mass = gh.Massarr[1];

          fread(&id, sizeof(int), 1, fp_3);

#ifdef VEL_PROFILE
          fread(vel, 3 * sizeof(float), 1, fp_4);
#endif
            
          update_data(pos, mass, id, vel);
        }
      
      fclose(fp_1);

      if(gh.Massarr[1] == 0) 
        fclose(fp_2);

      fclose(fp_3);

      printf("Done\n");
    }

  calculate_density();

  return 1;
}                                                              /* end process_single_file */

void calculate_density()
{
  double r_left, r_right, r_mid, rho;
#ifdef DECAY_DARK_MATTER
  double rho_decay;
#endif  
  int i;

#ifdef LOG_SCALE
  double log_dis = (log10(DET_RADIUS) - log10(min_cut_r)) / (double) BIN_NUM;
#endif

  for(i = 0; i < BIN_NUM; i++)
    {
      if(i == 0)
        r_left = min_cut_r;
      else
        r_left = bin_rt_bd[i - 1];

      r_right = bin_rt_bd[i];

      r_mid = (r_left + r_right) / 2.0;
      rho = profile[i][1] / (4.0 * PI / 3.0 * (r_right * r_right * r_right - r_left * r_left * r_left));
      
#ifdef DECAY_DARK_MATTER
      rho_decay = profile[i][4] / (4.0 * PI / 3.0 * (r_right * r_right * r_right - r_left * r_left * r_left));
      decay_fraction = decay_mass / total_mass;
#endif      
        
#ifdef LOG_SCALE
      profile[i][0] = log10(r_right) - (log_dis / 2.0);
      profile[i][1] = log10(rho); 
#else
      profile[i][0] = r_mid;
      profile[i][1] = rho;     
#endif

#ifdef DECAY_DARK_MATTER
#ifdef LOG_SCALE
      profile[i][4] = log10(rho_decay);
#else
      profile[i][4] = rho_decay;
#endif
      profile[i][5] = rho_decay / rho;
#endif

    }
}                                                                /* end calculate_density */

void update_data(float *pos, float mass, unsigned int id, float *vel)
{
  int i;
  double r = 0;
#ifdef VEL_PROFILE
  double v = 0;
  double v_interval;
  int j;
#endif

#ifdef DECAY_DARK_MATTER
  total_mass += mass;  
  if(id >= gh.Old_num)    /* this is a new particle */
    decay_mass += mass;
#endif

  for(i = 0; i < 3; i++)
    {
      r += pos[i] * pos[i];
#ifdef VEL_PROFILE
      v += vel[i] * vel[i];
#endif      
    }
  r = sqrt(r);
#ifdef VEL_PROFILE
  v = sqrt(v);
  v_interval = (VEL_MAX - VEL_MIN) / (double) VEL_BIN_NUM;
#endif

  if(r < DET_RADIUS && r > min_cut_r)  
    {
      for(i = 0; r > bin_rt_bd[i] && i < BIN_NUM; i++);
      profile[i][1] += mass;
      profile[i][2]++;
#ifdef DECAY_DARK_MATTER
      if(id >= gh.Old_num)    /* this is a new particle */
        {
        profile[i][3]++;
        profile[i][4] += mass;        
	    }
#endif        

#ifdef VEL_PROFILE
      i = i - VEL_BIN_START + 1;
      if(i >= 0 && i < VEL_BIN_END - VEL_BIN_START + 1)
        {
          if(v < VEL_MAX && v > VEL_MIN)
            {
              j = (v - VEL_MIN) / v_interval;
              if(j > VEL_BIN_NUM)
                j = VEL_BIN_NUM;
              
              vel_profile[i][j] += mass;
            }
        }
#endif
    }
}                                                                /* end update_data */

void read_the_head()
{
  unsigned int block_size;
  rewind(fp_1);

  fread(&block_size, sizeof(int), 1, fp_1);
  fread(&gh, sizeof(gadget_head), 1, fp_1);
  fread(&block_size, sizeof(int), 1, fp_1);

  rewind(fp_1);

  if(block_size != 256)
    printf("Error detected in reading the header.\n");
}                                                                /* end read_the_head */

void locate_pos_block()
{
  unsigned int block_size_1, block_size_2;
  long int disp;
  int i;

  for(i = 0; i < 1; i++)
    {
      fread(&block_size_1, sizeof(int), 1, fp_1);

      disp = block_size_1;

      fseek(fp_1, disp, 1);

      fread(&block_size_2, sizeof(int), 1, fp_1);     

      if(block_size_2 != block_size_1)
        printf("Error detected in locate POS block.\n");
    }

  fread(&block_size_2, sizeof(int), 1, fp_1); 
}                                                                /* end locate_pos_block */


void locate_mass_block()
{
  unsigned int block_size_1, block_size_2;
  long int disp;
  int i;

  for(i = 0; i < 4; i++)
    {
      fread(&block_size_1, sizeof(int), 1, fp_2);

      disp = block_size_1;

      fseek(fp_2, disp, 1);

      fread(&block_size_2, sizeof(int), 1, fp_2);     

      if(block_size_2 != block_size_1)
        printf("Error detected in locate MASS block.\n");
    }

  fread(&block_size_2, sizeof(int), 1, fp_2);  
}                                                               /* end locate_mass_block */

void locate_id_block()
{
  unsigned int block_size_1, block_size_2;
  long int disp;
  int i;

  for(i = 0; i < 3; i++)
    {
      fread(&block_size_1, sizeof(int), 1, fp_3);

      disp = block_size_1;

      fseek(fp_3, disp, 1);

      fread(&block_size_2, sizeof(int), 1, fp_3);     

      if(block_size_2 != block_size_1)
        printf("Error detected in locate ID block.\n");
    }

  fread(&block_size_2, sizeof(int), 1, fp_3); 
}                                                                /* end locate_id_block */

#ifdef VEL_PROFILE
void locate_vel_block()
{
  unsigned int block_size_1, block_size_2;
  long int disp;
  int i;

  for(i = 0; i < 2; i++)
    {
      fread(&block_size_1, sizeof(int), 1, fp_4);

      disp = block_size_1;

      fseek(fp_4, disp, 1);

      fread(&block_size_2, sizeof(int), 1, fp_4);     

      if(block_size_2 != block_size_1)
        printf("Error detected in locate VEL block.\n");
    }

  fread(&block_size_2, sizeof(int), 1, fp_4); 
}                                                                /* end locate_vel_block */
#endif
