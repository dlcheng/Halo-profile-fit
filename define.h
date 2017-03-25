/* Constants */
#define PI                3.14159265357
#define H0                1            /* hubble constant = 100 * H0 km/s/Mpc  */

/* Folders and files */
#define INPUT_FOLDER      "/users/s0902248/Install/Gadget2/Decay-dark-matter-halos-new/Result/sample-3/"   
#define OUTPUT_FOLDER     "./Result/sample-3/"

#define INPUT_FILE_BASE   "snapshot_"
#define OUTPUT_FILE_BASE  "halo_profile_"

/* Reading rules */
#define FILE_START_NUM    0            /* start form 0 */
#define FILE_DIS          1            /* file number increse step */
#define TOTAL_FILE_NUM    11           /* total file number */

/* Parameters */
#define DET_RADIUS        206.2        /* probe profile to this radius, in unit of kpc/h */
#define BIN_NUM           25           /* total bin numbers  */
#define RESO_NUM          0            /* only output bin density if its particle number > RESO_NUM */
#define ABANDON_R         (0.015 * DET_RADIUS)  
                                       /* particles inside this radius are discarded */

/* Not important */ 
#ifndef DECAY_DARAK_MATTER
#define SOFTEN_LENGTH     0.825        /* if not define decay dark matter */
#endif
#define FT_BIN_FACTOR     1            /* the first bin radius = Softening length * FT_BIN_FACTOR  */

#ifdef VEL_PROFILE
#define VEL_MAX           1000         /* consider particles with velocity less VEL_MAX for velocity statistic, in unit of km/s */
#define VEL_MIN           0
#define VEL_BIN_NUM       10           /* how many velocity bins */
#define VEL_BIN_START     2            /* start from the radius bin number VEL_BIN_START, if start from the first bin, set it to 1 */
#define VEL_BIN_END       3            /* end with the radius bin number VEL_BIN_END */
#endif

/* GUESS FIT PARAMETERS */
#define GUESS_LOG_RHOS   log10(6.0895824e-4)
#define GUESS_LOG_RS     log10(20)
#define GUESS_ALPHA      1
