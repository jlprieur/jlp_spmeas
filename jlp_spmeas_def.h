/********************************************************
* JLP
* Version 13/03/2020 
*********************************************************/
#ifndef __jlp_spmeas_def
#define __jlp_spmeas_def

/* See also M_PI in math.h ... */
#ifndef PI 
#define PI 3.14159265358979323846 
#endif

#ifndef MAXI
#define MAXI(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MINI
#define MINI(x,y) ((x) > (y) ? (y) : (x))
#endif
#ifndef ABS 
#define ABS(x) ((x) > 0. ? (x) : (-(x)))
#endif
#ifndef SQUARE 
#define SQUARE(x) ((x) * (x))
#endif
#ifndef NINT 
#define NINT(x) (int)((x) + 0.5)
#endif

typedef struct {
double error_theta10_deg_max;
double rel_error_rho10_max;
double rho10_min;
double n_left_over_right_min;
double n_left_over_right_max;
double neg_max_for_patch;
} BIN_PARAM;

#endif
