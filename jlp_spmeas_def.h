/********************************************************
* JLP
* Version 13/03/2020 
*********************************************************/
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
