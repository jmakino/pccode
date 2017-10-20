#ifndef PARMS
#define PARMS
#define MAXSTEPS 15 /* the highest  stepnumber of the integrator*/
#define MAXSTEPSP MAXSTEPS+1
#define REAL real*8
#define INTEGER integer*4
#define NDIM 3
#ifdef LARGEN
#   define MAXBODY 3000
#else
#   define MAXBODY 500
#endif
#define TIMEARG real*8
#endif
