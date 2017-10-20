c-----------------------------------------------------------------------
c
c  integrate.F
c
c  Integrator of Krogh-scheme type, Individual timestep PC scheme     
c
c  J. Makino
c
c  89-03-02 ver 0.00
c
c-----------------------------------------------------------------------
#include "pcheader.h"
      subroutine getdx(dmat, tvec, psteps, polow, csteps, 
     $                 xold, vold, xc, vc, tnew, xp, xpcoef)
c-----------------------------------------------------------------------     
c  dmat		in	divided difference table
c  tvec		in	time table (current time = 0)
c  psteps	in	steps of predictor
c  polow        in      order of lower-order predictor
c  csteps       in      steps of corrector
c  xold         in      position at t=t_-1
c  vold         in      velocity at t=t_-1
c  xc           out     position at t=t_0
c  vc           out     velocity at t=t_0
c  tnew         out     next time (t_0 =0)
c  xp           out     predicted position
c  xpcoef       out     predictor coefficients (including xc, vc)
c
c  Outline of the algorithm:
c
c  step1) form the power series expression 
c-----------------------------------------------------------------------     
      REAL dmat(NDIM,0:MAXSTEPS)
      REAL tvec(MAXSTEPS)
      INTEGER psteps, csteps, polow
      REAL xold(NDIM),vold(NDIM), xc(NDIM), vc(NDIM)
      REAL tnew, xp(NDIM), xpcoef(NDIM, 0:MAXSTEPS)
      
