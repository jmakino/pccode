c-----------------------------------------------------------------------
c  intgrtb.F : Krogh/Hermite scheme general order N-body integrator
c
c   Jun Makino  
c
c  89-03-27 Ver 0.1 created
c  89-07-06 Ver 0.2 RKG type round-off compensation is included in
c                   routine correc
c
c-----------------------------------------------------------------------
#include "pcheader.h"
      subroutine correb(pcoef, dtoldm, order, x, v, xq, vq)
c-----------------------------------------------------------------------      
c   correc : apply corrector term to x, v
c-----------------------------------------------------------------------
c   pcoef	in	power-series coefficients of acc ar t=t0
c   dtoldm	in	t_-1 - t0
c   order       in      order of polynominal
c   x		in/out	position 
c   v		in/out	velocity 
c   xq		in/out	position round-off compensation data
c   vq		in/out	velocity round-off compensation data
c-----------------------------------------------------------------------
      REAL pcoef(NDIM, 0:MAXSTEPS), dtoldm, x(NDIM), v(NDIM)
      INTEGER order
      REAL xnew(NDIM), vnew(NDIM)
      REAL xq(NDIM), vq(NDIM)
      INTEGER k 
      call corrb0(pcoef, dtoldm, x, v, order, xnew, vnew, xq, vq)
      do 10 k = 1, NDIM
         x(k)=xnew(k)
         v(k)=vnew(k)
10    continue
      end
c      
      subroutine corrb0(pcoef, dtoldm, xold, vold, order, xnew, vnew,
     $                  xq, vq)
c-----------------------------------------------------------------------      
c   correc : obtain corrector term for x, v
c-----------------------------------------------------------------------
c   pcoef	in	power-series coefficients of acc ar t=t0
c   dtoldm	in	t_-1 - t0
c   xold	in	position at t=t_-1
c   vold	in	velocity at t=t_-1
c   order       in      order of polynominal
c   xnew	out	new position at t=t0
c   vnew	out	new velosity at t=t0    
c   xq		in/out	position round-off compensation data
c   vq		in/out	velocity round-off compensation data
c-----------------------------------------------------------------------
      REAL pcoef(NDIM, 0:MAXSTEPS), dtoldm, xold(NDIM), vold(NDIM)
      INTEGER order
      REAL xnew(NDIM), vnew(NDIM)
      REAL xq(NDIM), vq(NDIM)
      REAL tpower, xfact(0:MAXSTEPS), vfact(0:MAXSTEPS)
      INTEGER i, k

      REAL xcoef(0:MAXSTEPS), vcoef(0:MAXSTEPS)
      common/rtabs/ xcoef, vcoef

      tpower = dtoldm*dtoldm
      xfact(0)=tpower*xcoef(0)
      vfact(0)=dtoldm
      do 20 i=1,order
         vfact(i) = tpower*vcoef(i)
         tpower = tpower * dtoldm
         xfact(i) = tpower*xcoef(i)
20    continue
      do 30 k = 1, NDIM
         vnew(k)=vold(k)
         xnew(k)=xold(k)
         do 40 i=0,order
c            xnew(k) = xnew(k) - pcoef(k,i)*xfact(i)
c            vnew(k) = vnew(k) - pcoef(k,i)*vfact(i)
            call sumval(xnew(k), - pcoef(k,i)*xfact(i), xq(k))
            call sumval(vnew(k), - pcoef(k,i)*vfact(i), vq(k))
40       continue
c         xnew(k) = xnew(k) - dtoldm * vnew(k)
          call sumval(xnew(k), - dtoldm * vnew(k), xq(k))
30    continue
      end
c
      subroutine sumval(sum, val, work)
c-----------------------------------------------------------------------      
c sumval : sum the value using round-off compemsating scheme
c-----------------------------------------------------------------------
c  sum		in/out	the accumulator
c  val		in	value to sum 
c  work		in/out	work valiable
c
c  
c-----------------------------------------------------------------------      
      REAL sum, val, work
      REAL r, w
      r = val - work
      w = sum
      sum = sum + r
      r = sum - w
      work = work + r - val
      end
