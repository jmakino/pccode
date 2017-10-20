c-----------------------------------------------------------------------
c intgrt3.f : General order Krogh scheme for first order O.D.Es
c-----------------------------------------------------------------------
#include "pcheader.h"
c
      subroutine fconvd(d, coef, order, ndim, pcoef)
c-----------------------------------------------------------------------
c   fconvd: transform divided difference to power series coefficient
c         using the table of the expanded newton interpolation
c         polynominal (by convolution)
c-----------------------------------------------------------------------
      INTEGER ndim
      REAL d(ndim,0:MAXSTEPS)
      REAL coef(MAXSTEPS, MAXSTEPS)
      INTEGER order
      REAL pcoef(ndim,0:MAXSTEPS)
      INTEGER j,k,l
      do 10 k=1,ndim
         pcoef(k,0)=d(k,0)
         do 20 j=1,order
            pcoef(k,j) = d(k,j)
            do 30 l=j+1,order
               pcoef(k,j) = pcoef(k,j) + coef(j,l)*d(k,l)
c               write(6,*) k,j,l,pcoef(k,j), coef(j,l), d(k,l)
30          continue
20       continue
10    continue
      end
c      
c                  
      subroutine fdtopo(d, t, order, ndim, pcoef)
c-----------------------------------------------------------------------
c   d		in	divided difference table
c   tlist	in 	time list for divided differences (t0=0)
c   order	in	order of the polynominal
c   ndim	in	number of equations
c   pcoef       out     power series coefficients at t0
c
c-----------------------------------------------------------------------
      INTEGER ndim
      REAL d(ndim,0:MAXSTEPS)
      REAL t(MAXSTEPS)
      INTEGER order
      REAL pcoef(ndim, 0:MAXSTEPS)
      REAL excoef(MAXSTEPS, MAXSTEPS)
c
      call expnwp(t, order, excoef)
      call fconvd(d, excoef, order, ndim, pcoef)
      end
c
      subroutine fupddt(ddt, t, t0, acc0, order, ndim, dtlist)
c-----------------------------------------------------------------------
c   upddt : update divided difference table and time list
c-----------------------------------------------------------------------
c   ddt		in/out	divided difference table
c   t		in/out	time list
c   t0		in 	new reference time
c   acc0	in	acceleration at t=t0
c   order	in	order of polynominal
c   ndim	in	number of equations
c   dtlist      out     list of (t_k - t_0)
c   
c   actual work of updating D is done in UPDDT1
c-----------------------------------------------------------------------
      INTEGER ndim
      REAL ddt(ndim,0:MAXSTEPS),t(0:MAXSTEPS), dtlist(MAXSTEPS)
      REAL t0, acc0(ndim)
      INTEGER order
      REAL dnew(NDIM*2,0:MAXSTEPS),tnew(0:MAXSTEPS)
      INTEGER i,k
      call fupdd0(ddt, t, t0, acc0, order, ndim, dnew, tnew)
      call fupdd2(dnew, tnew, order, ndim, ddt, t, dtlist)
      end
c      
      subroutine fupdd2(dnew, tnew, order, ndim, ddt, t, dtlist)
c-----------------------------------------------------------------------
c   fupdd2 : copy new divided difference table and time list 
c-----------------------------------------------------------------------
c   dnew	in	divided difference table
c   tnew	in	time list
c   order	in	order of polynominal
c   ddt		out	divided difference table
c   t		out	time list
c   dtlist      out     list of (t_k - t_0)
c   
c-----------------------------------------------------------------------
      INTEGER ndim
      REAL ddt(ndim,0:MAXSTEPS),t(0:MAXSTEPS), dtlist(MAXSTEPS)
      INTEGER order
      REAL dnew(ndim,0:MAXSTEPS),tnew(0:MAXSTEPS)
      INTEGER i,k
      do 20 k=1,ndim
         do 10 i=0,order
            ddt(k,i) = dnew(k,i)
10       continue
20    continue
      do 30 i=0,order
         t(i) = tnew(i)
30    continue
      do 40 i = 1, order
         dtlist(i) = t(i)-t(0)
40    continue
      end
c      
      subroutine setdtl(t,dt,order)
c-----------------------------------------------------------------------
c  setdtl : set the dt list from time list
c-----------------------------------------------------------------------
c  t		in	time list
c  dt		out	dt list
c  order	in	order of integrator
c-----------------------------------------------------------------------
      REAL t(0:MAXSTEPS), dt(MAXSTEPS)
      INTEGER order, i
      do 40 i = 1, order
         dt(i) = t(i)-t(0)
40    continue
      end
      
      subroutine fupdd0(dold, told, t0, acc0, order, ndim, dnew, tnew)
c-----------------------------------------------------------------------
c   upddt0 : update divided difference table and time list
c-----------------------------------------------------------------------
c   dold	in	old divided difference table
c   told	in	old time list
c   t0		in 	new reference time
c   acc0	in	acceleration at t=t0
c   order	in	order of polynominal
c   dnew	out	new divided difference table
c   tnew	out	new time list
c
c   D is updated using the following relation:
c   dnew_0 = acc0
c   dnew_k = (dnew_k-1 - dold_k-1)/(t_0-t_k)
c-----------------------------------------------------------------------
      INTEGER ndim
      REAL dold(ndim,0:MAXSTEPS),told(0:MAXSTEPS)
      REAL t0, acc0(ndim)
      INTEGER order
      REAL dnew(ndim,0:MAXSTEPS),tnew(0:MAXSTEPS)
      INTEGER i,k, im1
      REAL dtinv
c-----------------------------------------------------------------------      
c     First set D0, t0
c-----------------------------------------------------------------------      
      do 10 k=1,ndim
         dnew(k,0) = acc0(k)
10    continue
      tnew(0) = t0
c-----------------------------------------------------------------------      
c     calculate D using the iteration formula
c-----------------------------------------------------------------------      
      do 20 i = 1, order
         im1=i-1
         tnew(i) = told(im1)
         dtinv = 1.0d0/(t0 - told(im1))
         do 30 k = 1, ndim
            dnew(k,i) = (dnew(k,im1)-dold(k,im1))*dtinv
30       continue
20    continue
      end
c
c      
      subroutine fcorre(pcoef, dtoldm, order, ndim, x)
c-----------------------------------------------------------------------      
c   correc : apply corrector term to x, v
c-----------------------------------------------------------------------
c   pcoef	in	power-series coefficients of acc ar t=t0
c   dtoldm	in	t_-1 - t0
c   order       in      order of polynominal
c   x		in	position 
c-----------------------------------------------------------------------
      INTEGER ndim
      REAL pcoef(ndim, 0:MAXSTEPS), dtoldm, x(ndim)
      INTEGER order
      REAL xnew(NDIM*2)
      INTEGER k
      call fcorr0(pcoef, dtoldm, x, order, ndim, xnew)
      do 10 k = 1, ndim
         x(k)=xnew(k)
10    continue
      end
c      
      subroutine fcorr0(pcoef, dtoldm, xold, order, ndim, xnew)
c-----------------------------------------------------------------------      
c   correc : obtain corrector term for x, v
c-----------------------------------------------------------------------
c   pcoef	in	power-series coefficients of acc ar t=t0
c   dtoldm	in	t_-1 - t0
c   xold	in	position at t=t_-1
c   order       in      order of polynominal
c   xnew	out	new position at t=t0
c-----------------------------------------------------------------------
      INTEGER ndim
      REAL pcoef(ndim, 0:MAXSTEPS), dtoldm, xold(ndim)
      INTEGER order
      REAL xnew(ndim)
      REAL tpower, xfact(0:MAXSTEPS), vfact(0:MAXSTEPS)
      INTEGER i, k

      REAL xcoef(0:MAXSTEPS), vcoef(0:MAXSTEPS)
      common/rtabs/ xcoef, vcoef

      tpower = dtoldm*dtoldm
      vfact(0)=dtoldm
      do 20 i=1,order
         vfact(i) = tpower*vcoef(i)
         tpower = tpower * dtoldm
20    continue
      do 30 k = 1, ndim
         xnew(k)=xold(k)
         do 40 i=0,order
            xnew(k) = xnew(k) - pcoef(k,i)*vfact(i)
40       continue
30    continue
      end
c
      subroutine fpredi(pcoef, dt, xold, order, ndim, xnew)
c-----------------------------------------------------------------------      
c   predic : obtain predictor for x, v
c-----------------------------------------------------------------------
c   pcoef	in	power-series coefficients of acc ar t=t0
c   dt		in	the timestep
c   xold	in	position at t=t_0
c   order       in      order of polynominal
c   xnew	out	new position at t=t0+dt
c-----------------------------------------------------------------------
      INTEGER ndim
      REAL pcoef(ndim, 0:MAXSTEPS), dt, xold(ndim)
      INTEGER order
      REAL xnew(ndim)
      REAL tpower, xfact(0:MAXSTEPS), vfact(0:MAXSTEPS)
      INTEGER i, k

      REAL xcoef(0:MAXSTEPS), vcoef(0:MAXSTEPS)
      common/rtabs/ xcoef, vcoef

      tpower = dt*dt
      vfact(0)=dt
      do 20 i=1,order
         vfact(i) = tpower*vcoef(i)
         tpower = tpower * dt
20    continue
      do 30 k = 1, ndim
         xnew(k)=xold(k)
         do 40 i=0,order
            xnew(k) = xnew(k) + pcoef(k,i)*vfact(i)
40       continue
30    continue
      end
