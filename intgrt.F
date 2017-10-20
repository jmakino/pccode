c-----------------------------------------------------------------------
c  intgrt.F : Krogh/Hermite scheme general order N-body integrator
c
c   Jun Makino  
c
c  89-03-27 Ver 0.1
c
c-----------------------------------------------------------------------
#include "pcheader.h"
c
      subroutine convd(d, coef, order, pcoef)
c-----------------------------------------------------------------------
c   comv: transform divided difference to power series coefficient
c         using the table of the expanded newton interpolation
c         polynominal (by convolution)
c-----------------------------------------------------------------------
      REAL d(NDIM,0:MAXSTEPS)
      REAL coef(MAXSTEPS, MAXSTEPS)
      INTEGER order
      REAL pcoef(NDIM,0:MAXSTEPS)
      INTEGER j,k,l
      do 10 k=1,NDIM
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
      subroutine expnwp(t, order, coef)
c-----------------------------------------------------------------------
c   expnwp : expand newton interporation polynominal
c-----------------------------------------------------------------------
c   t		in	time list for the divided difference table
c   order	in	order of fitting polynominal
c   coef	out	coefficient table for expanded polynominal
c
c   for the polynominal
c   Qk(t) = t(t-t_1)....(t-t_k-1)
c   (Q0(t)=1)
c
c   coef gives the expanded form:
c   Qk(t) = coef(0,k)+coef(1,k)*t+.....coef(k,k)*t^k
c
c   the fomula for coef is
c                 1                           (i=k)
c    coef(i,k) =  coef(i-1,k-1)-t_k*coef(i,k-1) (1<i<k)
c                 -t_k*coef(1,k-1)              (i=1)
c-----------------------------------------------------------------------
      REAL t(MAXSTEPS)
      INTEGER order
      REAL coef(MAXSTEPS, MAXSTEPS)
      INTEGER i,k
      REAL tkm1
      coef(1,1)=1.0D0
c      write(6,*)' Coef table'
      do 10 k=2,order
         tkm1 = t(k-1)
         coef(k,k)=1.0d0
         coef(1,k)=-tkm1*coef(1,k-1)
         do 20 i=2,k-1
            coef(i,k)=coef(i-1,k-1)-tkm1*coef(i,k-1)
20       continue
c         write(6,600)k,(coef(i,k),i=1,k), tkm1
c600      format(i5, (4g15.8))
10    continue
      end
c                  
      subroutine dtopow(d, t, order, pcoef)
c-----------------------------------------------------------------------
c   d		in	divided difference table
c   tlist	in 	time list for divided differences (t0=0)
c   order	in	order of the polynominal
c   pcoef       out     power series coefficients at t0
c
c-----------------------------------------------------------------------
      REAL d(NDIM,0:MAXSTEPS)
      REAL t(MAXSTEPS)
      INTEGER order
      REAL pcoef(NDIM, 0:MAXSTEPS)
      REAL excoef(MAXSTEPS, MAXSTEPS)
c
      call expnwp(t, order, excoef)
      call convd(d, excoef, order, pcoef)
      end
c
      subroutine upddt(ddt, t, t0, acc0, order, dtlist)
c-----------------------------------------------------------------------
c   upddt : update divided difference table and time list
c-----------------------------------------------------------------------
c   ddt		in/out	divided difference table
c   t		in/out	time list
c   t0		in 	new reference time
c   acc0	in	acceleration at t=t0
c   order	in	order of polynominal
c   dtlist      out     list of (t_k - t_0)
c   
c   actual work of updating D is done in UPDDT1
c-----------------------------------------------------------------------
      REAL ddt(NDIM,0:MAXSTEPS),t(0:MAXSTEPS), dtlist(MAXSTEPS)
      REAL t0, acc0(NDIM)
      INTEGER order
      REAL dnew(NDIM,0:MAXSTEPS),tnew(0:MAXSTEPS)
      INTEGER i,k
      call upddt0(ddt, t, t0, acc0, order, dnew, tnew)
      call upddt2(dnew, tnew, order, ddt, t, dtlist)
      end
c      
      subroutine upddt2(dnew, tnew, order, ddt, t, dtlist)
c-----------------------------------------------------------------------
c   upddt2 : copy new divided difference table and time list 
c-----------------------------------------------------------------------
c   dnew	in	divided difference table
c   tnew	in	time list
c   order	in	order of polynominal
c   ddt		out	divided difference table
c   t		out	time list
c   dtlist      out     list of (t_k - t_0)
c   
c-----------------------------------------------------------------------
      REAL ddt(NDIM,0:MAXSTEPS),t(0:MAXSTEPS), dtlist(MAXSTEPS)
      INTEGER order
      REAL dnew(NDIM,0:MAXSTEPS),tnew(0:MAXSTEPS)
      INTEGER i,k
      do 20 k=1,NDIM
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
      subroutine upddt0(dold, told, t0, acc0, order, dnew, tnew)
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
      REAL dold(NDIM,0:MAXSTEPS),told(0:MAXSTEPS)
      REAL t0, acc0(NDIM)
      INTEGER order
      REAL dnew(NDIM,0:MAXSTEPS),tnew(0:MAXSTEPS)
      INTEGER i,k, im1
      REAL dtinv
c-----------------------------------------------------------------------      
c     First set D0, t0
c-----------------------------------------------------------------------      
      do 10 k=1,NDIM
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
         do 30 k = 1, NDIM
            dnew(k,i) = (dnew(k,im1)-dold(k,im1))*dtinv
30       continue
20    continue
      end
c
      subroutine strtab
c-----------------------------------------------------------------------      
c  strtab : set the value of valiables in the common block RTABS
c-----------------------------------------------------------------------
c  content of rtabs:
c  xcoef : 1, 1/(1*2), 1/(2*3),,,,,
c  vcoef : 1, 1/2, 1/3,,,,,,,
c-----------------------------------------------------------------------      

      REAL xcoef(0:MAXSTEPS), vcoef(0:MAXSTEPS)
      common/rtabs/ xcoef, vcoef

      integer i
      do 10 i = 0, MAXSTEPS
         vcoef(i) = 1.0d0/(i+1.0d0)
         xcoef(i) = 1.0d0/((i+1.0d0)*(i+2.0d0))
10    continue
      end
c      
      subroutine correc(pcoef, dtoldm, order, x, v)
c-----------------------------------------------------------------------      
c   correc : apply corrector term to x, v
c-----------------------------------------------------------------------
c   pcoef	in	power-series coefficients of acc ar t=t0
c   dtoldm	in	t_-1 - t0
c   order       in      order of polynominal
c   x		in	position 
c   v		in	velocity 
c-----------------------------------------------------------------------
      REAL pcoef(NDIM, 0:MAXSTEPS), dtoldm, x(NDIM), v(NDIM)
      INTEGER order
      REAL xnew(NDIM), vnew(NDIM)
      INTEGER k
      call corre0(pcoef, dtoldm, x, v, order, xnew, vnew)
      do 10 k = 1, NDIM
         x(k)=xnew(k)
         v(k)=vnew(k)
10    continue
      end
c      
      subroutine corre0(pcoef, dtoldm, xold, vold, order, xnew, vnew)
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
c-----------------------------------------------------------------------
      REAL pcoef(NDIM, 0:MAXSTEPS), dtoldm, xold(NDIM), vold(NDIM)
      INTEGER order
      REAL xnew(NDIM), vnew(NDIM)
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
            xnew(k) = xnew(k) - pcoef(k,i)*xfact(i)
            vnew(k) = vnew(k) - pcoef(k,i)*vfact(i)
40       continue
         xnew(k) = xnew(k) - dtoldm * vnew(k)
30    continue
      end
c
      subroutine predic(pcoef, dt, xold, vold, order, xnew, vnew)
c-----------------------------------------------------------------------      
c   predic : obtain predictor for x, v
c-----------------------------------------------------------------------
c   pcoef	in	power-series coefficients of acc ar t=t0
c   dt		in	the timestep
c   xold	in	position at t=t_0
c   vold	in	velocity at t=t_0
c   order       in      order of polynominal
c   xnew	out	new position at t=t0+dt
c   vnew	out	new velosity at t=t0+dt 
c-----------------------------------------------------------------------
      REAL pcoef(NDIM, 0:MAXSTEPS), dt, xold(NDIM), vold(NDIM)
      INTEGER order
      REAL xnew(NDIM), vnew(NDIM)
      REAL tpower, xfact(0:MAXSTEPS), vfact(0:MAXSTEPS)
      INTEGER i, k

      REAL xcoef(0:MAXSTEPS), vcoef(0:MAXSTEPS)
      common/rtabs/ xcoef, vcoef

      tpower = dt*dt
      xfact(0)=tpower*xcoef(0)
      vfact(0)=dt
      do 20 i=1,order
         vfact(i) = tpower*vcoef(i)
         tpower = tpower * dt
         xfact(i) = tpower*xcoef(i)
20    continue
      do 30 k = 1, NDIM
         vnew(k)=vold(k)
         xnew(k)=xold(k) + dt*vold(k)
         do 40 i=0,order
            xnew(k) = xnew(k) + pcoef(k,i)*xfact(i)
            vnew(k) = vnew(k) + pcoef(k,i)*vfact(i)
40       continue
30    continue
      end
c
      subroutine setpre(pcoef, order, ibody, mxptcl, mxord, 
     $                  xpcoef, vpcoef)
c-----------------------------------------------------------------------      
c  setpre : calculate the coefficient of predictor
c-----------------------------------------------------------------------
c  pcoef	in	power series coefficients of acc
c  order	in	order of pcoef
c  ibody	in	particle ID
c  xpcoef	out	power series for x (start from acc)
c  vpcoef	out	power series for v (start from acc)
c-----------------------------------------------------------------------
c  Algorithm:
c   xpcoef(ibody, k,i) = pcoef(k,i)/((i+1)*(i+2))
c   vpcoef(ibody, k,i) = pcoef(k,i)/(i+1)
c-----------------------------------------------------------------------      

      REAL pcoef(NDIM, 0:MAXSTEPS)
      INTEGER order, mxptcl, mxord, ibody
      REAL xpcoef(mxptcl, NDIM, 0:mxord)
      REAL vpcoef(mxptcl, NDIM, 0:mxord)

      REAL xcoef(0:MAXSTEPS), vcoef(0:MAXSTEPS)
      common/rtabs/ xcoef, vcoef

      INTEGER i,k
      
      do 20 k=1,NDIM
         do 10 i=0,order
            xpcoef(ibody,k,i)=pcoef(k,i)*xcoef(i)
            vpcoef(ibody,k,i)=pcoef(k,i)*vcoef(i)
10       continue
20    continue
      end
