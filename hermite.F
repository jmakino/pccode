#include "pcheader.h"
c-----------------------------------------------------------------------
c  hermite.f : worker routines for hermite polynominal
c              integrator scheme
c  J.Makino
c  Ver 0.00 89-03-17
c-----------------------------------------------------------------------
c
      subroutine upddth(ddt, t, t0, acc0, adot0, order, dtlist)
c-----------------------------------------------------------------------
c   upddth : update divided difference table and time list
c            for hermite interpolation
c-----------------------------------------------------------------------
c   ddt		in/out	divided difference table
c   t		in/out	time list
c   t0		in 	new reference time
c   acc0	in	acceleration at t=t0
c   adot0	in	da/dt
c   order	in	order of polynominal
c   dtlist      out     list of (t_k - t_0)
c   
c   actual work of updating D is done in UPDDT
c-----------------------------------------------------------------------
      REAL ddt(NDIM,0:MAXSTEPS),t(0:MAXSTEPS), dtlist(MAXSTEPS)
      REAL t0, acc0(NDIM), adot0(NDIM)
      INTEGER order
      REAL dnew(NDIM,0:MAXSTEPS),tnew(0:MAXSTEPS)
      INTEGER i,k
      call upddtw(ddt, t, t0, acc0, adot0, order, dnew, tnew)
      call upddt2(dnew, tnew, order, ddt, t, dtlist)
      end
c      
      subroutine upddtw(dold, told, t0, acc0, adot0, order, dnew, tnew)
c-----------------------------------------------------------------------
c   upddt0 : update divided difference table and time list
c            for the hermite interpolation
c-----------------------------------------------------------------------
c   dold	in	old divided difference table
c   told	in	old time list
c   t0		in 	new reference time
c   acc0	in	acceleration at t=t0
c   adot0       in      da/dt at t=t0
c   order	in	order of polynominal
c   dnew	out	new divided difference table
c   tnew	out	new time list
c
c   D is updated using the following relation:
c   t1 = t0
c   dnew_0 = acc0
c   dnew_1 = adot0
c   otherwise use the following relation
c   dmed_k = (dmed_k-1 - dold_k-1)/(t_1-t_k) (k=1,order)
c   dnew_k = (dnew_k-1 - dmed_k-1)/(t_0-t_k) (k=2,order)
c-----------------------------------------------------------------------
      REAL dold(NDIM,0:MAXSTEPS),told(0:MAXSTEPS)
      REAL t0, acc0(NDIM), adot0(NDIM)
      INTEGER order
      REAL dnew(NDIM,0:MAXSTEPS),tnew(0:MAXSTEPS)
      REAL dmed(NDIM,0:MAXSTEPS)
      INTEGER i,k, im1, im2
      REAL dtinv
c-----------------------------------------------------------------------      
c     set  Dmed0, Dnew0, Dnew1
c-----------------------------------------------------------------------      
      do 10 k=1,NDIM
         dnew(k,0) = acc0(k)
         dmed(k,0) = acc0(k)
         dnew(k,1) = adot0(k)
10    continue
c-----------------------------------------------------------------------      
c     calculate Dmed using the iteration formula
c-----------------------------------------------------------------------      
      do 20 i = 1, order-1
         im1=i-1
         dtinv = 1.0d0/(t0 - told(im1))
         do 30 k = 1, NDIM
            dmed(k,i) = (dmed(k,im1)-dold(k,im1))*dtinv
30       continue
20    continue
c-----------------------------------------------------------------------      
c     set t0, t1
c-----------------------------------------------------------------------      
      tnew(0)=t0
      tnew(1)=t0
c-----------------------------------------------------------------------      
c     calculate Dnew using the iteration formula
c-----------------------------------------------------------------------      
      do 40 i = 2, order
         im1=i-1
         im2=i-2
         tnew(i) = told(im2)
         dtinv = 1.0d0/(t0 - told(im2))
         do 50 k = 1, NDIM
            dnew(k,i) = (dnew(k,im1)-dmed(k,im1))*dtinv
50       continue
40    continue
      end
c
      subroutine fupddh(ddt, t, t0, acc0, adot0, order, ndim, dtlist)
c-----------------------------------------------------------------------
c   upddth : update divided difference table and time list
c            for hermite interpolation
c-----------------------------------------------------------------------
c   ddt		in/out	divided difference table
c   t		in/out	time list
c   t0		in 	new reference time
c   acc0	in	acceleration at t=t0
c   adot0	in	da/dt
c   order	in	order of polynominal
c   dtlist      out     list of (t_k - t_0)
c   
c   actual work of updating D is done in UPDDT
c-----------------------------------------------------------------------
      INTEGER ndim
      REAL ddt(ndim,0:MAXSTEPS),t(0:MAXSTEPS), dtlist(MAXSTEPS)
      REAL t0, acc0(ndim), adot0(ndim)
      INTEGER order
      REAL dnew(NDIM*2,0:MAXSTEPS),tnew(0:MAXSTEPS)
      INTEGER i,k
      call fupddw(ddt, t, t0, acc0, adot0, order, ndim, dnew, tnew)
      call fupdd2(dnew, tnew, order, ndim, ddt, t, dtlist)
      end
c      
      subroutine fupddw(dold, told, t0, acc0, adot0, order, ndim, 
     $                  dnew, tnew)
c-----------------------------------------------------------------------
c   upddt0 : update divided difference table and time list
c            for the hermite interpolation
c-----------------------------------------------------------------------
c   dold	in	old divided difference table
c   told	in	old time list
c   t0		in 	new reference time
c   acc0	in	acceleration at t=t0
c   adot0       in      da/dt at t=t0
c   order	in	order of polynominal
c   dnew	out	new divided difference table
c   tnew	out	new time list
c
c   D is updated using the following relation:
c   t1 = t0
c   dnew_0 = acc0
c   dnew_1 = adot0
c   otherwise use the following relation
c   dmed_k = (dmed_k-1 - dold_k-1)/(t_1-t_k) (k=1,order)
c   dnew_k = (dnew_k-1 - dmed_k-1)/(t_0-t_k) (k=2,order)
c-----------------------------------------------------------------------
      INTEGER ndim
      REAL dold(ndim,0:MAXSTEPS),told(0:MAXSTEPS)
      REAL t0, acc0(ndim), adot0(ndim)
      INTEGER order
      REAL dnew(ndim,0:MAXSTEPS),tnew(0:MAXSTEPS)
      REAL dmed(NDIM*2,0:MAXSTEPS)
      INTEGER i,k, im1, im2
      REAL dtinv
c-----------------------------------------------------------------------      
c     set  Dmed0, Dnew0, Dnew1
c-----------------------------------------------------------------------      
      do 10 k=1,ndim
         dnew(k,0) = acc0(k)
         dmed(k,0) = acc0(k)
         dnew(k,1) = adot0(k)
10    continue
c-----------------------------------------------------------------------      
c     calculate Dmed using the iteration formula
c-----------------------------------------------------------------------      
      do 20 i = 1, order-1
         im1=i-1
         dtinv = 1.0d0/(t0 - told(im1))
         do 30 k = 1, ndim
            dmed(k,i) = (dmed(k,im1)-dold(k,im1))*dtinv
30       continue
20    continue
c-----------------------------------------------------------------------      
c     set t0, t1
c-----------------------------------------------------------------------      
      tnew(0)=t0
      tnew(1)=t0
c-----------------------------------------------------------------------      
c     calculate Dnew using the iteration formula
c-----------------------------------------------------------------------      
      do 40 i = 2, order
         im1=i-1
         im2=i-2
         tnew(i) = told(im2)
         dtinv = 1.0d0/(t0 - told(im2))
         do 50 k = 1, ndim
            dnew(k,i) = (dnew(k,im1)-dmed(k,im1))*dtinv
50       continue
40    continue
      end
#ifdef TEST
      subroutine ftest(t,f,fdot)
c-----------------------------------------------------------------------
c  ftest : supply trial function for debug of hermite interpolation
c-----------------------------------------------------------------------
c  t		in	Independent variable
c  f		out	Value of function 
c  fdot		out	Value of df/dt
c-----------------------------------------------------------------------      
      REAL t, f(NDIM), fdot(NDIM)
      f(1) = t**2
      f(2) = t**3
      f(3) = t**5
      fdot(1) = 2d0*t
      fdot(2) = 3d0*t**2
      fdot(3) = 5d0*t**4
      end
      program tsther
c-----------------------------------------------------------------------      
c tsther : DEBUG driver for hermite interporation module
c-----------------------------------------------------------------------      
      REAL d(NDIM, 0:MAXSTEPS), t(0:MAXSTEPS), dt(MAXSTEPS)
      REAL a(NDIM), adot(NDIM)
      REAL t0
      INTEGER i, j, k
      t0 = 0d0
      do 10 i = 1, 7, 2
         call ftest(t0, a, adot)
         call upddth(d, t, t0, a, adot, i, dt)
         write(6,600)' a:   ', a
         write(6,600)' adot:', adot
600      format(a5, 3g18.9)
         do 20 j=0,i
            write(6,601)j,t(j),dt(j),(d(k,j),k=1,NDIM)
601         format(i3, 5g15.7)
20       continue
         t0 = t0 + 1
10    continue
      end
#endif         
