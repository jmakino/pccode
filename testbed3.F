#include "pcheader.h"

      program testmn
c-----------------------------------------------------------------------
c testmn : test main program to solve
c kepler one-body problem
c
c Version 0.1 : Both of Lagrange and Hermite interpolation schemes are
c		implemented.  
c		89-03-20
c-----------------------------------------------------------------------

      REAL t(0:MAXSTEPS), tnew(0:MAXSTEPS)
      REAL coef(MAXSTEPS,MAXSTEPS)
      REAL pcoef(NDIM, 0:MAXSTEPS)
      INTEGER order, ordmax
      INTEGER i, k, l, istep, nstep
      REAL d(NDIM, 0:MAXSTEPS), x(NDIM), v(NDIM), a(NDIM)
      REAL aold(NDIM, 0:MAXSTEPS), dtmp(NDIM,0:MAXSTEPS)
      REAL dt, dtl, tnow, xp(NDIM), vp(NDIM), dtlist(1:MAXSTEPS)
      REAL xc(NDIM), vc(NDIM)
      INTEGER dttype, gmode, ordflg, herflg, nrot
      REAL dtnew, eta, tolva, tolvr, tol1, tol2, exfact, dtmax
      REAL adot(NDIM), epsini, ecc, lz0, r2, rfact, tinit
      REAL xq(NDIM), vq(NDIM)
      REAL eerr, lerr
      call strtab
      write(6,*)'Enter order, hermit_flag, nrot, dtmax:'
      read(5,*) ordmax, herflg, nrot, dtmax
      write(6,612) ordmax,herflg,nrot,dtmax
612   format(' Order, Hermite, Nrot, dt :', 3i3, g15.5)
      write(6,*) 'Enter dttype(1-aarseth,   2-abs, lmode, 3-abs, gmode'
      write(6,*)'              4-rel, lmode,5-rel, gmode'
      read(5,*)dttype, tol1, tol2, exfact
      write(6,610)'dt expand limit = ', exfact
      if(dttype .eq. 1) then
         write(6,610) ' Aarseth dt, eta=', tol1
610      format(a, g15.7, a, i5)         
         eta = tol1
      else if (dttype .le. 3) then
         tolva = tol1
         gmode = dttype - 2
         write(6,610)' Vabs dt, tol = ', tolva, ' Mode=', gmode
      else if (dttype. le. 5) then
         tolva = tol1
         tolvr = tol2
         gmode = dttype - 4
         write(6,611)' Vrel dt, tol = ', tolva,tolvr, ' Mode=', gmode
611      format(a, 2g15.7, a, i5) 
      endif
c-----------------------------------------------------------------------      
c     set initial condition 
c-----------------------------------------------------------------------
      write(6,*)' Enter eccentricity:'
      read(5,*)ecc
      write(6,*)' eccentricity=',ecc
      lz0 = sqrt(1.0d0-ecc*ecc)
      v(2) = lz0
      v(1) = + ecc
      x(1)=1d0
      x(2)=0d0
      x(3)=0d0
      v(3)=0d0
      t(0)=0d0
600   format(' ',a, i5, (3g21.14))      

      dtl = 0.0

      do 20 k=1,NDIM
         xq(k) = 0d0
         vq(k) = 0d0
         aold(k,0)=-x(k)
         adot(k)=-v(k) -  3.0d0*x(1)*v(1)*aold(k,0)
20    continue

      call fset(x,v,a, adot)

      if(dttype .eq. 0) call terror('const dt cannot self start')

      write(6,*)'Enter epsinit:'
      read(5,*) epsini
      dt = epsini/sqrt(adot(1)**2+adot(2)**2+adot(3)**2)
#ifdef FACTOR2DT
      dtnew=1.0d0
900   if(dtnew .gt. dt) then
         dtnew=dtnew*0.5
         goto 900
      endif
      dt = dtnew
#endif      
      if(herflg .eq. 0) then
         order = 0
         ordflg = 0
         d(1,0)=aold(1,0)
         d(2,0)=aold(2,0)
         d(3,0)=aold(3,0)
         do 200 k = 1, 3
           d(k,1) = adot(k)
200      continue
         dtlist(1) = -dt
      else
         order=1
         ordflg=0
         tinit = t(0)
         call upddth(d, t, tinit, a, adot, order, dtlist)
      endif
      call dtopow(d, dtlist, order, pcoef)
      call predic(pcoef, dt, x, v, order, xp, vp)
      write(6,600)'xp ', 0, (xp(k), k=1,NDIM)
      write(6,600)'vp ', 0, (vp(k), k=1,NDIM)
c-----------------------------------------------------------------------
c   push the particle several steps
c-----------------------------------------------------------------------
      do 110 istep=1, 1000000
c---     NEW ACCELERATION      
         call fset(xp,vp,a,adot)
c---     CONTROL THE ORDER OF INTEGRATOR
         call adjord(ordmax, ordflg, herflg, order)
         if(herflg. eq. 0) then
c---        UPDATE DIVIDED DIFFERNCE  AND TIME LIST
            call upddt(d, t, t(0)+dt, a, order, dtlist)
c---        OBTAIN DERIVATIVES
            call dtopow(d, dtlist, order, pcoef)
c---        APPLY CORRECTOR
#ifdef ROUND_OFF_COMP
            call correb(pcoef, dtlist(1), order, x, v, xq, vq)
#else
            call correc(pcoef, dtlist(1), order, x, v)
#endif
         else
            call upddth(d, t, t(0)+dt, a, adot, order, dtlist)
            call dtopow(d, dtlist, order, pcoef)
#ifdef ROUND_OFF_COMP
            call correb(pcoef, dtlist(2), order, x, v, xq, vq)
#else
            call correc(pcoef, dtlist(2), order, x, v)
#endif
         endif
#ifdef STEPOUT
         write(6,600)'x ', 0, (x(k), k=1,NDIM)
         write(6,600)'v ', 0, (v(k), k=1,NDIM)
         call checke(x, v, -0.5d0 ,lz0, 0, 1, eerr, lerr)
#endif
#ifdef STEPOUT
         do 1010 i=0,order
            write(6,600)'pcoef ', i, (pcoef(k,i), k=1,NDIM)
1010     continue 
#endif
c-----------------------------------------------------------------------
c     determine new DT
c-----------------------------------------------------------------------
         if(dttype .gt. 0 .and. order .gt. ordmax - 3) then
            call dtgen(pcoef, eta, vp,v, a, tolva, tolvr, 
     $                 dt, order, exfact, 
     $                 (dttype+2)/2, mod(dttype,2), dtnew)
#ifdef STEPOUT
            write(6,600) 'dtold, new ', 0, dt, dtnew
#endif

#ifdef FACTOR2DT
            if(dtnew .ge. dt*2.0d0) then
               dtnew = 2.0d0*dt
            else if(dtnew .lt. dt) then
               dtnew = 0.5d0*dt
            else
               dtnew = dt
            endif
#endif            
            dt = min(dtmax, dtnew)
         endif

c---     OBTAIN NEXT PREDECTED POSITION
         call predic(pcoef, dt, x, v, order, xp, vp)
#ifdef STEPOUT
         write(6,600)'steps, tnow : ', istep,t(0)
#else
         if(mod(istep,50).eq.0) then
            write(6,600)'steps, tnow : ', istep,t(0)
            write(6,600)'x ', 0, (x(k), k=1,NDIM)
            write(6,600)'v ', 0, (v(k), k=1,NDIM)
            call checke(x, v, -0.5d0, lz0, 0, 1, eerr, lerr)
         else
            call checke(x, v, -0.5d0, lz0, 0, 0, eerr, lerr)
         endif
#endif
#ifdef STEPOUT
         write(6,600)'vp ', 0, (vp(k), k=1,NDIM)

#endif
601      format(4g18.9)               
         nstep = istep
         if(0.5*(t(0)+t(1)) .gt. 3.141592*2*nrot) goto 999
110   continue                  
999   write(6,600)'steps, tnow : ', nstep,t(0)
      call checke(x, v, -0.5d0, lz0, 1, 1, eerr, lerr)
      open(8,file='testbed.out')
      write(8,613)ordmax, herflg, nrot, dtmax, dttype, tol1, tol2,
     $            exfact, epsini, ecc, nstep, eerr, lerr
613   format(3i3,1x,1pg12.4,1x,i2,5g15.7,1x,i6,2g15.7)
      close(8)
      end
c
      subroutine checke(x, v, e0, lz0, maxout, preout, eerr, lerr)
      REAL x(NDIM), v(NDIM), e0, lz0, eerr, lerr
      INTEGER maxout, preout
      REAL e, lz
      REAL eemax, elmax
      data eemax,elmax /0.0,0.0/
      e= 0.5*(v(1)**2+v(2)**2+v(3)**2)
     $    -1.0d0/sqrt(x(1)**2+x(2)**2+x(3)**2)
      lz = v(2)*x(1)-v(1)*x(2)
      eemax = max(eemax, abs(e-e0))
      eerr = eemax
      elmax = max(elmax, abs(lz-lz0))
      lerr = elmax
      if (preout .ne. 0) then
         write(6,600) e, e-e0, lz, lz-lz0
600      format(' E, Eerr, L, Lerr : ', 1p4g15.7)
      endif
      if(maxout .ne. 0) then
          write(6,601) eemax, elmax
601       format(' EerrMAX, LerrMAX = ', 1p2g15.7)
      endif
      end
c
      subroutine fset(x, v, a, adot)
c-----------------------------------------------------------------------      
c  fset : calculate force amd fdot
c-----------------------------------------------------------------------
c  x		in	position
c  v		in	velocity
c  a		out	acceleration
c  adot		out	da/dt
c-----------------------------------------------------------------------      
      REAL x(NDIM), v(NDIM), a(NDIM), adot(NDIM)
      REAL r2, r3, afact
      INTEGER k
      r2 = x(1)**2 + x(2)**2 + x(3)**2
      r3 = sqrt(r2)*r2
      afact = (x(1)*v(1)+x(2)*v(2)+x(3)*v(3))/r2
      do 10 k = 1, NDIM
         a(k) = -x(k)/r3
         adot(k) = -v(k)/r3 - 3.0d0*afact*a(k)
10    continue
      end
