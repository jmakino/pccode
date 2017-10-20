/*#define STEPOUT*/
#include "pcheader.h"

      program testmn
c-----------------------------------------------------------------------
c testmn : test main program to solve
c kepler one-body problem
c
c Version 0.1 : Both of Lagrange and Hermite interpolation schemes are
c		implemented.  
c		89-03-20
c Version 0.1a: test main for Symmetric schemes (first order solvers)
c-----------------------------------------------------------------------
      INTEGER neqs
      PARAMETER (neqs = NDIM*2)
      REAL t(0:MAXSTEPS), tnew(0:MAXSTEPS)
      REAL coef(MAXSTEPS,MAXSTEPS)
      REAL pcoef(neqs, 0:MAXSTEPS)
      INTEGER order, ordmax
      INTEGER i, k, l, istep, nstep
      REAL d(neqs, 0:MAXSTEPS), x(NDIM), v(NDIM), a(NDIM)
      REAL dt, dtl, tnow, xp(NDIM), vp(NDIM), dtlist(1:MAXSTEPS)
      REAL dnew(neqs, 0:MAXSTEPS), dtlnew(1:MAXSTEPS)
      REAL xc(NDIM), vc(NDIM)
      REAL xa(neqs), xap(neqs), xac(neqs), aa(neqs), adota(neqs)
      INTEGER dttype, gmode, ordflg, herflg, nrot
      REAL dtnew, eta, tolva, tolvr, tol1, tol2, exfact, dtmax
      REAL adot(NDIM), epsini, ecc, lz0, r2, rfact, tinit
      REAL eerr, lerr
      INTEGER imode
      EQUIVALENCE (a,aa)
      EQUIVALENCE (adot,adota)
      EQUIVALENCE (x,xa)
      EQUIVALENCE (v, xa(4))
      EQUIVALENCE (vp, xap(4))
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
      do 1010 k = 1, NDIM
        xa(k) = x(k)
        xa(k+NDIM)=v(k)
1010  continue        
      t(0)=0d0
600   format(' ',a, i5, (3g21.14))      

      dtl = 0.0

      call ffset(xa,a, adot)

      if(dttype .eq. 0) call terror('const dt cannot self start')

      write(6,*)'Enter epsinit:'
      read(5,*) epsini
      dt = epsini/sqrt(adot(1)**2+adot(2)**2+adot(3)**2)
      if(herflg .eq. 0) then
         order = 0
         ordflg = 0
         do 1020 k =1, NDIM
            d(k,0)=v(k)
            d(k+NDIM, 0) = a(k)
1020     continue
         dtlist(1) = -dt
      else
         order=1
         ordflg=0
         tinit = t(0)
         call fupddh(d, t, tinit, a, adot, order, neqs, dtlist)
      endif
      call fdtopo(d, dtlist, order, neqs, pcoef)
      call fpredi(pcoef, dt, xa, order, neqs, xap)
      write(6,600)'xp ', 0, (xap(k), k=1,NDIM)
      write(6,600)'vp ', 0, (xap(k+NDIM), k=1,NDIM)
c-----------------------------------------------------------------------
c   push the particle several steps
c-----------------------------------------------------------------------
      do 110 istep=1, 10000
c---     TEST P(EC)^3 MODE      
c---        INITIAL GUESS FOR CORRECTOR : PREDICTOR
         do 1060 k = 1, neqs
            xac(k) = xap(k)
1060     continue
         do 1030 imode=1,10
c---        NEW ACCELERATION      
            call ffset(xac,a,adot)
#ifdef STEPOUT
            write(6,600)'xc ', 0, (xac(k), k=1,NDIM)
            write(6,600)'vc ', 0, (xac(k+NDIM), k=1,NDIM)
            write(6,600)'a ', 0, (a(k), k=1,NDIM)
            write(6,600)'da ', 0, (adot(k), k=1,NDIM)
#endif
c---        RESTORE OLD POSITION into CORRECTOR            
            do 1040 k = 1, neqs
               xac(k) = xa(k)
1040        continue
c---        CONTROL THE ORDER OF INTEGRATOR
            if(imode .eq. 1) then
               call adjord(ordmax, ordflg, herflg, order)
            endif
            if(herflg. eq. 0) then
c---           UPDATE DIVIDED DIFFERNCE  AND TIME LIST
               call fupdd0(d, t, t(0)+dt, a, order, neqs, dnew, tnew)
               call setdtl(tnew, dtlnew, order)
c---           OBTAIN DERIVATIVES
               call fdtopo(dnew, dtlnew, order, neqs,  pcoef)
c---           APPLY CORRECTOR
               call fcorre(pcoef, dtlnew(1), order, neqs, xac)
            else
               call fupddw(d, t, t(0)+dt, a, adot, order, neqs, dnew, 
     $                     tnew)
               call setdtl(tnew, dtlnew, order)
               call fdtopo(dnew, dtlnew, order, neqs, pcoef)
               call fcorre(pcoef, dtlnew(2), order, neqs, xac)
            endif
1030     continue            
c---     UPDATE DDT, POSITION, TIMELIST
         call fupdd2(dnew, tnew, order, neqs, d, t, dtlist)
         do 1050 k = 1, neqs
            xa(k) = xac(k)
1050     continue
#ifdef STEPOUT
         write(6,600)'x ', 0, (xa(k), k=1,NDIM)
         write(6,600)'v ', 0, (xa(k+NDIM), k=1,NDIM)
         call checke(x, x(NDIM+1), -0.5d0 ,lz0, 0, 1, eerr, lerr)
#endif
#ifdef STEPOUT
         do 200 i=0,order
            write(6,600)'pcoef ', i, (pcoef(k,i), k=1,neqs)
200      continue 
#endif
c-----------------------------------------------------------------------
c     determine new DT (NOW CONSTANT dt only)
c-----------------------------------------------------------------------
         if(dttype .gt. 0 .and. order .gt. ordmax - 3) then
            call dtgen(pcoef, eta, vp,v, a, tolva, tolvr, 
     $                 dt, order, exfact, 
     $                 (dttype+2)/2, mod(dttype,2), dtnew)
#ifdef STEPOUT
            write(6,600) 'dtold, new ', 0, dt, dtnew
#endif
#define FACTOR2DT
#ifdef FACTOR2DT
            if(dtnew .ge. dt*1.999) then
               dtnew = 2.0d0*dt
            else if(dtnew .lt. 0.999*dt) then
               dtnew = 0.5d0*dt
            endif
#endif            
            dt = min(dtmax, dtnew)
         endif

c---     OBTAIN NEXT PREDECTED POSITION
         call fpredi(pcoef, dt, xa, order, neqs, xap)
#ifdef STEPOUT
         write(6,600)'steps, tnow : ', istep,t(0)
#else
         if(mod(istep,50).eq.0) then
            write(6,600)'steps, tnow : ', istep,t(0)
            write(6,600)'x ', 0, (xa(k), k=1,NDIM)
            write(6,600)'v ', 0, (xa(k+NDIM), k=1,NDIM)
            call checke(xa, xa(4), -0.5d0, lz0, 0, 1, eerr, lerr)
         else
            call checke(xa, xa(4), -0.5d0, lz0, 0, 0, eerr, lerr)
         endif
#endif
601      format(4g18.9)               
         nstep = istep
         if(0.5*(t(0)+t(1)) .gt. 3.141592*2*nrot) goto 999
110   continue                  
999   write(6,600)'steps, tnow : ', nstep,t(0)
      call checke(xa, xa(4), -0.5d0, lz0, 1, 1, eerr, lerr)
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
      subroutine ffset(xv, f, fdot)
c-----------------------------------------------------------------------      
c  fset : calculate force amd fdot
c-----------------------------------------------------------------------
c  xv		in	position & velocity
c  f		out	dxv/dt
c  fdot		out	da/dt
c-----------------------------------------------------------------------      
      REAL xv(6), f(6), fdot(6)
      REAL x(NDIM), v(NDIM), a(NDIM), adot(NDIM)
      REAL r2, r3, afact
      INTEGER k
      do 30 k = 1, NDIM
         x(k)=xv(k)
         v(k)=xv(k+NDIM)
30    continue
      r2 = x(1)**2 + x(2)**2 + x(3)**2
      r3 = sqrt(r2)*r2
      afact = (x(1)*v(1)+x(2)*v(2)+x(3)*v(3))/r2
      do 10 k = 1, NDIM
         a(k) = -x(k)/r3
         adot(k) = -v(k)/r3 - 3.0d0*afact*a(k)
10    continue
      do 20 k=1,NDIM
         f(k)=v(k)
         f(k+NDIM)=a(k)
         fdot(k)=a(k)
         fdot(k+NDIM)=adot(k)
20    continue
      end
