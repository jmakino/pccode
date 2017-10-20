#include "pcheader.h"

      program testmn
c-----------------------------------------------------------------------
c testmn : test main program to solve
c y'' = alpha(alpha-1)y^{(alpha-2)/alpha}
c
c this equation gives the solution : y = x^alpha
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
      INTEGER dttype, gmode, ordflg
      REAL dtnew, eta, tolva, tolvr, tol1, tol2, exfact
      REAL adot(NDIM), epsini, alpha, xtrue, vtrue
      call strtab
      write(6,*)'Enter order, dt:'
      read(5,*) ordmax, dt
      write(6,*) ' Order = ', ordmax, '  DT =',dt
      write(6,*)'Enter alpha:'
      read(5,*) alpha
      write(6,*)' Alpha = ', alpha
#ifdef VARIABLEDT
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
#endif
c-----------------------------------------------------------------------      
c     set initial condition 
c-----------------------------------------------------------------------
      x(1)=1d0
      v(1)=alpha
      x(2)=1d0
      v(2)=alpha
      x(3)=cos(1d0)
      v(3)=-sin(1d0)

600   format(' ',a, i5, (3g21.14))      
      order = ordmax
      t(0)=1d0
      do 20 i=0, order
         dtl = i * dt
         t(i)=1.0d0-dtl
         aold(1,i)=alpha*(alpha-1d0)*(t(i)**(alpha-2d0))
         aold(2,i)= aold(1,i)
         aold(3,i)=-cos(t(i))

         write(6,600)'aold', i, (aold(k,i), k=1,NDIM), t(i)
20    continue         

c-----------------------------------------------------------------------      
c     form divided difference
c    
c     
c-----------------------------------------------------------------------
c
      do 30 i=order,0,-1
c---     SET D_0 term      
         do 40 k=1,NDIM
            d(k,0) = aold(k,i)
40       continue
c---     SET D_1 through D_order by bootstrapping
         do 50 l = 1,order-i
            dtl = t(i)-t(i+l)
            do 60 k = 1, NDIM
               d(k,l) = (d(k,l-1)-dtmp(k,l-1))/dtl
60          continue
50       continue
c         write(6,600) 'd init ', i, (d(1,l),l=0,order-i)
         do 70 l=0,order-i
            do 80 k=1,NDIM
               dtmp(k,l)=d(k,l)
80          continue
70       continue
30    continue
c---
      do 100 i=0,order
         write(6,600)'d ', i, (d(k,i), k=1,NDIM)
100   continue 
c-----------------------------------------------------------------------
c   obtain next predicted position
c-----------------------------------------------------------------------       
      do 102 i = 1, order
         dtlist(i) = t(i)-t(0)
102   continue
      call dtopow(d, dtlist, order, pcoef)
#ifdef STEPOUT
       do 1000 i=0,order
         write(6,600)'pcoef ', i, (pcoef(k,i), k=1,NDIM)
1000   continue 
#endif
      call predic(pcoef, dt, x, v, order, xp, vp)
      write(6,600)'xp ', 0, (xp(k), k=1,NDIM)
      write(6,600)'vp ', 0, (vp(k), k=1,NDIM)
c-----------------------------------------------------------------------
c   push the particle several steps
c-----------------------------------------------------------------------
      do 110 istep=1, 100000000
c---     NEW ACCELERATION      
         a(1) = alpha*(alpha-1d0)*(xp(1)**((alpha-2d0)/alpha))
         a(2) = a(1)
         a(3) = -xp(3)
#ifdef STEPOUT
         write(6,600)'a ', 0, a
#endif
#ifdef SELF_START
c---     CONTROL THE ORDER OF INTEGRATOR
         call adjord(ordmax, ordflg, order)
#endif
c---     UPDATE DIVIDED DIFFERNCE  AND TIME LIST
         call upddt(d, t, t(0)+dt, a, order, dtlist)
c---     OBTAIN DERIVATIVES
         call dtopow(d, dtlist, order, pcoef)
c---     APPLY CORRECTOR
         call correc(pcoef, dtlist(1), order, x, v)
#ifdef STEPOUT
         write(6,600)'v ', 0, (v(k), k=1,NDIM)
#endif
#ifdef STEPOUT
         do 1010 i=0,order
            write(6,600)'pcoef ', i, (pcoef(k,i), k=1,NDIM)
1010     continue 
#endif
#ifdef VARIABLEDT
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
            dt = dtnew
         endif
#endif
c---     OBTAIN NEXT PREDECTED POSITION
         call predic(pcoef, dt, x, v, order, xp, vp)

         xtrue = t(0)**alpha
         vtrue = alpha*t(0)**(alpha-1d0)
#ifdef STEPOUT
         write(6,600)'steps, tnow : ', istep,t(0)
         write(6,601) x(1), xtrue, v(1), vtrue
         write(6,601) x(1)-xtrue, v(1)-vtrue
#else
         if(mod(istep,int(5.0/dt/10.0)).eq.0) then
            write(6,600)'steps, tnow : ', istep,t(0)
            write(6,601) x(1), v(1),  x(1)-xtrue, v(1)-vtrue
            write(6,601) x(3), v(3), x(3)-cos(t(0)), v(3)+sin(t(0))
         endif
#endif
601      format(1p4g18.9)               
         nstep = istep
         if(t(0) .gt. 5) goto 999
110   continue                  
999   write(6,600)'steps, tnow : ', nstep,t(0)
      write(6,601) x(1), v(1),  x(1)-xtrue, v(1)-vtrue
      write(6,601) x(3), v(3), x(3)-cos(t(0)), v(3)+sin(t(0))
      end

