#define DBGWRITE(x) WRITE(6,*)'x :', x
#include "pcheader.h"
c-----------------------------------------------------------------------
c   intgrt2.F : integrator driver/utility routines
c
c   Modified 89-03-18 : herflg added to args of adjord
c            89-03-27 : N-body driver added
c-----------------------------------------------------------------------
c
      subroutine adjord(ordmax, ordflg, herflg, order)
c-----------------------------------------------------------------------
c   adjord : adjust the order of the integrator
c-----------------------------------------------------------------------
c   ordmax	in	maximum order specified
c   ordflg	in/out	flag to control the increment of order
c   herflg	in	flag for Hermite interpolation (0:Lagrange, 
c							1:Hermite)
c   order	in/out	current/new order
c
c   Strategy to change the order:
c   After the order is incremented to p, it will not be changed for next
c   p steps to prevent the integrator becomming unstable
c
c   Thus, after the order is incremented, the ordflg is set to order.
c   At each call of adjord, ordflg is decremented. 
c   when ordflg is 0, order is allowed to change
c
c   89-03-18 herflg introduced
c   if herflg <> 0, increment order by 2 (not 1)
c-----------------------------------------------------------------------
      INTEGER ordmax, ordflg, order, herflg
      if(order .lt. ordmax) then
         if(ordflg .ne. 0) then
            ordflg = ordflg - 1
         else
            if(herflg .eq. 0) then 
               order = order + 1
            else
               order = order + 2
            endif
            ordflg = order
         endif
      endif
      end
c
      subroutine selnxt(tnext, nbody, inext)
c-----------------------------------------------------------------------      
c  selnxt : Select Next Particle to Push
c-----------------------------------------------------------------------
c   tnext	in	list of next time of particles
c   nbody	in	number of particles
c   inext	out	particle to move next
c
c   Return the particle with smallest tnext(i)
c
c   This version is impremented as simple scan-all (O(N)) routine
c-----------------------------------------------------------------------      
      REAL tnext(MAXBODY)
      INTEGER nbody, inext
      REAL tmin
      INTEGER i, imin
      tmin = tnext(1)
      imin =1
      do 10 i = 2, NBODY
         if(tmin .gt. tnext(i)) then
            imin = i
            tmin = tnext(i)
         endif
10    continue
      inext = imin
      end
c
      subroutine accel(xdata, vdata, xpre, vpre, mass, tnow, 
     $                    order, nbody, posi, veli, ti, i, eps2, 
     $                    acci, adoti)
c-----------------------------------------------------------------------      
c  accel : calculate the acceleration and da/dt of particle i
c-----------------------------------------------------------------------
c   xdata	in	position of particles
c   vdata	in	velocity of particles
c   xpre	in	predictor polynominal for position (order of
c                       acceleration and upper)
c   vpre	in	predictor polynominal for velocity (order of
c                       acceleration and upper)
c   mass 	in	mass of particles
c   tnow	in	time of particles
c   order	in	order of predictor polynominal
c   nbody	in	number of particles
c   posi	in	position of particle i (predicted)
c   veli	in	velocity of particle i (predicted)
c   ti		in	time of particle i     (for predicted position)
c   i		in	index of particle to calculate the force on it
c   eps2 	in	softening parameter squared
c   acci	out	acceleration
c   adoti	out	da/dt
c
c   This routine calculate BOTH a and adot.
c   The routine that calculate a only is provided with different name
c-----------------------------------------------------------------------      
      REAL xdata(MAXBODY, NDIM), vdata(MAXBODY, NDIM)
      REAL xpre(MAXBODY, NDIM, 0:MAXSTEPS)
      REAL vpre(MAXBODY, NDIM, 0:MAXSTEPS)
      REAL mass(MAXBODY), tnow(MAXBODY)
      INTEGER order, nbody
      REAL posi(NDIM), veli(NDIM), ti
      INTEGER i
      REAL eps2, acci(NDIM), adoti(NDIM)
      REAL dx(NDIM), posj(NDIM), velj(NDIM), dv(NDIM)
      REAL rfact, dtj, r2, r2inv, rdfact
      INTEGER j,k,l
#ifdef VECTOR
      REAL dxv(MAXBODY,NDIM), posjv(MAXBODY,NDIM)
      REAL veljv(MAXBODY,NDIM), dvv(MAXBODY,NDIM)
      REAL dtjv(MAXBODY)
      INTEGER j1,j2
#   ifdef FAST_SQRT
      REAL r, rinv, eps
#   endif            
#endif
      if(NDIM .ne. 3) call terror(' NDIM <>3 unsupported in sumfrc')      
      do 10 k = 1, NDIM
         acci(k) = 0.0d0
         adoti(k) = 0.0d0
10    continue
#ifndef VECTOR
      do 20 j = 1, nbody
         if(j .eq. i) goto 20
         dtj = ti - tnow(j)
         do 30 k = 1, NDIM
            posj(k) = xpre(j, k, order)
#ifndef STANDARD
            velj(k) = vpre(j, k, order)
#endif
            do 40  l=order-1, 0, -1
               posj(k) = posj(k)*dtj + xpre(j,k,l)
#ifndef STANDARD
               velj(k) = velj(k)*dtj + vpre(j,k,l)
#endif
40          continue
            dx(k) = (posj(k)*dtj+vdata(j,k))*dtj+xdata(j,k)-posi(k)
#ifdef STANDARD
            dv(k) = velj(k)*dtj+vdata(j,k)-veli(k)
#endif
30       continue

         r2 = dx(1)*dx(1)+dx(2)*dx(2)+dx(3)*dx(3)+eps2
         r2inv = 1.0d0/r2
         rfact = mass(j)*r2inv*sqrt(r2inv)
#ifndef STANDARD
         rdfact = -3d0*(dx(1)*dv(1)+dx(2)*dv(2)+dx(3)*dv(3))*rfact*r2inv
#endif
         
         do 50 k = 1, NDIM
            acci(k) = acci(k) + dx(k)*rfact
#ifndef STANDARD
            adoti(k) = adoti(k) + dx(k)*rdfact + dv(k)*rfact
#endif
50       continue
20    continue
#else
#   ifdef FAST_SQRT
      eps = sqrt(eps2)
#   endif      
      do 20 j = 1, nbody
         dtjv(j)=ti-tnow(j)
20    continue
cvd$ noconcur
      do 30 k=1,NDIM
         do 40 j=1,nbody
            posjv(j,k) = xpre(j, k, order)
            veljv(j,k) = vpre(j, k, order)
40       continue
         do 50  l=order-1, 0, -1
            do 60 j=1,nbody
               posjv(j,k) = posjv(j,k)*dtjv(j) + xpre(j,k,l)
               veljv(j,k) = veljv(j,k)*dtjv(j) + vpre(j,k,l)
60          continue
50       continue
         do 70 j=1,nbody    
            dxv(j,k) = (posjv(j,k)*dtjv(j)+vdata(j,k))*dtjv(j)
     $               +xdata(j,k)-posi(k)
            dvv(j,k) = veljv(j,k)*dtjv(j)+vdata(j,k)-veli(k)
70       continue
30    continue
cvd$ noconcur
      do 90 l=1,2
         if(l .eq. 1) then
            j1=1
            j2=i-1
         else
            j1=i+1
            j2=nbody
         endif
         do 80 j=j1,j2
            r2 = dxv(j,1)*dxv(j,1)+dxv(j,2)*dxv(j,2)
     $              +dxv(j,3)*dxv(j,3)+eps2
#ifndef FAST_SQRT
            r2inv = 1.0d0/r2
            rfact = mass(j)*r2inv*sqrt(r2inv)
#else
            r = 0.5*(abs(dxv(j,1)) + abs(dxv(j,2))
     $               + abs(dxv(j,3)) + eps)
            rinv = 2.0*r/(r*r + r2)
            rinv = 0.5*rinv*(3.0-rinv**2*r2)
            rinv = 0.5*rinv*(3.0-rinv**2*r2)
            rinv = 0.5*rinv*(3.0-rinv**2*r2)
            rinv = 0.5*rinv*(3.0-rinv**2*r2)
            rinv = 0.5*rinv*(3.0-rinv**2*r2)
            r2inv = rinv*rinv
            rfact = mass(j)*r2inv*rinv
#endif            
            
     
            rdfact = -3d0*(dxv(j,1)*dvv(j,1)
     $            +dxv(j,2)*dvv(j,2)+dxv(j,3)*dvv(j,3))*rfact*r2inv
            acci(1) = acci(1) + dxv(j,1)*rfact
            acci(2) = acci(2) + dxv(j,2)*rfact
            acci(3) = acci(3) + dxv(j,3)*rfact
            adoti(1) = adoti(1) + dxv(j,1)*rdfact + dvv(j,1)*rfact 
            adoti(2) = adoti(2) + dxv(j,2)*rdfact + dvv(j,2)*rfact 
            adoti(3) = adoti(3) + dxv(j,3)*rdfact + dvv(j,3)*rfact 

80       continue
90    continue
#endif
      end

      subroutine stptcl(ibody, ordi, xi, vi, xip, vip, ddti, tli, dti)
c-----------------------------------------------------------------------      
c  stptcl : store particle data to common arrays
c-----------------------------------------------------------------------
c  ibody	in	particle id to get data 
c  ordi		in	order of integrator for particle i
c  xi		in	particle position
c  vi		in	velocity
c  xip		in	predicted position at next timestep
c  vip		in     predicted velocity at next timestep
c  ddti		in	divided difference table
c  tli		in	old time list
c  dti  	in	current time step 
c-----------------------------------------------------------------------      
#include "pcinc.h"
      INTEGER ibody
      REAL xi(NDIM), vi(NDIM), xip(NDIM), vip(NDIM)
      REAL ddti(NDIM, 0:MAXSTEPS), dti, tli(0:MAXSTEPS)
      INTEGER ordi
      INTEGER j, k
      
      do 10 k = 1, NDIM
         post(ibody, k) = xi(k)
         velt(ibody, k) = vi(k)
         posnxt(k, ibody) = xip(k)
         velnxt(k, ibody) = vip(k)
10    continue
      dt(ibody) = dti
      tnow(ibody) = tli(0)
      tnext(ibody) = dti + tli(0)
      ordnow(ibody) = ordi
      do 20 j = 0, ordi
         do 30 k = 1, NDIM
            ddt(k, j, ibody) = ddti(k,j)
30       continue
         tlist(j, ibody) = tli(j)
20    continue
      end

      subroutine gtptcl(ibody, xi,vi,xip,vip,ddti,tli, dti,  ordi)
c-----------------------------------------------------------------------      
c  gtptcl : get particle data and save it to local variables
c-----------------------------------------------------------------------
c  ibody	in	particle id to get data 
c  xi		out	particle position
c  vi		out	velocity
c  xip		out	predicted position at next timestep
c  vip		out     predicted velocity at next timestep
c  ddti		out	divided difference table
c  tli		out	old time list
c  dti  	out	current time step 
c  ordi		out	order of integrator for particle i
c-----------------------------------------------------------------------      
#include "pcinc.h"
      INTEGER ibody
      REAL xi(NDIM), vi(NDIM), xip(NDIM), vip(NDIM)
      REAL ddti(NDIM, 0:MAXSTEPS), dti, tli(0:MAXSTEPS)
      INTEGER ordi
      INTEGER j, k
      
      do 10 k = 1, NDIM
         xi(k) = post(ibody, k)
         vi(k) = velt(ibody, k)
         xip(k) = posnxt(k, ibody)
         vip(k) = velnxt(k, ibody)
10    continue
      dti = dt(ibody)
      ordi = ordnow(ibody)
      do 20 j = 0, ordi
         do 30 k = 1, NDIM
            ddti(k,j) = ddt(k, j, ibody)
30       continue
         tli(j) = tlist(j, ibody)
20    continue
      end
            
      subroutine intgrt(tout)
c-----------------------------------------------------------------------      
c  intgrt : update one particle
c-----------------------------------------------------------------------
#include "pcinc.h"
      REAL tout
      INTEGER ibody
      REAL ti, xi(NDIM), vi(NDIM), xip(NDIM), vip(NDIM), acci(NDIM)
      REAL adoti(NDIM), dti, dtnewi, ddti(NDIM, 0:MAXSTEPS)
      REAL pcoefi(NDIM, 0:MAXSTEPS), tli(0:MAXSTEPS), dtli(MAXSTEPS)
      INTEGER ordi, k, j
c---  SELECT PARTICLE TO MOVE
#ifdef DEBUG
      write(6,*)'Enter intgrt'
#endif
      call selnxt(tnext, nbody, ibody)
      ti = tnext(ibody)
      tsys = ti
      if(tsys .gt. tout) then
         tsys = tout
         return
      endif
#ifdef DEBUG
      DBGWRITE(ti)
      DBGWRITE(ibody)
#endif
c--   EXTRACT PARTICLE DATA INTO LOCAL VARIABLES      
      call gtptcl(ibody, xi,vi,xip,vip,ddti,tli, dti, ordi)
#ifdef DEBUG
      DBGWRITE(xi)
      DBGWRITE(vi)
      DBGWRITE(xip)
      DBGWRITE(vip)
      do 9000 j=0,ordi
9000     DBGWRITE((ddti(k,j),k=1,NDIM))
      DBGWRITE((tli(j),j=0,ordi))
      DBGWRITE(dti)
      DBGWRITE(ordi)
#endif
c---  CALCULATE FORCE ETC.
      call accel(post, velt, pospre, velpre, mass, tnow, loword, 
     $           nbody, xip, vip, ti, ibody, eps2, acci, adoti)
#ifdef DEBUG
      DBGWRITE(acci)
      DBGWRITE(adoti)
#endif      
c---  CONTROL THE ORDER OF THE INTEGRATOR
      call adjord(ordmax, ordflg(ibody), herflg, ordi)      
#ifdef DEBUG
      DBGWRITE(ordflg(ibody))
      DBGWRITE(ordi)
#endif      
c---  UPDATE DIVIDED DIFFERENCE AND TIME LIST
      if(herflg .eq. 0) then
         call upddt(ddti, tli, ti, acci, ordi, dtli)
      else
         call upddth(ddti, tli, ti, acci, adoti, ordi, dtli)
      endif
#ifdef DEBUG
      do 9010 j=0,ordi
9010     DBGWRITE((ddti(k,j),k=1,NDIM))
      DBGWRITE((tli(j),j=0,ordi))
#endif      
c---  OBTAIN DERIVATIVES
      call dtopow(ddti, dtli, ordi, pcoefi)
#ifdef DEBUG
      DBGWRITE(pcoefi)
#endif      
c---  APPLY CORRECTOR
      if(herflg .eq. 0) then
         call correc(pcoefi, dtli(1), ordi, xi, vi)
      else
         call correc(pcoefi, dtli(2), ordi, xi, vi)
      endif
#ifdef DEBUG
      DBGWRITE(xi)
      DBGWRITE(vi)
#endif      
c---  DETERMINE NEW TIMESTEP
      if(dttyp1 .gt. 0 .and. ordi .gt. ordmax - 3) then
      call dtgen(pcoefi, eta, vip,vi, acci, tolva, tolvr, 
     $                 dti, ordi, exfact, 
     $                 dttyp1, dttyp2, dtnewi)
#ifdef DEBUG
         write(6,600) 'dtold, new ', 0, dti, dtnewi
600      format(' ',a,i5,4g15.8)
#endif
         dti = min(dtmax, dtnewi)
      endif

c---  OBTAIN NEXT PREDECTED POSITION
      call predic(pcoefi, dti, xi, vi, ordi, xip, vip)
c---  STORE UPDATED VARIABLES INTO COMMON BLOCK
      call stptcl(ibody, ordi, xi, vi, xip, vip, ddti, tli, dti)
c---  SET LOWER ORDER PREDICTOR
      call setpre(pcoefi, ordi, ibody, MAXBODY, MAXSTEPS, pospre, velpre)
      end

