#define DBGWRITE(x) WRITE(6,*)'x :', x
#include "pcheader.h"
c-----------------------------------------------------------------------
c pccode.F : main routine for the Krogh/Hermite Integrator
c
c Jun Makino
c
c Ver 0.0 89-03-27
c 
c-----------------------------------------------------------------------
      program PCCODE
#include "pcinc.h"
      call iniclk
#ifndef FOOLISH_OS
      open(5,file='indata',status='old')
#endif
      call gtparm
      call getsys(11)
#ifndef FOOLISH_OS
      close(5)
      if(snapou) open(12,file=outnam,status='new')
#endif      
      call iniout(6)
      call start
      call push
      call endout
      end
c
      subroutine gtparm
c-----------------------------------------------------------------------
c  gtparm : Read in parameters
c-----------------------------------------------------------------------
#include "pcinc.h"
699   format(' ',a)
      read(5,500)header
      read(5,500)inname		/* initial condition file */
      read(5,500)outnam		/* snapshot output file */
 500  format(a)
      if(outnam(1:1).eq.' ') then /* if outnam is not file name, does not */
        snapou = .false.	  /* produce snapshot output */
      else
        snapou = .true.
      endif
      read(5,*) nbody, tstop, dtout, eps, 
     $          ordmax, loword, herflg, dttyp1, dttyp2, dtmax, eta, 
     $          tolva, tolvr, exfact, epsini
c----------------------------------------------------------------------
c     nbody:  number of particles
c     TSTOP:  time to stop integration
c     DTOUT:  output interval
c     eps:    The softening parameter
c     ordmax: order of the integrator
c     loword: order of predictor for force calculation (0 -> ordmax)
c     herflg: flag to use Hermite integrator (1: Hermite (upto fdot)
c					      0: Lagrange ) 
c     ddtyp1:  Timestep adjust scheme  	0:Constant
c         				1:Aarseth
c					2:Absolute velocity correction
c					3:Composite Velocity Correction
c     dttyp2: Minor mode for timestep scheme 2 & 3
c             0: Local mode (error < tol)
c	      1: GLobal mode (error < tol*dt)
c     dtmax:  Maximum allowable timestep
c     eta:    Aarseth DT parameter
c     tolva:  Absolute velocity error tolerance
c     tolvr:  Relative velocity error tolerance
c     exfact: Maximum expansion factor of timestep
c     epsini: parameter for initial timestep
c----------------------------------------------------------------------
      eps2=eps**2
      if (loword .eq. 0) loword = ordmax
      if(herflg .eq. 1 .and. ordmax .eq. 1) then
         write(6,*)' herflg = ', herflg, ' ordmax=', ordmax, ' unsupp'
         stop
      endif
      end

      subroutine getsys(inunit)
c-----------------------------------------------------------------------      
c  getsys : setup an NBODY system
c-----------------------------------------------------------------------
c  inunit	in	unit # to read in the nbody system
c
c
c  if NBODY in the input parameter file is 0, read in the system from
c  file with name INNAME (atos format).
c  if NBODY > 0, create plummer model.
c
c  This routine set values in variables/arrays
c  NBODY
c  POST
c  VELT
c  MASS
c  TSYS
c-----------------------------------------------------------------------      
#include "pcinc.h"
      INTEGER inunit
      REAL mtot, etot
      INTEGER k
#ifdef DEBUG_2BODY
      REAL ecc
#endif
      if(nbody .eq. 0) then
#ifndef FOOLISH_OS
         open(inunit,file=inname, status='old')
#endif         
         call inbody(inunit, nbody, tsys, mass, post, velt, 
     $			MAXBODY, MAXBODY, NDIM)
#ifndef FOOLISH_OS
         close(inunit)
#endif                  
      else
#ifdef DEBUG_2BODY
         read(5,*) ecc
         write(6,*) ' 2-body test ecc = ', ecc
         nbody=2
         mass(1)=0.5
         mass(2)=0.5
         post(1,1)=1.0
         post(1,2)=0.0
         post(1,3)=0.0
         velt(1,1)=sqrt(0.125d0)*ecc
         velt(1,2)=sqrt(0.125d0)*sqrt(1D0-ecc**2)
         velt(1,3)=0.0
         do 10 k=1,3
           post(2,k)=-post(1,k)
           velt(2,k)=-velt(1,k)
 10   continue
#else
         mtot = 1D0
         etot = 0.25D0
         call setplm(nbody, mtot, etot, mass, post, velt, 
     $			MAXBODY, MAXBODY, NDIM)
         tsys = 0D0
#endif
      endif
      end

      subroutine inbody (inunit, nbody, t0, mass, post, velt, 
     $			 npary, nvary, ndim)
      INTEGER inunit, nbody, npary, nvary, ndim
      REAL t0
      REAL mass(npary), post(npary, ndim),  velt(nvary, ndim)
      INTEGER inndim, i, k
      
      read(inunit, *) nbody, inndim, t0
      if(inndim .ne. ndim) then
         call terror(' ERROR in INBODY: input file format error')
      endif
      read(inunit,*)(mass(i), i=1, nbody), 
     $              ((post(i,k), k=1, ndim), i=1, nbody),
     $              ((velt(i,k), k=1, ndim), i=1, nbody)
      end
      
      subroutine iniout(outuni)
c-----------------------------------------------------------------------      
c  iniout: program start diagnostics
c-----------------------------------------------------------------------
c  outuni	in	unit to write
c-----------------------------------------------------------------------      
#include "pcinc.h"
      INTEGER outuni
      write(outuni, 600)
600   format(' -------- PCCODE INITIAL LOGGING --------')
      write(outuni,601) header
601   format(' ', a)
      write(outuni,602)nbody, tsys, tstop, dtout, eps
602   format(' N=', i5, ' Tstart=', f10.5, ' Tstop=', f10.5, 
     $       ' DTout=', f9.6, ' eps=', f9.6)
      write(outuni,603)ordmax, loword, herflg, dttyp1*2+dttyp2, eta, 
     $                 tolva, tolvr, exfact, epsini
603   format(' Order=', i2, ' Preord=', i2,' Hermite=', i1,
     $       ' DTtype=', i1, 
     $       ' Aarseth eta = ', f8.5, ' Vtols=', 1p2g12.5,/,
     $       ' DT expansion factor = ', g9.1, 
     $       ' Initial dt eps = ', 1pg12.5)
      end
c
      subroutine start
c-----------------------------------------------------------------------
c  start : initialize the integrator
c-----------------------------------------------------------------------
#include "pcinc.h"
      REAL a(MAXBODY, NDIM), adot(MAXBODY, NDIM)
      REAL ai(NDIM), adoti(NDIM)
      REAL xi(NDIM), vi(NDIM), xip(NDIM), vip(NDIM)
      REAL ddti(NDIM, 0:MAXSTEPS), dti, tli(0:MAXSTEPS)
      REAL pcoefi(NDIM, 0:MAXSTEPS), dtli(MAXSTEPS)
      INTEGER ordi
      INTEGER i,j, k
c--   SET UP TABLE OF 1,1/2,,,, etc
      call strtab
c--   CLEAR timestep counter
      steps = 0
c--   CLEAR PREDICTOR ARRAY
      do 40 j=0,ordmax
         do 50 k=1,NDIM
            do 60 i=1, nbody
               pospre(i,k,j) = 0D0
               velpre(i,k,j) = 0D0
60          continue
50       continue
40    continue
c--   CHOOSE APPROPRIATE INITIAL ORDER     
      if(herflg .eq. 0) then
         ordi = 0
      else
         ordi = 1
      endif
#ifdef DEBUG
      DBGWRITE(ordi)
#endif
      do 10 i = 1, nbody
c--      Set the initial order of the integrator
         ordnow(i) = ordi
         ordflg(i) = 0
c--      Set time of particles
         tnow(i) = tsys
10    continue                  
c--   CALCULATE acc and adot      
      call initav(nbody, MAXBODY, NDIM, post, velt, mass, eps2, a, adot)
      do 20 i=1, nbody
c--      calculate initial timestep
         do 30 k = 1, NDIM
            xi(k) = post(i,k)
            vi(k) = velt(i,k)
            ai(k) = a(i,k)
            adoti(k) = adot(i,k)
30       continue
         call dtinit(vi, ai, adoti, epsini, dti)
         if(herflg .eq. 1) then
            call upddth(ddti, tli, tsys, ai, adoti, ordi, dtli)
         else
            call upddt(ddti, tli, tsys, ai, ordi, dtli)
         endif
         call dtopow(ddti, dtli, ordi, pcoefi)
         call predic(pcoefi, dti, xi, vi, ordi, xip, vip)
         call stptcl(i, ordi, xi, vi, xip, vip, ddti, tli, dti)
c-- -    SET LOWER ORDER PREDICTOR
         call setpre(pcoefi, ordi, i, MAXBODY, MAXSTEPS, pospre, velpre)
#ifdef DEBUG
         DBGWRITE(i)
         DBGWRITE(xi)
         DBGWRITE(vi)
         DBGWRITE(xip)
         DBGWRITE(vip)
         do 9000 j=0,ordi
            DBGWRITE((pcoefi(k,j),k=1,NDIM))
9000     continue
         DBGWRITE(ai)
         DBGWRITE(adoti)
         DBGWRITE(dt(i))
         DBGWRITE(tnext(i))
#endif
20    continue
      end
c
      subroutine endout
#include "pcinc.h"
      write(7,600) ordmax, herflg, steps - steps0,
     $    demax, real(steps-steps0)/nbody
600   format(2i4, i10, 1p2e15.6)
      end
c      
      subroutine push
c-----------------------------------------------------------------------      
c  push : push system to time tstop
c-----------------------------------------------------------------------      
#include "pcinc.h"
      REAL tnout
      tnout = tsys+dtout
      call diag(0)
      call outbod(12) 
10    if(tnout .le. tstop) then
20       if(tsys .lt. tnout) then
            call intgrt(tnout)
            steps = steps + 1
            goto 20
         endif
         call diag(1)
         call outbod(12)
         tnout = tnout + dtout
         goto 10
      endif
      end 
