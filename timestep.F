#define DBGWRITE(x) WRITE(6,*)'x :', x
#include "pcheader.h"
c-----------------------------------------------------------------------
c   timestep.F : timestep control routines for Krogh
c 		 scheme LLM for second-order ODE.
c
c   89-03-04   J.Makino
c-----------------------------------------------------------------------

      subroutine dtvabs(vp,vc, tolva, dtold, order, exfact, gmode,newdt)
c-----------------------------------------------------------------------      
c  dtvabs : Gives timestep using the absolute error of velocity
c-----------------------------------------------------------------------
c   vp		in	predicted velocity
c   vc		in	corrected velocity
c   tolva	in	tolerance for the velocity absolute error
c   dtold	in	previous timestep
c   order       in      order of power-series
c   exfact	in	maximum expand factor for timestep
c   gmode 	in	flag to limit global error
c   newdt       out 	new timestep
c
c   Timesteps is given by the following formula:
c
c   estimated error = |vc-vp|
c   tolerance   = tolva 	(gmode = 0)
c                 tolva*dtold 	(gmode = 1)
c   dt =  dtold * (tolerance/estimated error)^(1.0/(order+2))
c   if dt > dtold*exfact then dt = dtold *exfact
c
c   Derivation:
c   for both of predictor and corrector, the formula is
c   vnew=vold+p0*dt+p1*dt^2+.....pk*dt^k+1 + O(dt^k+2)
c
c   Thus, the error term is O(dt^k+2). Since the coefficient for the
c   error term of predictor is by far larger than that for the corrector,
c   |vp-vc| gives good estimate for the error of predictor. For the
c   next timestep, we choose dt so that the error term of the predictor
c   will be roughly equal to the given tolerance.
c
c   If the timestep is constant, error_predictor ~= C eror_corrector
c   (C >>1). So the actual truncation error is also roughly limited.
c 
c   It is possible to apply the difference of predictor and corrector
c   as the "modifier", to eliminate O(dt^k+2) term. This procedure is
c   not adopted in this code to reduce round-off sensitivity.
c-----------------------------------------------------------------------      
      REAL vp(NDIM), vc(NDIM)
      REAL tolva, dtold, exfact
      INTEGER order, gmode
      REAL newdt
      REAL dvabs2, tol
      
      if(NDIM .ne. 3) call terror('DIM not supported in dtvabs')
      if(gmode .eq. 0) then
         tol = tolva
      else
         tol = tolva*dtold
      endif
      dvabs2 = (vp(1)-vc(1))**2+(vp(2)-vc(2))**2+(vp(3)-vc(3))**2
c      write(6,*)vp
c      write(6,*)vc
      if(dvabs2 .ne. 0d0) then
         newdt = dtold * min(exfact,
     $                     (tol**2/dvabs2)**(0.5/(order+2)))
#ifdef DEBUG_TIMESTEP     
         write(6,600)' tol2, dvabs2 = ', tol**2, dvabs2, newdt
600      format(a, 3g18.9)
#endif
      else
         newdt = dtold*exfact
      endif
      end
c
      subroutine dtvcmp(vp,vc, a, tolva, tolvr, dtold, order, exfact,
     $                gmode, newdt)
c-----------------------------------------------------------------------      
c  dtvcmp : Gives timestep using composite error of velocity
c-----------------------------------------------------------------------
c   vp		in	predicted velocity
c   vc		in	corrected velocity
c   a           in      acceleration
c   tolva	in	tolerance for velocity absolute error
c   tolvr       in      tolerance for velocity relative error 
c   dtold	in	previous timestep
c   order       in      order of power-series
c   exfact	in	maximum expand factor for timestep
c   gmode	in	flag for global/local error estimate
c   newdt	out	new timestep
c   Timesteps is given by the following formula:
c
c   estimated error = |vc-vp|
c   dt =  dtold * (tolerance/estimated error)^(1.0/(order+2))
c   if dt > dtold*exfact then dt = dtold *exfact
c
c   tolerance = tol_abs + a*tol_rel     (gmode = 0)
c   tolerance = dt(tol_abs + a*tol_rel) (gmode = 1)
c-----------------------------------------------------------------------      
      REAL vp(NDIM), vc(NDIM), a(NDIM)
      REAL tolva, tolvr, dtold, exfact
      INTEGER order, gmode
      REAL newdt
      REAL dvabs2, tol2
      
      if(NDIM .ne. 3) call terror('DIM not supported in dtvabs')
      tol2 = tolva**2 + (tolvr**2)*(a(1)**2+a(2)**2+a(3)**2)
      if(gmode .ne. 0) tol2 = tol2*dtold**2
#ifdef DEBUG
      DBGWRITE(dtold)
      DBGWRITE(a(1)**2+a(2)**2+a(3)**2)
      DBGWRITE(gmode)
      DBGWRITE(tolva)
      DBGWRITE(tolvr)
      DBGWRITE(tol2)
#endif
      dvabs2 = (vp(1)-vc(1))**2+(vp(2)-vc(2))**2+(vp(3)-vc(3))**2
      if(dvabs2 .ne. 0d0) then
         newdt = dtold * min(exfact,
     $                     (tol2/dvabs2)**(0.5/(order+2)))
#ifdef DEBUG_TIMESTEP     
         write(6,600)' tol2, dvabs2 = ', tol2, dvabs2, newdt
600      format(a, 3g18.9)
#endif
      else
         newdt = dtold*exfact
      endif
      end
c
      subroutine dtvcmx(vp,vc, a, tolva, tolvr, dtold, order, exfact,
     $                gmode, newdt)
c-----------------------------------------------------------------------      
c  dtvcmx : Gives timestep using composite error of velocity
c-----------------------------------------------------------------------
c   vp		in	predicted velocity
c   vc		in	corrected velocity
c   a           in      acceleration
c   tolva	in	tolerance for velocity absolute error
c   tolvr       in      tolerance for velocity relative error 
c   dtold	in	previous timestep
c   order       in      order of power-series
c   exfact	in	maximum expand factor for timestep
c   gmode	in	flag for global/local error estimate
c   newdt	out	new timestep
c   Timesteps is given by the following formula:
c
c   estimated error = |vc-vp|
c   dt =  dtold * (tolerance/estimated error)^(1.0/(order+2))
c   if dt > dtold*exfact then dt = dtold *exfact
c
c   tolerance = 1/(1/tol_abs + 1/a*tol_rel)     (gmode = 0)
c   tolerance = dt/(1/tol_abs + 1/a*tol_rel) (gmode = 1)
c-----------------------------------------------------------------------      
      REAL vp(NDIM), vc(NDIM), a(NDIM)
      REAL tolva, tolvr, dtold, exfact
      INTEGER order, gmode
      REAL newdt
      REAL dvabs2, tol2
      
      if(NDIM .ne. 3) call terror('DIM not supported in dtvabs')
      tol2 = 1.0/(1.0/tolva**2 
     $      + 1.0/((tolvr**2)*(a(1)**2+a(2)**2+a(3)**2)))
      if(gmode .ne. 0) tol2 = tol2*dtold**2
#ifdef DEBUG
      DBGWRITE(dtold)
      DBGWRITE(a(1)**2+a(2)**2+a(3)**2)
      DBGWRITE(gmode)
      DBGWRITE(tolva)
      DBGWRITE(tolvr)
      DBGWRITE(tol2)
#endif
      
      dvabs2 = (vp(1)-vc(1))**2+(vp(2)-vc(2))**2+(vp(3)-vc(3))**2
      if(dvabs2 .ne. 0d0) then
         newdt = dtold * min(exfact,
     $                     (tol2/dvabs2)**(0.5/(order+2)))
#ifdef DEBUG_TIMESTEP     
         write(6,600)' tol2, dvabs2 = ', tol2, dvabs2, newdt
600      format(a, 3g18.9)
#endif
      else
         newdt = dtold*exfact
      endif
      end
c
      subroutine dtaars(pcoef, eta, dtold, order, exfact, newdt)
c-----------------------------------------------------------------------      
c  dtaars : Gives Aarseth Timesteps
c-----------------------------------------------------------------------
c   pcoef	in	power-series coefficients for acceleration
c   eta		in	accuracy parameter (Aarseth definition)
c   dtold	in	previous timestep
c   order       in      order of power-series
c   exfact	in	maximum expand factor for timestep
c   newdt	out	new timestep
c
c   Timesteps is given by the following formula:
c
c   dt = sqrt[eta(f*f2+f1^2)/(f1*f3+f2^2)]
c   if dt > dtold*exfact then dt = dtold *exfact
c-----------------------------------------------------------------------      
      REAL pcoef(NDIM,0:MAXSTEPS)
      REAL eta, dtold, exfact
      INTEGER order
      REAL newdt
      
      REAL fabs(0:MAXSTEPS)
      INTEGER i,k
      if(order .lt. 3) then
         newdt = dtold
         return
      endif
      DO 10 i=0, 3
         fabs(i)=0.0d0
         do 20 k=1, NDIM
            fabs(i) = fabs(i)+pcoef(k,i)**2
20       continue
         fabs(i)=sqrt(fabs(i))
c         write(6,*)' fabs ', i, fabs(i)
10    continue
      fabs(2)=fabs(2)*2
      fabs(3)=fabs(3)*6
c      write(6,*) ' Fabs:', (fabs(i),i=0,3)
      newdt = min(dtold*exfact,
     $ sqrt(eta*(fabs(0)*fabs(2)+fabs(1)**2)
     $         /(fabs(1)*fabs(3)+fabs(2)**2)))
      end
c      
      subroutine dtgen(pcoef, eta, vp,vc, a, tolva, tolvr, 
     $                 dtold, order, exfact, method, gmode, newdt)
c-----------------------------------------------------------------------      
c  dtvabs : Gives dt by specified scheme
c-----------------------------------------------------------------------
c   pcoef	in	power-series coefficients for acceleration
c   eta		in	accuracy parameter (Aarseth definition)
c   vp		in	predicted velocity
c   vc		in	corrected velocity
c   a		in	acceleration
c   tolva	in	tolerance for the velocity absolute error
c   tolvr       in      tolerance for velocity relative error 
c   dtold	in	previous timestep
c   order       in      order of power-series
c   exfact	in	maximum expand factor for timestep
c   method	in	scheme to use (1: aarseth, 2:Vabs, 3:Vrel)
c   gmode 	in	flag to limit global error(1:tol is for global
c							error)
c   newdt       out 	new timestep
c-----------------------------------------------------------------------
      REAL pcoef(NDIM,0:MAXSTEPS), eta
      REAL vp(NDIM), vc(NDIM), a(NDIM)
      REAL tolva, tolvr, dtold, exfact
      INTEGER order, gmode, method
      REAL newdt
      if(method .eq. 1) then
         call dtaars(pcoef, eta, dtold, order, exfact, newdt)
      else if(method .eq. 2) then
         call dtvabs(vp, vc, tolva, dtold, order, exfact, gmode, newdt)
      else if(method .eq. 3) then
         call dtvcmp(vp, vc, a, tolva, tolvr, dtold, 
     $               order, exfact, gmode, newdt)
      else if(method .eq. 4) then
         call dtvcmx(vp, vc, a, tolva, tolvr, dtold, 
     $               order, exfact, gmode, newdt)
      else
         newdt = dtold
      endif
      end

