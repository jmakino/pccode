#include "pcheader.h"
      subroutine outbod(outuni)
c-----------------------------------------------------------------------
c  outbod : Write snapshot to file (atos format)
c-----------------------------------------------------------------------
c  outuni	in	unit # to write
c-----------------------------------------------------------------------

#include "pcinc.h"
      INTEGER outuni
      integer i,k, ndim
      ndim = NDIM
      if ( .not. snapou) return
      write(outuni,*)nbody
      write(outuni,*)ndim
      write(outuni,600) tsys
 600  format(1p3e23.15)
      do 10 i=1,nbody
        write(outuni,600)mass(i)
 10   continue
      do 20 i=1,nbody
        write(outuni,600)(posout(i,k),k=1,ndim)
 20   continue
      do 30 i=1,nbody
        write(outuni,600)(velout(i,k),k=1,ndim)
 30   continue
#ifdef UNIX
      call flush(outuni)
#endif      
      end
      
      subroutine iniclk
c-----------------------------------------------------------------------      
c iniclk : Initialize the clock (HITAC ONLY)
c-----------------------------------------------------------------------      
#ifdef HITAC
      call clock
#endif
      end            

      subroutine getsec(t) 
c-----------------------------------------------------------------------
c getsec : timer routine interface
c-----------------------------------------------------------------------
c  t		out	consumed CPU time
c
c  note: MSDOS version returns 0
c-----------------------------------------------------------------------
      TIMEARG t
#ifdef CYBER
      t=second()
#endif
#ifdef MSDOS
      t=0
#endif
#ifdef UNIX
      real etime
      real tarray(2), tmp
      tmp=etime(tarray)
      t=tarray(1)
#endif
#ifdef FACOM
      real * 8 cputim
      call clock(cputim, 0, 2)
      t = cputim
#endif
#ifdef SX2
      real * 8 cputim
      call clock(cputim)
      t = cputim
#endif
#ifdef HITAC
      real * 4 time4
      call clock(time4, 5)
      t = time4
#endif
      end
c      
      subroutine prnttm	
c-----------------------------------------------------------------------
c prnttm : Print out current CPU time to stdout (unit 6)
c-----------------------------------------------------------------------
      TIMEARG t
      call getsec(t)
      write(6,600)t/60.0
 600  format(' CPU min = ',1pg17.7)
      end
c      
      subroutine sync
c-----------------------------------------------------------------------      
c  sync : synchronize all particles  to tsys
c
c  apply predictor to all particles and set the result to arrays
c  posout & velout
c-----------------------------------------------------------------------
#include "pcinc.h"
      INTEGER i, j, k
      REAL dti
      do 10 i=1, nbody
         dti= tsys - tnow(i)
         do 20 k = 1, NDIM
            posout(i,k)=pospre(i,k,ordmax)
            velout(i,k)=velpre(i,k,ordmax)
            do 30 j=ordmax-1,0,-1
               posout(i,k)=posout(i,k)*dti+pospre(i,k,j)            
               velout(i,k)=velout(i,k)*dti+velpre(i,k,j)            
30          continue
            posout(i,k)=(posout(i,k)*dti+velt(i,k))*dti+post(i,k)
            velout(i,k)=velout(i,k)*dti+velt(i,k)
20       continue
10    continue
      end
c      
      subroutine diag(eflag)
c-----------------------------------------------------------------------
c  diag : Diagnostics (E, L etc)
c-----------------------------------------------------------------------
c eflag		in	Flag to indicate the first call 
c
c EFLAG=1 means that the initial total energy is not obtained yet
c-----------------------------------------------------------------------
#include "pcinc.h"
      integer eflag
      integer i,j,k
      integer stpflg
      REAL cm(NDIM), mtot, cmv(NDIM), ke, pe, etot, de, htot
      call sync
      do 10 k = 1,NDIM
        cm(k)=0.0
        cmv(k)=0.0
 10   continue
      mtot=0.0
      do 20 k=1,NDIM		/* center of mass and its velocity */
        do 30 i=1,nbody
          cm(k)=cm(k)+mass(i)*posout(i,k)
          cmv(k)=cmv(k)+mass(i)*velout(i,k)
 30     continue
 20   continue
      ke=0.0
      pe=0.0
      do 40 i=1,nbody		/* kinetic energy & total mass */
        mtot=mtot+mass(i)
        ke=ke+mass(i)
     $       *(velout(i,1)**2+velout(i,2)**2+velout(i,3)**2)
40    continue
      do 70 k=1,NDIM		
        cm(k)=cm(k)/mtot
        cmv(k)=cmv(k)/mtot
 70   continue
      ke=ke*0.5
      do 50 i=1,nbody-1		/* potential energy */
        do 60 j=i+1,nbody
          pe=pe-mass(i)*mass(j)
     $         /sqrt((posout(i,1)-posout(j,1))**2
     $         +(posout(i,2)-posout(j,2))**2
     $         +(posout(i,3)-posout(j,3))**2+eps2)
 60     continue
 50   continue
      etot = pe+ke
      if(eflag .eq. 0) then	/* first call. record einit */
         einit=etot
         de=0.0
         demax = 0.0
         stpflg = 0
      else			/* not first call. calculate error */
         if(stpflg .eq. 0) then
            steps0 = steps
            stpflg = 1
         endif
         enow=etot
         de=(etot-einit)/(einit)
         demax = max(demax, abs(de))
      endif
      write(6,600) tsys, etot, de, pe, -ke/pe
 600  format(' T: ',f9.2,' E: ',f11.7, ' DE: ',e10.3, ' PE: ',f8.4,
     $     ' Q: ',f5.2)
      write(6,*) ' steps:', steps

      write(6,601)'CM ', cm
      write(6,601)'CMV', cmv
 601  format(' ',a3, ' : ', 1p3g17.7)
      call prnttm		/* print CPU time */
#ifdef UNIX
      call flush(6)
#endif      
      end
