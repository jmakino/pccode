#include "pcheader.h"
c----------------------------------------------------------------------
c
c     randam number in (0,1)
c
c----------------------------------------------------------------------
      function rand()
      REAL   rand
      REAL   rnow, rmod
      integer irdum
      common /rndcom/ rnow, rmod
      rnow = mod(rnow*48828125, rmod)
      irdum = rnow
      rnow = irdum
      rand   = rnow*0.4656613e-9
      end
c
c----------------------------------------------------------------------
c
c     random number generator initializer
c
c     should be called with some large prime number
c
c----------------------------------------------------------------------
      subroutine srand(seed)
      integer * 4 seed
      REAL   rnow, rmod
      common /rndcom/ rnow, rmod
      rmod = 2.0**31
      rnow = seed
c     write(6,*) 'rmod = ', rmod
      end
c
c----------------------------------------------------------------------
c
c     random number in (-1,1)
c
c----------------------------------------------------------------------
      function rand2()
      REAL rand2
      REAL   rand
      rand2 = rand()*2.0 - 1.0
      end
c
      subroutine wrvect(scale)
      REAL   rand2
      REAL   scale
      write(8,600) rand2()*scale, rand2()*scale, rand2()*scale
600   format(3g18.9)
      end
c
c----------------------------------------------------------------------
c
c     random vector in 3-D sphere with radius SCALE
c
c----------------------------------------------------------------------
      subroutine grvect(scale,x1,x2,x3)
      REAL   rand2
      REAL   scale, x1, x2, x3
10    continue
         x1 =  rand2()*scale
         x2 =  rand2()*scale
         x3 =  rand2()*scale
      if(x1*x1+x2*x2+x3*x3 .gt. scale*scale) goto 10
      end
c
      subroutine ranout
c---------------------------------------------------------------------------
c
c     random initial condition file generator
c
c---------------------------------------------------------------------------
      integer i, n
      REAL mass
      read(5,*) n
      write(8,*) n
      write(8,*) 3
      write(8,*) 0.0
      mass = 1.0/n
      call srand(584287)
      do 10 i = 1, n
         write(8,600) mass
600      format(g18.9)
10    continue
      do 20 i = 1, n
         call wrvect(1.0d0)
20    continue
      do 30 i = 1, n
         call wrvect(0.707d0)
30    continue
      end
c
c---------------------------------------------------------------------------
c
c     random initial condition array generator
c
c---------------------------------------------------------------------------
      subroutine srunif(n, m, e, mass, x, v, nary, nvary,dim)
      integer n, nary, nvary, dim
      REAL m, e, x(nary,dim), v(nvary, dim), mass(nary)
      REAL   masseq, scale
      INTEGER   i
      call srand(584287)
c     call srand(184287)
      masseq = m/n
      do 10 i = 1, n
         mass(i) = masseq
10    continue
      scale = 0.3*m*m/e
      do 20 i = 1, n
         call grvect(scale, x(i,1), x(i,2), x(i,3))
20    continue
      scale = sqrt(e/m/0.3)
      do 30 i = 1, n
         call grvect(scale, v(i,1), v(i,2), v(i,3))
c        write(6,600) v(i,1), v(i,2), v(i,3), i
600      format(3f15.5, i10)
30    continue
      end
c
      subroutine set3dv(r,x,y,z)
      REAL   r,x,y,z
      REAL   rand
      REAL   rxy,theta
      z=(1-2*rand())*r
      theta=3.1415926535*2*rand()
      rxy=sqrt(r*r-z*z)
      x=rxy*sin(theta)
      y=rxy*cos(theta)
      end
c
c---------------------------------------------------------------------------
c
c     plummer model configulation (one particle)
c
c     For details see Aarseth et. al., Astron. & Astrophys. 37, 183
c     (1974).
c
c---------------------------------------------------------------------------
      subroutine stplcn(r, x, y, z)
      REAL   x, y, z, r
      REAL   rand
c
      r=(rand()**(-0.666666666666)-1)**(-0.5)
      call set3dv(r,x,y,z)
      end
c
      function qofve()
      REAL   qofve
      REAL   rand
      REAL   q, g, limit
10    continue
         limit=rand()*0.1
         q = rand()
         g=q*q*(1-q*q)**(3.5)
         if (limit .ge. g)goto 10
      qofve=q
c     write(6,*) 'q=',qofve
      end
c
c---------------------------------------------------------------------------
c
c     plummer model velosity      (one particle)
c
c     For details see Aarseth et. al., Astron. & Astrophys. 37, 183
c     (1974).
c
c---------------------------------------------------------------------------
      subroutine stplv (r,vx,vy,vz)
      REAL   r,vx,vy,vz
      REAL   rand,qofve
      REAL   vabs, vesc
c
      vesc=1.41421356*(1+r*r)**(-0.25)
      vabs = vesc*qofve()
      call set3dv(vabs,vx,vy,vz)
c     write(6,*) r,vabs, vx, vy, vz
      end
c
c---------------------------------------------------------------------------
c
c     plummer model (M=R=1)
c
c---------------------------------------------------------------------------
      subroutine setpl0(n, mass, x, v, nary, nvary,dim)
      integer n, nary, nvary, dim
      REAL m, e, x(nary,dim), v(nvary, dim), mass(nary)
      REAL   masseq, rabs
      INTEGER   i
      call srand(584287)
c     call srand(184287)
      masseq = 1.0/n
      do 10 i = 1, n
         mass(i) = masseq
10    continue
      do 20 i = 1, n
         call stplcn(rabs, x(i,1), x(i,2), x(i,3))
         call stplv (rabs, v(i,1), v(i,2), v(i,3))
20    continue
      end
c
c---------------------------------------------------------------------------
c
c     adjust vectors so that the sum=0
c
c---------------------------------------------------------------------------
      subroutine adjvec(n,mass,vector,nvec,dim)
      integer n,nvec,dim
      REAL   mass(nvec), vector(nvec, dim)
      integer i,j
      REAL   sum, offset,mtot
      do 30 j=1,dim
         sum=0.0
         mtot = 0.0
         do 10 i=1,n
            sum=sum+mass(i)*vector(i,j)
            mtot = mtot+mass(i)
10       continue
         offset = sum/mtot
         write(6,*)i, j, offset
         do 20 i=1,n
            vector(i,j)=vector(i,j)-offset
20       continue
30    continue
      end
      
      subroutine adjene(n,x,v,mass,nary, nvary, dim, e)     
      integer n, nary, nvary, dim
      REAL e, x(nary,dim), v(nvary, dim), mass(nary)
      REAL   spot, skin, sx, sv, etot
      integer i,j,k
      spot = 0
      skin = 0
      do 10 i=1, n
         spot = spot + mass(i)*(v(i,1)**2+v(i,2)**2+v(i,3)**2)*0.5
         do 20 j=1,i-1
            spot = spot - mass(i)*mass(j)
     $                   /sqrt((x(i,1)-x(j,1))**2
     $                        +(x(i,2)-x(j,2))**2
     $                        +(x(i,3)-x(j,3))**2)
20       continue
10    continue
      etot = spot + skin
      sx = abs(etot/e)
      sv = sqrt(1D0/sx)
      do 30 i=1,n
         do 40 k=1,3
            x(i,k)=x(i,k)*sx
            v(i,k)=v(i,k)*sv
40       continue
30    continue
      end
c
c---------------------------------------------------------------------------
c
c     plummer model (scaled)
c
c---------------------------------------------------------------------------
      subroutine setplm (n, m, e, mass, x, v, nary, nvary,dim)
      integer n, nary, nvary, dim
      REAL m, e, x(nary,dim), v(nvary, dim), mass(nary)
      REAL   sm, sx, sv, s0
      integer i,j
      call setpl0(n,mass, x, v, nary, nvary, dim)
      sm = m
      s0=3.1415926535*3/64.0
      sx=s0*m*m/e
      sv=sqrt(e/(m*s0))
      do 10 j=1,dim
         do 20 i=1,n
            x(i,j)=x(i,j)*sx
            v(i,j)=v(i,j)*sv
20       continue
10    continue
      do 30 i=1,n
         mass(i)=mass(i)*sm
30    continue
      call adjvec(n,mass,x,nary,dim)
      call adjvec(n,mass,v,nvary,dim)
      call adjene(n,x,v,mass,nary, nvary, dim, e)
      end
