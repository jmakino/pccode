c----------------------------------------------------------------------
c    file startup.f 
c
c    start-up module
c    for general n-body integrator using Krogh variable step integrator
c
c    ver 0.00 1989-03-18
c
c    by jun makino
c----------------------------------------------------------------------
c
#include "pcheader.h"
c
      function random(i)
c-----------------------------------------------------------------------
c random : return the pseudo uniform random number in (0,1)
c-----------------------------------------------------------------------
c  i		in/out	seed
c-----------------------------------------------------------------------
      REAL   random
      INTEGER   i
      REAL   tmp
      tmp = i
      tmp = mod(tmp*56789,1d9)
      random = tmp*1d-9
      i = int(tmp)
      end
c
      subroutine initav(nbody, nmax, ndim, x, v, mass, eps2, a, adot)
c----------------------------------------------------------------------
c initav : get the initial value of acceleration & da/dt
c          This routine is vectorized. Array argument has the form
c          array(nmax, ndim) 
c          G = 1 is assumed
c-----------------------------------------------------------------------
c  nbody	in	number of particles
c  nmax		in	1st dimension of array (max number of particles)
c  ndim		in	number of space dimensions
c  x		in	position of particles
c  v		in	velocity of particles
c  mass		in	mass of particles
c  eps2  	in	softening parameter squared
c  a		out	acceleration
c  adot		out	da/dt
c-----------------------------------------------------------------------
      INTEGER nbody, nmax, ndim
      REAL x(nmax, ndim), v(nmax, ndim), mass(nmax), eps2
      REAL a(nmax, ndim), adot(nmax, ndim)
      REAL   dx,dy,dz,r2,r,r3inv,dvx,dvy,dvz
      REAL   fix,fiy,fiz,fjx,fjy,fjz,fi,fj,fdi,fdj,ai,aj,aterm
      INTEGER   i,j,k
      do 10 k = 1, ndim
         do 20 i = 1, nbody
            a(i,k) = 0.0d0
            adot(i,k) = 0.0d0
20       continue
10    continue
      do 30 i = 1, nbody - 1
         do 40 j = i+1, nbody
            dx = x(i,1) - x(j,1)
            dy = x(i,2) - x(j,2)
            dz = x(i,3) - x(j,3)
            dvx = v(i,1) - v(j,1)
            dvy = v(i,2) - v(j,2)
            dvz = v(i,3) - v(j,3)
            r2 =  dx*dx+dy*dy+dz*dz+eps2
            r3inv = 1.0d0/(r2*sqrt(r2))
            fi = -Mass(j)*r3inv
            fj = Mass(i)*r3inv
            fix = fi*dx
            fiy = fi*dy
            fiz = fi*dz
            fjx = fj*dx
            fjy = fj*dy
            fjz = fj*dz
            aterm = -3d0*(dx*dvx + dy*dvy + dz*dvz)/r2
            a(i,1) = a(i,1) + fix
            a(i,2) = a(i,2) + fiy
            a(i,3) = a(i,3) + fiz
            a(j,1) = a(j,1) + fjx
            a(j,2) = a(j,2) + fjy
            a(j,3) = a(j,3) + fjz
            adot(i,1) = adot(i,1) +fi*dvx + aterm*fix
            adot(i,2) = adot(i,2) +fi*dvy + aterm*fiy
            adot(i,3) = adot(i,3) +fi*dvz + aterm*fiz
            adot(j,1) = adot(j,1) +fj*dvx + aterm*fjx
            adot(j,2) = adot(j,2) +fj*dvy + aterm*fjy
            adot(j,3) = adot(j,3) +fj*dvz + aterm*fjz
40       continue
30    continue
      end
c
      subroutine dtinit(v, a, adot, epsini, dt0)
c-----------------------------------------------------------------------
c dtinit :  determine the size of the first time step
c-----------------------------------------------------------------------
c  v		in	velocity of a particle
c  a		in	acceleration of a particle
c  adot		in	da/dt of a particle
c  epsini	in	tolerance parameter for the initial timestep
c  dt0		out	the initial timestep
c-----------------------------------------------------------------------
      REAL v(NDIM), a(NDIM), adot(NDIM), epsini
      REAL dt0
      INTEGER   i,j,k
      REAL   vabs,aabs,adabs
      REAL   dtfact, f
      parameter (dtfact = 1d0)
      vabs = 0d0
      aabs = 0d0
      adabs =0d0
      do 10 k = 1, NDIM
         vabs = vabs + v(k)**2
         aabs = aabs + a(k)**2
         adabs = adabs + adot(k)**2
10    continue
      dt0 = epsini * sqrt(max(vabs/aabs,aabs/adabs))
#ifdef DEBUG
      write(6,*)'DTINIT vabs, aabs, adabs, epsini, dt0'
      write(6,*)vabs, aabs, adabs, epsini, dt0
#endif
      end
