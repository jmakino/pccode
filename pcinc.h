#include "pcheader.h"
c---  NUMBER OF PARTICLES
      INTEGER nbody
c---  Lower order predictor array for pos & vel
      REAL pospre(MAXBODY, NDIM, 0:MAXSTEPS)
      REAL velpre(MAXBODY, NDIM, 0:MAXSTEPS)
c--   Divided Difference Table      
      REAL ddt(NDIM, 0:MAXSTEPS, MAXBODY)
c--   Time Table
      REAL tlist(0:MAXSTEPS, MAXBODY)
c--   particle data
      REAL pos(NDIM, MAXBODY), vel(NDIM, MAXBODY), mass(MAXBODY)
c--   PARTICLE DATA TRANSPOSED (FOR VECTORIZED CALCULATION)
      REAL post(MAXBODY, NDIM), velt(MAXBODY, NDIM)
c--   PARTICLE DATA FOR OUTPUT/DIAGNOSTICS
      REAL posout(MAXBODY, NDIM), velout(MAXBODY, NDIM)
c--   current timestep
      REAL dt(MAXBODY)
c--   current time of particles
      REAL tnow(MAXBODY)
c--   Next time of particles
      REAL tnext(MAXBODY)
c--   Next Predicted position/velocity
      REAL posnxt(NDIM,MAXBODY), velnxt(NDIM,MAXBODY)
c--   Integrator order control variables/flags
      INTEGER ordnow(MAXBODY), ordflg(MAXBODY), ordmax, loword
c--   Hermite integrator flag
      INTEGER herflg
c--   Softening parameter
      REAL eps2, eps
c--   TIMESTEP SCHEME INDICATOR
      INTEGER dttyp1, dttyp2
c--   TIME STEP PARAMETERS
      REAL dtmax, eta, tolva, tolvr, exfact, epsini
c--
      CHARACTER * 60 header, inname, outnam
      LOGICAL snapou
c--   
      REAL tstop, dtout, tsys
c--   
      REAL einit, enow, demax
c--   
      INTEGER steps, steps0
      common /rdata/pospre, velpre, ddt, tlist, pos, vel, mass, 
     $              dt, tnow, tnext, posnxt, velnxt, post, velt,
     $              posout, velout
      common /idata/ ordnow, ordflg
      common /iparms/ nbody, ordmax, herflg, dttyp1, dttyp2, steps,
     $                loword, steps0
      common /rparms/ eps, eps2, dtmax, eta, tolva, tolvr, exfact,
     $                epsini, tstop, dtout, tsys,  einit, enow, demax
      common /cparms/ header, inname, outnam
      common /lparms/ snapou
