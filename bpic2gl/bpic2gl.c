/*--------------------------------------------------------------------*/
/* Skeleton 2-1/2D Electromagnetic Gridless PIC code */
/* written by Viktor K. Decyk, UCLA */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <sys/time.h>
#include "bpush2gl.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

int main(int argc, char *argv[]) {
/* indx/indy = exponent which determines grid points in x/y direction: */
/* nx = 2**indx, ny = 2**indy. */
   int indx =   7, indy =   7;
/* npx/npy = number of electrons distributed in x/y direction. */
   int npx =  768, npy =   768;
/* ndim = number of velocity coordinates = 3 */
   int ndim = 3;
/* tend = time at end of simulation, in units of plasma frequency. */
/* dt = time interval between successive calculations. */
/* qme = charge on electron, in units of e. */
   double tend = 10.0, dt = 0.04, qme = -1.0;
/* vtx/vty = thermal velocity of electrons in x/y direction */
/* vx0/vy0 = drift velocity of electrons in x/y direction. */
   double vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0;
/* vtx/vz0 = thermal/drift velocity of electrons in z direction */
   double vtz = 1.0, vz0 = 0.0;
/* ax/ay = smoothed particle size in x/y direction */
/* ci = reciprocal of velocity of light. */
   double ax = 0.4, ay = 0.4, ci = 0.1;
/* idimp = number of particle coordinates = 5 */
/* ipbc = particle boundary condition: 1 = periodic */
/* relativity = (no,yes) = (0,1) = relativity is used */
   int idimp = 5, ipbc = 1, relativity = 1;
/* wke/we = particle kinetic/electrostatic field energy */
/* wf/wm/wt = magnetic field/transverse electric field/total energy */
   double wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0;
/* declare scalars for standard code */
   int j;
   int np, nx, ny, nxh, nyh, nxe, nye, nxeh, nxhm, nyhm;
   int ntime, nloop, isign;
   double qbme, affp, at, dth;

/* declare arrays for standard code: */
/* part = particle array */
   double *part = NULL;
/* qe/qi = electron/ion charge density */
   double complex *qe = NULL, *qi = NULL;
/* cue = electron current density */
/* fxyze/bxyze = smoothed electric/magnetic field */
   double complex *cue = NULL, *fxyze = NULL, *bxyze = NULL;
/* exyz/bxyz = transverse electric/magnetic field in fourier space */
/* dcu = time derivative of electron current density */
   double complex *exyz = NULL, *bxyz = NULL, *dcu = NULL;
/* ffc = form factor array for poisson solver */
   double complex *ffc = NULL;
/* sctx = scratch array for sines and cosines */
   double complex *sctx = NULL;

/* declare and initialize timing data */
   double time;
   struct timeval itime;
   double tdpost = 0.0, tfield = 0.0, tdjpost = 0.0, tpush = 0.0;
   double dtime;

/* initialize scalars for standard code */
/* np = total number of particles in simulation */
/* nx/ny = number of grid points in x/y direction */
   np = npx*npy; nx = 1L<<indx; ny = 1L<<indy;
   nxh = nx/2; nyh = 1 > ny/2 ? 1 : ny/2;
   nxe = nx + 2; nye = ny + 1;
/* nloop = number of time steps in simulation */
/* ntime = current time step */
   nloop = tend/dt + .0001; ntime = 0;
   qbme = qme;
   affp = (double) (nx*ny)/(double ) np;
/* nxhm/nyhm = number of fourier modes kept in x/y */
   nxhm = nxh; nyhm = nyh;
/* nxhm = nxh+1; nyhm = nyh+1 */
   nxeh = nxe/2 > nxhm ? nxe/2 : nxhm;
   nye = nye > 2*nyhm ? nye : 2*nyhm;
   dth = 0.0;

/* allocate data for standard code */
   part = (double *) malloc(idimp*np*sizeof(double));
   qe = (double complex *) malloc(nxeh*nye*sizeof(double complex));
   qi = (double complex *) malloc(nxeh*nye*sizeof(double complex));
   fxyze = (double complex *) malloc(ndim*nxeh*nye*sizeof(double complex));
   cue = (double complex *) malloc(ndim*nxeh*nye*sizeof(double complex));
   bxyze = (double complex *) malloc(ndim*nxeh*nye*sizeof(double complex));
   exyz = (double complex *) malloc(ndim*nxeh*nye*sizeof(double complex));
   bxyz = (double complex *) malloc(ndim*nxeh*nye*sizeof(double complex));
   ffc = (double complex *) malloc(nxhm*nyhm*sizeof(double complex));
   sctx = (double complex *) malloc(nxeh*sizeof(double complex));

/* calculate form factors for gridless code: updates ffc */
   cffc2initgl(ffc,ax,ay,affp,nx,ny,nxhm,nyhm);
/* calculate form factors for conventional PIC code: updates ffc */
/* isign = 0; */
/* cpois23gl(qe,fxyze,isign,ffc,ax,ay,affp,&we,nx,ny,nxhm,nyhm,nxeh, */
/*           nye);                                                   */
/* initialize electrons */
   cdistr2h(part,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,idimp,np,nx,ny,ipbc);
/* initialize ion background */
   for (j = 0; j < nxeh*nye; j++) {
      qi[j] = 0.0 + 0.0*_Complex_I;
   }
   cdpost2gl(part,qi,sctx,-qme,np,idimp,nx,ny,nxhm,nyhm,nxeh,nye);

/* initialize transverse electromagnetic fields */
   for (j = 0; j < ndim*nxeh*nye; j++) {
      exyz[j] = 0.0 + 0.0*_Complex_I;
      bxyz[j] = 0.0 + 0.0*_Complex_I;
    }
 
   at = sqrt(pow(((double) nxhm)/((double) nxh),2) + 
             pow(((double) nyhm)/((double) nyh),2));
   if (dt > (0.63*ci/at)) {
      printf("Warning: Courant condition may be exceeded!\n");
   }

/* * * * start main iteration loop * * * */

L500: if (nloop <= ntime)
         goto L2000;
/*    printf("ntime = %i\n",ntime); */

/* deposit current with standard procedure: updates part, cue */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < ndim*nxeh*nye; j++) {
         cue[j] = 0.0 + 0.0*_Complex_I;
      }
      if (relativity==1)
         crdjpost2gl(part,cue,sctx,qme,dth,ci,np,idimp,nx,ny,nxhm,nyhm,
                     nxeh,nye,ipbc);
      else
         cdjpost2gl(part,cue,sctx,qme,dth,np,idimp,nx,ny,nxhm,nyhm,nxeh,
                    nye,ipbc);
      dtimer(&dtime,&itime,1);
      time = dtime;
      tdjpost += time;

/* deposit charge with standard procedure: updates qe */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < nxeh*nye; j++) {
         qe[j] = 0.0 + 0.0*_Complex_I;
      }
      cdpost2gl(part,qe,sctx,qme,np,idimp,nx,ny,nxhm,nyhm,nxeh,nye);
      for (j = 0; j < nxeh*nye; j++) {
         qe[j] += qi[j];
      }
      dtimer(&dtime,&itime,1);
      time = dtime;
      tdpost += time;

/* take transverse part of current with standard procedure: updates cue */
      dtimer(&dtime,&itime,-1);
      ccuperp2gl(cue,nx,ny,nxhm,nyhm,nxeh,nye);
      dtimer(&dtime,&itime,1);
      time = dtime;
      tfield += time;

/* calculate electromagnetic fields in fourier space with standard */
/* procedure: updates exyz, bxyz, wf, wm */
      dtimer(&dtime,&itime,-1);
/* initialization assumes initial acceleration of particles is zero */
      if (ntime==0) {
         dcu = (double complex *) malloc(ndim*nxeh*nye*sizeof(double complex));
         for (j = 0; j < ndim*nxeh*nye; j++) {
            dcu[j] = 0.0 + 0.0*_Complex_I;
         }
/* initialize with darwin electric and magnetic fields */
/* deposit time derivative of current: updates dcu */
         if (relativity==1)
            crd2jpost2gl(part,dcu,sctx,qme,ci,np,idimp,nx,ny,nxhm,nyhm,
                         nxeh,nye);
         else
            cd2jpost2gl(part,dcu,sctx,qme,np,idimp,nx,ny,nxhm,nyhm,nxeh,
                        nye);
/* initialize electromagnetic fields from free-streaming particles */
/* deposit scaled current and time derivative of current: */
/* updates cue, dcu */
/*       for (j = 0; j < ndim*nxeh*nye; j++) { */
/*          cue[j] = 0.0 + 0.0*_Complex_I;     */
/*       }                                     */
/*       if (relativity==1) */
/*          crdsjpost2gl(part,cue,dcu,sctx,qme,ci,np,idimp,nx,ny,nxhm, */
/*                       nyhm,nxeh,nye);                               */
/*       else                                                          */
/*          cdsjpost2gl(part,cue,dcu,sctx,qme,ci,np,idimp,nx,ny,nxhm,  */
/*                      nyhm,nxeh,nye);                                */
/*       ccuperp2gl(cue,nx,ny,nxhm,nyhm,nxeh,nye); */
/* initialize magnetic field: updates bxyz, wm */
         cibpois23gl(cue,bxyz,ffc,ci,&wm,nx,ny,nxhm,nyhm,nxeh,nye);
         wf = 0.0;
/* calculates transverse part of the derivative of the current */
/* updates: dcu */
         ccuperp2gl(dcu,nx,ny,nxhm,nyhm,nxeh,nye);
/* calculates transverse electric field: updates exyz, wf */
         isign = 1;
         cepois23gl(dcu,exyz,isign,ffc,ax,ay,affp,0.0,ci,&wf,nx,ny,nxhm,
                    nyhm,nxeh,nye);
         free(dcu);
         dth = 0.5*dt;
      }
/* calculates transverse electromagnetic fields: */
/* updates exyz, bxyz, wf, wm */
      else {
         camaxwel2gl(exyz,bxyz,cue,ffc,ci,dt,&wf,&wm,nx,ny,nxhm,nyhm,
                     nxeh,nye);
      }
      dtimer(&dtime,&itime,1);
      time = dtime;
      tfield += time;

/* calculate force/charge in fourier space with standard procedure: */
/* updates fxyze, we */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cpois23gl(qe,fxyze,isign,ffc,ax,ay,affp,&we,nx,ny,nxhm,nyhm,nxeh,
                nye);
      dtimer(&dtime,&itime,1);
      time = dtime;
      tfield += time;

/* add longitudinal and transverse electric fields with standard */
/* procedure: updates fxyze */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cemfield2gl(fxyze,exyz,ffc,isign,nxhm,nyhm,nxeh,nye);
/* copy magnetic field with standard procedure: updates bxyze */
      isign = -1;
      cemfield2gl(bxyze,bxyz,ffc,isign,nxhm,nyhm,nxeh,nye);
      dtimer(&dtime,&itime,1);
      time = dtime;
      tfield += time;

/* push particles with standard procedure: updates part, wke */
      wke = 0.0;
      dtimer(&dtime,&itime,-1);
      if (relativity==1)
/* analytic Boris mover, gamma constant during time step */
         carbpush23gl(part,fxyze,bxyze,sctx,qbme,dt,dth,ci,&wke,idimp,
                      np,nx,ny,nxhm,nyhm,nxeh,nye,ipbc);
/* exact analytic Boris mover, gamma varies during time step */
/*       cearbpush23gl(part,fxyze,bxyze,sctx,qbme,dt,dth,ci,&wke,idimp, */
/*                     np,nx,ny,nxhm,nyhm,nxeh,nye,ipbc);               */
      else
         cabpush23gl(part,fxyze,bxyze,sctx,qbme,dt,dth,&wke,idimp,np,nx,
                     ny,nxhm,nyhm,nxeh,nye,ipbc);
      dtimer(&dtime,&itime,1);
      time = dtime;
      tpush += time;

      if (ntime==0) {
         wt = we + wf + wm;
         printf("Initial Total Field, Kinetic and Total Energies:\n");
         printf("%e %e %e\n",wt,wke,wke+wt);
         printf("Initial Electrostatic, Transverse Electric and Magnetic \
Field Energies:\n");
         printf("%e %e %e\n",we,wf,wm);
      }
      ntime += 1;
      goto L500;
L2000:

/* * * * end main iteration loop * * * */

   printf("ntime, relativity = %i,%i\n",ntime,relativity);
   printf("nxhm, nyhm = %i,%i\n",nxhm,nyhm);
   wt = we + wf + wm;
   printf("Final Total Field, Kinetic and Total Energies:\n");
   printf("%e %e %e\n",wt,wke,wke+wt);
   printf("Final Electrostatic, Transverse Electric and Magnetic Field \
Energies:\n");
   printf("%e %e %e\n",we,wf,wm);

   printf("\n");
   printf("deposit time = %f\n",tdpost);
   printf("current deposit time = %f\n",tdjpost);
   tdpost += tdjpost;
   printf("total deposit time = %f\n",tdpost);
   printf("solver time = %f\n",tfield);
   printf("push time = %f\n",tpush);
   time = tdpost + tpush;
   printf("total particle time = %f\n",time);
   wt = time + tfield;
   printf("total time = %f\n",wt);
   printf("\n");

   wt = 1.0e+09/(((float) nloop)*((float) np));
   printf("Push Time (nsec) = %f\n",tpush*wt);
   printf("Deposit Time (nsec) = %f\n",tdpost*wt);
   printf("Total Particle Time (nsec) = %f\n",time*wt);

   return 0;
}
