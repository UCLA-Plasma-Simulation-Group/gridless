/*---------------------------------------------------------------------*/
/* Skeleton 1D Electrostatic Gridless PIC code */
/* written by Viktor K. Decyk, UCLA */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <sys/time.h>
#include "push1gl.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

int main(int argc, char *argv[]) {
/* indx = exponent which determines grid points in x direction: */
/* nx = 2**indx */
   int indx =   9;
/* npx = number of electrons distributed in x direction */
   int npx =  18432;
/* tend = time at end of simulation, in units of plasma frequency */
/* dt = time interval between successive calculations */
/* qme = charge on electron, in units of e */
   double tend = 10.0, dt = 0.1, qme = -1.0;
/* vtx = thermal velocity of electrons in x direction */
/* vx0 = drift velocity of electrons in x direction */
   double vtx = 1.0, vx0 = 0.0;
/* ax = smoothed particle size in x direction */
   double ax = 0.4;
/* idimp = number of particle coordinates = 2 */
/* ipbc = particle boundary condition: 1 = periodic */
   int idimp = 2, ipbc = 1;
/* wke/we/wt = particle kinetic/electric field/total energy */
   double wke = 0.0, we = 0.0, wt = 0.0;
/* declare scalars for standard code */
   int j;
   int np, nx, nxh, nxe, nxeh, nxhm;
   int ntime, nloop, isign;
   double qbme, affp;

/* declare array for standard code: */
/* part = particle array */
   double *part = NULL;
/* qe/qi = electron/ion charge density with guard cells */
   double complex *qe = NULL, *qi = NULL;
/* fxe = smoothed electric field with guard cells */
   double complex *fxe = NULL;
/* ffc = form factor array for poisson solver */
   double complex *ffc = NULL;
/* sctx = scratch array for sines and cosines */
   double complex *sctx = NULL;
  
/* declare and initialize timing data */
   double time;
   struct timeval itime;
   double tdpost = 0.0, tfield = 0.0, tpush = 0.0;
   double dtime;

/* initialize scalars for standard code */
/* np = total number of particles in simulation */
/* nx = number of grid points in x direction */
   np = npx; nx = 1L<<indx; nxh = nx/2;
   nxe = nx + 2;
/* nloop = number of time steps in simulation */
/* ntime = current time step */
   nloop = tend/dt + .0001; ntime = 0;
   qbme = qme;
   affp = (double) nx/(double) np;
/* nxhm = number of fourier modes kept */
   nxhm = nxh; nxeh = nxe/2 > nxhm ? nxe/2 : nxhm;

/* allocate data for standard code */
   part = (double *) malloc(idimp*np*sizeof(double));
   qe = (double complex *) malloc(nxeh*sizeof(double complex));
   qi = (double complex *) malloc(nxeh*sizeof(double complex));
   fxe = (double complex *) malloc(nxeh*sizeof(double complex));
   ffc = (double complex *) malloc(nxhm*sizeof(double complex));
   sctx = (double complex *) malloc(nxeh*sizeof(double complex));

/* calculate form factors for gridless code: updates ffc */
   cffc1initgl(ffc,ax,affp,nx,nxhm);
/* calculate form factors for conventional PIC code: updates ffc */
/* isign = 0; */
/* cpois1gl(qe,fxe,isign,ffc,ax,affp,&we,nx,nxhm,nxeh); */
/* initialize electrons */
   cdistr1(part,vtx,vx0,npx,idimp,np,nx,ipbc);
/* initialize ion background */
   for (j = 0; j < nxeh; j++) {
      qi[j] = 0.0 + 0.0*_Complex_I;
   }
   cdpost1gl(part,qi,sctx,-qme,np,idimp,nx,nxhm,nxeh);

/* * * * start main iteration loop * * * */
 
L500: if (nloop <= ntime)
         goto L2000;
/*    printf("ntime = %i\n",ntime); */
 
/* deposit charge and add ion background with standard procedure: */
/* updates qe */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < nxeh; j++) {
         qe[j] = 0.0 + 0.0*_Complex_I;
      }
      cdpost1gl(part,qe,sctx,qme,np,idimp,nx,nxhm,nxeh);
      for (j = 0; j < nxeh; j++) {
         qe[j] += qi[j];
      }
      dtimer(&dtime,&itime,1);
      time = dtime;
      tdpost += time;

/* calculate force/charge in fourier space with standard procedure: */
/* updates fxe, we */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cpois1gl(qe,fxe,isign,ffc,ax,affp,&we,nx,nxhm,nxeh);
      dtimer(&dtime,&itime,1);
      time = dtime;
      tfield += time;

/* push particles with standard procedure: updates part, wke */
      wke = 0.0;
      dtimer(&dtime,&itime,-1);
      cpush1gl(part,fxe,sctx,qbme,dt,&wke,idimp,np,nx,nxhm,nxeh);
      dtimer(&dtime,&itime,1);
      time = dtime;
      tpush += time;

      if (ntime==0) {
         printf("Initial Field, Kinetic and Total Energies:\n");
         printf("%e %e %e\n",we,wke,wke+we);
      }
      ntime += 1;
      goto L500;
L2000:

/* * * * end main iteration loop * * * */

   printf("ntime = %i\n",ntime);
   printf("nxhm = %i\n",nxhm);
   printf("Final Field, Kinetic and Total Energies:\n");
   printf("%e %e %e\n",we,wke,wke+we);

   printf("\n");
   printf("deposit time = %f\n",tdpost);
   printf("solver time = %f\n",tfield);
   printf("push time = %f\n",tpush);
   time = tdpost + tpush;
   printf("total particle time = %f\n",time);
   wt = time + tfield;
   printf("total time = %f\n",wt);
   printf("\n");

   wt = 1.0e+09/(((double) nloop)*((double) np));
   printf("Push Time (nsec) = %f\n",tpush*wt);
   printf("Deposit Time (nsec) = %f\n",tdpost*wt);
   printf("Total Particle Time (nsec) = %f\n",time*wt);

   return 0;
}
