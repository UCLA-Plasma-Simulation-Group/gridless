/*--------------------------------------------------------------------*/
/* Skeleton 2D Electrostatic Gridless PIC code */
/* written by Viktor K. Decyk, UCLA */
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <sys/time.h>
#include "push2gl.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

int main(int argc, char *argv[]) {
/* indx/indy = exponent which determines grid points in x/y direction: */
/* nx = 2**indx, ny = 2**indy. */
   int indx =   7, indy =   7;
/* npx/npy = number of electrons distributed in x/y direction. */
   int  npx =  768, npy =   768;
/* ndim = number of velocity coordinates = 2 */
   int ndim = 2;
/* tend = time at end of simulation, in units of plasma frequency. */
/* dt = time interval between successive calculations. */
/* qme = charge on electron, in units of e. */
   double tend = 10.0, dt = 0.1, qme = -1.0;
/* vtx/vty = thermal velocity of electrons in x/y direction */
/* vx0/vy0 = drift velocity of electrons in x/y direction. */
   double vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0;
/* ax/ay = smoothed particle size in x/y direction */
   double ax = 0.4, ay = 0.4;
/* idimp = number of particle coordinates = 4 */
/* ipbc = particle boundary condition: 1 = periodic */
   int idimp = 4, ipbc = 1;
/* wke/we/wt = particle kinetic/electric field/total energy */
   double wke = 0.0, we = 0.0, wt = 0.0;
/* declare scalars for standard code */
   int j;
   int np, nx, ny, nxh, nyh, nxe, nye, nxeh, nxhm, nyhm;
   int ntime, nloop, isign;
   double qbme, affp;

/* declare arrays for standard code: */
/* part = particle array */
   double *part = NULL;
/* qe/qi = electron/ion charge density with guard cells */
   double complex *qe = NULL, *qi = NULL;
/* fxye = smoothed electric field with guard cells */
   double complex *fxye = NULL;
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
/* nx/ny = number of grid points in x/y direction */
   np = npx*npy; nx = 1L<<indx; ny = 1L<<indy;
   nxh = nx/2; nyh = 1 > ny/2 ? 1 : ny/2;
   nxe = nx + 2; nye = ny + 1; nxeh = nxe/2;
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

/* allocate data for standard code */
   part = (double *) malloc(idimp*np*sizeof(double));
   qe = (double complex *) malloc(nxeh*nye*sizeof(double complex));
   qi = (double complex *) malloc(nxeh*nye*sizeof(double complex));
   fxye = (double complex *) malloc(ndim*nxeh*nye*sizeof(double complex));
   ffc = (double complex *) malloc(nxhm*nyhm*sizeof(double complex));
   sctx = (double complex *) malloc(nxeh*sizeof(double complex));

/* calculate form factors for gridless code: updates ffc */
   cffc2initgl(ffc,ax,ay,affp,nx,ny,nxhm,nyhm);
/* calculate form factors for conventional PIC code: updates ffc */
/* isign = 0; */
/* cpois22gl(qe,fxye,isign,ffc,ax,ay,affp,&we,nx,ny,nxhm,nyhm,nxeh,nye); */
/* initialize electrons */
   cdistr2(part,vtx,vty,vx0,vy0,npx,npy,idimp,np,nx,ny,ipbc);
/* initialize ion background */
   for (j = 0; j < nxeh*nye; j++) {
      qi[j] = 0.0 + 0.0*_Complex_I;
   }
   cdpost2gl(part,qi,sctx,-qme,np,idimp,nx,ny,nxhm,nyhm,nxeh,nye);

/* * * * start main iteration loop * * * */

L500: if (nloop <= ntime)
         goto L2000;
/*    printf("ntime = %i\n",ntime); */

/* deposit charge and add ion background with standard procedure: */
/* updates qe*/
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

/* calculate force/charge in fourier space with standard procedure: */
/* updates fxye, we */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cpois22gl(qe,fxye,isign,ffc,ax,ay,affp,&we,nx,ny,nxhm,nyhm,nxeh,
                nye)  ; 
      dtimer(&dtime,&itime,1);
      time = dtime;
      tfield += time;

/* push particles with standard procedure: updates part, wke */
      wke = 0.0;
      dtimer(&dtime,&itime,-1);
      cpush2gl(part,fxye,sctx,qbme,dt,&wke,idimp,np,nx,ny,nxhm,nyhm,
               nxeh,nye,ipbc);
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
   printf("nxhm, nyhm = %i,%i\n",nxhm,nyhm);
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

   wt = 1.0e+09/(((float) nloop)*((float) np));
   printf("Push Time (nsec) = %f\n",tpush*wt);
   printf("Deposit Time (nsec) = %f\n",tdpost*wt);
   printf("Total Particle Time (nsec) = %f\n",time*wt);

   return 0;
}
