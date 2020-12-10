/*--------------------------------------------------------------------*/
/* Skeleton 1-2/2D Electromagnetic Gridless OpenMP PIC code */
/* written by Viktor K. Decyk, UCLA */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <sys/time.h>
#include "mbpush1gl.h"
#include "omplib.h"

void dtimer(double *time, struct timeval *itime, int icntrl);

int main(int argc, char *argv[]) {
/* indx = exponent which determines grid points in x direction: */
/* nx = 2**indx. */
   int indx =   9;
/* npx = number of electrons distributed in x direction. */
   int  npx =  18432;
/* tend = time at end of simulation, in units of plasma frequency. */
/* dt = time interval between successive calculations. */
/* qme = charge on electron, in units of e. */
   double tend = 10.0, dt = 0.05, qme = -1.0;
/* vtx/vty = thermal velocity of electrons in x/y direction */
/* vx0/vy0 = drift velocity of electrons in x/y direction. */
   double vtx = 1.0, vty = 1.0, vx0 = 0.0, vy0 = 0.0;
/* vtx/vz0 = thermal/drift velocity of electrons in z direction */
   double vtz = 1.0, vz0 = 0.0;
/* omx = magnetic field electron cyclotron frequency in x */
   double omx = 0.0;
/* ax = smoothed particle size in x direction */
/* ci = reciprocal of velocity of light. */
   double ax = 0.4, ci = 0.1;
/* idimp = number of particle coordinates = 4 */
/* ipbc = particle boundary condition: 1 = periodic */
/* relativity = (no,yes) = (0,1) = relativity is used */
   int  idimp = 4, ipbc = 1, relativity = 1;
/* wke/we = particle kinetic/electrostatic field energy */
/* wf/wm/wt = magnetic field/transverse electric field/total energy */
   double wke = 0.0, we = 0.0, wf = 0.0, wm = 0.0, wt = 0.0;
/* declare scalars for standard code */
   int j;
   int  np, nx, nxh, nxe, nxeh, nxhm;
   int  ntime, nloop, isign;
   double qbme, affp, dth, cux0;

/* declare scalars for OpenMP code */
   int irc, nvp;

/* declare arrays for standard code: */
/* part = particle array */
   double *part = NULL;
/* qe/qi = electron charge density with guard cells */
/* fxe = smoothed longitudinal electric field with guard cells */
   double complex *qe = NULL, *qi = NULL, *fxe = NULL;
/* cue = electron current density with guard cells */
/* fxyze/byze = smoothed electric/magnetic field with guard cells */
   double complex *cue = NULL, *fxyze = NULL, *byze = NULL;
/* eyz/byz = transverse electric/magnetic field in fourier space */
/* dcu = time derivative of electron current density */
   double complex *eyz = NULL, *byz = NULL, *dcu = NULL;
/* ffc = form factor array for poisson solver */
   double complex *ffc = NULL;

/* declare and initialize timing data */
   double time;
   struct timeval itime;
   double tdpost = 0.0, tfield = 0.0, tdjpost = 0.0, tpush = 0.0;
   double dtime;

   irc = 0;
/* nvp = number of shared memory nodes  (0=default) */
   nvp = 0;
/* printf("enter number of nodes:\n"); */
/* scanf("%i",&nvp);                   */
/* initialize for shared memory parallel processing */
   cinit_omp(nvp);

/* initialize scalars for standard code */
/* np = total number of particles in simulation */
/* nx = number of grid points in x direction */
   np = npx; nx = 1L<<indx; nxh = nx/2;
   nxe = nx + 2;
/* nloop = number of time steps in simulation */
/* ntime = current time step */
   nloop = tend/dt + .0001; ntime = 0;
   qbme = qme;
   affp = (double) nx/(double ) np;
/* nxhm = number of fourier modes kept */
   nxhm = nxh; nxeh = nxe/2 > nxhm ? nxe/2 : nxhm;
   dth = 0.0;

/*allocate data for standard code */
   part = (double *) malloc(idimp*np*sizeof(double));
   qe = (double complex *) malloc(nxeh*sizeof(double complex));
   qi = (double complex *) malloc(nxeh*sizeof(double complex));
   fxe = (double complex *) malloc(nxeh*sizeof(double complex));
   fxyze = (double complex *) malloc(3*nxeh*sizeof(double complex));
   cue = (double complex *) malloc(2*nxeh*sizeof(double complex));
   byze = (double complex *) malloc(2*nxeh*sizeof(double complex));
   eyz = (double complex *) malloc(2*nxeh*sizeof(double complex));
   byz = (double complex *) malloc(2*nxeh*sizeof(double complex));
   ffc = (double complex *) malloc(nxhm*sizeof(double complex));

/* calculate form factors for gridless code: updates ffc */
   cffc1initgl(ffc,ax,affp,nx,nxhm);
/* calculate form factors for conventional PIC code: updates ffc */
/* isign = 0; */
/* cpois1gl(qe,fxe,isign,ffc,ax,affp,&we,nx,nxhm,nxeh); */
/* initialize electrons */
   cdistr1h(part,vtx,vty,vtz,vx0,vy0,vz0,npx,idimp,np,nx,ipbc);
/* initialize ion background */
   for (j = 0; j < nxeh; j++) {
      qi[j] = 0.0 + 0.0*_Complex_I;
   }
   cmdpost1gl(part,qi,-qme,np,idimp,nx,nxhm,nxeh);

/* initialize transverse electromagnetic fields */
   for (j = 0; j < 2*nxeh; j++) {
      eyz[j] = 0.0 + 0.0*_Complex_I;
      byz[j] = 0.0 + 0.0*_Complex_I;
   }
   cux0 = 0.0;

   if (dt > (0.63*ci*((double) nxh)/((double) nxhm))) {
      printf("Info: Courant condition may be exceeded!\n");
   }

/* * * * start main iteration loop * * * */

L500: if (nloop <= ntime)
         goto L2000;
/*    printf("ntime = %i\n",ntime); */

/* deposit current with standard procedure: updates part, cue, cux0 */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < 2*nxeh; j++) {
         cue[j] = 0.0 + 0.0*_Complex_I;
      }
      if (relativity==1)
         cmrdjpost1gl(part,cue,&cux0,qme,dth,ci,np,idimp,nx,nxhm,nxeh);
      else
         cmdjpost1gl(part,cue,&cux0,qme,dth,np,idimp,nx,nxhm,nxeh);
      dtimer(&dtime,&itime,1);
      time = dtime;
      tdjpost += time;

/* deposit charge with standard procedure: updates qe */
      dtimer(&dtime,&itime,-1);
      for (j = 0; j < nxeh; j++) {
         qe[j] = 0.0 + 0.0*_Complex_I;
      }
      cmdpost1gl(part,qe,qme,np,idimp,nx,nxhm,nxeh);
      for (j = 0; j < nxeh; j++) {
         qe[j] += qi[j];
      }
      dtimer(&dtime,&itime,1);
      time = dtime;
      tdpost += time;

/* calculate electromagnetic fields in fourier space with standard */
/* procedure: updates eyz, byz */
      dtimer(&dtime,&itime,-1);
/* initialization assumes initial acceleration of particles is zero */
      if (ntime==0) {
         dcu = (double complex *) malloc(2*nxeh*sizeof(double complex));
         for (j = 0; j < 2*nxeh; j++) {
            dcu[j] = 0.0 + 0.0*_Complex_I;
         }
/* initialize with darwin electric and magnetic fields */
/* deposit time derivative of current: updates dcu */
         if (relativity==1)
            cmrd2jpost1gl(part,dcu,qme,ci,np,idimp,nx,nxhm,nxeh);
         else
            cmd2jpost1gl(part,dcu,qme,np,idimp,nx,nxhm,nxeh);
/* initialize electromagnetic fields from free-streaming particles */
/* deposit scaled current and time derivative of current: */
/* updates cue, dcu */
/*       for (j = 0; j < 2*nxeh; j++) { */
/*          cue[j] = 0.0 + 0.0*_Complex_I; */
/*       } */
/*       if (relativity==1) */
/*          cmrdsjpost1gl(part,cue,dcu,qme,ci,np,idimp,nx,nxhm,nxeh); */
/*                                                                    */
/*       else                                                     */
/*          cmdsjpost1gl(part,cue,dcu,qme,ci,np,idimp,nx,nxhm,nxeh); */
/*                                                                   */
/* initialize magnetic field: updates byz, wm */
         cibpois13gl(cue,byz,ffc,ci,&wm,nx,nxhm,nxeh);
         wf = 0.0;
/* initialize transverse electric field: updates eyz */
         isign = 1;
         cepois13gl(dcu,eyz,isign,ffc,ax,affp,0.0,ci,&wf,nx,nxhm,nxeh);
         free(dcu);
         dth = 0.5*dt;
      }
      else {
         camaxwel1gl(eyz,byz,cue,ffc,ci,dt,&wf,&wm,nx,nxhm,nxeh);
      }
      dtimer(&dtime,&itime,1);
      time = dtime;
      tfield += time;

/* calculate force/charge in fourier space with standard procedure: */
/* updates fxe, we */
      dtimer(&dtime,&itime,-1);
      isign = -1;
      cpois1gl(qe,fxe,isign,ffc,ax,affp,&we,nx,nxhm,nxeh);
      dtimer(&dtime,&itime,1);
      time = dtime;
      tfield += time;

/* add longitudinal and transverse electric fields with standard */
/* procedure: updates fxyze */
      dtimer(&dtime,&itime,-1);
      isign = 1;
      cemfield1gl(fxyze,fxe,eyz,ffc,nxhm,nxeh);
/* copy magnetic field with standard procedure: updates byze */
      isign = -1;
      cbmfield1gl(byze,byz,ffc,nxhm,nxeh);
      dtimer(&dtime,&itime,1);
      time = dtime;
      tfield += time;

/* push particles with standard procedure: updates part, wke */
      wke = 0.0;
      dtimer(&dtime,&itime,-1);
      if (relativity==1)
/* analytic Boris mover, gamma constant during time step */
         cmarbpush13gl(part,fxyze,byze,omx,qbme,dt,dth,ci,&wke,idimp,
                       np,nx,nxhm,nxeh);
/* exact analytic Boris mover, gamma varies during time step */
/*       cmearbpush13gl(part,fxyze,byze,omx,qbme,dt,dth,ci,&wke,idimp, */
/*                      np,nx,nxhm,nxeh);                              */
      else
         cmabpush13gl(part,fxyze,byze,omx,qbme,dt,dth,&wke,idimp,np,nx,
                      nxhm,nxeh);
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
   printf("nxhm = %i\n",nxhm);
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

