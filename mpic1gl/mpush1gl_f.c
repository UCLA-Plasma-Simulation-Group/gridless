/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

double randum_();

void distr1_(double *part, double *vtx, double *vdx, int *npx, int *idimp,
             int *nop, int *nx, int *ipbc);

void mpush1gl_(double *part, double complex *fx, double *qbm,
               double *dt, double *ek, int *idimp, int *nop, int *nx,
               int *nxhm, int *nxvh);

void mdpost1gl_(double *part, double complex *q, double *qm, int *nop,
                int *idimp, int *nx, int *nxhm, int *nxvh);

void ffc1initgl_(double complex *ffc, double *ax, double *affp, int *nx,
                 int *nxhm);

void pois1gl_(double complex *q, double complex *fx, int *isign,
              double complex *ffc, double *ax, double *affp, double *we,
              int *nx, int *nxhm, int *nxvh);

void mefield1gl_(double complex *fx, double *gx, double *xp, int *nx,
                 int *nxhm, int *nxvh);

/* Interfaces to C */

double ranorm() {
  return ranorm_();
}

double randum() {
   return randum_();
}

/*--------------------------------------------------------------------*/
void cdistr1(double part[], double vtx, double vdx, int npx, int idimp,
             int nop, int nx, int ipbc) {
   distr1_(part,&vtx,&vdx,&npx,&idimp,&nop,&nx,&ipbc);
   return;
}

/*--------------------------------------------------------------------*/
void cmpush1gl(double part[], double complex fx[], double qbm,
               double dt, double *ek, int idimp, int nop, int nx,
               int nxhm, int nxvh) {
   mpush1gl_(part,fx,&qbm,&dt,ek,&idimp,&nop,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cmdpost1gl(double part[], double complex q[], double qm, int nop,
                int idimp, int nx, int nxhm, int nxvh) {
   mdpost1gl_(part,q,&qm,&nop,&idimp,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cffc1initgl(double complex ffc[], double ax, double affp, int nx,
                 int nxhm) {
   ffc1initgl_(ffc,&ax,&affp,&nx,&nxhm);
   return;
}

/*--------------------------------------------------------------------*/
void cpois1gl(double complex q[], double complex fx[], int isign,
              double complex ffc[], double ax, double affp, double *we,
              int nx, int nxhm, int nxvh) {
   pois1gl_(q,fx,&isign,ffc,&ax,&affp,we,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cmefield1gl(double complex fx[], double *gx, double xp, int nx,
                 int nxhm, int nxvh) {
   mefield1gl_(fx,gx,&xp,&nx,&nxhm,&nxvh);
   return;
}
