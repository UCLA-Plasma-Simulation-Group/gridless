/* Wrappers for calling the Fortran routines from a C main program */

#include <complex.h>

double ranorm_();

double randum_();

void distr1_(double *part, double *vtx, double *vdx, int *npx, int *idimp,
             int *nop, int *nx, int *ipbc);

void push1gl_(double *part, double complex *fx, double complex *sctx,
              double *qbm, double *dt, double *ek, int *idimp, int *nop,
              int *nx, int *nxhm, int *nxvh);

void dpost1gl_(double *part, double complex *q, double complex *sctx,
               double *qm, int *nop, int *idimp, int *nx, int *nxhm,
               int *nxvh);

void ffc1initgl_(double complex *ffc, double *ax, double *affp, int *nx,
                 int *nxhm);

void pois1gl_(double complex *q, double complex *fx, int *isign,
              double complex *ffc, double *ax, double *affp, double *we,
              int *nx, int *nxhm, int *nxvh);

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
void cpush1gl(double part[], double complex fx[], double complex sctx[],
              double qbm, double dt, double *ek, int idimp, int nop,
              int nx, int nxhm, int nxvh) {
   push1gl_(part,fx,sctx,&qbm,&dt,ek,&idimp,&nop,&nx,&nxhm,&nxvh);
   return;
}

/*--------------------------------------------------------------------*/
void cdpost1gl(double part[], double complex q[], double complex sctx[],
               double qm, int nop, int idimp, int nx, int nxhm,
               int nxvh) {
   dpost1gl_(part,q,sctx,&qm,&nop,&idimp,&nx,&nxhm,&nxvh);
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
